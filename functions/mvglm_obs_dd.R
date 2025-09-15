mvglm_obs_dd <- function(x, add_var = NULL, block = NULL, n_permutations, xi.vec = seq(1.2, 1.8, length.out=5), p.anova=0.01) {
  
  # Check dataset for all columns
  cols <- c("uid", "biomass", "species", "observed", add_var, block)
  if(length(setdiff(cols, colnames(x))) > 0 ) {
    stop(paste0(
      "Column(s) ", paste0(setdiff(cols, colnames(x)), collapse = ", "),
      " were missing in the dataset" ))
  }
  
  
  if(!is.numeric(xi.vec)) stop("'xi.vec' needs to be a numeric vector to run the tweedie profile!")
  
  if(!length(unique(x$species))>1) stop("need at least 2 species")
  
  if(!length(unique(x$observed))>1) stop("need both observed and unobserved records")
  
  
  
  # Specify the model formula here
  profile_call <- "biomass ~ species + observed + 1"
  if(!is.null(block)) {
    
    if(length(unique(x[[block]])) > 1) {
        profile_call <- paste(profile_call, "+", block) 
        } else {
        #cat(paste0("'", block, "' had only 1 group and was therefore dropped."))
        block <- NULL
        }
  }
  if(!is.null(add_var)) profile_call <- paste(profile_call, "+", add_var)
  # Prepare the datsaet
  profile_dat <- subset(x, select = cols)
  if(is.null(block)) profile_dat$block <- "block"
  
  ## Determine proper family and link: -------
  t_profile <- tweedie::tweedie.profile(as.formula(profile_call), 
                                        data = profile_dat, 
                                        xi.vec = xi.vec,
                                        do.plot = FALSE, verbose = FALSE, method = "series",
                                        do.ci=FALSE, do.smooth = TRUE, fit.glm = FALSE,
                                        control = list(maxit = 1000, epsilon = 1e-8)) 
  
  
  ## Examine whether a main effect or interaction term for days.fished is needed ----
  main_model <- glm(data = profile_dat,
                    formula = as.formula(profile_call),
                    family = statmod::tweedie(var.power = as.numeric(t_profile$xi.max), link.power = 0))
  
  # If both block and add_var are specified, evaluate the model with an interaction between the two
  if(!is.null(block) & !is.null(add_var)) {
    
    int_call <- gsub(paste("\\+", add_var), paste("*", add_var), profile_call )
    int_model <- glm(data = profile_dat,
                     formula = as.formula(int_call),
                     family = statmod::tweedie(var.power = as.numeric(t_profile$xi.max), link.power = 0))#)
    aov_out <- anova(main_model, int_model, test = "F")
    aov_p <- aov_out$`Pr(>F)`[2]
    
    ### Assign model call to whatever form is depicted by the F test above ----
    
    if(aov_p < p.anova | is.na(aov_p)) {                         
      model_call <- int_call
    } else model_call <- profile_call
    
  } else {
    model_call <- profile_call
    int_model <- NULL
    aov_out <- NULL
    aov_p <- NULL
  }
  
  
  # Prepare the inputs for the mvabund package functions
  
  # Replace 'biomass' with `abund` and remove `species` from the call
  model_call <- gsub("\\bbiomass\\b", "abund", model_call)
  model_call <- gsub("species \\+ ", "", model_call)
  
  # Make a wide-format version of the dataset. There is where 'uid' is important to retain.
  model_dat <- profile_dat %>%
    tidyr::pivot_wider(names_from = species, values_from = biomass, values_fill = 0)
  if(!is.null(block)) { 
    model_dat <- model_dat %>% arrange(BLOCK)
  }
  abund <- mvabund::mvabund(subset(model_dat, select = unique(profile_dat$species)))
  
  # Full Model
  model_formula.full <- as.formula(model_call)
  model_dat.full <- model.frame(model_formula.full, data = model_dat)
  mv_out_full <- manyany_obs(fn = "glm", data = model_dat.full,
                             formula = model_formula.full, #interactions term allows direction to change
                             family = statmod::tweedie(var.power = as.numeric(t_profile$xi.max), link.power = 0),
                             var.power = as.numeric(t_profile$xi.max),
                             control = list(maxit = 10000), get.what = "models")
  
  # Reduced Model,  Removing observed flag
  model_formula.reduced <-  as.formula(gsub("observed \\+", "", model_call))
  #model_formula.reduced <- as.formula("abund ~ 1")
  model_dat.reduced <- model.frame(model_formula.reduced, data = model_dat)
  mv_out_reduced <- manyany_obs(fn = "glm", data = model_dat.reduced, 
                                formula = model_formula.reduced , #Note here the call is updated
                                family = statmod::tweedie(var.power = as.numeric(t_profile$xi.max), link.power = 0), 
                                var.power = as.numeric(t_profile$xi.max),
                                control = list(maxit = 10000), get.what = "models")
  
  #Prepare for anova
  if(!is.null(block)) {
    CTRL <- permute::how(blocks = model_dat[[block]], within = permute::Within(type = "free")) 
  } else CTRL <- permute::how(blocks = model_dat$block, within = permute::Within(type = "free")) 
  perm <- permute::shuffleSet(nset = n_permutations, n = nrow(model_dat), control = CTRL) #one minute to do 5, 13 minutes to do 100.
  
  anova_out <- anova_obs(mv_out_reduced, mv_out_full, bootID = perm, p.uni = "unadjusted")
  
  # Outputs ---- 
  
  return(
    results = list(
      ps = data.frame(t(anova_out$uni.p), p=anova_out$p),
      xi.max = t_profile$xi.max,
      aov_p = aov_p
    )
  )
  
  # return(
  #   list(
  #     tweedie_profile = list(xi.max = t_profile$xi.max,
  #                            L.max = t_profile$L.max,
  #                            x = t_profile$x, 
  #                            y = t_profile$y
  #     ),
  #     model_aov = list(
  #       formula.main = main_model$formula,
  #       formula.int = int_model$formula,
  #       aov = aov_out
  #     ),
  #     model_dat = model_dat,
  #     full_model = list(
  #       formula = mv_out_full$formula,
  #       residuals = mv_out_full$residuals,
  #       linear.predictor = mv_out_full$linear.predictor,
  #       fitted.values = mv_out_full$fitted.values,
  #       fits = mv_out_full$fits
  #     ),
  #     reduced_model = list(
  #       formula = mv_out_reduced$formula,
  #       residuals = mv_out_reduced$residuals,
  #       linear.predictor = mv_out_reduced$linear.predictor,
  #       fitted.values = mv_out_reduced$fitted.values,
  #       fits = mv_out_reduced$fits
  #     ),
  #     results = anova_out,
  #     perm = perm
  #   )
  # )
}

