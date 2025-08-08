mvglm_obs_dd <- function(dat, n_permutations, xi.vec = seq(1.2, 1.8, length.out=5)) {
  
  # Specify the model formula here
  model_call <- "biomass ~ species + observed + BLOCK+days.fished"
  
  ## Determine proper family and link: -------
  t_profile <- tweedie::tweedie.profile(as.formula(model_call), #Formula from fit2 below using cpglm
                                        data = dat, 
                                        #xi.vec = xi.vec,
                                        do.plot = FALSE, verbose = FALSE, method = "series",
                                        do.ci=FALSE, do.smooth = TRUE, fit.glm = FALSE,
                                        control = list(maxit = 1000, epsilon = 1e-8)) 
  #cat(paste("xi.max:", t_profile$xi.max, "\n"))
  
  # Prepare the inputs for the mvabund package functions
  
  # Replace 'biomass' with `abund` and remove `species` from the call
  model_call <- gsub("\\bbiomass\\b", "abund", model_call)
  model_call <- gsub("species \\+ ", "", model_call)
  
  # Make a wide-format version of the dataset. There is where 'uid' is important to retain.
  model_dat <- dat %>%
    tidyr::pivot_wider(names_from = species, values_from = biomass, values_fill = 0) %>% arrange(BLOCK)
  abund <- mvabund::mvabund(subset(model_dat, select = unique(dat$species)))
  
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
  model_dat.reduced <- model.frame(model_formula.reduced, data = model_dat)
  mv_out_reduced <- manyany_obs(fn = "glm", data = model_dat.reduced, 
                                formula = model_formula.reduced , #Note here the call is updated
                                family = statmod::tweedie(var.power = as.numeric(t_profile$xi.max), link.power = 0), 
                                var.power = as.numeric(t_profile$xi.max),
                                control = list(maxit = 10000), get.what = "models")
  
  #Prepare for anova
  CTRL <- permute::how(blocks = model_dat$BLOCK, within = permute::Within(type = "free")) 
  perm <- permute::shuffleSet(nset = n_permutations, n = nrow(model_dat), control = CTRL, quietly=TRUE) 
  
  anova_out <- anova_obs(mv_out_reduced, mv_out_full, bootID = perm, p.uni = "unadjusted")
  
  # Outputs ---- 
  
  return(
    results = list(
      uni.p = data.frame(t(anova_out$uni.p)),
      p = anova_out$p
    )
  )
}