#' Blurb of what this function does
#' 
#' @param bootRows description.
#' @param ... other arguments to pass onward
#' @return what it returns.
#' @export
#' @examples
#' \dontrun{
#' out <- bootAnova_obs(
#'   bootRows = 1:nBoot, 
#'   bootID = bootID,
#'   block = block, 
#'   blockIDs = blockIDs, 
#'   n.rows = n.rows,
#'   n.vars = n.vars,
#'   replace = replace,
#'   unlistIDs = unlistIDs,
#'   n.levels = n.levels,
#'   object1 = object1,
#'   object2 = object2,
#'   qfn = qfn,
#'   nCores = 1)
#' }
bootAnova_obs_dd <- function(bootRows, ...){
  
  if(F){
    # For testing bootAnova_obs
    bootRows <- 1:nBoot; argList = list(
      bootID = bootID, block = block, blockIDs = blockIDs,
      n.rows = n.rows, n.vars = n.vars,
      replace = replace, unlistIDs = unlistIDs, n.levels = n.levels,
      object1 = object1, object2 = object2, qfn = qfn, nCores = 1
    )
  }
  
  nBooti = length(bootRows)
  
  # Don't run this in testing
  dots = list(...)
  if ( any(names(dots) == "nCores") ) { # if nCores=1, take ... and call it argList, to match parLapply call
    if(dots$nCores == 1) argList <- list(...)
  }
  #initialise parameters for bootstrapping
  
  #require(mvabund)
  yMat = matrix(NA, argList$n.rows, argList$n.vars)
  #next two lines no longer correct, requires matrix input to use formula
  #  if(argList$object1$family[[1]]$family=="ordinal") 
  #    yMat=data.frame(yMat)
  argList$object1$call$get.what <- "none" #to avoid wasting time computing residuals etc when resampling
  argList$object2$call$get.what <- "none" #to avoid wasting time computing residuals etc when resampling
  
  if(is.null(argList$bootID)) boot.Resamp <- rep(NA,argList$n.rows)
  
  # find where in object1$call and object2$call the response matrix is so it can be replaced with bootstrapped version
  mf              <- argList$object2$model
  nameOfResponse  <- as.character(argList$object2$formula[[2]])
  whichIsResponse <- which(names(mf) == nameOfResponse)
  
  #suggestions from chatgpt
  ## ---------------------------------------------
  ## Precompute invariants OUTSIDE the bootstrap
  ## ---------------------------------------------
  
  # Precompute pnorm(residuals)
  pnorm_res <- pnorm(argList$object1$residuals)
  
  # Prebuild z-templates with q removed
  z_templates <- lapply(argList$object1$params, function(z) {
    z$q <- NULL
    z
  })
  
  # Shortcuts
  qfun_list <- argList$qfn
  n_species  <- length(z_templates)
  
  # Prebuild model environment for eval()
  model_env <- new.env(parent = environment(argList$object1$formula))
  
  model_name1 <- deparse(argList$object1$call$data)
  model_name2 <- deparse(argList$object2$call$data)
  
  # Formulas (used but never change)
  model_formula.reduced <- argList$object1$formula
  model_formula.full    <- argList$object2$formula
  
  # t_profile is stored in the formula environment
  t_profile <- get("t_profile",
                   envir = attributes(argList$object1$formula)$.Environment)
  
  ## ---------------------------------------------
  ## Bootstrap loop
  ## ---------------------------------------------
  
  
  print(paste0(now(), ": starting loop"))
  out.lst <- lapply(seq_len(nrow(argList$bootID)), function(i) {
    
    x <- argList$bootID[i, ]   # bootstrap row (indices)
    
    # Build list-of-vectors of pnorm residuals at bootstrap indices
    # This replaces apply(..., simplify=FALSE)
    y_list <- lapply(seq_len(ncol(pnorm_res)), function(j) pnorm_res[x, j])
    
    # Compute yMat species-by-species with qtweedie_obs
    yMat_list <- vector("list", n_species)
    
    for (j in seq_len(n_species)) {
      
      z_temp <- z_templates[[j]]   # template already has q=NULL
      z_temp$p <- y_list[[j]]      # insert bootstrap p-values
      
      # Call qtweedie_obs or whatever function lives in qfn
      yMat_list[[j]] <- do.call(qfun_list[[j]], z_temp)
      #yMat_list[[j]] <- do.call("qtweedie_obs", z_temp)
      #print(paste0(now(), ": qfun ", i, ", ", j, ""))
    }
    
    # cbind the list → matrix
    yMat <- do.call(cbind, yMat_list)
    
    ## ---------------------------------
    ## Zeroton detection
    ## ---------------------------------
    if (argList$object1$family[[1]]$family == "ordinal") {
      is.zeroton <- apply(yMat, 2, function(col) length(table(col)) == 1)
    } else {
      is.zeroton <- colSums(yMat, na.rm = TRUE) == 0
    }
    
    # Update mf with filtered yMat
    mf[[whichIsResponse]] <- yMat[, !is.zeroton, drop = FALSE]
    
    ## ---------------------------------
    ## Inject updated data into model env
    ## ---------------------------------
    model_env[[model_name1]] <- mf
    model_env[[model_name2]] <- mf
    
    
    ## ---------------------------------
    ## Refit reduced & full models
    ## ---------------------------------
    stati.i  <- NA
    statj.ii <- rep(NA, times = argList$n.vars)
    
    if (sum(!is.zeroton) > 0) {
      
      ft.1i <- tryCatch(eval(argList$object1$call, envir = model_env),
                        error = function(e) NA)
      
      ft.2i <- tryCatch(eval(argList$object2$call, envir = model_env),
                        error = function(e) NA)
      
      if (length(ft.1i) > 1 && length(ft.2i) > 1) {
        
        # Compute per-species LRTs
        statj.ii[!is.zeroton] <-
          2 * (logLik(ft.2i) - logLik(ft.1i))
        
        stati.i <- sum(statj.ii, na.rm = TRUE)
      }
      
    } else {
      stati.i <- 0
    }
    
    list(
      statj.ii = statj.ii,
      stati.i  = stati.i
    )
  })
  print(paste0(now(), ": ending loop"))
  
  # Create final list of outputs
  list(
    stati.i = sapply(out.lst, "[[", "stati.i"),
    statj.ii = sapply(out.lst, "[[", "statj.ii")
  )
  
} 
