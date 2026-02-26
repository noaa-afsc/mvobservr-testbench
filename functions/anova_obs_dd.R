#' Blurb of what this function does
#' 
#' @param object description
#' @param ... other arguments to pass onward
#' @param nBoot description
#' @param p.uni description
#' @param block description
#' @param nCores description
#' @param bootID description
#' @param replace description
#' @return what it returns
#' @export
#' @examples
#' \dontrun{
#' out <- anova_obs(bootRows = 1:nBoot, 
#'                  bootID = bootID,
#'                  block = block, 
#'                  blockIDs = blockIDs, 
#'                  n.rows = n.rows,
#'                  n.vars = n.vars,
#'                  replace = replace,
#'                  unlistIDs = unlistIDs,
#'                  n.levels = n.levels,
#'                  object1 = object1,
#'                  object2 = object2,
#'                  qfn = qfn,
#'                  nCores = 1)
#' }

anova_obs_dd = function(object, ..., nBoot = 99, p.uni="none", block = object1$block, 
                     nCores = 1, bootID = NULL, replace=TRUE) {
  
  # Modified from mavbund::manyany()
  
  if(F){
    # For testing for anova_obs
    object <- mv_out_reduced
    object2 <- mv_out_full; 
    bootID <- perm; nBoot <- 99; p.uni <- "unadjusted"; nCores <- 1; replace <- TRUE
    dots <- list(full_model = mv_out_full); ndots <- 1
  }
  
  # set default na.action to exclude in order to not change dimensions of anything when NAs are present
  naOptInit <- getOption("na.action")
  options(na.action <- "na.exclude")
  
  object1 <- object 
  
  # get object 2
  # Dont run this in testing
  dots <- list(...)
  ndots <- length(dots)
  
  fndObj2 <- FALSE
  if (ndots == 0) {
    stop("missing a second manyany object")
  } else {
    if (ndots > 1) warning("current version only compares two manyany objects")
    for (i in 1:ndots) {
      if (any(class(dots[[i]]) == "manyany")) {
        object2 <- dots[[i]]
        fndObj2 <- TRUE
        break
      }
    }
    if (!fndObj2) stop("cannot find object 2")
  }
  
  if(any(names(object1$call) == "composition")) {
    if(object1$call$composition == TRUE) { #recode so that it fits compositional models as univariate, to save time and fuss
      object1$call$formula <- object1$formula
      object2$call$formula <- object2$formula
      object1$call$data <- object1$model
      object2$call$data <- object2$model
      object1$residuals <- as.matrix(c(object1$residuals))
      object1$call$composition <- FALSE
      object2$call$composition <- FALSE
    }
  }
  
  #DW, 18/1/18: check for same composition arguments in each call
  if(all(names(object1$call) != "composition")) object1$call$composition <- FALSE
  if(all(names(object2$call) != "composition")) object2$call$composition <- FALSE
  
  if(object1$call$composition != object2$call$composition) stop("Sorry, either both manyany objects will need to be compositional, or neither")
  
  n.rows <- dim(object1$resid)[1]
  n.vars <- dim(object1$resid)[2]
  
  qfn <- rep(NA,n.vars)
  for(i.var in 1:n.vars){
    if(grepl("egative",object1$family[[i.var]]$family) || object1$family[[i.var]]$family == "negbinomial") qfn[i.var] <- "qnbinom"
    if(object1$family[[i.var]]$family == "poisson") qfn[i.var] <- "qpois"
    if(object1$family[[i.var]]$family == "binomial"){
      qfn[i.var] <- "qbinom"
      warning("The binomial option of manyany currently assumes you have presence/absence data")
    } 
    if(object1$family[[i.var]]$family == "gaussian") qfn[i.var] <- "qnorm"  
    if(object1$family[[i.var]]$family == "Tweedie") qfn[i.var] <- "qtweedie_obs_dd"   # *[NOTE] Using custom qtweedie_obs*
    if(object1$family[[i.var]]$family == "ordinal") qfn[i.var] <- "qordinal"
  }
  
  if(is.null(bootID) == FALSE) {
    bootID <- as.matrix(bootID)
    if(dim(bootID)[2] != n.rows) stop("Number of columns of bootID must match number of rows in data")
    nBoot = dim(bootID)[1] #overwriting nBoot with value implied by user-entered ID matrix
    block = NULL #overwriting previous value for block
    #print("User-entered bootID matrix will be used to generate bootstrap samples")
  }  
  n.levels <- n.rows; unlistIDs = NULL #objects needed for block resampling otherwise ignorable 
  blockIDs <- NULL
  
  if(is.null(block) == FALSE) {
    tb <- table(block)
    n.levels <- length(tb)
    if(any(tb != n.rows/n.levels)) {   
      print(tb) 
      stop("Sorry, block needs to be a balanced factor - same number of rows for each level")
    } else {
      blockIDs <- vector("list",n.levels)
      for(i.level in 1:n.levels) {
        blockIDs[[i.level]] <- which(block == names(tb)[i.level])
      }
      unlistIDs <- unlist(blockIDs) #needed to match each resampled observation with its correct location
    }
  }
  #get observed test stat
  #ft.1i=eval(object1$call) #this call not needed but good to check that eval is working, compare logLik to logLik(object1)
  statj <- 2 * ( logLik(object2) - logLik(object1) )
  stat <- sum(statj)
  
  if(nCores > 1) {
    nBooti <- ceiling(nBoot/nCores)
    
    # construct a list which says which rows of bootID to use in which cluster: only needed when bootID provided
    bootRows <- vector(length = nCores, mode = "list")
    for(iCore in 1:nCores) {
      bootRows[[iCore]] <- 1:nBooti + nBooti * (iCore - 1)
      bootRows[[nCores]] <- pmin( bootRows[[nCores]], nBoot )
    }
    
    # set up clusters, pass through arguments
    cl <- parallel::makeCluster(nCores)
    argList <- list(bootID = bootID, 
                    block = block, 
                    blockIDs = blockIDs, 
                    n.rows = n.rows, 
                    n.vars = n.vars, 
                    replace = replace, 
                    unlistIDs = unlistIDs, 
                    n.levels = n.levels, 
                    object1 = object1, 
                    object2 = object2, 
                    qfn = qfn)
    
    parallel::clusterExport(cl,"argList", envir = environment())
    #clusterExport(cl,c("nBooti","bootID","block","n.rows","n.vars","replace","unlistIDs","n.levels","object1","object2","qfn"), envir=environment())
    out <- parallel::parLapply(cl, bootRows, bootAnova_obs_dd)
    #why not clusterapply??
    parallel::stopCluster(cl)
    
    # store results in vectors/matrices not lists
    stat.i  <- rep(NA, nBooti * nCores)
    statj.i <- matrix(NA, n.vars, nBooti * nCores)
    for(i.core in 1:nCores) {
      stat.i[(i.core-1) * nBooti + 1:nBooti] <- out[[i.core]]$stati.i
      statj.i[,(i.core-1) * nBooti + 1:nBooti] <- out[[i.core]]$statj.ii
    }  
    stat.i <- stat.i[1:nBoot]
    statj.i <- statj.i[,1:nBoot]
  } else {
    # bootAnova_obs runs here
    print(paste0(now(), ": starting bootAnova"))
    out <- bootAnova_obs_dd(bootRows = 1:nBoot, 
                         bootID = bootID,
                         block = block, 
                         blockIDs = blockIDs, 
                         n.rows = n.rows,
                         n.vars = n.vars,
                         replace = replace,
                         unlistIDs = unlistIDs,
                         n.levels = n.levels,
                         object1 = object1,
                         object2 = object2,
                         qfn = qfn,
                         nCores = 1)
    print(paste0(now(), ": ending bootAnova"))
    stat.i <- out$stati.i
    statj.i <- out$statj.ii
  }
  if(n.vars > 1) {
    dimnames(statj.i)[[1]] <- dimnames(object1$residuals)[[2]]
    # CF Suggested correction, and apply below
    # If length of stat.i = nBoot, my correction will work
    # nBoot is the number of rows in bootID (l. 81)
    # stati is from stati.i which is a rep(NA of nBooti). 
    # nBooti is length(bootRows).
    # bootRows = 1:nBoot.  Yes, they are equal.
    nBootSuccess <- as.numeric(length(stat.i[!is.na(stat.i)])) #Number of successful bootstraps
    
    # p = (1 + sum(stat.i > stat - 1.e-8) ) / (nBoot + 1)  OLD.  This seems strange at low nBoot numbers.  
    #                                                           Consider nBoot = 1, p = 0.5!  Need at least 100 nBoot!
    # The old way will no longer work because for failed model bootstraps we have NAs in stat.i.  
    # So although we did two bootstraps, only one has a stat.i
    p <- (1 + sum(stat.i > (stat - 1.e-8), na.rm = TRUE) ) / (nBootSuccess + 1) 
    # The na.rm = T above converts stat.i to successful model runs.
    
    if(length(statj) > 1) {
      # pj = (1 + apply(statj.i > statj-1.e-8, 1, sum) ) / (nBoot + 1)
      pj <- (1 + apply(statj.i > (statj - 1.e-8), 1, sum, na.rm = TRUE) ) / (nBootSuccess + 1) 
      # Same as p applied to each species (row).
    }
  } else {
    # pj = (1 + sum(statj.i > statj-1.e-8) ) / (nBoot + 1) OLD
    pj <- (1 + sum(statj.i > (statj - 1.e-8), na.rm = TRUE) ) / (nBootSuccess + 1)
  }
  
  class(stat.i) <- "numeric"
  
  if(p.uni == "unadjusted") {
    result <- list(stat = stat,
                   p = p,
                   uni.test = statj,
                   uni.p = pj,
                   stat.i = stat.i,
                   statj.i = statj.i,
                   p.uni = p.uni,
                   nBoot = nBoot,
                   nBootSuccess = nBootSuccess) # Adding Successes to output.
  }
  
  if(p.uni == "none") {
    result <- list(stat = stat,
                   p = p,
                   stat.i = stat.i,
                   p.uni = p.uni,
                   nBoot = nBoot,
                   nBootSuccess = nBootSuccess) # Adding Successes to output.
  }
  
  options(na.action = naOptInit) #restore previous default for na.action
  
  class(result) <- "anova_obs"
  result
  
}
