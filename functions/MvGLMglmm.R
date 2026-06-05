MvGLMglmm <- function(dat) {
  start_time <- Sys.time()
  
  tryCatch({
    # 1. Pivot and format data
    long <- dat %>% 
      pivot_longer(cols = starts_with("sp_"), names_to = "species", values_to = "biomass") %>%
      mutate(species = factor(species), obs = factor(obs))
    
    # 2. Fit models 
    # https://www.jstatsoft.org/article/view/v112i01 for source doc and example on how to use with boot and nested models
    # https://journal.r-project.org/articles/RJ-2017-066/index.html is source for original package article and has run times.
    m_null <- glmmTMB(biomass ~ 0 + species, data = long, family = tweedie(link = "log")) #default is to set these up serially.
    m_full <- glmmTMB(biomass ~ 0 + species + species:obs, data = long, family = tweedie(link = "log"))
    
    # 3. Tripwire: Force failure if the optimizer didn't converge cleanly
    if (m_full$fit$convergence != 0) stop() 
    
    # 4. If successful, calculate runtime and return metrics
    data.frame(
      p_glmm       = anova(m_null, m_full)$`Pr(>Chisq)`[2],
     # pwr_glmm     = unname(family_params(m_full)), #tweedie power
      runtime_secs_glmm = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
    
  }, error = function(e) {
    # If any error or convergence failure occurs, return entirely blank rows
    data.frame(p_glmm = NA_real_, #pwr_glmm = NA_real_, 
               runtime_secs_glmm = NA_real_)
  })
}

