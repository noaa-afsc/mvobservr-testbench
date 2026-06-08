MvGLMgllvm <- function(dat) {
  start_time <- Sys.time()
  
  tryCatch({
    # 1. Isolate the response matrix (Y) and predictor (X)
    Y <- dat %>% select(starts_with("sp_"))
    X <- dat %>% select(obs)
    
    # 2. Create silenced versions of both gllvm and anova.  Quirk of the package compared to glmm.
    quiet_gllvm <- purrr::quietly(gllvm)
    quiet_anova <- purrr::quietly(anova)
    
    # 3. Fit the Null Model silently
    m_null_out <- quiet_gllvm(
      y = Y, 
      family = "tweedie", 
      num.lv = 0, 
      Power = NULL,
      control = list(trace = 0) #another silencing tool via gemini.
    )
    m_null <- m_null_out$result 
    
    # Fit the Full Model silently
    m_full_out <- quiet_gllvm(
      y = Y, 
      X = X, 
      formula = ~ obs,
      family = "tweedie", 
      num.lv = 0, 
      Power = NULL,
      control = list(trace = 0)
    )
    m_full <- m_full_out$result
    
    # 4. Tripwire: Force failure if the optimizer didn't converge cleanly
    if (m_full$convergence != TRUE) stop() 
    
    # 5. Run the anova test silently to protect the progress bar
    anova_out <- quiet_anova(m_null, m_full)
    anova_res <- anova_out$result
    
    # 6. Calculate runtime and return metrics
    data.frame(
      p_gllvm            = as.numeric(anova_res$`P.value`[2]),
      runtime_secs_gllvm = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
    
  }, error = function(e) {
    # If any error or convergence failure occurs, return entirely blank rows
    data.frame(p_gllvm = NA_real_, runtime_secs_gllvm = NA_real_)
  })
}


