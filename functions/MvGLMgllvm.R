MvGLMgllvm <- function(dat, n_lv = 0) {
  start_time <- Sys.time()
  
  tryCatch({ 
    R.utils::withTimeout({
      
      # 1. Isolate the response matrix (Y) and predictor (X)
      Y <- dat %>% select(starts_with("sp_"))
      X <- dat %>% select(obs)
      
      # 2. Create silenced versions of both gllvm and anova
      quiet_gllvm <- purrr::quietly(gllvm::gllvm)
      quiet_anova <- purrr::quietly(anova)
      
      # 3. Conditional logic for Tweedie power estimation
      # If no latent variables, estimate it (NULL). Otherwise, lock it at 1.6.
      tweedie_power_arg <- if (n_lv == 0) NULL else 1.6
      
      # 4. Fit the Null Model silently
      m_null_out <- quiet_gllvm(
        y = Y, 
        family = "tweedie", 
        num.lv = n_lv, 
        Power = tweedie_power_arg,             
        control = list(trace = 0, maxit = 1000) 
      )
      m_null <- m_null_out$result 
      
      # Fit the Full Model silently
      m_full_out <- quiet_gllvm(
        y = Y, 
        X = X, 
        formula = ~ obs,
        family = "tweedie", 
        num.lv = n_lv, 
        Power = tweedie_power_arg,             
        control = list(trace = 0, maxit = 1000) 
      )
      m_full <- m_full_out$result
      
      # 5. Tripwire: Force failure if the optimizer didn't converge cleanly
      if (m_full$convergence != TRUE) stop() 
      
      # 6. Run the anova test silently
      anova_out <- quiet_anova(m_null, m_full)
      anova_res <- anova_out$result
      
      # 7. Calculate runtime and conditionally return text or numeric for power
      data.frame(
        p_gllvm            = as.numeric(anova_res$`P.value`[2]),
        pwr_gllvm          = if (n_lv == 0) as.numeric(m_full$Power) else "1.6 (fixed)",
        runtime_secs_gllvm = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      )
    }, timeout = 60, onTimeout = "error") 
    
  }, error = function(e) { 
    # The insurance policy: conditionally matching the NA type to the success type
    data.frame(
      p_gllvm = NA_real_, 
      pwr_gllvm = if (n_lv == 0) NA_real_ else NA_character_, 
      runtime_secs_gllvm = NA_real_
    ) 
  }) 
}
