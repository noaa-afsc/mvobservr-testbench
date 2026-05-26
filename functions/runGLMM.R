## using mgcv ----
# https://library.virginia.edu/data/articles/getting-started-tweedie-models-0
runGLMM <- function(dat, metric) { 
  
    dat[,"met"] <- dat[,metric] 
    dat$metz <- dat$met
    
    glm_full_mgcv <- gam(metz ~ factor(obs), #gam without a smoother is a glm
                       family=tw(link="log"),
                       data=dat)   
    
    theta_est <- glm_full_mgcv$family$getTheta(TRUE) #gets the power from the full model
  
    glm_null_mgcv <- gam(metz ~ 1, 
                    family=tw(theta = theta_est, link="log"), #ensuring we are using the correct and same power from tweedie
                    data=dat) 
  
  # 2. Check for convergence failure immediately
  if (!glm_null_mgcv$converged || !glm_full_mgcv$converged) {
    return(data.frame(glm_mgcv_p = NA)) #did not converge
  }
  
  anova_res <- anova(glm_null_mgcv, glm_full_mgcv, test = "Chisq")
  dev_diff <-  anova_res$Deviance[2]
  if (!is.na(dev_diff) && dev_diff <= 0) {
    return(data.frame(glm_mgcv_p = 1.0)) #numerical noise
  }
  # 4. If it passes the checks, run the ANOVA safely
  data.frame(glm_mgcv_p = anova_res$`Pr(>Chi)`[2]) #p-value

#check residuals
#par(mfrow = c(2,2))
#gam.check(glm_fullmgcv)
}