runGLMM <- function(dat, metric) {
  dat[,"met"] <- dat[,metric]
  dat$metz <- dat$met + ifelse(dat$met==0,jitter(.001),0) #remove 0s for log link
  
  glm1 <- gam(metz ~ 1+factor(obs), #gam without a smoother is a glm
              family=tw(link="log"),
              data=dat) 
  
  glm1u <- update(glm1, metz~1) #remove observer effect, intercept only
  AIC1 <-  AIC(glm1)
  AIC2 <-  AIC(glm1u)
  AICd <- AIC2 - AIC1
  
  glmconv <- ifelse(glm1$converged== TRUE,1,0)
  
  data.frame(cbind(AICd, glmconv))
}
