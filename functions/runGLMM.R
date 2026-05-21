# runGLMM <- function(dat, metric) {
#   dat[,"met"] <- dat[,metric]
#   dat$metz <- dat$met
#   
#   glm1 <- gam(metz ~ 1+factor(obs), #gam without a smoother is a glm
#               family=tw(link="log"),
#               data=dat) 
#   
#   glm1u <- update(glm1, metz~1) #remove observer effect, intercept only
#   AIC1 <-  AIC(glm1)
#   AIC2 <-  AIC(glm1u)
#   AICd <- AIC2 - AIC1
#   
#   glmconv <- ifelse(glm1$converged== TRUE,1,0)
#   
#   data.frame(cbind(AICd, glmconv))
# }

## using mgcv ----
# https://library.virginia.edu/data/articles/getting-started-tweedie-models-0
runGLMM <- function(dat, metric) { #Copied from existing
  dat[,"met"] <- dat[,metric] #TODO - Why all the assignment here if we just use 'obs' later?
  dat$metz <- dat$met
  
  glm_null_mgcv <- gam(metz ~ 1, #gam without a smoother is a glm
                    family=tw(link="log"),
                    data=trips) 
  
  glm_full_mgcv <- update(glm_nullmgcv,  .~. + factor(obs), #TODO - I dont see why we use obs explicitly here.
                       family=tw(link="log"),
                       data=trips) 
  
  data.frame(glm_mgcv_p = anova(glm_null_mgcv, glm_full_mgcv, test = "F")$
                             `Pr(>Chi)`[2]) #p-value
#check residuals
#par(mfrow = c(2,2))
#gam.check(glm_fullmgcv)
}