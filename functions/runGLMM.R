## using mgcv ----
# https://library.virginia.edu/data/articles/getting-started-tweedie-models-0
runGLMM <- function(dat, metric) { 
  dat[,"met"] <- dat[,metric] 
  dat$metz <- dat$met
  
  glm_null_mgcv <- gam(metz ~ 1, #gam without a smoother is a glm
                    family=tw(link="log"),
                    data=trips) 
  
  glm_full_mgcv <- update(glm_null_mgcv,  .~. + factor(obs), 
                       family=tw(link="log"),
                       data=trips) 
  
  data.frame(glm_mgcv_p = anova(glm_null_mgcv, glm_full_mgcv, test = "F")$
                             `Pr(>Chi)`[2]) #p-value
#check residuals
#par(mfrow = c(2,2))
#gam.check(glm_fullmgcv)
}