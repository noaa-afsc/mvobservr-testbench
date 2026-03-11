runGLMM <- function(dat, metric) {
  dat[,"met"] <- dat[,metric]
  dat$metz <- dat$met + ifelse(dat$met==0,jitter(.001),0) #remove 0s for log link
  
  #filter out vessels with fewer than 2 observed or 2 unobserved trips
  dat2 <- dat %>%
    group_by(vessnum) %>%
    mutate(nobs=sum(obs==1), nunobs=sum(obs==0)) %>%
    filter(nobs>=2 & nunobs>=2)
  dat2$vessnum <- factor(dat2$vessnum)

  glmm1 <- lmer(metz~1+factor(obs) + (1|vessnum),
          data=transform(dat2, metz = log(metz)),
          REML = FALSE)
  glmm1.s <- summary(glmm1)
  glmm1u <- update(glmm1, metz~1 + (1|vessnum)) #remove observer effect
  glmm1u.s <- summary(glmm1u)
  AIC1 <-  glmm1.s$AICtab[1]
  AIC2 <-  glmm1u.s$AICtab[1]
  AICd <- AIC2 - AIC1
  factest <- glmm1.s$coef[2,1]
  factse <- glmm1.s$coef[2,2]
  glmmt <- glmm1.s$coef[2,3]
  glmmnves <- length(levels(dat2$vessnum))
  glmmntrp <- nrow(dat2)
  glmmconv <- ifelse(glmm1.s$optinfo$conv$opt == 0,1,0) #opt=0 means OK = converged
  
  bhate = exp(factest) - 1
  bhatelo = exp(factest-1.96*factse) - 1
  bhatehi = exp(factest+1.96*factse) - 1
  
  data.frame(cbind(AICd, glmmnves, glmmntrp, glmmconv, bhate, bhatelo, bhatehi))
}
