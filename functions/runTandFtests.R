runTandFtests <- function(dat, metric){
    dat[,"met"] <- dat[,metric]
   
    ttest <- t.test(met~factor(obs), data=dat)
    tstat <- ttest$statistic
    tmu <- ttest$estimate[1] #unobs
    tmo <- ttest$estimate[2] #obs
    tse <- ttest$stderr #CI = tval +/- 1.96*tse
    tp <- ttest$p.value

    Ftest <- anova(lm(met~factor(obs), data=dat))
    Fstat <- Ftest$`F value`[1]
    Fp <- Ftest$`Pr(>F)`[1]
    
    bhatt = tmo/tmu - 1
    bhattlo = (tmo-1.96*tse)/tmu - 1
    bhatthi = (tmo+1.96*tse)/tmu - 1
    
    data.frame(cbind(tstat, tmu, tmo, tse, tp, bhatt, bhattlo, bhatthi, Fstat, Fp))
}
