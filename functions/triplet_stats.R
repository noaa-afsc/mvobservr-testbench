#returns a tibble of results by metric(col)
triplet_stats <- function(pairslongbymetric, bootstrap_reps){
  tryCatch({
    OU <- pairslongbymetric %>% 
      filter(triplet_type=="010") %>% 
      dplyr::select(val) %>%
      unlist() %>%
      unname() %>%
      na.omit()
    UU <- pairslongbymetric %>% 
      filter(triplet_type=="000") %>% 
      dplyr::select(val) %>%
      unlist() %>%
      unname() %>%
      na.omit()
    
    if(length(OU)>0 & length(UU)>0) {
      KSD <- suppressWarnings(ks.test(OU,UU)$statistic)
      KSp <- suppressWarnings(ks.test(OU,UU)$p.value)
    } else {
      KSD <- NaN
      KSp <- NaN
    }
    median_diff <- median(UU, na.rm = TRUE) - median(OU, na.rm = TRUE)
    
    bs <- replicate(bootstrap_reps, {
      OU_bs <- sample(OU,length(OU),replace=TRUE)
      UU_bs <- sample(UU,length(UU),replace=TRUE)
      median(UU_bs, na.rm=TRUE) - median(OU_bs, na.rm=TRUE)
    }
    )
    ci_lo <- quantile(bs,.025, na.rm = TRUE)
    ci_hi <- quantile(bs,.975, na.rm = TRUE)
    
    bhatm = -median_diff
    bhatmlo = -ci_hi
    bhatmhi = -ci_lo
    
    data.frame(cbind(KSD, KSp, median_diff, ci_lo, ci_hi, bhatm, bhatmlo, bhatmhi))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})  
}
