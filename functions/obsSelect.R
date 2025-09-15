obsSelect <- function(trips, pcov, preference = NULL, prefby = NULL,  ...) {
  if(!is.null(preference) & !is.null(prefby)) {
    #add mu proxy if it doesn't exist (i.e. for real trips)
    if(!"mu" %in% colnames(trips)){
      trips <- trips %>%
        group_by(vessnum) %>%
        mutate(mu = mean(!!sym(prefby))) %>%
        ungroup()
    }
    
    if(preference == "high") {
      rnk <- rank(trips$mu)/nrow(trips) #rank 1 = lowest
    } else { #preference must be low
      rnk <- rank(desc(trips$mu))/nrow(trips) #rank 1 = highest
    }
    rnk <- ifelse(rnk*pcov*2 > 1, 1, rnk*pcov*2)
    #multiply target cov by 2 since rnk averages to 0.5
    #don't let it go over 100% coverage
    trips$obs = rbinom(nrow(trips), 1, rnk)
    
  } else {
    #no deployment effect, simple random selection
    trips$obs <- rbinom(nrow(trips), 1, pcov)
  }
  return(trips)
}
