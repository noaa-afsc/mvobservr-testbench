MakeTrips <- function(nvess, tripsize, tripprob, n_species, mumean, musd, tweedie_power, tweedie_phi){
  # Tweedie power parameter (lambda). 1 < p < 2 is typical for biomass.
  # Tweedie_phi = dispersion parameter for the Tweedie distribution.
  # Higher values create more variance (more rare/dominant species).
  # musd = how similar vessels are to each other
  
  #select the number of trips by vessel
  ntrips <- rnbinom(nvess, size=tripsize, prob=tripprob) 
  
  #select the tweedie mu for catches by vessel & species
  mus <- matrix(nrow=nvess, ncol=n_species)
  for(s in 1:n_species){
    mus[,s] <- rlnorm(nvess, mean=mumean+rnorm(1, sd=mumean*0.3), sd=musd)
    #add a small amount of variability (CV = 0.3) so species are not identical
  }
  mus <- t(apply(mus, 1, sort, decreasing=TRUE)) #sort columns within each row by descending species weights
  
  #make list of vessels with parameters
  #**species
  vesparams <- data.frame(vessnum = 1:nvess) %>% 
    cbind(mus) %>%
    rename_with(.cols=-vessnum, .fn=~paste0("sp_",.x))
  
  trips <- data.frame(vessnum = rep(1:nvess, times=ntrips)) %>%
    inner_join(vesparams, by="vessnum") %>%
    rowwise() %>%
    mutate(across(c(starts_with("sp_")), ~rtweedie(1, power=tweedie_power, mu=.x, phi=tweedie_phi))) %>%
    ungroup() %>%
    mutate(uid = row_number())
  
  return(trips)
}

