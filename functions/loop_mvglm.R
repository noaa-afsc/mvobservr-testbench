loop_mvglm <- function(dat, covrate, nperm, bias=NA, add_var = "days.fished", block = "BLOCK"){
  idat <- dat %>% 
    mutate(observed = ifelse(runif(nrow(.)) <= covrate, 'Y', 'N')) %>% 
    pivot_longer(cols=-c(uid, days.fished, observed, BLOCK),
                 names_to = 'species', values_to = 'biomass') 
  
  if(!is.na(bias)) {
    #add bias to all species equally for observed trips only
    idat <- idat %>% 
      mutate(biomass =  biomass * ifelse(observed == 'Y',(1+bias), 1))
  }

  mvglm_obs_dd(idat, add_var = add_var, block = "BLOCK", n_permutations=nperm) 
}
