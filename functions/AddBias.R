AddBias <- function(trips, metric, amount){
  trips <- trips %>%
    mutate_at(metric, function(x) ifelse(.$obs==1, x*(1+amount), x)) %>%
    rowwise() %>%
    mutate(biomass_total = sum(c_across(starts_with("sp_")))) %>%
    ungroup() 
  return(trips)
}