AddBias <- function(trips, metric, amount){
  trips <- trips %>%
    mutate_at(metric, function(x) ifelse(.$obs==1, x*(1+amount), x))
  return(trips)
}