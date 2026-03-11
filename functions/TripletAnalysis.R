TripletAnalysis <- function(trips, metrics, bootstrap_reps=500){
  triplets <- trips %>% 
    dplyr::select(vessnum, obs, all_of(metrics)) %>%
    arrange(vessnum) %>% #and date if had it
    group_by(vessnum) %>%
    mutate(tripnum = row_number()) %>%
    ungroup() %>%
    mutate_all(list(~lag(., 1), #trip before
                    ~lead(., 1)) #trip after
    ) %>%
    #make sure next 2 trips are comparable to this one but not duplicates
    filter(vessnum==vessnum_lag & vessnum==vessnum_lead)
  
  #select just sequence types we want
  triplets <- triplets %>%
    mutate(triplet_type = paste(obs_lag, obs, obs_lead, sep="")) %>%
    filter(triplet_type %in% c("000", "010")) #UUU and UOU
  
  triplets$rand <- runif(nrow(triplets))
  pairs <- bind_rows( 
    triplets %>% filter(rand<=.5) %>% # if <=0.5 then choose trip before
      dplyr::select(-ends_with("_lead")) %>%
      rename_all(list(~str_replace(.,"_lag",".2")))
    ,
    triplets %>% filter(rand>.5) %>% # if >0.5 then choose trip after
      dplyr::select(-ends_with("_lag")) %>%
      rename_all(list(~str_replace(.,"_lead",".2")))
  )
  
  #remove if used again
  pairs <- pairs %>%
    mutate(trip1 = ifelse(tripnum < tripnum.2, tripnum, tripnum.2),
           trip2 = ifelse(tripnum > tripnum.2, tripnum, tripnum.2)) %>%
    arrange(vessnum, trip1, trip2, rand) %>% #sort by rand so if A-B and B-A the lower rand will be prioritized
    group_by(vessnum) %>%
    mutate(prevtrip2 = lag(trip2,1, order_by = trip2)) %>% #what was the later trip used in the previous pair? 
    filter(is.na(prevtrip2) | trip1 > prevtrip2) %>%
    ungroup()
  
  #create new table for vessel means for each metric in "metrics"
  trips.means <- trips %>%
    filter(obs == 0) %>%
    group_by(vessnum) %>%
    summarize_at(metrics,list(m = ~ mean(.)), .groups="drop") 
  #returns col_m if metrics > 1, just m if metrics == 1
  if(length(metrics)==1){
    colnames(trips.means)[colnames(trips.means)=="m"] <- paste(metrics,"_m",sep="")
  }
  
  pairs.means <- dplyr::inner_join(pairs, trips.means, by="vessnum")
  for (metric in metrics) {
    metric_delta <- paste(metric,"_delta",sep="")
    metric2 <- paste(metric,".2",sep="")
    metric_m <- paste(metric,"_m",sep="")
    pairs.means[,metric_delta] <- (pairs.means[,metric] -
                                  pairs.means[,metric2])/pairs.means[,metric_m]
  }
  
  #drop vessels that don't have at least 3 sequences
  under <- pairs.means %>% 
    group_by(vessnum) %>% 
    dplyr::summarize(n=n(), .groups="drop") %>% 
    filter(n<3) %>% 
    ungroup()
  if(nrow(under) > 0) {
    pairs.means <- pairs.means %>% 
      filter(!vessnum %in% under$vessnum)
  }
  
  vessels_retained = pairs.means %>%
    dplyr::summarize(ves_ret = n_distinct(vessnum)) 
  
  ##stats
  pairs_long <- pairs.means %>%
    dplyr::select(triplet_type, ends_with("_delta")) %>%
    pivot_longer(cols=ends_with("_delta"),
                 names_to = "metric",
                 values_to = "val",
                 names_pattern = "(.*?)_delta"
    ) %>%
    filter(metric %in% metrics)
  
  metric_stats <- pairs_long %>%
    split(.$metric) %>% 
    map(~triplet_stats(.x, bootstrap_reps)) %>%
    list_rbind(names_to="metric")
  
  if (nrow(metric_stats) < 1) {
    metric_stats <- data.frame(KSD = NaN, KSp = NaN, median_diff = NaN, ci_lo = NaN, ci_hi = NaN)
  }
  
  return(cbind(metric_stats, vessels_retained))
  
} #end function TripletAnalysis
