#takes a dataframe of trips and returns a one row tibble with test statistics and significance values
ObserverEffectStats <- function(trips, metrics, nperm_mvglm, nboot_triplet, nperm_adon){
  results_list <- setNames(lapply(metrics, function(x) data.frame()), metrics)
  
  for (metric in metrics) {
    #loop over functions that can't take more than one metric at a time
    r_desc <- getDescriptiveStats(trips, metric)
    r_tandf <- runTandFtests(trips, metric)
    r_glmm <- runGLMM(trips, metric)
    results_list[[metric]] <- remove_rownames(cbind(r_desc, r_tandf, r_glmm))
  }
  
  res <- results_list %>%
    list_rbind(names_to = "metric")
  
  triplet_res <- TripletAnalysis(trips, metrics, bootstrap_reps = nboot_triplet)
  
  mvglmobs_res <- trips %>%
    mutate(BLOCK = "block") %>% #needed for mvglm_obs_dd
    mutate(observed = ifelse(obs==1, 'Y', 'N')) %>%
    pivot_longer(cols=starts_with("sp_"),
                 names_to = 'species', values_to = 'biomass') %>%
    mvglm_obs_dd(block = "BLOCK", n_permutations=nperm_mvglm) %>%
    pluck("ps") %>%
    rename(biomass_total = p) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "mvglm_p") 
  
  #permanova, one value for all species vs. obs factor
  Y <- trips %>% dplyr::select(starts_with("sp_"))
  adon <- vegan::adonis2(Y ~ factor(obs), data = test_trips, permutations = nperm_adon, method="euclidean")
  permanova <- data.frame(metric = "biomass_total", perma_p = adon$`Pr(>F)`[1])
  
  
  res %>%
    left_join(triplet_res, by="metric") %>%
    left_join(mvglmobs_res, by="metric") %>%
    left_join(permanova, by="metric")
}
