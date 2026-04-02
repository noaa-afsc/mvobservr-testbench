#' Clean Permutation Function
#' @param data Input dataframe
#' @param metric_var String name of the column to test (e.g., "biomass_total")
#' @param obs_var String name of the observation flag (assumed 0/1)
#' @param n_rep Number of permutations
perm_fxn <- function(data, metric_var, obs_var, n_rep = 1000) {
  
  # 1. Standardize the data
  target <- data[[metric_var]]
  is_obs <- data[[obs_var]]
  
  # 2. Calculate Actual Observed Statistics
  mean_obs <- mean(target[is_obs == 1], na.rm = TRUE)
  mean_unobs <- mean(target[is_obs == 0], na.rm = TRUE)
  observed_diff <- mean_obs - mean_unobs
  
  # 3. Resampling logic
  perm_dist <- replicate(n_rep, {
    shuffled_labels <- sample(is_obs)
    m1 <- mean(target[shuffled_labels == 1], na.rm = TRUE)
    m0 <- mean(target[shuffled_labels == 0], na.rm = TRUE)
    m1 - m0
  })
  
  # 4. Return results as a single-row tibble
  results <- tibble(
    mean_obs = mean_obs,
    mean_unobs = mean_unobs,
    observed_diff = observed_diff,
    p_val = sum(abs(perm_dist) >= abs(observed_diff)) / n_rep
  )
  
  return(results)
}