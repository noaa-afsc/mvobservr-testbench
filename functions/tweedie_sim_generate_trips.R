tweedie_sim_generate_trips <- function(n_species, n_samples_per_level, tolerance, tweedie_power, phi, fixed_total_biomass, nvess, tripsize, tripprob, mumean, musd, target_bp_levels, max_attempts, seed = NA){
  
# for reproducibility
  set.seed(seed)

# --- Data Storage ---

# Create a list to act as "buckets" for each target level.
results_list <- setNames(lapply(target_bp_levels, function(x) data.frame()), target_bp_levels)
# Create a counter to track how many samples we've found for each level.
sample_counts <- setNames(rep(0, length(target_bp_levels)), target_bp_levels)

# --- Simulation Loop ---

message("Generating populations and sorting into target levels...")
attempts <- 0

# Calculate the total number of samples needed for the progress bar.
total_needed <- n_samples_per_level * length(target_bp_levels)
# Initialize the progress bar.
pb <- txtProgressBar(min = 0, max = total_needed, style = 3)

# Loop until all levels have the desired number of samples.
while(any(sample_counts < n_samples_per_level) && attempts < max_attempts) {
  attempts <- attempts + 1
  
  # Step 1: Generate a single realistic community.
  #mu_values <- rlnorm(n_species, meanlog = 2, sdlog = 1.5)
  #biomass <- rtweedie(n_species, p = tweedie_power, mu = mu_values, phi = phi)
  
  gen_dat <- MakeTrips(nvess = nvess, tripsize=tripsize, 
                                   tripprob=tripprob, n_species = n_species, mumean=mumean, 
                                   musd=musd, tweedie_power = tweedie_power, tweedie_phi = phi) %>%
    rowwise() %>%
    mutate(biomass_total = sum(c_across(starts_with("sp_")))) %>%
    ungroup() 

  biomass <- gen_dat %>% 
    summarize(across(starts_with("sp_"), sum))
  
  if (sum(biomass) == 0) next
  
  # Step 2: Calculate its Berger-Parker index.
  proportions <- biomass / sum(biomass)
  bp_index <- max(proportions)
  
  # Step 3: Find the closest target level.
  closest_target_index <- which.min(abs(target_bp_levels - bp_index))
  closest_target_bp <- target_bp_levels[closest_target_index]
  bp_key <- as.character(closest_target_bp)
  
  # Step 4: Check if it's a match and if that bucket still needs samples.
  if (abs(bp_index - closest_target_bp) <= tolerance && sample_counts[bp_key] < n_samples_per_level) {
    
    # If it's a match, calculate other indices.
    shannon_entropy <- -sum(proportions * log(proportions + 1e-9))
    calculated_n1 <- exp(shannon_entropy)
    richness <- sum(biomass > 0)
    if (richness <= 1) next
    calculated_J <- shannon_entropy / log(richness)
    
    # Normalize biomass across all trip records
    normalized_biomass_scalar <- fixed_total_biomass / sum(biomass)
    gen_dat <- gen_dat %>%
      mutate(across(c(starts_with("sp"), biomass_total), function(x) x*normalized_biomass_scalar)) 
    
    new_dat <- gen_dat %>%
      add_column(PopulationID = sum(sample_counts) + 1, .before=1) #PopulationID is unique across all lists
    
    # Add the new row to the correct bucket in the list.
    results_list[[bp_key]] <- rbind(results_list[[bp_key]], new_dat)
    
    # Increment the counter for that level.
    sample_counts[bp_key] <- sample_counts[bp_key] + 1
    
    # Update the progress bar with the current total number of samples found.
    setTxtProgressBar(pb, sum(sample_counts))
  }
}

# Close the progress bar.
close(pb)

message(paste("\nSimulation complete after", attempts, "attempts."))
return(list(results_list = results_list, sample_counts = sample_counts))
}
