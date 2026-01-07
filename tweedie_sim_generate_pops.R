# Craig.Faunce@noaa.gov
# January 2026

# The goal of this script is to generate simulated populations of different species that have a pre-specified
# Berger parker index (dominance) and number of species for use in testing mvobservr.
# This version was passed to Debra Duarte and was eventually incorporated into the main script.

# Load necessary libraries for plotting and the Tweedie distribution.
# If you don't have them installed, run:
# install.packages("ggplot2")
# install.packages("tweedie")
library(ggplot2)
library(tweedie)

# --- Simulation Setup ---

# Set the number of species to a more realistic value.
n_species <- 15
# Set how many valid communities we want to find for each target level.
n_samples_per_level <- 500
# Set the tolerance for matching the target Berger-Parker index.
# This needs to be small since our steps are small.
tolerance <- 0.005
# Set the Tweedie power parameter (lambda). 1 < p < 2 is typical for biomass.
tweedie_power <- 1.6
# Set the dispersion parameter for the Tweedie distribution.
# Higher values create more variance (more rare/dominant species).
phi <- 6
# Set the desired total biomass for every saved population.
fixed_total_biomass <- 10000
# for reproducibility
set.seed(123)

# --- Define Target Berger-Parker Levels ---

# Set target Berger-Parker index (dominance) levels.
target_bp_levels <- seq(0.5, 0.95, by = 0.05)

# --- Data Storage ---

# Create a list to act as "buckets" for each target level.
results_list <- setNames(lapply(target_bp_levels, function(x) data.frame()), target_bp_levels)
# Create a counter to track how many samples we've found for each level.
sample_counts <- setNames(rep(0, length(target_bp_levels)), target_bp_levels)

# --- Simulation Loop ---

message("Generating populations and sorting into target levels...")
attempts <- 0
max_attempts <- 3000000 # Set a generous ceiling for the total search.

# Calculate the total number of samples needed for the progress bar.
total_needed <- n_samples_per_level * length(target_bp_levels)
# Initialize the progress bar.
pb <- txtProgressBar(min = 0, max = total_needed, style = 3)

# Loop until all levels have the desired number of samples.
while(any(sample_counts < n_samples_per_level) && attempts < max_attempts) {
  attempts <- attempts + 1
  
  # Step 1: Generate a single realistic community.
  mu_values <- rlnorm(n_species, meanlog = 2, sdlog = 1.5)
  biomass <- rtweedie(n_species, p = tweedie_power, mu = mu_values, phi = phi)
  
  if (sum(biomass) == 0) next
  
  # Step 2: Calculate its Berger-Parker index.
  proportions <- biomass / sum(biomass)
  bp_index <- max(proportions)
  
  # Step 3: Find the closest target level.
  closest_target_index <- which.min(abs(target_bp_levels - bp_index))
  closest_target_bp <- target_bp_levels[closest_target_index]
  
  # Step 4: Check if it's a match and if that bucket still needs samples.
  if (abs(bp_index - closest_target_bp) <= tolerance && sample_counts[as.character(closest_target_bp)] < n_samples_per_level) {
    
    # If it's a match, calculate other indices.
    shannon_entropy <- -sum(proportions * log(proportions + 1e-9))
    calculated_n1 <- exp(shannon_entropy)
    richness <- sum(biomass > 0)
    if (richness <= 1) next
    calculated_J <- shannon_entropy / log(richness)
    
    # Normalize biomass and create the data row.
    normalized_biomass <- proportions * fixed_total_biomass
    biomass_list <- as.list(setNames(normalized_biomass, paste0("sp", 1:n_species)))
    
    new_row <- data.frame(
      PopulationID = sample_counts[as.character(closest_target_bp)] + 1,
      BergerParker = bp_index,
      N1 = calculated_n1,
      Pielou_J = calculated_J,
      TargetBP = as.factor(closest_target_bp),
      biomass_list
    )
    
    # Add the new row to the correct bucket in the list.
    results_list[[as.character(closest_target_bp)]] <- rbind(results_list[[as.character(closest_target_bp)]], new_row)
    
    # Increment the counter for that level.
    sample_counts[as.character(closest_target_bp)] <- sample_counts[as.character(closest_target_bp)] + 1
    
    # Update the progress bar with the current total number of samples found.
    setTxtProgressBar(pb, sum(sample_counts))
  }
}

# Close the progress bar.
close(pb)

message(paste("\nSimulation complete after", attempts, "attempts."))

# --- Data Reordering ---
message("Reordering species abundances in the final results...")

# Loop through each data frame in the results list (one for each target BP level)
for (i in 1:length(results_list)) {
  df <- results_list[[i]]
  
  # Skip if the data frame is empty (no samples were found for this level)
  if (nrow(df) == 0) next
  
  # Identify the columns that contain the biomass data
  biomass_cols_indices <- grep("sp", names(df))
  biomass_data <- df[, biomass_cols_indices]
  
  # Sort each row (population) in descending order
  # t() is used because apply returns a matrix with columns as populations
  sorted_biomass <- t(apply(biomass_data, 1, sort, decreasing = TRUE))
  
  # Convert back to a data frame and assign original names
  sorted_biomass_df <- as.data.frame(sorted_biomass)
  names(sorted_biomass_df) <- names(biomass_data)
  
  # Replace the original biomass columns with the sorted ones
  df[, biomass_cols_indices] <- sorted_biomass_df
  
  # Put the modified data frame back into the list
  results_list[[i]] <- df
}

# --- Data Preparation and Plotting ---

# Combine the list of data frames into one single data frame for plotting.
final_results_df <- do.call(rbind, results_list)

# You can view the final data frame with the full biomass data by running:
# View(final_results_df)

# Generate the scatter plot using ggplot2.
ggplot(final_results_df, aes(x = BergerParker, y = N1, color = Pielou_J)) +
  # Add points with some transparency.
  geom_point(alpha = 0.8, size = 2.5) +
  # Set titles and axis labels.
  labs(
    title = "N1 vs. Dominance in Tweedie-Distributed Communities",
    subtitle = paste("Based on simulations of a", n_species, "-species community"),
    x = "Berger-Parker Index (Dominance)",
    y = "Hill's N1 (Effective Number of Species)",
    color = "Evenness (J')"
  ) +
  # Use a clean theme and a continuous color scale suitable for gradients.
  theme_minimal(base_size = 14) +
  scale_color_viridis_c(direction = -1)


