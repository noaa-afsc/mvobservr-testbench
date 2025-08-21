# Load necessary libraries for plotting and the Tweedie distribution.
# If you don't have them installed, run:
# install.packages("ggplot2")
# install.packages("tweedie")
library(ggplot2)
library(tweedie)

# --- Simulation Setup ---

# Set the number of species to a more realistic value.
n_species <- 25
# Set how many valid communities we want to find for each target level.
n_samples_per_level <- 50
# Set the tolerance for matching the target Berger-Parker index.
# This needs to be small since our steps are small.
tolerance <- 0.005
# Set the Tweedie power parameter (lambda). 1 < p < 2 is typical for biomass.
tweedie_power <- 1.6
# Set the dispersion parameter for the Tweedie distribution.
# Higher values create more variance (more rare/dominant species).
phi <- 5
# for reproducibility
set.seed(123)

# --- Define Target Berger-Parker Levels ---

# Set target Berger-Parker index (dominance) levels.
target_bp_levels <- seq(0.5, 0.95, by = 0.05)

# --- Data Storage ---

# We'll store the results for each level in a list and combine them at the end.
results_list <- list()

# --- Simulation Loop ---

# Loop over each of our target Berger-Parker values using an index 'j'.
for (j in 1:length(target_bp_levels)) {
  
  # Get the target value for this iteration.
  target_bp <- target_bp_levels[j]
  
  # Create an empty data frame to store results for the current target level.
  level_results <- data.frame()
  
  # Print a message to the console to show progress.
  message(paste("Searching for communities with Berger-Parker ≈", target_bp))
  
  # Keep track of attempts to prevent an infinite loop.
  attempts <- 0
  max_attempts <- 1000000 # Increased attempts as criteria are stricter.
  
  # Loop until we have found the desired number of matching communities.
  while (nrow(level_results) < n_samples_per_level && attempts < max_attempts) {
    attempts <- attempts + 1
    
    # Step 1: Generate a realistic community using the Tweedie distribution.
    mu_values <- rlnorm(n_species, meanlog = 2, sdlog = 1.5)
    biomass <- rtweedie(n_species, p = tweedie_power, mu = mu_values, phi = phi)
    
    # Skip if total biomass is zero.
    if (sum(biomass) == 0) next
    
    # Step 2: Calculate diversity indices.
    proportions <- biomass / sum(biomass)
    bp_index <- max(proportions)
    
    # Step 3: Check if the community's dominance matches our target.
    if (abs(bp_index - target_bp) <= tolerance) {
      
      # Step 4: If it matches, calculate other indices and store everything.
      # We add a very small number to proportions to avoid log(0) issues.
      shannon_entropy <- -sum(proportions * log(proportions + 1e-9))
      calculated_n1 <- exp(shannon_entropy)
      
      # Calculate Pielou's J'.
      richness <- sum(biomass > 0)
      if (richness <= 1) next # Evenness is undefined for 1 or 0 species.
      calculated_J <- shannon_entropy / log(richness)
      
      # Store the results.
      new_row <- data.frame(
        BergerParker = bp_index,
        N1 = calculated_n1,
        Pielou_J = calculated_J,
        TargetBP = as.factor(target_bp)
      )
      
      level_results <- rbind(level_results, new_row)
    }
  }
  
  # Add the results for this level to our main list using the correct index 'j'.
  results_list[[j]] <- level_results
  message(paste("   ...Found", nrow(level_results), "samples in", attempts, "attempts."))
  
  # Warn the user if the loop timed out.
  if (attempts >= max_attempts) {
    warning(paste("Could not find enough samples for target Berger-Parker =", target_bp))
  }
}

# --- Data Preparation and Plotting ---

# Combine the list of data frames into one single data frame for plotting.
final_results_df <- do.call(rbind, results_list)

# You can view the final data frame with the abundances by running:
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


