library(dplyr)
library(furrr)
library(purrr)
library(tidyr)
library(ggplot2)
library(gdrive)       # devtools::install_github("noaa-afsc/gdrive")
library(mvobservr)    # devtools::install_github("noaa-afsc/mvobservr")
library(tweedie)
library(progressr)

mvobservr_dribble <- gdrive_set_dribble(folder_id = "1Wh-ZQlJ3AIVaQZTWk4QNuyiMfoVECQgt")

set.seed(123) 

set_skip_prompt <- T

#setup parameters
bias <- c(0, -0.25) #-0.25, -0.10, -0.40
# Set target coverage rate
trip_coverage <- 0.25

## Test Parameters ----

# Set the Tweedie power parameter (lambda). 1 < p < 2 is typical for biomass.
tweedie_power <- 1.6
# Set the dispersion parameter for the Tweedie distribution.
# Higher values create more variance (further from target BP).
phi <- 3
# Set the desired total biomass for every saved population.
fixed_total_biomass <- 1000000
# Set the number of trips 
ntrips <- 500
# Set target Berger-Parker index (dominance) levels.
target_bp_level <- 0.6
# Set scalar for Tweedie distributions
mu_scalar <- 100

make_trips <- function() {
  #create a single set of trips
  
  # 1. Generate data until we have enough valid rows
  # We generate a bit extra (e.g., 20% more) to reduce the chance of needing a second loop
  valid_catches <- data.frame()
  
  while(nrow(valid_catches) < ntrips) {
    # Determine how many more rows we need
    needed <- ntrips - nrow(valid_catches)
    # Oversampling slightly (1.5x) to ensure we hit the target in one go
    batch_size <- ceiling(needed * 1.5) 
    
    #set up data frame for trips
    temp_catches <- data.frame(id = 1:batch_size)
    
    #create catches
    temp_catches$sp_1 <- rtweedie(batch_size, p = tweedie_power, mu = mu_scalar*target_bp_level, phi = phi)
    temp_catches$sp_2 <- rtweedie(batch_size, p = tweedie_power, mu = mu_scalar*(1-target_bp_level), phi = phi)
    # the 100 acts as a raising factor to get away from a lot of 0s
    
    # Keep only rows where at least one species is > 0
    temp_catches <- temp_catches %>% filter(sp_1 > 0 | sp_2 > 0)
    
    valid_catches <- bind_rows(valid_catches, temp_catches)
  }
  
  # 2. Trim to exactly ntrips and add UID
  catches <- valid_catches %>%
    slice(1:ntrips) %>%
    select(-id) %>% #drop temp id
    mutate(uid = 1:ntrips)
  
  
  #standardize catches
  catches <- catches %>%
    mutate(total_biomass = sum(sp_1 + sp_2)) %>%
    mutate(scalar = fixed_total_biomass/total_biomass) %>%
    mutate(across(c(sp_1, sp_2), function(x) x*scalar)) %>%
    select(-scalar, -total_biomass) %>%
    mutate(biomass_total = sp_1 + sp_2) %>% #this is at the trip level, previously was at the fleet level
    mutate(bp_level = target_bp_level,
           bp_true = sum(sp_1)/sum(biomass_total))
  
  
  #select obs and add bias
  catches <- catches %>%
    mutate(obs = rbinom(nrow(.), 1, trip_coverage)) %>%
    mutate(sp_1 = ifelse(obs == 1, sp_1 * (1+bias[1]), sp_1),
           sp_2 = ifelse(obs == 1, sp_2 * (1+bias[2]), sp_2),
           biomass_total = sp_1 + sp_2
    )
  
  return(catches)
}

# loop over different permutation sizes

perm_levels <- rep(c(15, 30, 50, 100, 250, 500, 1000, 2500, 5000), each = 100)# Define permutations (avoid nperm < cores error)
perm_levels_shuffled <- sample(perm_levels)# SHUFFLE the order randomly
pvals_by_perm <- map(perm_levels_shuffled, ~{# Feed the shuffled vector into map
#Run model 
  df <- make_trips()
  df_formatted <- df %>%
    mutate(observed = ifelse(obs == 1, 'Y', 'N')) %>%
    pivot_longer(cols = starts_with("sp_"),
                 names_to = 'species', values_to = 'biomass')
#run silently 
  capture.output({
    suppressMessages({
      model_out <- mvglm_obs(df_formatted, block = NULL, add_var = NULL, n_permutations = .x, nCores = parallel::detectCores()-2)})
  })
  res <- data.frame(nperm = .x, pval = model_out$results$p)
  gc(verbose = FALSE) 
  return(res)
}, #add custom progress bar
.progress = list(
  name = "Running Models",
  format = "{cli::pb_name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
)) %>% #put it all together
  list_rbind() %>% 
  arrange(nperm) #Sort the final dataframe back in order!


#save data to gdrive
save(pvals_by_perm, file="output_data/SUPPL_pvals_by_perm.Rdat")

gdrive_upload(local_path = "output_data/SUPPL_pvals_by_perm.Rdat",
              gdrive_dribble = mvobservr_dribble,
              skip_prompt = set_skip_prompt)

# # Quick load and plot -----------------------------------------------------
 load(gdrive_download(local_path = "output_data/SUPPL_pvals_by_perm.Rdat", gdrive_dribble = mvobservr_dribble))

suppl_fig_pval <-
  pvals_by_perm  %>%
  ggplot(aes(x=nperm, y=pval, group = nperm)) +  
 # geom_violin() +
  geom_boxplot(fill = "grey", color = "black", outlier.shape = NA) +
  geom_point(alpha = 0.5, pch = 21) +
 # geom_hline(yintercept = 0.01) +
  stat_summary(fun = "median", geom = "point", color = "black", fill = "red", size = 3, pch = 21) +
 #scale_x_log10() + 
  geom_hline(yintercept = 0.001) +
  theme_bw() + 
  scale_y_continuous(limits = c(0, 0.1)) +
  labs(x="Number of permutations", y="p-values")
suppl_fig_pval

save(pvals_by_perm, suppl_fig_pval, file="output_data/SUPPL_pvals_by_perm.Rdat")

gdrive_upload(local_path = "output_data/SUPPL_pvals_by_perm.Rdat", 
              gdrive_dribble = mvobservr_dribble, 
              skip_prompt = set_skip_prompt)
