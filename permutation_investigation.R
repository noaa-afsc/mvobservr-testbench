library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(gdrive)       # devtools::install_github("noaa-afsc/gdrive")
library(mvobservr)    # devtools::install_github("noaa-afsc/mvobservr")
mvobservr_dribble <- gdrive_set_dribble(folder_id = "1xQTE9ap6GBnz4ErSrULbEqvPUtQzpHt_")
library(tweedie)


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


#loop over different permutation sizes
pvals_by_perm <- map(rep(c(5, 50, 100, 500, 1000), each=50), ~{
    make_trips() %>%
    mutate(observed = ifelse(obs==1, 'Y', 'N')) %>%
    pivot_longer(cols = starts_with("sp_"),
                 names_to = 'species', values_to = 'biomass') %>%
    mvglm_obs(block = NULL, add_var = NULL, n_permutations = .x, nCores = T) %>%
    pluck("results") %>%
    {data.frame(nperm = .x, pval = .$p)} 
}, .progress=TRUE) %>%
  list_rbind()


#save data to gdrive
save(pvals_by_perm, file="output_data/pvals_by_perm.Rdat")
gdrive_upload(local_path = "output_data/pvals_by_perm.Rdat", gdrive_dribble = mvobservr_dribble)



#load & plot
load(gdrive_download(local_path = "output_data/pvals_by_perm.Rdat", gdrive_dribble = mvobservr_dribble))

pvals_by_perm  %>%
  ggplot(aes(x=nperm, y=pval)) +
  geom_point(alpha = 0.2) +
  geom_smooth() + 
  scale_x_log10()

