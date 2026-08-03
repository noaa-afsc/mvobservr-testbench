#' Problem statement: species do not always have constant observer effects - relative biomass can mask effects on 
#' non-dominant species.
#' 
#' Goal: test a new method to do multivariate analysis on species catches and compare these results to those derived 
#' from alternative multivariate and single species or univariate approaches.  We do this over multiple degrees of 
#' species interactions (dominance) and degrees of bias (observer effects). 

# Packages and Functions ----
library(tidyverse)
library(gdrive)       # devtools::install_github("noaa-afsc/gdrive")
library(mvobservr)    # devtools::install_github("noaa-afsc/mvobservr")
library(mgcv) #univariate glm
library(glmmTMB) #MvGLM   
library(lme4)
library(gllvm)
library(tweedie)

source("functions/triplet_stats.R")
source("functions/TripletAnalysis.R")
source("functions/getDescriptiveStats.R")
source("functions/runTandFtests.R")
source("functions/runGLMM.R")
source("functions/MvGLMglmm.R")
source("functions/MvGLMgllvm.R")
source("functions/assign_sig.R")
source("functions/ObserverEffectStats.R")
source("functions/perm_fxn.R")

#' `===================================================================================================================`

#' *Set whether to automatically upload to the Gdrive. Manually change to 'TRUE' when performing full runs!*
set_skip_prompt <- F

# Parameter Setup ----

#' Set the destination folder for data outputs on the Google shared drive
output_dribble <- gdrive_set_dribble(folder_id = "1Wh-ZQlJ3AIVaQZTWk4QNuyiMfoVECQgt")

## Batch Setup ----

#' Currently building 6 different trip populations to run on separate google cloud workstations. `set` determines which
#' random seed will be used and `bias` specifies the magnitude of the observer effect

# Set number (results in a separate seed)
set_number <- 1  # 1 or 2
# Set bias levels for species (change on observed trips; 0 = no bias, -0.25 = 25% reduction)
bias <- c(0, -0.10) #-0.25, -0.10, -0.40
# Set target coverage rate
trip_coverage <- 0.25

## Test Parameters ----

# For tests using permutations or bootstraps, set the number of iterations
nperm <- 1000 
# Set how many populations per scenario to generate
n_samples_per_level <- 100
# Set the Tweedie power parameter (lambda). 1 < p < 2 is typical for biomass.
tweedie_power <- 1.6
# Set the dispersion parameter for the Tweedie distribution.
# Higher values create more variance (further from target BP).
phi <- 3
# Set the desired total biomass for every saved population.
fixed_total_biomass <- 1000000
# Set the number of vessels in the fleet
# Vessel needs to be defined for some of the analyses
nvess <- 1
# Set the number of trips 
ntrips <- 500
# Set target Berger-Parker index (dominance) levels.
target_bp_levels <- seq(0.5, 0.9, by = 0.1)
# Set scalar for Tweedie distributions
mu_scalar <- 100

#'`====================================================================================================================`

# Generate Trip Populations ============================================================================================

if(set_number == 1) {
  set.seed(123)
} else if(set_number == 2) {
  set.seed(456)
} else stop("'set_number' needs to be specified as '1' or '2'!")

trip_sets <- map(rep(1:length(target_bp_levels), n_samples_per_level), ~{
  
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
    temp_catches$sp_1 <- rtweedie(batch_size, p = tweedie_power, mu = mu_scalar*target_bp_levels[.x], phi = phi)
    temp_catches$sp_2 <- rtweedie(batch_size, p = tweedie_power, mu = mu_scalar*(1-target_bp_levels[.x]), phi = phi)
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
    mutate(bp_level = target_bp_levels[.x],
           bp_true = sum(sp_1)/sum(biomass_total))
  
  #add vessel randomly
  catches <- catches %>%
    mutate(vessnum = sample(nvess, size=nrow(.), replace=TRUE)) 
  
}, .progress = TRUE)

## Add Bias ------------------------------------------------------------------------------------------------------------

trip_sets_adj <- map(trip_sets, ~{.x %>%
    mutate(obs = rbinom(nrow(.), 1, trip_coverage)) %>%
    mutate(sp_1 = ifelse(obs == 1, sp_1 * (1+bias[1]), sp_1),
           sp_2 = ifelse(obs == 1, sp_2 * (1+bias[2]), sp_2),
           biomass_total = sp_1 + sp_2
    )
}, .progress = TRUE)


## Check number of trips with 0 biomass
map(trip_sets_adj, ~{
  .x %>%
    summarize(propdouble0s = mean(sp_1 == 0 & sp_2 == 0))
}) %>%
  list_rbind() %>%
  count(propdouble0s)


## Check distribution of overall bias
 ### bias applied at trip level

map(trip_sets_adj, ~{
  .x %>%
    summarize(sp2 = mean(sp_2), sp2unobs = mean(ifelse(obs==0, sp_2, NA), na.rm=TRUE),
              sp2obs = mean(ifelse(obs==1, sp_2, NA), na.rm=TRUE))
}) %>%
  list_rbind() %>%
  mutate(b = (sp2obs-sp2unobs)/sp2unobs) %>%
  ggplot(aes(x=b)) + geom_histogram()

#'`====================================================================================================================`

# Statistical Tests ====================================================================================================

## T and F -------------------------------------------------------------------------------------------------------------

res_t <- map(trip_sets_adj, ~runTandFtests(.x, "biomass_total"))
res_t <- list_rbind(res_t, names_to = "set")

res_t %>% 
  mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mean(tp < 0.05), mean(Fp < 0.05))


## Univariate Permutation ----------------------------------------------------------------------------------------------

res_p_list <- map(trip_sets_adj, ~perm_fxn(data = .x, 
                                           metric_var = "biomass_total", 
                                           obs_var = "obs", 
                                           n_rep = nperm), 
                  .progress = "Running Permutations")
# Combine results into a 500-row table
res_permute <- list_rbind(res_p_list, names_to = "set")

res_permute %>% 
  mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mn = mean(p_val<0.05, na.rm = TRUE))

## Univariate GLM -----------------------------------------------------------------------------------------------------

res_g <- map(trip_sets_adj, ~ suppressMessages(runGLMM(.x, "biomass_total")), .progress = TRUE)
res_g <- list_rbind(res_g, names_to = "set")

res_g %>% 
  mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mn = mean(glm_mgcv_p<0.05, na.rm = TRUE), 
            n_success_of_100 = sum(!is.na(glm_mgcv_p)))

gc(verbose = FALSE)

## Triplet Analysis ----------------------------------------------------------------------------------------------------

res_tt <- map(trip_sets_adj, ~TripletAnalysis(.x, "biomass_total", bootstrap_reps = nperm), .progress = TRUE)
res_tt <- list_rbind(res_tt, names_to = "set")

res_tt %>% 
  mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mean(KSp < 0.05), mean(ci_lo > 0 | ci_hi < 0))

## Permanova -----------------------------------------------------------------------------------------------------------

adon_formula = "Y~factor(obs)"
res_p <- map(trip_sets_adj, ~{
  Y <- .x %>% dplyr::select(starts_with("sp_"))
  adon <- vegan::adonis2(as.formula(adon_formula), data = .x, permutations = nperm, method="bray")
  data.frame(metric = "biomass_total", perma_p = adon$`Pr(>F)`[1])
}, .progress=TRUE) #45min for 500
res_p <- list_rbind(res_p, names_to = "set")

res_p %>% 
  mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mean(perma_p<0.05))

# MvGLM with glmmTMB ---------------------------------------------------------------------------------------------------

res_glmm <- trip_sets_adj %>% map(MvGLMglmm, lv = 0, .progress = TRUE) %>% 
  list_rbind(names_to = "set") %>% mutate(metric = "biomass_total")

res_glmm %>% mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mean(p_glmm < 0.05))

# MvGLM with GLLVM ----------------------------------------------------------------------------------------------------------------

res_gllvm <- map(trip_sets_adj, MvGLMgllvm, .progress = TRUE) %>% 
  list_rbind(names_to = "set") %>% mutate(metric = "biomass_total")

res_gllvm %>% mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mean(p_gllvm < 0.05))

## *Save all but MvGLMperm --------------------------------------------------------------------------------------------------

allbutmv_name <- paste0("output_data/allbutmv_b", max(abs(bias))*100, "_set", set_number, ".Rdata")
save(trip_sets, trip_sets_adj, res_g, res_p, res_permute, res_t, 
     res_tt, res_glmm, res_gllvm, file = allbutmv_name)
# Upload to the Google Shared Drive
gdrive_upload(allbutmv_name, output_dribble, skip_prompt = set_skip_prompt)

## MVGLMperm ---------------------------------------------------------------------------------------------------------------

res_m <- imap(1:length(trip_sets_adj), ~{
  print(paste0("***set ", .x, " ", now()))
  trip_sets_adj[[.x]] %>%
    mutate(observed = ifelse(obs==1, 'Y', 'N')) %>%
    pivot_longer(cols = starts_with("sp_"),
                 names_to = 'species', values_to = 'biomass') %>%
    mvglm_obs(block = NULL, add_var = NULL, n_permutations = nperm, nCores = TRUE) %>%
    pluck("results") %>%
    {data.frame(t(c(uni.p = .$uni.p, biomass_total = .$p)))} %>%
    pivot_longer(everything(), names_to = "metric", values_to = "mvglm_p") 
})
res_m <- list_rbind(res_m, names_to = "set") %>%
  filter(metric == "biomass_total")

res_m %>% 
  mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mean(mvglm_p<0.05))

## *Save and upload mvobservr results --------------------------------------------------------------------------------------

mvglm_name <- paste0("output_data/mvglmloop_b", max(abs(bias))*100, "_set", set_number, ".Rdata")
save(res_m, file = mvglm_name)
gdrive_upload(mvglm_name, output_dribble, skip_prompt = set_skip_prompt)

#'`====================================================================================================================`

# Save Batch Results ===================================================================================================
if(set_skip_prompt != F){
  load(gdrive_download(mvglm_name, output_dribble))
  load(gdrive_download(allbutmv_name, output_dribble))
}

res_comb <- map(trip_sets_adj, ~{
  .x %>%
    summarize(bp_level = max(bp_level), bp_true = max(bp_true), cov = mean(obs), sp1 = sum(sp_1), sp2 = sum(sp_2))
}) %>%
  list_rbind(names_to = "set") %>%
  left_join(res_t, by = "set") %>%
  left_join(res_g, by = "set") %>%
  left_join(res_tt, by = "set") %>%
  left_join(res_p, by = "set") %>%
  left_join(res_m, by = "set") %>%
  left_join(res_glmm , by = "set") %>%
  left_join(res_gllvm , by = "set") %>%
  left_join(res_permute, by = "set")
res_comb

alltests_name <- paste0(
  "output_data/alltests_2spp_oe", max(abs(bias))*100, "_cov", trip_coverage*100, "_set", set_number, ".Rdata"
)
save(res_comb, file = alltests_name)
gdrive_upload(alltests_name, output_dribble, skip_prompt = set_skip_prompt)
