# Packages and Functions ----
library(tidyverse)
library(gdrive)       # devtools::install_github("noaa-afsc/gdrive")
library(mvobservr)    # devtools::install_github("noaa-afsc/mvobservr")
library(mgcv)         # for gam
library(lme4)
library(tweedie)

source("functions/triplet_stats.R")
source("functions/TripletAnalysis.R")
source("functions/getDescriptiveStats.R")
source("functions/runTandFtests.R")
source("functions/runGLMM.R")
source("functions/assign_sig.R")
source("functions/ObserverEffectStats.R")
source("functions/perm_fxn.R")

#' `===================================================================================================================`

#' *Set whether to automatically upload to the Gdrive. Manually change to 'TRUE' when performing full runs!*
set_skip_prompt <- T

#' Set the destination folder for data outputs on the Google shared drive
output_dribble <- gdrive_set_dribble(folder_id = "1Wh-ZQlJ3AIVaQZTWk4QNuyiMfoVECQgt")

## Test Parameters ----

# For tests using permutations or bootstraps, set the number of iterations
nperm <- 1000 
# Set target Berger-Parker index (dominance) levels.
target_bp_levels <- 0.5                                   #' *For false positive test, only testing bp level of 0.5*
# Set bias levels for species (change on observed trips; 0 = no bias, -0.25 = 25% reduction)
bias <- c(0, 0)                                           #' *For false positive test, there will be no bias.*
# Set target coverage rate
trip_coverage <- 0.25

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
# Set how many populations per scenario to generate
n_samples_per_level <- 500                                #' *Increased from 100 to 500*
# Set the number of trips 
ntrips <- 500

#'`====================================================================================================================`

# Generate Trip Populations ============================================================================================

set.seed(123)

trip_sets <- map(rep(1:length(target_bp_levels), n_samples_per_level), ~{
  #set up data frame for trips
  catches <- data.frame(uid = 1:ntrips)
  
  #create catches
  catches$sp_1 <- rtweedie(ntrips, p = tweedie_power, mu = 100*target_bp_levels[.x], phi = phi)
  catches$sp_2 <- rtweedie(ntrips, p = tweedie_power, mu = 100*(1-target_bp_levels[.x]), phi = phi)
  # the 100 acts as a raising factor to get away from a lot of 0s
  
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

# Using the same machinery as in final_comparisons_simulations.R, but we aren't adding any
# bias this time around
map(trip_sets_adj, ~{
  .x %>%
    summarize(sp2 = mean(sp_2), sp2unobs = mean(ifelse(obs==0, sp_2, NA), na.rm=TRUE),
              sp2obs = mean(ifelse(obs==1, sp_2, NA), na.rm=TRUE))
}) %>%
  list_rbind() %>%
  mutate(b = (sp2obs-sp2unobs)/sp2unobs) %>%
  ggplot(aes(x=b)) + geom_histogram()

if(!all(mapply(function(x, y) all(x$sp_2 == y$sp_2), x = trip_sets, y = trip_sets_adj)) ) {
  stop("Somehow sp_2 got modified and bias got added when it shouldn't have been!")
}

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

## Univariate GLMM -----------------------------------------------------------------------------------------------------

res_g <- map(trip_sets_adj, ~ suppressMessages(runGLMM(.x, "biomass_total")), .progress = TRUE)
res_g <- list_rbind(res_g, names_to = "set")

res_g %>% 
  mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mean(AICd>2))

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
  adon <- vegan::adonis2(as.formula(adon_formula), data = .x, permutations = nperm, method="euclidean")
  data.frame(metric = "biomass_total", perma_p = adon$`Pr(>F)`[1])
}, .progress=TRUE) 
res_p <- list_rbind(res_p, names_to = "set")

res_p %>% 
  mutate(bp_level = rep(target_bp_levels, n_samples_per_level)) %>% 
  group_by(bp_level) %>% 
  summarize(mean(perma_p<0.05))

## *Save all but mvglm --------------------------------------------------------------------------------------------------

allbutmv_name <- paste0("output_data/allbutmv_falsepos.Rdata")
save(trip_sets, trip_sets_adj, res_g, res_p, res_permute, res_t, res_tt, file = allbutmv_name)
# Upload to the Google Shared Drive
gdrive_upload(allbutmv_name, output_dribble, skip_prompt = set_skip_prompt)

## MVGLM ---------------------------------------------------------------------------------------------------------------

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

## *Save and upload mvglm results --------------------------------------------------------------------------------------

mvglm_name <- paste0("output_data/mvglmloop_falsepos.Rdata")
save(res_m, file = mvglm_name)
gdrive_upload(mvglm_name, output_dribble, skip_prompt = set_skip_prompt)

#'`====================================================================================================================`

# Save Batch Results ===================================================================================================

load(gdrive_download(mvglm_name, output_dribble))
load(gdrive_download(allbutmv_name, output_dribble))

res_comb <- map(trip_sets_adj, ~{
  .x %>%
    summarize(
      bp_level = max(bp_level), 
      bp_true = max(bp_true), 
      cov = mean(obs), 
      sp1 = sum(sp_1), 
      sp2 = sum(sp_2))
}) %>%
  list_rbind(names_to = "set") %>%
  left_join(res_t, by = "set") %>%
  left_join(res_g, by = "set") %>%
  left_join(res_tt, by = "set") %>%
  left_join(res_p, by = "set") %>%
  left_join(res_m, by = "set") %>%
  left_join(res_permute, by = "set")
res_comb

alltests_name <- paste0("output_data/alltests_falsepos.Rdata")
save(res_comb, file = alltests_name)
gdrive_upload(alltests_name, output_dribble, skip_prompt = set_skip_prompt)


##figure: boxplot for false positives

n_bootstraps <- 1000
set.seed(123)
res_fp <- res_comb %>%
  mutate(t_test_sig = ifelse(tp < 0.05, 1, 0),
         f_test_sig = ifelse(Fp < 0.05, 1, 0),
         dAIC_sig = ifelse(glmconv==1, ifelse(AICd > 2, 1, 0), NA),
         KS_test_sig = ifelse(KSp < 0.05, 1, 0),
         median_test_sig = ifelse(ci_lo > 0 | ci_hi < 0, 1, 0),
         mvglm_sig = ifelse(mvglm_p < 0.05, 1, 0),
         perma_sig = ifelse(perma_p < 0.05, 1, 0)
  ) %>%
  select(ends_with("sig")) %>%
  pivot_longer(cols = ends_with("sig"), names_to = "test", values_to = "sig") %>%
  mutate(test = factor(test, levels = c("t_test_sig", "f_test_sig", "dAIC_sig","KS_test_sig","median_test_sig","mvglm_sig", "perma_sig"),
                       ordered = TRUE,
                       labels = c("t-test", "F test", "uni. GLM", "K-S test", "med. diff.", "MV GLM", "permanova")))
map(1:n_bootstraps, ~{slice_sample(res_fp, n=n_samples_per_level, by=test, replace=TRUE) %>%
    mutate(boot = .x)
}) %>%
  list_rbind() %>%
  group_by(boot, test) %>%
  summarize(pctsig = mean(sig), .groups="drop") %>%
  ggplot(aes(x=test, y=pctsig, fill=test)) +
  geom_boxplot() +
  scale_color_brewer(palette = 'Set2') +
  theme_bw() +
  labs(x=NA,
       y="% of simulations with positive (significant) result",
       title = "false positive rate of tests on total biomass",
       subtitle = paste0(round(trip_coverage*100,0), "% coverage")
  ) +
  scale_y_continuous(labels=scales::label_percent()) +
  geom_hline(yintercept = c(0.1, 0.05), linetype='dashed') +
  theme(legend.position = "none")
