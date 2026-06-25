# Packages and Functions ----
library(tidyverse)
library(gdrive)       # pak::pak("noaa-afsc/gdrive")
library(mvobservr)    # pak::pak("noaa-afsc/mvobservr")
library(mgcv)         # for gam
library(glmmTMB)
library(gllvm)
library(lme4)
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
set_skip_prompt <- T

#' Set the destination folder for data outputs on the Google shared drive
output_dribble <- gdrive_set_dribble(folder_id = "1Wh-ZQlJ3AIVaQZTWk4QNuyiMfoVECQgt")

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
n_samples_per_level <- 2        # 
# Set the number of trips 
ntrips <- 500
# Set scalar for Tweedie distributions
mu_scalar <- 100

#'`====================================================================================================================`
for(i in 1:4){
  set_number <- i # Set this
  
  # Generate Trip Populations ============================================================================================
  seed_max <- set_number*3
  seed_seq <- seq(seed_max - 2, seed_max, 1)
  seed_num <- as.numeric(paste(seed_seq, collapse = ""))
  set.seed(seed_num)
  
  
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
  
  if(!all(mapply(function(x, y) all(x$sp_2 == y$sp_2), x = trip_sets, y = trip_sets_adj)) ) {
    stop("Somehow sp_2 got modified and bias got added when it shouldn't have been!")
  }
  
  #'`====================================================================================================================`
  
  # Statistical Tests ====================================================================================================
  
  ## T and F -------------------------------------------------------------------------------------------------------------
  res_t <- map(trip_sets_adj, ~runTandFtests(.x, "biomass_total"))
  res_t <- list_rbind(res_t, names_to = "set")
  
  ## Univariate Permutation ----------------------------------------------------------------------------------------------
  res_p_list <- map(trip_sets_adj, ~perm_fxn(data = .x, 
                                             metric_var = "biomass_total", 
                                             obs_var = "obs", 
                                             n_rep = nperm), 
                    .progress = "Running Permutations")
  res_permute <- list_rbind(res_p_list, names_to = "set")
  
  ## Univariate GLMM -----------------------------------------------------------------------------------------------------
  res_g <- map(trip_sets_adj, ~ suppressMessages(runGLMM(.x, "biomass_total")), .progress = TRUE)
  res_g <- list_rbind(res_g, names_to = "set")
  
  gc(verbose = FALSE)
  
  ## Triplet Analysis ----------------------------------------------------------------------------------------------------
  res_tt <- map(trip_sets_adj, ~TripletAnalysis(.x, "biomass_total", bootstrap_reps = nperm), .progress = TRUE)
  res_tt <- list_rbind(res_tt, names_to = "set")
  
  # MvGLM with glmmTMB ---------------------------------------------------------------------------------------------------
  res_glmm <- trip_sets_adj %>% map(MvGLMglmm, lv = 0, .progress = TRUE)
  res_glmm  <- list_rbind(res_glmm, names_to = "set")
  
  # MvGLM with glmTMB and lv ---------------------------------------------------------------------------------------------
  res_glmm_lv1 <- trip_sets_adj %>% map(MvGLMglmm, lv = 1, .progress = TRUE)
  res_glmm_lv1  <- list_rbind(res_glmm_lv1, names_to = "set") %>% 
    rename(p_glmm1 = p_glmm,
           pwr_glmm1 = pwr_glmm,
           runtime_secs_glmm1 = runtime_secs_glmm)
  
  # MvGLM with GLLVM -----------------------------------------------------------------------------------------------------
  res_gllvm <- map(trip_sets_adj, MvGLMgllvm, .progress = TRUE)
  res_gllvm  <- list_rbind(res_gllvm, names_to = "set")
  
  # MvGLM with GLLVM and 1 latent variables ------------------------------------------------------------------------------
  res_gllvm_lv1 <- map(trip_sets_adj, ~MvGLMgllvm(.x, n_lv = 1), .progress = TRUE)
  res_gllvm_lv1  <- list_rbind(res_gllvm_lv1, names_to = "set") %>% 
    rename(p_gllvm1 = p_gllvm,
           runtime_secs_gllvm1 = runtime_secs_gllvm,
           pwr_gllvm1 = pwr_gllvm)
  
  ## Permanova -----------------------------------------------------------------------------------------------------------
  adon_formula = "Y~factor(obs)"
  res_p <- map(trip_sets_adj, ~{
    Y <- .x %>% dplyr::select(starts_with("sp_"))
    adon <- vegan::adonis2(as.formula(adon_formula), data = .x, permutations = nperm, method="bray")
    data.frame(metric = "biomass_total", perma_p = adon$`Pr(>F)`[1])
  }, .progress=TRUE) 
  res_p <- list_rbind(res_p, names_to = "set")
  
  ## *Save all but mvglm --------------------------------------------------------------------------------------------------
  allbutmv_name <- paste0("output_data/allbutmv_falsepos_set", set_number, ".Rdata")
  save(trip_sets, trip_sets_adj, res_g, res_p, res_permute, res_t, res_tt, res_glmm, res_gllvm, 
       res_glmm_lv1, res_gllvm_lv1, file = allbutmv_name)
  
  # Upload to the Google Shared Drive
  gdrive_upload(allbutmv_name, output_dribble, skip_prompt = set_skip_prompt)
  
  ## Clean memory explicitly before heavy parallel execution
  rm(res_p_list)
  gc(verbose = FALSE)
  Sys.sleep(10) # Extended break for VM memory to settle
  
  ## MvGLM permutation (RESTRUCTURED FOR MEMORY) -------------------------------------------------------------------------
  # Using a standard loop here instead of imap so we can force garbage collection 
  # after every single parallel processing round.
  
  res_m_list <- vector("list", length(trip_sets_adj))
  
  for(j in seq_along(trip_sets_adj)) {
    start_time <- Sys.time()
    
    # 1. Prepare data
    processed_df <- trip_sets_adj[[j]] %>%
      mutate(observed = ifelse(obs == 1, 'Y', 'N')) %>%
      pivot_longer(cols = starts_with("sp_"),
                   names_to = 'species', values_to = 'biomass')
    
    # 2. Run model and capture output
    capture.output(suppressMessages(
      out <- mvglm_obs(processed_df, block = NULL, add_var = NULL, n_permutations = nperm, nCores = TRUE)
    ))
    
    elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # 3. Build clean, wide 1-row data frame
    res_m_list[[j]] <- data.frame(
      set = j,
      p_mvobs    = out$results$p,
      pwr_mvobs   = out$tweedie_profile$xi.max,
      runtime_secs_mvobs = elapsed_time
    )
    
    # 4. Strictly clear iteration-specific heavy objects and force garbage collection
    rm(processed_df, out)
    closeAllConnections() # Closes any hanging parallel socket connections
    gc(verbose = FALSE)
    Sys.sleep(1) # Tiny micro-pause between cores
  }
  
  # Combine all list elements into a single data frame
  res_m <- list_rbind(res_m_list)
  
  ## *Save and upload mvglm results --------------------------------------------------------------------------------------
  mvglm_name <- paste0("output_data/mvglmloop_falsepos_set", set_number, ".Rdata")
  save(res_m, file = mvglm_name)
  gdrive_upload(mvglm_name, output_dribble, skip_prompt = set_skip_prompt)
  
  ## Clean memory before combining
  rm(res_m_list)
  gc(verbose = FALSE)
  Sys.sleep(5)
  
  # Save Batch Results ===================================================================================================
  
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
    left_join(res_glmm, by = "set") %>%
    left_join(res_gllvm, by = "set") %>%
    left_join(res_glmm_lv1, by = "set") %>%
    left_join(res_gllvm_lv1, by = "set") %>%
    left_join(res_permute, by = "set")
  
  alltests_name <- paste0("output_data/alltests_falsepos_set", set_number, ".Rdata")
  save(res_comb, n_samples_per_level, file = alltests_name)
  gdrive_upload(alltests_name, output_dribble, skip_prompt = set_skip_prompt)
  
  # Final cleanup for this set iteration
  rm(trip_sets, trip_sets_adj, res_t, res_g, res_tt, res_p, res_m, res_glmm, 
     res_gllvm, res_glmm_lv1, res_gllvm_lv1, res_permute, res_comb)
  gc(verbose = FALSE)
  Sys.sleep(10) # 10 second breather before the next major set begins
  
} # End of the main for(i in 1:4) loop!

# go to false_positive_test_tables_figures.r
