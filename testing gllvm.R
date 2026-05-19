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
library(mgcv)         # for gam
library(lme4)
library(tweedie)
library(gllvm)

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

#'`====================================================================================================================`

# Generate Trip Populations ============================================================================================

if(set_number == 1) {
  set.seed(123)
} else if(set_number == 2) {
  set.seed(456)
} else stop("'set_number' needs to be specified as '1' or '2'!")

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

map(trip_sets_adj, ~{
  .x %>%
    summarize(sp2 = mean(sp_2), sp2unobs = mean(ifelse(obs==0, sp_2, NA), na.rm=TRUE),
              sp2obs = mean(ifelse(obs==1, sp_2, NA), na.rm=TRUE))
}) %>%
  list_rbind() %>%
  mutate(b = (sp2obs-sp2unobs)/sp2unobs) %>%
  ggplot(aes(x=b)) + geom_histogram()

#'`====================================================================================================================`

trips <- as.data.frame(trip_sets_adj[2])

Y <- as.matrix(trips %>% select(starts_with("sp_"))) #Species
X_TR <- as.matrix(trips %>% select("obs"))
#X_TR$obs <- as.factor(X_TR$obs)

#Setup test.

# 1. PREPARE DATA STRUCTURES
# Y: Your species matrix (sites as rows, species as columns)
# X_env: Dataframe containing your fixed observer factor


# 2. FIT THE NULL MODEL
# Controls for block and builds the ordination space, ignoring the observer effect
system.time({ # 194 seconds for 1 set.
  model_null <- gllvm(par = c(1.3, 1.8),
  y = Y, 
 # X = NULL, 
family = "tweedie" # Use "tweedie" if your data is biomass/weight!
)

# 3. FIT THE FULL MODEL
# Introduces the ObserverStatus as an environmental predictor
model_full <- gllvm(
  y = Y, 
  X = X_TR, 
  family = "tweedie"
)

# 4. TEST FOR SIGNIFICANCE (The P-Value Step)
# Computes a Likelihood Ratio Test between the nested models
anova(model_null, model_full)
})

#From instructions- these diagnostic plots look like shit
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
plot(model_null, var.colors = 1)

par(mfrow = c(1, 2))
plot(model_null, which = 1:2, var.colors = 1)


coefplot(model_full) #shows the 10% effect.

## New data - example of testing with sites, locations, env parameters, etc.
# http://jenniniku.github.io/gllvm/articles/vignette2.html


