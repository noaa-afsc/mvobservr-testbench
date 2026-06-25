#setup ----
library(tidyverse)
library(tweedie)
library(lme4)
library(glmmTMB)
source("functions/MvGLMglmm.R")
library(gllvm)
source("functions/MvGLMgllvm.R")
library(vegan)
library(gdrive)
library(mvobservr) # devtools::install_github("noaa-afsc/mvobservr")

set <- 1  # 1 or 2
output_name <- paste0("output_data/param_plot_set", set, ".Rdata")

# Constants ----
# Set how many populations per scenario to generate
# There are 30 combinations, so 1 n_samples_per_level is 30 populations to test.
n_samples_per_level <- 1 #3000 populations at 100

# Set bias levels for species (change on observed trips; 0 = no bias, -0.25 = 25% reduction)
bias <- c(0, -0.25)
# Set target Berger-Parker index (dominance) levels.
bp_level <- 0.6
# Set number of permutations
n_permutations <- 1000

# Defaults for varying params ----
# Set the Tweedie power parameter (lambda). 1 < p < 2 is typical for biomass.
tweedie_power <- 1.6
# Set the dispersion parameter for the Tweedie distribution.
# Higher values create more variance (further from target BP).
phi <- 3
# Set the desired total biomass for every saved population.
fixed_total_biomass <- 6 #log scale
# Set the number of trips 
ntrips <- 500
# Set target coverage rate
trip_coverage <- 0.25
# Set the mu scalar
mu_scalar <- 100

# Set up varying parameters and their levels ----
params_to_test <- data.frame(
  param = c("tweedie_power", "phi", "fixed_total_biomass", "ntrips", "trip_coverage", "mu_scalar"),
  param_display_name = c("Tweedie power parameter (poisson-gamma)",
                         "Tweedie phi parameter (dispersion)",
                         "Total Biomass (log scale)",
                         "Total trips",
                         paste0("Number of observed trips (out of ",ntrips," total)"),
                         "Tweedie mu parameter (mean)"
  ),
  param_values = c("1.1, 1.3, 1.5, 1.7, 1.9", #power
                   "1, 3, 5, 7, 9", #phi
                   "3, 4, 5, 6, 7", #log total biomass
                   "50, 100, 250, 500, 1000", #total trips
                   "0.1, 0.3, 0.5, 0.7, 0.9", #coverage rate
                   "1, 10, 100, 1000, 2000" #mu scalar
  ),
  defaults = c(tweedie_power, phi, fixed_total_biomass, ntrips, trip_coverage, mu_scalar)
)
# Need to run to here for figure to work with quick pull down.

param_values_df <- params_to_test %>% 
  select(param, param_values) %>%
  separate_wider_delim(
    cols = param_values, 
    delim = ",", 
    names = c("val1", "val2", "val3", "val4", "val5")
  ) %>%
  pivot_longer(cols=starts_with("val"), values_to = "param_val") %>%
  mutate(param_val = as.numeric(param_val)) %>%
  select(-name)

# Generate data ----
# for reproducibility
if(set == 1) {
  set.seed(123)
} else if (set == 2) {
  set.seed(456)
} else stop("'set' needs to be either 1 or 2!")

trip_sets_adj <- map(rep(1:nrow(param_values_df), n_samples_per_level), ~{
  #replace the variable parameter with the loop value, all other values remain t default
  assign(param_values_df$param[.x], param_values_df$param_val[.x])
  #print(paste(bp_level, ntrips, tweedie_power, mu_scalar, phi, fixed_total_biomass, trip_coverage))
  
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
    temp_catches$sp_1 <- rtweedie(batch_size, p = tweedie_power,mu = mu_scalar*bp_level, phi = phi)
    temp_catches$sp_2 <- rtweedie(batch_size, p = tweedie_power, mu = mu_scalar*(1-bp_level), phi = phi)
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
    mutate(scalar = (10^fixed_total_biomass)/total_biomass) %>% #fixed_total_biomass on log scale
    mutate(across(c(sp_1, sp_2), function(x) x*scalar)) %>%
    select(-scalar, -total_biomass) %>%
    mutate(biomass_total = sp_1 + sp_2) %>%
    mutate(bp_level = bp_level,
           bp_true = sum(sp_1)/sum(biomass_total),
           ntrips = ntrips,
           tweedie_power = tweedie_power,
           mu_scalar = mu_scalar,
           phi = phi,
           fixed_total_biomass = fixed_total_biomass,
           trip_coverage = trip_coverage) %>%
    mutate(param_mod = param_values_df$param[.x],
           param_value = param_values_df$param_val[.x],
           .before=1)
  
  #select observed and add bias
  catches <- catches %>%
    mutate(obs = rbinom(nrow(.), 1, trip_coverage)) %>%
    mutate(sp_1 = ifelse(obs == 1, sp_1 * (1+bias[1]), sp_1),
           sp_2 = ifelse(obs == 1, sp_2 * (1+bias[2]), sp_2),
           biomass_total = sp_1 + sp_2
    )
  
  return(catches)
}, .progress=TRUE)

# Begin tests ----
##permanova
adon_formula = "Y~factor(obs)"

res_p <- map(trip_sets_adj, ~{
  Y <- .x %>% dplyr::select(starts_with("sp_"))
  adon <- vegan::adonis2(as.formula(adon_formula), data = .x, permutations = n_permutations, method="bray")
  data.frame(metric = "biomass_total", perma_p = adon$`Pr(>F)`[1])
}, .progress=TRUE) #1 min for 60 sets
res_p <- list_rbind(res_p, names_to = "set")

##mvglm
#with glmmTMB
res_glmm <- trip_sets_adj %>% map(MvGLMglmm, lv = 0, .progress = TRUE) %>% 
  list_rbind(names_to = "set") %>% mutate(metric = "biomass_total")

#with gllvm
res_gllvm <- map(trip_sets_adj, MvGLMgllvm, .progress = TRUE) %>% 
  list_rbind(names_to = "set") %>% mutate(metric = "biomass_total")

#with mvobservr
res_m <- imap(1:length(trip_sets_adj), ~{
  print(paste0("***set ", .x, " ", now()))
  trip_sets_adj[[.x]] %>%
    mutate(observed = ifelse(obs==1, 'Y', 'N')) %>%
    pivot_longer(cols=starts_with("sp_"),
                 names_to = 'species', values_to = 'biomass') %>%
    mvglm_obs(block = NULL, add_var = NULL, n_permutations= n_permutations, nCores = TRUE) %>%
    pluck("results") %>%
    {data.frame(t(c(uni.p = .$uni.p, biomass_total = .$p)))} %>%
    pivot_longer(everything(), names_to = "metric", values_to = "mvglm_p") 
})
res_m <- list_rbind(res_m, names_to = "set") %>%
  filter(metric == "biomass_total")

#combine output
res_comb <- map(trip_sets_adj, ~ .x[1, c("param_mod", "param_value")]) %>%
  list_rbind(names_to = "set") %>%
  left_join(res_p, by="set") %>%
  left_join(res_glmm, by=c("set","metric")) %>%
  left_join(res_gllvm, by=c("set","metric")) %>%
  left_join(res_m, by=c("set","metric"))
res_comb

save(res_comb, file = output_name)

#upload files to gdrive only once when transferring from data ran on cloud workstation
mvobservr_dribble <- gdrive_set_dribble(folder_id = "1xQTE9ap6GBnz4ErSrULbEqvPUtQzpHt_")
gdrive_upload(local_path = output_name, gdrive_dribble = mvobservr_dribble)

# Quick pulldown and figure generation ----
#pull down and load data files from gdrive
mvobservr_dribble <- gdrive_set_dribble(folder_id = "1xQTE9ap6GBnz4ErSrULbEqvPUtQzpHt_")
load(gdrive_download(local_path = "output_data/param_plot_set1.Rdata", gdrive_dribble = mvobservr_dribble))
res_comb_tbl <- res_comb
load(gdrive_download(local_path = "output_data/param_plot_set2.Rdata", gdrive_dribble = mvobservr_dribble))
res_comb <- res_comb |> mutate(set = set + max(res_comb_tbl$set))
res_comb_tbl <- rbind(res_comb_tbl, res_comb)
rm(res_comb)

#make plot
perf_plot_dat <- 
  res_comb_tbl %>%
  pivot_longer(cols=c(perma_p, mvglm_p, p_glmm, p_gllvm), names_to = "Test", values_to = "p") %>%
  group_by(param_mod, param_value, Test) %>%
  summarize(
    psig = mean(p < 0.05, na.rm = TRUE), 
    num  = sum(!is.na(p)), # Counts only non-NA p-values for accurate SE calculation
    se   = sqrt((psig * (1 - psig)) / num), 
    .groups = "drop"
  ) %>%
  inner_join(params_to_test, by=c("param_mod" = "param")) %>%
  mutate(Test = factor(Test, levels = c("perma_p", "mvglm_p", "p_glmm", "p_gllvm"),
                       labels = c("PERMANOVA", "mvobservr", "glmm", "gllvm")),
         #fix panel titles
         param_display_name = str_to_title(param_display_name) %>%
           str_replace_all(" Of ", " of ")
  ) %>%
  mutate(param_value = case_when(
    param_mod == "trip_coverage" ~ param_value*ntrips, 
    #convert to number of trips
    TRUE ~ param_value)) %>%
  mutate(defaults = case_when(
    param_mod == "trip_coverage" ~ defaults*ntrips, 
    #convert to number of trips
    TRUE ~ defaults))

pd <- position_dodge(width = 0.1) # Now 0.4 means "40% of the gap between points"

perf_plot <- ggplot(perf_plot_dat, 
                    aes(x = factor(param_value), # Treat X as labels, not numbers
                        y = psig, 
                        shape = Test,
                        linetype = Test,
                        fill = Test,
                        group = Test,
                        color = Test)) +         # Group keeps the lines connected
  geom_line(position = pd) +
  geom_errorbar(aes(ymin = psig - (1.96 * se), ymax = psig + (1.96 * se)), 
                width = 0.2,                     # Width is now relative to the discrete gap
                position = pd) + 
  geom_point(size = 3, color = "black", position = pd) + 
  scale_color_manual(values = c("#377EB8",  # PERMANOVA (Blue)
                                "#E41A1C", # MvGLM (Red)
                                "cyan3", #glmmTMB MvGLM
                                "violet" #gllvm)
  )) +
  scale_shape_manual(values = c(21, 21, 21, 21)) +
  scale_fill_manual(values = c("#377EB8",  # PERMANOVA (Blue)
                               "#E41A1C", # MvGLM (Red)
                               "cyan3", #glmmTMB MvGLM
                               "violet" #gllvm)
  )) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
  facet_wrap(~param_display_name, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Parameter Value", y = "Power")

perf_plot

ggsave("figures/power_analysis_plot.png", perf_plot, 
       width = 10, height = 7, dpi = 300, units = "in")