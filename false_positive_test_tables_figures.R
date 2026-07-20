# Stitch sets  ===================================================================================================
# loads data previously run from false_positive_test.r
library(tidyverse)
library(broom)
library(gdrive)
library(flextable)
library(officer)

mvobservr_dribble <- gdrive_set_dribble(folder_id = "1Wh-ZQlJ3AIVaQZTWk4QNuyiMfoVECQgt")

#list all files in google drive folder that match your naming pattern
# The 'pattern' argument acts as a search filter
drive_files <- gdrive_ls(mvobservr_dribble)
drive_files <- drive_files %>% filter(grepl(pattern = "^alltests_falsepos_set_", x = name))

file_names <- drive_files$name

# Initialize our empty master dataframe
res_comb_tbl <- tibble()

# 4. Loop over the EXACT file names we just pulled from the cloud
for(current_name in file_names) {
  # Reconstruct the path exactly how your gdrive_download function expects it
  current_file <- paste0("output_data/", current_name)
  cat("Downloading and processing:", current_file, "\n")
  tryCatch({
    load(gdrive_download(local_path = current_file, gdrive_dribble = mvobservr_dribble))
    if(nrow(res_comb_tbl) == 0) {
      # For the very first file, assign it and add a tracking column
      res_comb_tbl <- res_comb |> 
        mutate(source_file = current_name)
    } else {
      # For all subsequent files, offset the set ID, add the tracking column, and bind
      res_comb <- res_comb |> 
        mutate(set = set + max(res_comb_tbl$set), source_file = current_name)
      res_comb_tbl <- bind_rows(res_comb_tbl, res_comb)
    }
    # Clean up the single batch
    rm(res_comb)
  }, error = function(e) {
    message("--> Skipped ", current_file, " (Error: ", e$message, ")")
  })
}

#So we have an extra copy here for set1_3k.  Remove it.
Set1_3k <- res_comb_tbl |> filter(source_file == "alltests_falsepos_set_Set1_3k.Rdata")
res_comb_tbl <- res_comb_tbl |> filter(source_file != "alltests_falsepos_set_Set1_3k.Rdata")

#Lets take a quick look at our runtimes - these will identify which sets were run with furrr.
res_comb_tbl |> group_by(source_file) |> summarize(M = mean(runtime_secs_gllvm1, na.rm = TRUE))
#Sets 1 & 2 were run serially, other sets run with furrr and have longer model run times.  
#This is important when we summarize model results.

## figure: boxplot for false positives ----
res_fp <- res_comb_tbl %>%
  mutate(
    # Standard p-value threshold is < 0.05
    t_test_sig = ifelse(tp <= 0.05, 1, 0),
    f_test_sig = ifelse(Fp <= 0.05, 1, 0),
    glm_mgcv_sig = ifelse(glm_mgcv_p <= 0.05, 1, 0), 
    KS_test_sig = ifelse(KSp <= 0.05, 1, 0),
    permu_sig = ifelse(p_val <= 0.05, 1, 0), # permutation
    median_test_sig = ifelse(ci_lo >= 0 | ci_hi <= 0, 1, 0),
    mvglm_sig = ifelse(p_mvobs <= 0.05, 1, 0),
    gllvm_sig = ifelse(p_gllvm <= 0.05, 1, 0),
    # dplyr version (strict and predictable)
    gllvm1_sig = dplyr::if_else(p_gllvm1 <= 0.05, 1, 0, missing = NA_real_), #gllvm1_sig = ifelse(p_gllvm1< 0.05, 1, 0),
    perma_sig = ifelse(perma_p <= 0.05, 1, 0), # permanova
    glmm_sig = ifelse(p_glmm <= 0.05, 1, 0),
    glmm1_sig = ifelse(p_glmm1 <= 0.05, 1, 0)
  ) %>%
  select(ends_with("sig")) %>%
  pivot_longer(
    cols      = ends_with("sig"), 
    names_to  = "test", 
    values_to = "sig"
  ) %>%
  mutate(
    # 1. Swap the raw column names for your clean plotting labels
    test = fct_recode(factor(test),
                      "t-test"      = "t_test_sig",
                      "GLM"         = "glm_mgcv_sig",
                      "F test"      = "f_test_sig",
                      "Permutation" = "permu_sig",
                      "mvobservr"   = "mvglm_sig",
                      "K-S test"    = "KS_test_sig",
                      "PERMANOVA"   = "perma_sig",
                      "Median"      = "median_test_sig",
                      "glmm"        = "glmm_sig",
                      "glmm1"       = "glmm1_sig",
                      "gllvm"       = "gllvm_sig",
                      "gllvm1"      = "gllvm1_sig"
    ),
    # 2. Reorder the factor levels based on the mean of the 'sig' column
    test = fct_reorder(test, sig, .fun = mean, .na_rm = TRUE)
  )

fpr_stats <- res_fp %>%
  group_by(test) %>%
  summarize(
    total_runs = sum(!is.na(sig)),
    false_positives = sum(sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Run an exact binomial test against the expected 0.05 rate
    binom_res = map2(false_positives, total_runs, ~binom.test(x = .x, n = .y, p = 0.05, alt = "two.sided")),
    # Extract the point estimate, p-value, and 95% CI into columns
    tidied = map(binom_res, tidy)
  ) %>%
  unnest(tidied) %>%
  select(test, estimate, conf.low, conf.high, p.value)
# For extracting tidy test results

#Prep bootstrap figure.
set.seed(123) # Ensure reproducibility for the bootstrap
n_bootstraps <- 1000

fpr_plot_dat <- map(1:n_bootstraps, ~{
res_fp %>%
     filter(!is.na(sig)) %>% # Removes crashed models independently per test
     slice_sample(n = n_samples_per_level, by = test, replace = TRUE) %>%
     mutate(boot = .x)
 }) %>%
   list_rbind() %>%
   group_by(boot, test) %>%
   summarize(pctsig = mean(sig, na.rm = TRUE), .groups = "drop")

# Define the tests that should be white
mv_tests <- c("glmm", "gllm1", "mvobservr", "gllvm", "gllvm1")

fp_boot_fig <-
ggplot(fpr_plot_dat, aes(x = test, y = pctsig)) +
  geom_hline(yintercept = 0.05, linetype = 'solid', linewidth = 1, color = "red") +
  geom_jitter(shape = 1, width = 0.1, alpha = 0.05) +
  # Map a logical condition directly to the fill and use the updated quantiles syntax
  geom_violin(aes(fill = test %in% mv_tests), alpha = 0.5, 
              quantiles = c(0.25, 0.5, 0.75), quantile.linetype = 1) +
  geom_point(data = fpr_stats, aes(y = estimate)) +
  # Map TRUE to white, FALSE to dodgerblue
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "dodgerblue")) + 
  theme_classic() +
  labs(x = NULL, y = "Percentage of Positive Tests") +
  scale_y_continuous(labels = scales::label_percent()) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
fp_boot_fig 

# Below is for making an exact binomial figure
## 1. Calculate Exact Binomial Confidence Intervals ----

## 2. Generate the Point-Range Plot ----  this is the results of a hypothesis test of "is the fp rate 0.05%?"
ci_plot <- fpr_stats %>%
  # Order the tests from highest false positive rate to lowest for clean reading
  mutate(test = fct_reorder(test, estimate)) %>%
  ggplot(aes(y = test, x = estimate)) +
  # Add the error bars (95% CI) and the point estimate using the modern geom_errorbar
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, linewidth = 0.8) +
  geom_point(shape = 21, size = 3, color = "black", aes(fill = test %in% mv_tests)) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "dodgerblue")) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  theme_classic() +
  labs(y = NULL, x = "False Positive Rate (95% Confidence Interval)") +
  # Format x-axis as percentages
  scale_x_continuous(labels = scales::percent_format()) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
ci_plot

ggsave("figures/false_positive_plot.png", plot = fpr_plot, 
       width = 5, height = 4, units = "in", dpi = 300) 


## Make two summary tables of runtimes and power from Tweedie
create_summary_table <- function(data, select_prefix, remove_pattern) {
  data %>%
    select(starts_with(select_prefix)) %>%
    rename_with(~str_remove(., remove_pattern)) %>%
    rename("mvobservr" = any_of("mvobs")) %>% 
    pivot_longer(everything(), names_to = "Model", values_to = "value") %>%
    group_by(Model) %>%
    summarize(
      Mean_num = mean(value, na.rm = TRUE),
      SD       = sd(value, na.rm = TRUE),
      N        = sum(!is.na(value)),
      .groups  = "drop"
    ) %>%
    mutate(
      SE       = SD / sqrt(N),
      Mean     = sprintf("%.2f", Mean_num),
      ci_lower = sprintf("%.2f", Mean_num - (1.96 * SE)),
      ci_upper = sprintf("%.2f", Mean_num + (1.96 * SE))
    ) %>%
    select(Model, N, Mean, ci_lower, ci_upper) %>%
    arrange(as.numeric(Mean))
}

# 2. MEE-compliant Flextable function
make_mee_word_table <- function(df, caption_text, value_label) {
  df %>%
    mutate(Estimate = sprintf("%s (%s, %s)", Mean, ci_lower, ci_upper)) %>%
    select(Model, N, !!sym(value_label) := Estimate) %>%
    
    flextable() %>%
    # BES Styling Guidelines
    theme_booktabs() %>% 
    set_caption(caption = caption_text, align_with_table = TRUE) %>% 
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 11, part = "all") %>%
    bold(part = "header") %>%
    
    # Text left-aligned, numbers/CIs centered
    align(align = "left", j = "Model", part = "all") %>%
    align(align = "center", j = c("N", value_label), part = "all") %>%
    autofit()
}

# 3. Generate DataFrames
runtimes_df <- create_summary_table(res_comb_tbl, "runtime", "runtime_secs_")
power_df    <- create_summary_table(res_comb_tbl, "pwr", "pwr_")

# 4. Generate MEE-formatted Word Tables
word_table1 <- runtimes_df %>% make_mee_word_table("Table 1. Computational runtimes for false positives.", "Mean runtime in seconds (95% CI)")
word_table2 <- power_df    %>% make_mee_word_table("Table 2. Tweedie power statistics for false positives.", "Mean Tweedie power (95% CI)")

# 5. Export directly to your MEE manuscript directory
save_as_docx(
  "Table 1" = word_table1, 
  "Table 2" = word_table2, 
  path = "output_data/MEE_manuscript_tables.docx"
)