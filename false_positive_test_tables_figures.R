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

#TODO - revert (decision - to the former bootstrap version of this figure - see false_positive_test_ccf5693.r)

## figure: boxplot for false positives ----
res_fp <- res_comb_tbl %>%
  mutate(
    # Standard p-value threshold is < 0.05
    t_test_sig = ifelse(tp < 0.05, 1, 0),
    f_test_sig = ifelse(Fp < 0.05, 1, 0),
    glm_mgcv_sig = ifelse(glm_mgcv_p <= 0.05, 1, 0), 
    KS_test_sig = ifelse(KSp < 0.05, 1, 0),
    permu_sig = ifelse(p_val < 0.05, 1, 0), # permutation
    median_test_sig = ifelse(ci_lo > 0 | ci_hi < 0, 1, 0),
    mvglm_sig = ifelse(p_mvobs < 0.05, 1, 0),
    gllvm_sig = ifelse(p_gllvm < 0.05, 1, 0),
    # dplyr version (strict and predictable)
    gllvm1_sig = dplyr::if_else(p_gllvm1 < 0.05, 1, 0, missing = NA_real_), #gllvm1_sig = ifelse(p_gllvm1< 0.05, 1, 0),
    perma_sig = ifelse(perma_p < 0.05, 1, 0), # permanova
    glmm_sig = ifelse(p_glmm < 0.05, 1, 0),
    glmm1_sig = ifelse(p_glmm1 < 0.05, 1, 0)
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

#Prep bootstrap figure.
set.seed(123) # Ensure reproducibility for the bootstrap
n_bootstraps <- 1000

fpr_plot <- map(1:n_bootstraps, ~{
  res_fp %>%
    filter(!is.na(sig)) %>% # Removes crashed models independently per test
    slice_sample(n = n_samples_per_level, by = test, replace = TRUE) %>%
    mutate(boot = .x)
}) %>%
  list_rbind() %>%
  group_by(boot, test) %>%
  summarize(pctsig = mean(sig, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = test, y = pctsig)) + #fill was = test
  geom_violin(alpha = 0.5, fill = 'black', quantiles = c(0.25, 0.5, 0.75), 
              quantile.linetype = 1, quantile.color = "white") +
  stat_summary(fun = "mean", geom = "point", size = 2, color = "white") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Percentage of Positive (Significant) Tests",
    # title = "False Positive Rate of Tests on Total Biomass",
    #subtitle = paste0(round(trip_coverage * 100, 0), "% coverage")
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  # Reference lines for standard alpha levels
  geom_hline(yintercept = c(0.05), linetype = 'dashed', linewidth = 1) +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
fpr_plot

# Below is for making an exact binomial figure
{## 1. Calculate Exact Binomial Confidence Intervals ----
# fpr_stats <- res_fp %>%
#   group_by(test) %>%
#   summarize(
#     total_runs = n(),
#     false_positives = sum(sig, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     # Run an exact binomial test against the expected 0.05 rate
#     binom_res = map2(false_positives, total_runs, ~binom.test(x = .x, n = .y, p = 0.05)),
#     # Extract the point estimate, p-value, and 95% CI into columns
#     tidied = map(binom_res, tidy)
#   ) %>%
#   unnest(tidied) %>%
#   select(test, estimate, conf.low, conf.high, p.value)
#  # For extracting tidy test results
# 
# ## 2. Generate the Point-Range Plot ----
# fpr_plot <- fpr_stats %>%
#   # Order the tests from highest false positive rate to lowest for clean reading
#   mutate(test = fct_reorder(test, estimate)) %>%
#   ggplot(aes(y = test, x = estimate)) +
#   
#   # Add the error bars (95% CI) and the point estimate using the modern geom_errorbar
#   geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, linewidth = 0.8) +
#   geom_point(size = 3, color = "black") +
#   geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
#   
#   theme_classic() +
#   labs(y = NULL, x = "False Positive Rate (95% Confidence Interval)") +
#   
#   # Format x-axis as percentages
#   scale_x_continuous(labels = scales::percent_format()) +
#   
#   theme(
#     axis.text = element_text(size = 11, color = "black"),
#     axis.title.x = element_text(margin = margin(t = 10))
#   )
# fpr_plot
  }

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