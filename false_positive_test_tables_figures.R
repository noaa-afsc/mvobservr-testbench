# Stitch sets  ===================================================================================================
# loads data previously run from false_positive_test.r

library(kableExtra)
library(ggplot2)
library(dplyr)

mvobservr_dribble <- gdrive_set_dribble(folder_id = "1Wh-ZQlJ3AIVaQZTWk4QNuyiMfoVECQgt")
load(gdrive_download(local_path = "output_data/alltests_falsepos_set1.Rdata", gdrive_dribble = mvobservr_dribble))
res_comb_tbl <- res_comb
load(gdrive_download(local_path = "output_data/param_plot_batch2.Rdata", gdrive_dribble = mvobservr_dribble))
#res_comb is now overwritten
res_comb <- res_comb |> mutate(set = set + max(res_comb_tbl$set))
res_comb_tbl <- rbind(res_comb_tbl, res_comb)
rm(res_comb)

##figure: boxplot for false positives ----
res_fp <- res_comb_tbl %>%
  mutate(
    # Standard p-value threshold is < 0.05
    t_test_sig = ifelse(tp < 0.05, 1, 0),
    f_test_sig = ifelse(Fp < 0.05, 1, 0),
    glm_mgcv_sig = ifelse(glm_mgcv_p <= 0.05, 1, 0), 
    KS_test_sig = ifelse(KSp < 0.05, 1, 0),
    permu_sig = ifelse(p_val < 0.05, 1, 0), #permutation
    median_test_sig = ifelse(ci_lo > 0 | ci_hi < 0, 1, 0),
    mvglm_sig = ifelse(p_mvobs < 0.05, 1, 0),
    gllvm_sig = ifelse(p_gllvm < 0.05, 1, 0),
    gllvm1_sig = ifelse(p_gllvm1< 0.05, 1, 0),
    perma_sig = ifelse(perma_p < 0.05, 1, 0), #permanova
    glmm_sig = ifelse(p_glmm < 0.05, 1, 0),
    glmm1_sig = ifelse(p_glmm1 < 0.05, 1, 0)
  ) %>%
  select(ends_with("sig")) %>%
  pivot_longer(
    cols      = ends_with("sig"), 
    names_to  = "test", 
    values_to = "sig"
  ) %>%
  mutate(test = factor(
    test, 
    levels  = c("t_test_sig", "glm_mgcv_sig", "f_test_sig", "permu_sig", "mvglm_sig",  
                "KS_test_sig", "perma_sig", "median_test_sig", "glmm_sig", "glmm1_sig", "gllvm_sig", "gllvm1_sig"),
    ordered = TRUE,
    labels  = c("t-test", "GLM", "F test", "Permutation","mvobservr",  
                "K-S test", "PERMANOVA", "Median", "glmm", "glmm1", "gllvm", "gllvm1")
  ))

set.seed(123) # Ensure reproducibility for the bootstrap
n_bootstraps <- 1000

fpr_plot <- map(1:n_bootstraps, ~{
  slice_sample(res_fp, n = n_samples_per_level, by = test, replace = TRUE) %>%
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
    y = "Percentage of Positive (Significant) Tests"#,
    # title = "False Positive Rate of Tests on Total Biomass",
    # subtitle = paste0(round(trip_coverage * 100, 0), "% coverage")
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  # Reference lines for standard alpha levels
  geom_hline(yintercept = c(0.05), linetype = 'dashed', linewidth = 1) +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
fpr_plot

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

# Function to turn those data frames into publication-ready kables
make_publication_kable <- function(df, caption_text, value_label = "") {
  
  # 1. Base table setup (keeping Mean and CI separated)
  kbl_out <- df %>% 
    mutate(CI = sprintf("(%s, %s)", ci_lower, ci_upper)) %>%
    select(Model, N, Mean, CI) %>%
    kbl(
      caption = caption_text, 
      booktabs = TRUE,
      align = c("l", "c", "r", "l"), # Model=l, N=c, Mean=r, CI=l
      col.names = c("Model", "N", "Mean", "(95% CI)")
    ) %>%
    kable_classic(full_width = FALSE, html_font = "Times New Roman")
  
  # 2. Only add the spanning header if value_label is NOT blank/null
  if (!is.null(value_label) && value_label != "") {
    spanning_header <- setNames(c(1, 1, 2), c(" ", " ", value_label))
    kbl_out <- kbl_out %>% add_header_above(spanning_header, bold = TRUE)
  }
  
  # 3. Clean up the final row headers
  kbl_out %>% row_spec(0, bold = TRUE)
}

# 3. Generate the underlying data
runtimes_df <- create_summary_table(res_comb_tbl, "runtime", "runtime_secs_")
power_df    <- create_summary_table(res_comb_tbl, "pwr", "pwr_")

## Table 1: Runtimes ----
runtimes_df %>% 
  make_publication_kable(
    caption_text = "Table 1: Computational Runtimes for False Positives", 
    value_label  = "" #Mean Runtime in Seconds (95% CI)"
  )

## Table 2: Tweedie Power ----
power_df %>% 
  make_publication_kable(
    caption_text = "Table 2: Tweedie Power Statistics for False Positives", 
    value_label  = "" #Mean Tweedie Power (95% CI)"
  )