assign_sig <- function(dat) {
  dat %>%
    mutate(t_test_sig = ifelse(tp < 0.05, 1, 0),
           f_test_sig = ifelse(Fp < 0.05, 1, 0),
           dAIC_sig = ifelse(glmmconv==1, ifelse(AICd > 2, 1, 0), NA),
           KS_test_sig = ifelse(KSp < 0.05, 1, 0),
           median_test_sig = ifelse(ci_lo > 0 | ci_hi < 0, 1, 0),
           mvglm_sig = ifelse(mvglm_p < 0.05, 1, 0),
           perma_sig = ifelse(perma_p < 0.05, 1, 0)
    ) 
}
