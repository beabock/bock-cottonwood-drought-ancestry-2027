# ==============================================================================
# Inoculum-vs-control contrasts and their multiplicity.
# Models as in r_code/08_sample_flow_extra_tests.R and 186-210: reduced LMM + myco_cat:treatment; contrasts
#   contrast(emmeans(m_int, ~ myco_cat | treatment, lmer.df = "satterthwaite"), "trt.vs.ctrl")   (08_sample_flow_extra_tests.R)
# For every contrast: unadjusted P, the reported adjustment (emmeans default for trt.vs.ctrl = "dunnettx", Hsu's
# approximation), the exact multivariate-t Dunnett ("mvt"), and Holm across all 18 contrasts.
# Run from the project root:  Rscript r_code/17_inoculation_contrasts.R
# ==============================================================================
suppressMessages({ source("r_code/helpers/utilities_data_prep.R"); library(tidyverse); library(lmerTest); library(emmeans) })
out <- "results/17_inoculation_contrasts"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
cat("INOCULATION CONTRASTS run:", "|", R.version.string, "| emmeans", as.character(packageVersion("emmeans")), "\n")
cat("emmeans default adjustment for 'trt.vs.ctrl':", emmeans:::trt.vs.ctrl.emmc(1:4)$.offset %||% "", attr(emmeans:::trt.vs.ctrl.emmc(letters[1:4]), "adjust"), "\n\n")
set.seed(1998)
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, geno_id = paste(pop, type, geno, sep = "_"))
bio_red <- lmerTest::lmer(log(dryshoot_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
al <- read_csv("results/06_allocation_models/allocation_data.csv", show_col_types = FALSE) %>% select(tree_idold, total_new, scored)
ta <- ds %>% inner_join(al, by = "tree_idold") %>% mutate(log_total_new = log(total_new)); stopifnot(nrow(ta) == 414)
tot_a <- lmerTest::lmer(log_total_new ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ta)
tot_b <- update(tot_a, data = filter(ta, !scored))
mods <- list(`New stem biomass (log)` = bio_red, `Total new biomass (log), (a) estimated root mass` = tot_a, `Total new biomass (log), (b) measured root mass only` = tot_b)
rows <- list()
for (k in names(mods)) {
  m_int <- update(mods[[k]], . ~ . + myco_cat:treatment)
  em <- emmeans(m_int, ~ myco_cat | treatment, lmer.df = "satterthwaite")
  ct <- contrast(em, "trt.vs.ctrl")
  cat(k, "| n =", nobs(m_int), "| adjustment attached by emmeans:", attr(summary(ct), "mesg")[grepl("adjust", attr(summary(ct), "mesg"))], "\n")
  s_none <- as.data.frame(summary(ct, adjust = "none", infer = c(TRUE, TRUE)))
  s_def <- as.data.frame(summary(ct, infer = c(TRUE, TRUE)))
  s_mvt <- as.data.frame(summary(ct, adjust = "mvt", infer = c(TRUE, TRUE)))
  rows[[k]] <- s_none %>% transmute(response = k, n = nobs(m_int), treatment, contrast, ratio = exp(estimate), t = t.ratio, df, P_unadjusted = p.value) %>%
    mutate(P_dunnettx_reported = s_def$p.value, ratio_lo_dunnettx = exp(s_def$lower.CL), ratio_hi_dunnettx = exp(s_def$upper.CL), P_dunnett_mvt = s_mvt$p.value)
}
tab <- bind_rows(rows) %>% mutate(P_holm_all18 = p.adjust(P_unadjusted, "holm"), P_bonferroni_all18 = p.adjust(P_unadjusted, "bonferroni"))
stopifnot(nrow(tab) == 18)
print(as.data.frame(tab %>% mutate(across(where(is.double), ~ signif(.x, 4)))), right = FALSE, row.names = FALSE)
cat("\nContrasts with Dunnett-adjusted P < 0.05:\n"); print(as.data.frame(tab %>% filter(P_dunnettx_reported < 0.05) %>% select(response, treatment, contrast, ratio, P_unadjusted, P_dunnettx_reported, P_dunnett_mvt, P_holm_all18)))
cat("\nContrasts with unadjusted P < 0.05:", sum(tab$P_unadjusted < 0.05), "| Holm across 18 < 0.05:", sum(tab$P_holm_all18 < 0.05), "\n")
write.csv(tab, file.path(out, "inoculation_contrasts.csv"), row.names = FALSE)
sink()
