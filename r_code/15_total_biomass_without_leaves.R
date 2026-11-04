# ==============================================================================
# Total new biomass without leaves (new stem + root), both allocation samples.
# Reported model: reduced LMM log(total_new) ~ pf50 + treatment + I(pf50^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id)
# (r_code/06_allocation_models.R; strings in results/07_allocation_appendix_S7/appendix_S7_source.csv). The same structure is
# fitted to log(new stem + root). Ancestry x drought is tested by LRT as in the reported analysis.
# Run from the project root:  Rscript r_code/15_total_biomass_without_leaves.R
# ==============================================================================
suppressMessages({ source("r_code/helpers/utilities_data_prep.R"); library(tidyverse); library(lmerTest) })
out <- "results/15_total_biomass_without_leaves"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
cat("TOTAL BIOMASS WITHOUT LEAVES run:", "|", R.version.string, "\n\n")
EN <- "–"; dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, pf50 = per_fc / 50, geno_id = paste(pop, type, geno, sep = "_"))
al <- read_csv("results/06_allocation_models/allocation_data.csv", show_col_types = FALSE) %>% select(tree_idold, scored, root_est, new_stem, leaf, total_new) %>%
  inner_join(ds %>% select(tree_idold, pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), by = "tree_idold") %>%
  mutate(stem_root = new_stem + root_est)
samp <- list(`(a) estimated root mass` = al, `(b) measured root mass only` = filter(al, !scored))
cov <- "og_drc_mm + myco_cat + soil"; re <- "(1 | pop) + (1 | pop:geno_id)"
red <- paste("pf50 + treatment + I(pf50^2) +", cov, "+", re); full <- paste("pf50 * treatment + I(pf50^2) * treatment +", cov, "+", re)
rows <- list()
for (s in names(samp)) for (y in c("log(total_new)", "log(stem_root)")) {
  d <- samp[[s]]; mr <- lmerTest::lmer(as.formula(paste(y, "~", red)), data = d); mf <- lmerTest::lmer(as.formula(paste(y, "~", full)), data = d)
  lr <- anova(mr, mf, refit = TRUE); cs <- coef(summary(mr))
  for (tm in c("pf50", "I(pf50^2)", "treatmentDrought")) {
    sc <- c(pf50 = 50, `I(pf50^2)` = 2500, treatmentDrought = 1)[[tm]]
    rows[[length(rows) + 1]] <- tibble(sample = s, response = ifelse(y == "log(total_new)", "Total new biomass (reported; stem + leaf + root)", "New stem + root (leaves excluded)"),
      n = nobs(mr), term = c(pf50 = "Linear ancestry (per 1%)", `I(pf50^2)` = "Quadratic ancestry", treatmentDrought = "Drought")[[tm]],
      estimate = cs[tm, 1] / sc, t_value = cs[tm, "t value"], df = cs[tm, "df"], P = cs[tm, "Pr(>|t|)"],
      string = sprintf("β = %s, t = %s, df = %s, %s", num(cs[tm, 1] / sc), dec(cs[tm, "t value"]), dff(cs[tm, "df"]), P(cs[tm, "Pr(>|t|)"])),
      back_transformed = ifelse(tm == "treatmentDrought", sprintf("drought × %.3f (%s%%)", exp(cs[tm, 1]), dec(100 * (exp(cs[tm, 1]) - 1), 1)), ""),
      lrt_ancestry_x_drought = sprintf("χ²(%d) = %s, %s", lr$Df[2], dec(lr$Chisq[2]), P(lr$`Pr(>Chisq)`[2])))
  } }
res <- bind_rows(rows)
print(as.data.frame(res %>% select(sample, response, n, term, string, back_transformed, lrt_ancestry_x_drought)), right = FALSE, row.names = FALSE)
write.csv(res, file.path(out, "total_without_leaves.csv"), row.names = FALSE)
sink()
