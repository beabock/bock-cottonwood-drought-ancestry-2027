# ==============================================================================
# Porometer (stomatal conductance) by leaf side and measurement day.
# Raw file gh_poromeasurements_allgroups.csv: per tree, leaves 1-3 x top/bottom,
# measured on a "wet" day (watering day, 9/22/2021) and the "dry" (driest) day
# (9/17/2021), all at 12 PM (metadata workbook). 01_merge_datasets.R reshapes this to
# l{1-3}_{t,b}_{wet,dry} and builds por_avg = mean(c(por_bottoms, por_tops, por_drys, por_wets)).
# Exploratory: same model structure as the log-scale stomatal model of 02a, fit per side.
# Run from the project root:  Rscript r_code/03_conductance_leaf_sides.R
# Writes results/03_conductance_leaf_sides/output.txt
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(lme4); library(lmerTest); library(tidyverse)
sink("results/03_conductance_leaf_sides/output.txt", split = TRUE)
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, geno_id = paste(pop, type, geno, sep = "_"))
cols <- as.vector(outer(outer(paste0("l", 1:3), c("t", "b"), paste, sep = "_"), c("wet", "dry"), paste, sep = "_"))
stopifnot(all(cols %in% names(ds)))
p <- ds %>% filter(!is.na(por_avg), por_avg > 0) %>% mutate(across(all_of(cols), as.numeric))
cat("Trees in the stomatal model of 02a (por_avg > 0):", nrow(p), "\n\n")

n_side_day <- sapply(c("t_wet", "t_dry", "b_wet", "b_dry"), function(s) rowSums(!is.na(p[, grep(paste0("_", s, "$"), cols, value = TRUE)])))
cat("Leaves measured per tree (0-3), by side and day:\n"); print(apply(n_side_day, 2, function(x) table(factor(x, levels = 0:3))))
cat("\nTrees with at least one TOP and one BOTTOM reading:", sum(rowSums(n_side_day[, c("t_wet", "t_dry")]) > 0 & rowSums(n_side_day[, c("b_wet", "b_dry")]) > 0), "\n")
cat("Trees with both days:", sum(rowSums(n_side_day[, c("t_wet", "b_wet")]) > 0 & rowSums(n_side_day[, c("t_dry", "b_dry")]) > 0), "\n")
cat("Readings per tree (of 12): "); print(summary(rowSums(n_side_day)))

p <- p %>% mutate(top = rowMeans(across(matches("^l[1-3]_t_(wet|dry)$")), na.rm = TRUE),
                  bottom = rowMeans(across(matches("^l[1-3]_b_(wet|dry)$")), na.rm = TRUE),
                  wet = rowMeans(across(matches("^l[1-3]_[tb]_wet$")), na.rm = TRUE),
                  dry = rowMeans(across(matches("^l[1-3]_[tb]_dry$")), na.rm = TRUE),
                  grand = rowMeans(across(all_of(cols)), na.rm = TRUE))
cat("\nMeans across trees: top", round(mean(p$top), 4), "| bottom", round(mean(p$bottom), 4),
    "| top/bottom ratio (median per tree)", round(median(p$top / p$bottom), 2), "| trees with top > bottom:", sum(p$top > p$bottom), "of", nrow(p), "\n")
cat("Wet day", round(mean(p$wet, na.rm = TRUE), 4), "| dry day", round(mean(p$dry, na.rm = TRUE), 4), "\n")
cat("Trees with a non-positive side mean (negative porometer readings): top", sum(p$top <= 0), "| bottom", sum(p$bottom <= 0), "\n")
pos <- p$top > 0 & p$bottom > 0
cat("Correlation top vs bottom (log, both positive, n =", sum(pos), "):", round(cor(log(p$top[pos]), log(p$bottom[pos])), 2), "\n")
cat("Share of the grand mean coming from the top surface (median):", round(median(p$top / (p$top + p$bottom)), 2), "\n")
cat("por_avg vs simple mean of all readings: max abs difference", signif(max(abs(p$por_avg - p$grand)), 3),
    "(identical here because missing readings are whole days, so por_avg is the plain mean of all readings)\n")

cat("\nDrought vs well-watered, side x day means:\n")
print(p %>% group_by(treatment) %>% summarise(top_wet = mean(rowMeans(across(matches("^l[1-3]_t_wet$")), na.rm = TRUE), na.rm = TRUE),
                                              top_dry = mean(rowMeans(across(matches("^l[1-3]_t_dry$")), na.rm = TRUE), na.rm = TRUE),
                                              bottom_wet = mean(rowMeans(across(matches("^l[1-3]_b_wet$")), na.rm = TRUE), na.rm = TRUE),
                                              bottom_dry = mean(rowMeans(across(matches("^l[1-3]_b_dry$")), na.rm = TRUE), na.rm = TRUE), .groups = "drop") %>%
        mutate(across(-treatment, ~ round(., 4))) %>% as.data.frame())

fit <- function(y) lmerTest::lmer(as.formula(paste0("log(", y, ") ~ per_fc * treatment + I(per_fc^2) + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id)")),
                                   data = p %>% filter(.data[[y]] > 0))
cat("\nSame model as the log-scale stomatal model of 02a, fit to each side separately:\n")
for (y in c("top", "bottom", "grand")) {
  m <- fit(y); cs <- coef(summary(m))
  cat(sprintf("\n%s (n = %d)\n", y, nobs(m)))
  print(round(cs[c("per_fc", "I(per_fc^2)", "treatmentDrought", "per_fc:treatmentDrought"), c("Estimate", "df", "t value", "Pr(>|t|)")], 5))
}
long <- p %>% select(tree_id, pop, geno_id, per_fc, treatment, og_drc_mm, soil, myco_cat, top, bottom) %>%
  pivot_longer(c(top, bottom), names_to = "side", values_to = "gsw") %>% filter(gsw > 0) %>% mutate(side = factor(side, levels = c("bottom", "top")))
ms <- lmerTest::lmer(log(gsw) ~ per_fc * side + I(per_fc^2) + per_fc * treatment + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id) + (1 | tree_id), data = long)
cat("\nBoth sides in one model, ancestry x side term (does the ancestry slope differ between surfaces?):\n")
print(round(coef(summary(ms))[c("per_fc", "sidetop", "per_fc:sidetop"), c("Estimate", "df", "t value", "Pr(>|t|)")], 5))
sink()
