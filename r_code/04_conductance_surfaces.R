# ==============================================================================
# CONDUCTANCE BY LEAF SURFACE, AND SURVIVAL GROUPING: Appendix S5 (aggregated-survival grouping sensitivity table) and
# surface-specific stomatal conductance (adaxial = top, abaxial = bottom; the
# stomatal density data use the same convention).
# ==============================================================================
# Run from the project root:  Rscript r_code/04_conductance_surfaces.R   (about 1 min)
# Writes: results/04_conductance_surfaces/output.txt,
#         results/04_conductance_surfaces/surface_sensitivity.csv, and appends a marked
#         section to results/02_main_models/formatted_results_main_models_df.md (replaced on every run).
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(lme4); library(lmerTest); library(tidyverse)
sink("results/04_conductance_surfaces/output.txt", split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("CONDUCTANCE SURFACES run:", "\n")
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, geno_id = paste(pop, type, geno, sep = "_"))

# AJB helpers (as in 02c_main_models_df.R)
EN <- "–"; BETA <- "β"; TIMES <- "×"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
tdf <- function(t, df, p) sprintf("t = %s, df = %s, %s", dec(t), dff(df), P(p))     # format requested for the new strings

# ------------------------------------------------------------------------------
# 1. APPENDIX S5: aggregated survival, grouping sensitivity
# ------------------------------------------------------------------------------
hdr("1. Appendix S5: aggregated survival grouping sensitivity")
build <- function(gvars) {
  cells <- ds %>% group_by(across(all_of(gvars))) %>%
    summarise(n_total = n(), survival_rate = sum(alive_01_10_2021, na.rm = TRUE) / n(),
              prop_tam = mean(soil == "Tam Legacy"), prop_AM = mean(myco_cat == "AM"), prop_EM = mean(myco_cat == "EM"),
              prop_Both = mean(myco_cat == "Both"), .groups = "drop") %>% filter(n_total >= 3)
  colz <- ds %>% filter(alive_01_10_2021 == 1) %>% mutate(sc = !is.na(amf) & !is.na(emf) & !is.na(dse)) %>%
    group_by(across(all_of(gvars))) %>%
    summarise(n_scored = sum(sc), mean_emf = mean(emf, na.rm = TRUE), mean_amf = mean(amf, na.rm = TRUE), mean_dse = mean(dse, na.rm = TRUE),
              mean_per_f = mean(per_f, na.rm = TRUE), mean_drc_mm = mean(og_drc_mm, na.rm = TRUE), .groups = "drop")
  inner_join(cells, colz, by = gvars) %>% mutate(mean_per_fc = mean_per_f - 50) %>%
    filter(is.finite(mean_emf), is.finite(mean_amf), is.finite(mean_dse))
}
fungi <- c(AMF = "scale(mean_amf):treatmentDrought", EMF = "scale(mean_emf):treatmentDrought", DSE = "scale(mean_dse):treatmentDrought")
fit_agg <- function(d, covs, label) {
  out <- tibble(grouping = label, n_groups = nrow(d), n_plants = sum(d$n_total),
                scored_median = if (nrow(d)) median(d$n_scored) else NA_real_,
                scored_min = if (nrow(d)) min(d$n_scored) else NA_real_, scored_max = if (nrow(d)) max(d$n_scored) else NA_real_,
                dispersion = NA_real_, estimable = FALSE, note = "")
  f <- as.formula(paste("survival_rate ~ (scale(mean_emf) + scale(mean_dse) + scale(mean_amf)) * treatment + mean_per_fc + I(mean_per_fc^2) + mean_drc_mm +", covs))
  npar <- if (nrow(d) && n_distinct(d$treatment) == 2) ncol(model.matrix(f, d)) else NA
  if (is.na(npar) || nrow(d) < npar + 5) {
    out$note <- if (nrow(d) == 0) "no groups" else if (is.na(npar)) sprintf("not estimable: %d groups, only one watering treatment", nrow(d)) else sprintf("not estimable: %d groups for %d parameters", nrow(d), npar)
    return(out)
  }
  m <- glm(f, data = d, weights = n_total, family = binomial); cs <- coef(summary(m))
  mq <- glm(f, data = d, weights = n_total, family = quasibinomial); cq <- coef(summary(mq))
  out$dispersion <- sum(residuals(m, type = "pearson")^2) / df.residual(m); out$estimable <- TRUE
  for (k in names(fungi)) {
    out[[paste0(k, "_beta")]] <- cs[fungi[k], "Estimate"]; out[[paste0(k, "_z")]] <- cs[fungi[k], "z value"]
    out[[paste0(k, "_P_wald")]] <- cs[fungi[k], "Pr(>|z|)"]; out[[paste0(k, "_P_quasi")]] <- cq[fungi[k], "Pr(>|t|)"]
  }
  out$quasi_df <- df.residual(mq)
  out
}
cur <- build(c("pop", "new_class", "treatment", "myco_cat", "soil")); coa <- build(c("pop", "new_class", "treatment"))
s4 <- bind_rows(
  fit_agg(cur, "soil + myco_cat", paste("Population", TIMES, "cross type", TIMES, "drought", TIMES, "inoculum", TIMES, "soil (reported grouping)")),
  fit_agg(filter(cur, n_scored >= 2), "soil + myco_cat", "Reported grouping, groups with ≥ 2 scored plants"),
  fit_agg(filter(cur, n_scored >= 3), "soil + myco_cat", "Reported grouping, groups with ≥ 3 scored plants"),
  fit_agg(coa, "prop_tam + prop_AM + prop_EM + prop_Both", paste("Population", TIMES, "cross type", TIMES, "drought (inoculum and soil as group proportions)")),
  fit_agg(filter(coa, n_scored >= 2), "prop_tam + prop_AM + prop_EM + prop_Both", "Coarser grouping, groups with ≥ 2 scored plants"),
  fit_agg(filter(coa, n_scored >= 3), "prop_tam + prop_AM + prop_EM + prop_Both", "Coarser grouping, groups with ≥ 3 scored plants"))
print(as.data.frame(s4 %>% mutate(across(where(is.numeric), ~ signif(., 3)))), right = FALSE)
# The grouping-sensitivity table is printed above only (not an appendix table).

# ------------------------------------------------------------------------------
# 2. SURFACE-SPECIFIC STOMATAL CONDUCTANCE
# ------------------------------------------------------------------------------
hdr("2. Surface-specific stomatal conductance (adaxial = top, abaxial = bottom)")
side_mean <- function(d, side, days) {
  cc <- as.vector(outer(paste0("l", 1:3, "_", side), days, paste, sep = "_"))
  m <- rowMeans(as.matrix(sapply(d[, cc], as.numeric)), na.rm = TRUE); m[is.nan(m)] <- NA; m }
has_day <- function(d, day) rowSums(!is.na(as.matrix(sapply(d[, as.vector(outer(paste0("l", 1:3), c("t", "b"), paste, sep = "_")) %>% paste0("_", day)], as.numeric)))) > 0
por <- ds %>% filter(!is.na(por_avg)) %>%
  mutate(adaxial = side_mean(., "t", c("wet", "dry")), abaxial = side_mean(., "b", c("wet", "dry")),
         adaxial_wet = side_mean(., "t", "wet"), abaxial_wet = side_mean(., "b", "wet"),
         both_days = has_day(., "wet") & has_day(., "dry"))
cat("Trees with porometer data:", nrow(por), "| with both days:", sum(por$both_days), "\n")
cat("Non-positive surface means (cannot be logged): adaxial", sum(por$adaxial <= 0, na.rm = TRUE), "| abaxial", sum(por$abaxial <= 0, na.rm = TRUE),
    "| adaxial wet-day", sum(por$adaxial_wet <= 0, na.rm = TRUE), "| abaxial wet-day", sum(por$abaxial_wet <= 0, na.rm = TRUE), "\n")

rhs <- "per_fc * treatment + I(per_fc^2) + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id)"
fit_side <- function(y, d) {
  d <- d %>% filter(!is.na(.data[[y]]), .data[[y]] > 0)
  lmerTest::lmer(as.formula(paste0("log(", y, ") ~ ", rhs)), data = d)
}
term_lab <- c("per_fc" = "Linear ancestry (per 1% Fremont)", "I(per_fc^2)" = "Quadratic ancestry", "treatmentDrought" = "Drought (at 50% ancestry)",
              "per_fc:treatmentDrought" = paste0("Ancestry ", TIMES, " drought"), "og_drc_mm" = "Initial cutting diameter (per mm)",
              "soilTam Legacy" = "Soil legacy (tamarisk vs. cottonwood)", "myco_catAM" = "Inoculum: AMF vs. control",
              "myco_catEM" = "Inoculum: EMF vs. control", "myco_catBoth" = "Inoculum: dual vs. control")
slopes <- function(m) {
  b <- fixef(m); V <- as.matrix(vcov(m))
  s <- c(ww = unname(b["per_fc"]), dr = unname(b["per_fc"] + b["per_fc:treatmentDrought"]))
  se <- c(ww = sqrt(V["per_fc", "per_fc"]), dr = sqrt(V["per_fc", "per_fc"] + V["per_fc:treatmentDrought", "per_fc:treatmentDrought"] + 2 * V["per_fc", "per_fc:treatmentDrought"]))
  tibble(treatment = c("Well-watered", "Drought"), slope_log = s, se = se,
         pct_per_10 = 100 * (exp(10 * s) - 1), lo = 100 * (exp(10 * (s - 1.96 * se)) - 1), hi = 100 * (exp(10 * (s + 1.96 * se)) - 1))
}
models <- list(adaxial = fit_side("adaxial", por), abaxial = fit_side("abaxial", por))
term_rows <- list(); slope_rows <- list()
for (s in names(models)) {
  m <- models[[s]]; cs <- coef(summary(m))
  cat(sprintf("\n--- %s: n = %d trees, %d genotypes, %d populations; singular = %s ---\n", s, nobs(m), ngrps(m)["pop:geno_id"], ngrps(m)["pop"], isSingular(m)))
  print(round(cs[names(term_lab), c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")], 5)); print(VarCorr(m))
  term_rows[[s]] <- tibble(surface = s, term = names(term_lab), label = term_lab, beta = cs[names(term_lab), "Estimate"],
                           t = cs[names(term_lab), "t value"], df = cs[names(term_lab), "df"], p = cs[names(term_lab), "Pr(>|t|)"])
  sl <- slopes(m); print(sl %>% mutate(across(where(is.numeric), ~ round(., 4))))
  slope_rows[[s]] <- sl %>% mutate(surface = s)
}
terms_tab <- bind_rows(term_rows); slope_tab <- bind_rows(slope_rows)

hdr("2b. Ancestry x surface model (both surfaces in one model)")
long <- por %>% select(tree_id, pop, geno_id, per_fc, treatment, og_drc_mm, soil, myco_cat, adaxial, abaxial) %>%
  pivot_longer(c(adaxial, abaxial), names_to = "surface", values_to = "gsw") %>%
  filter(!is.na(gsw), gsw > 0) %>% mutate(surface = factor(surface, levels = c("abaxial", "adaxial")))
joint <- lmerTest::lmer(log(gsw) ~ surface * (per_fc * treatment + I(per_fc^2)) + og_drc_mm + soil + myco_cat +
                          (1 | pop) + (1 | pop:geno_id) + (1 | tree_id), data = long)
jc <- coef(summary(joint))
cat("Formula: log(gsw) ~ surface * (per_fc * treatment + I(per_fc^2)) + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id) + (1 | tree_id)\n")
mf <- model.frame(joint)
cat(sprintf("n = %d surface-level observations from %d trees (%d with both surfaces), %d genotypes, %d populations; singular = %s\n",
            nobs(joint), n_distinct(long$tree_id[as.integer(rownames(mf))]), sum(table(long$tree_id) == 2), ngrps(joint)["pop:geno_id"], ngrps(joint)["pop"], isSingular(joint)))
print(round(jc, 5)); print(VarCorr(joint))

hdr("2c. Sensitivity: watering-day readings only; plants with both days only")
sens <- list(
  `All readings (main)` = list(adaxial = models$adaxial, abaxial = models$abaxial),
  `Watering-day readings only` = list(adaxial = fit_side("adaxial_wet", por), abaxial = fit_side("abaxial_wet", por)),
  `Plants with both days only` = list(adaxial = fit_side("adaxial", filter(por, both_days)), abaxial = fit_side("abaxial", filter(por, both_days))))
key <- c("per_fc", "I(per_fc^2)", "treatmentDrought", "per_fc:treatmentDrought")
sens_tab <- map_dfr(names(sens), function(a) map_dfr(c("adaxial", "abaxial"), function(s) {
  m <- sens[[a]][[s]]; cs <- coef(summary(m))
  tibble(analysis = a, surface = s, n = nobs(m), term = key, beta = cs[key, "Estimate"], t = cs[key, "t value"], df = cs[key, "df"], p = cs[key, "Pr(>|t|)"])
})) %>% group_by(surface, term) %>%
  mutate(sig = p < 0.05, same_sign_as_main = sign(beta) == sign(beta[analysis == "All readings (main)"]),
         same_significance_as_main = sig == sig[analysis == "All readings (main)"]) %>% ungroup()
print(as.data.frame(sens_tab %>% mutate(across(c(beta, t, df, p), ~ signif(., 3)))), right = FALSE)
cat("\nConclusion changes (sign or P < 0.05 status differs from the main model):\n")
chg <- sens_tab %>% filter(!same_sign_as_main | !same_significance_as_main)
if (nrow(chg)) print(as.data.frame(chg %>% select(analysis, surface, term, beta, p) %>% mutate(across(c(beta, p), ~ signif(., 3)))), right = FALSE) else cat("none\n")
write.csv(sens_tab, "results/04_conductance_surfaces/surface_sensitivity.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Strings appended to formatted_results_main_models_df.md (marked block, replaced each run)
# ------------------------------------------------------------------------------
blk <- c("<!-- BEGIN surface-specific conductance (generated by r_code/04_conductance_surfaces.R) -->", "",
         "## Results: Stomatal conductance by leaf surface (adaxial = top, abaxial = bottom)", "",
         "Format for these rows: t = x, df = y, P = z (Satterthwaite df). Models: log surface-mean conductance (mean of three leaves on both measurement days) ~ linear and quadratic Fremont ancestry (centered at 50%), drought, linear ancestry × drought, initial cutting diameter, soil legacy, inoculum category; random intercepts for population and maternal genotype within population.", "",
         "| Surface | Term | Paste-ready string | Flag / note |", "|---|---|---|---|")
for (i in seq_len(nrow(terms_tab))) {
  r <- terms_tab[i, ]
  blk <- c(blk, sprintf("| %s (n = %d) | %s | %s = %s, %s | %s |", ifelse(r$surface == "adaxial", "Adaxial", "Abaxial"), nobs(models[[r$surface]]), r$label,
                        BETA, num(r$beta), tdf(r$t, r$df, r$p),
                        if (r$term == "per_fc" && r$surface == "abaxial") "opposite sign to adaxial and to the pooled model" else ""))
}
for (i in seq_len(nrow(slope_tab))) {
  r <- slope_tab[i, ]
  blk <- c(blk, sprintf("| %s | Back-transformed slope, %s | %s%% per 10 percentage points of Fremont ancestry (95%% CI %s to %s%%) | average slope from 0 to 100%% ancestry; negative = decline |",
                        ifelse(r$surface == "adaxial", "Adaxial", "Abaxial"), tolower(r$treatment), dec(r$pct_per_10, 1), dec(r$lo, 1), dec(r$hi, 1)))
}
jr <- function(t) sprintf("%s = %s, %s", BETA, num(jc[t, "Estimate"]), tdf(jc[t, "t value"], jc[t, "df"], jc[t, "Pr(>|t|)"]))
blk <- c(blk, sprintf("| Both | Ancestry %s surface (joint model) | %s | difference in linear ancestry slope, adaxial minus abaxial, well-watered |", TIMES, jr("surfaceadaxial:per_fc")),
         sprintf("| Both | Ancestry %s drought %s surface (joint model) | %s | |", TIMES, TIMES, jr("surfaceadaxial:per_fc:treatmentDrought")), "",
         "<!-- END surface-specific conductance -->")
f <- "results/02_main_models/formatted_results_main_models_df.md"
old <- readLines(f, encoding = "UTF-8", warn = FALSE)
b0 <- grep("BEGIN surface-specific", old); b1 <- grep("END surface-specific", old)
if (length(b0) && length(b1)) old <- old[-(b0:b1)]
while (length(old) && old[length(old)] == "") old <- old[-length(old)]
con <- file(f, "w", encoding = "UTF-8"); writeLines(enc2utf8(c(old, "", blk, "")), con); close(con)
cat("\nAppended", length(blk), "lines to", f, "\n")
sink()
