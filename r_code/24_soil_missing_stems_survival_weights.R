# ==============================================================================
# Soil legacy, survivors without new stem mass, and survivorship weighting.
#   6a  soil legacy term in the reduced models of record: survival GLMM (r_code/09_allocation_survival_followup.R), new stem LMM (line 161),
#       total new biomass (both root-mass samples; as r_code/17_inoculation_contrasts.R) and leaf drop GLMM (09_allocation_survival_followup.R).
#   6b  the 15 survivors with no recorded new stem mass, by watering, inoculum, soil, population and ancestry tercile;
#       Fisher exact tests of missingness among the 441 survivors. Nothing is imputed.
#   6c  inverse-probability-of-survival weights from the survival GLMM: w = mean(survival) / p_i, with p_i the fitted probability
#       (genotype random effect included; a second version uses population-level probabilities). Weights enter lmer as
#       precision weights, so this is an approximate sensitivity check (no robust SEs).
# Run from the project root:  Rscript r_code/24_soil_missing_stems_survival_weights.R
# Writes results/24_soil_missing_stems_survival_weights/: output.txt, strings.csv, missing_stem_table.csv, ipw.csv
# ==============================================================================
source("r_code/helpers/model_helpers.R")
OUT <- "results/24_soil_missing_stems_survival_weights"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "output.txt"), split = TRUE); run_stamp("Soil, missing stems, survival weights")
src <- "results/24_soil_missing_stems_survival_weights/strings.csv"
cat("Soil levels (first = reference):", paste(levels(factor(ds$soil)), collapse = ", "), "\n")
red <- paste("per_fc + treatment + I(per_fc^2) +", cov, "+", re)
surv_red <- glmer(alive_01_10_2021 ~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat + (1 | geno_id), data = ds, family = binomial, control = ctl)
bio_red <- lmerTest::lmer(as.formula(paste("log(dryshoot_new) ~", red)), data = ds)
ta <- al %>% mutate(log_total_new = log(total_new))
tot_a <- lmerTest::lmer(as.formula(paste("log_total_new ~", red)), data = ta); tot_b <- update(tot_a, data = filter(ta, !scored))
lz <- ds %>% filter(alive_01_10_2021 %in% 1, !is.na(dryshoot_new)) %>% mutate(zero_leaf = as.integer(coalesce(dryleaf_weight, 0) == 0)) %>%
  filter(if_all(c(pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), ~ !is.na(.x)))
leaf_m <- glmer(zero_leaf ~ pf50 + treatment + I(pf50^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = lz, family = binomial, control = ctl)
stopifnot(nobs(surv_red) == 563, nobs(bio_red) == 425, nobs(tot_a) == 414, nobs(tot_b) == 242, nobs(leaf_m) == 425)
st <- grep("^soil", names(fixef(bio_red)), value = TRUE); stopifnot(length(st) == 1)

hdr("6a. Soil legacy (tamarisk-legacy vs. control soil) in the reduced models")
for (k in list(list("Survival (binomial GLMM)", surv_red, "z", "odds ratio"), list("New stem biomass (log)", bio_red, "t", "×"),
               list("Total new biomass (log), (a) estimated root mass", tot_a, "t", "×"), list("Total new biomass (log), (b) measured root mass only", tot_b, "t", "×"),
               list("Leafless at harvest (binomial GLMM)", leaf_m, "z", "odds ratio"))) {
  r <- if (k[[3]] == "z") row_z(k[[2]], st) else row_t(k[[2]], st)
  rec(6, "a", k[[1]], "reduced model of record", nobs(k[[2]]), sprintf("Soil legacy (%s)", sub("^soil", "", st)),
      sprintf("%s; %s %.2f", r$s, ifelse(k[[4]] == "odds ratio", "odds ratio", "ratio"), exp(r$b)), r$p, src) }

hdr("6b. Survivors without a recorded new stem mass")
sv <- ds %>% filter(alive_01_10_2021 %in% 1) %>% mutate(missing = is.na(dryshoot_new),
  tercile = cut(per_f, quantile(per_f, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("lower third", "middle third", "upper third")))
cat("Survivors:", nrow(sv), "| missing new stem mass:", sum(sv$missing), "| ancestry tercile cut points:", paste(round(quantile(sv$per_f, c(1/3, 2/3)), 1), collapse = ", "), "% Fremont\n")
tabs <- map_dfr(c(watering = "treatment", inoculum = "myco_cat", soil = "soil", population = "pop", ancestry_tercile = "tercile"), function(v) {
  t <- table(sv[[v]], sv$missing); p <- fisher.test(t, simulate.p.value = nrow(t) > 2, B = 1e5)$p.value
  tibble(level = rownames(t), survivors = rowSums(t), missing = t[, "TRUE"], pct_missing = round(100 * t[, "TRUE"] / rowSums(t), 1), fisher_P = p) }, .id = "factor")
print(as.data.frame(tabs), row.names = FALSE); readr::write_excel_csv(tabs, file.path(OUT, "missing_stem_table.csv"))
for (f in unique(tabs$factor)) { x <- filter(tabs, factor == f)
  rec(6, "b", "Survivors without new stem mass", "441 survivors", 441, paste("Missingness by", gsub("_", " ", f)),
      sprintf("%s; Fisher exact %s", paste(sprintf("%s %d/%d", x$level, x$missing, x$survivors), collapse = ", "), P(x$fisher_P[1])), x$fisher_P[1], "results/24_soil_missing_stems_survival_weights/missing_stem_table.csv") }
gm <- glm(missing ~ per_f + treatment + soil + myco_cat, data = sv, family = binomial)
cat("\nLogistic regression of missingness (all factors together):\n"); print(round(coef(summary(gm)), 4))
cat("Notes recorded for these plants are in results/09_allocation_survival_followup/table_unharvested_survivors.csv.\n")

hdr("6c. Inverse-probability-of-survival weights")
sp <- ds %>% filter(if_all(c(alive_01_10_2021, pf50, treatment, pop, og_drc_mm, soil, myco_cat, geno_id), ~ !is.na(.x)))
sp$p_cond <- fitted(surv_red); sp$p_marg <- predict(surv_red, re.form = NA, type = "response")
pbar <- mean(sp$alive_01_10_2021)
wt <- sp %>% transmute(tree_idold, w_cond = pbar / p_cond, w_marg = pbar / p_marg)
cat(sprintf("Stabilized weights among survivors: conditional %.2f–%.2f, population-level %.2f–%.2f\n",
            min(wt$w_cond[sp$alive_01_10_2021 == 1]), max(wt$w_cond[sp$alive_01_10_2021 == 1]), min(wt$w_marg[sp$alive_01_10_2021 == 1]), max(wt$w_marg[sp$alive_01_10_2021 == 1])))
ipw <- list()
for (k in list(list("New stem biomass (log)", "log(dryshoot_new)", ds), list("Total new biomass (log), (a) estimated root mass", "log_total_new", ta))) {
  d <- k[[3]] %>% inner_join(wt, by = "tree_idold")
  for (w in c("none (reported)", "w_cond", "w_marg")) {
    dd <- d; dd$w <- if (w == "none (reported)") 1 else dd[[w]] / mean(dd[[w]][!is.na(dd[[sub("log\\((.*)\\)", "\\1", k[[2]])]])], na.rm = TRUE)
    m <- lmerTest::lmer(as.formula(paste(k[[2]], "~", red)), data = dd, weights = w)
    wl <- c("none (reported)" = "unweighted (reported)", w_cond = "IPW, fitted survival probability (with genotype effect)", w_marg = "IPW, population-level survival probability")[[w]]
    for (tm in c("per_fc", "I(per_fc^2)", "treatmentDrought")) { r <- row_t(m, tm); lab <- c(per_fc = "Linear ancestry (per 1%)", "I(per_fc^2)" = "Quadratic ancestry", treatmentDrought = "Drought")[[tm]]
      ipw[[length(ipw) + 1]] <- tibble(response = k[[1]], weights = wl, n = nobs(m), term = lab, string = r$s, P = r$p)
      rec(6, "c", k[[1]], wl, nobs(m), lab, paste0(r$s, ifelse(tm == "treatmentDrought", sprintf(" (%s%%)", pc(100 * (exp(r$b) - 1))), "")), r$p, "results/24_soil_missing_stems_survival_weights/ipw.csv") } } }
readr::write_excel_csv(bind_rows(ipw), file.path(OUT, "ipw.csv"))
write_rec("strings.csv"); cat("\nDone.\n"); sink()
