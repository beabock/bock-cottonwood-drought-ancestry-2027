# ==============================================================================
# Allocation framing.
#   4a  ancestry and drought effects on the three pairwise log ratios (root:new stem, leaf:new stem, leaf:root) among plants
#       with leaves; the first two reproduce results/14_allocation_fraction_checks/alr_models.csv (14_allocation_fraction_checks.R section 3b).
#   4b  predicted new stem, leaf and root dry mass at 0% and 100% Fremont ancestry (reduced models; both root-mass samples).
#   4d  why the dry:fresh ratio comes from 116 plants.
#   4e  total new biomass (and new stem + root) with one pooled dry:fresh ratio.
# Reduced model (as 14_allocation_fraction_checks.R): y ~ pf50 + treatment + I(pf50^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id);
# ancestry fitted per 50 points (pf50) and reported per 1% Fremont.
# Run from the project root:  Rscript r_code/22_allocation_log_ratios.R
# Writes results/22_allocation_log_ratios/: output.txt, strings.csv, log_ratios.csv, predicted_mass.csv, pooled_ratio.csv,
#        root_mass_counts.csv. (Figure 2 is drawn by r_code/28_figures.R.)
# ==============================================================================
source("r_code/helpers/model_helpers.R")
OUT <- "results/22_allocation_log_ratios"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "output.txt"), split = TRUE); run_stamp("Allocation log ratios")
set.seed(1998); src <- "results/22_allocation_log_ratios/strings.csv"
red <- paste("pf50 + treatment + I(pf50^2) +", cov, "+", re); full <- paste("pf50 * treatment + I(pf50^2) * treatment +", cov, "+", re)
samp <- list(`(a) estimated root mass` = al, `(b) measured root mass only` = filter(al, !scored))
fit_red <- function(y, d) { mr <- lmerTest::lmer(as.formula(paste(y, "~", red)), data = d); mf <- lmerTest::lmer(as.formula(paste(y, "~", full)), data = d)
  list(m = mr, lrt = anova(mr, mf, refit = TRUE)) }
lrt_s <- function(l) sprintf("%s(%d) = %s, %s", CHI, l$Df[2], dec(l$Chisq[2]), P(l$`Pr(>Chisq)`[2]))

hdr("4a. Pairwise log ratios among plants with leaves")
lr <- list()
for (s in names(samp)) {
  d <- samp[[s]] %>% filter(leaf > 0, root_est > 0, new_stem > 0) %>% mutate(root_stem = log(root_est / new_stem), leaf_stem = log(leaf / new_stem), leaf_root = log(leaf / root_est))
  for (y in c(root_stem = "log(root:new stem)", leaf_stem = "log(leaf:new stem)", leaf_root = "log(leaf:root)") %>% names()) {
    lab <- c(root_stem = "log(root:new stem)", leaf_stem = "log(leaf:new stem)", leaf_root = "log(leaf:root)")[[y]]
    f <- fit_red(y, d)
    for (tm in c("pf50", "treatmentDrought")) { r <- row_t(f$m, tm, ifelse(tm == "pf50", 50, 1))
      chg <- if (tm == "pf50") sprintf("%s%% per 10 points", pc(100 * (exp(10 * r$b) - 1))) else sprintf("drought × %.3f", exp(r$b))
      lr[[length(lr) + 1]] <- tibble(sample = s, n = nobs(f$m), response = lab, term = ifelse(tm == "pf50", "Linear ancestry (per 1%)", "Drought"), estimate = r$b, t = r$t, df = r$df, P = r$p,
                                     string = r$s, change = chg, lrt_ancestry_x_drought = lrt_s(f$lrt))
      rec(4, "a", paste(lab, "(plants with leaves)"), s, nobs(f$m), ifelse(tm == "pf50", "Linear ancestry (per 1%)", "Drought"), paste0(r$s, "; ", chg), r$p, "results/22_allocation_log_ratios/log_ratios.csv") }
    rq <- row_t(f$m, "I(pf50^2)", 2500); rec(4, "a", paste(lab, "(plants with leaves)"), s, nobs(f$m), "Quadratic ancestry", rq$s, rq$p, "results/22_allocation_log_ratios/log_ratios.csv")
    rec(4, "a", paste(lab, "(plants with leaves)"), s, nobs(f$m), paste("Ancestry", TIMES, "drought (LRT, 2 df)"), lrt_s(f$lrt), f$lrt$`Pr(>Chisq)`[2], "results/22_allocation_log_ratios/log_ratios.csv") } }
lr <- bind_rows(lr); readr::write_excel_csv(lr, file.path(OUT, "log_ratios.csv"))
old <- read_csv("results/14_allocation_fraction_checks/alr_models.csv", show_col_types = FALSE)
chk <- lr %>% filter(response != "log(leaf:root)") %>% mutate(response = ifelse(response == "log(root:new stem)", "alr_root", "alr_leaf")) %>%
  inner_join(old %>% select(sample, response, term, string_v55 = string), by = c("sample", "response", "term"))
cat("Check against alr_models.csv: ", sum(chk$string == chk$string_v55), "of", nrow(chk), "strings identical\n"); stopifnot(all(chk$string == chk$string_v55))

hdr("4b. Predicted dry mass (g) at 0% and 100% Fremont ancestry")
marg <- function(m, d, g, trt) { b <- fixef(m); nd <- d; nd$pf50 <- (g - 50) / 50; nd$treatment <- factor(trt, levels = levels(ds$treatment))
  fe <- delete.response(terms(as.formula(paste("y ~ pf50 + treatment + I(pf50^2) +", cov)))); x <- colMeans(model.matrix(fe, nd)[, names(b)])
  c(fit = sum(x * b), se = sqrt(as.numeric(t(x) %*% as.matrix(vcov(m)) %*% x))) }
pr <- list()
for (s in names(samp)) {
  d <- samp[[s]]; dl <- filter(d, leaf > 0) %>% mutate(log_leaf = log(leaf))
  ms <- lmerTest::lmer(as.formula(paste("log(new_stem) ~", red)), data = d); mr <- lmerTest::lmer(as.formula(paste("log(root_est) ~", red)), data = d)
  mlf <- lmerTest::lmer(as.formula(paste("log_leaf ~", red)), data = dl)
  mz <- glmer(as.formula(paste("I(leaf == 0) ~", red)), data = d, family = binomial, control = ctl)
  for (g in c(0, 100)) for (trt in c(levels(ds$treatment), "averaged")) {
    tr <- if (trt == "averaged") levels(ds$treatment) else trt
    mm <- function(m, dd) { v <- sapply(tr, function(z) marg(m, dd, g, z)); c(fit = mean(v["fit", ]), se = mean(v["se", ])) }   # se: mean of the two (approximate for the average)
    st <- mm(ms, d); rt <- mm(mr, d); lf <- mm(mlf, dl); pz <- mean(sapply(tr, function(z) plogis(marg(mz, d, g, z)["fit"])))
    pr[[length(pr) + 1]] <- tibble(sample = s, n = nrow(d), n_leafy = nrow(dl), ancestry = paste0(g, "% Fremont"), watering = trt,
      new_stem_g = exp(st["fit"]), new_stem_lo = exp(st["fit"] - 1.96 * st["se"]), new_stem_hi = exp(st["fit"] + 1.96 * st["se"]),
      root_g = exp(rt["fit"]), root_lo = exp(rt["fit"] - 1.96 * rt["se"]), root_hi = exp(rt["fit"] + 1.96 * rt["se"]),
      leaf_if_leafy_g = exp(lf["fit"]), leaf_lo = exp(lf["fit"] - 1.96 * lf["se"]), leaf_hi = exp(lf["fit"] + 1.96 * lf["se"]),
      p_leafless = pz, leaf_expected_g = (1 - pz) * exp(lf["fit"])) } }
pr <- bind_rows(pr) %>% mutate(stem_plus_leaf_plus_root_g = new_stem_g + leaf_expected_g + root_g)
print(as.data.frame(pr %>% mutate(across(where(is.double), ~ round(.x, 3)))), row.names = FALSE)
readr::write_excel_csv(pr, file.path(OUT, "predicted_mass.csv"))
for (s in names(samp)) { x <- pr %>% filter(sample == s, watering == "averaged"); a0 <- x[x$ancestry == "0% Fremont", ]; a1 <- x[x$ancestry == "100% Fremont", ]
  fmt <- function(v, lo, hi) sprintf("%.3f g (95%% CI %.3f–%.3f)", v, lo, hi)
  rec(4, "b", "Predicted dry mass, averaged over watering (geometric-mean scale)", s, x$n[1], "New stem: 0% → 100% Fremont", paste(fmt(a0$new_stem_g, a0$new_stem_lo, a0$new_stem_hi), "→", fmt(a1$new_stem_g, a1$new_stem_lo, a1$new_stem_hi)), NA, "results/22_allocation_log_ratios/predicted_mass.csv")
  rec(4, "b", "Predicted dry mass, averaged over watering (geometric-mean scale)", s, x$n_leafy[1], "Leaf, plants with leaves: 0% → 100% Fremont", paste(fmt(a0$leaf_if_leafy_g, a0$leaf_lo, a0$leaf_hi), "→", fmt(a1$leaf_if_leafy_g, a1$leaf_lo, a1$leaf_hi)), NA, "results/22_allocation_log_ratios/predicted_mass.csv")
  rec(4, "b", "Predicted dry mass, averaged over watering (geometric-mean scale)", s, x$n[1], "Leaf × P(has leaves): 0% → 100% Fremont", sprintf("%.3f g (P leafless %.2f) → %.3f g (P leafless %.2f)", a0$leaf_expected_g, a0$p_leafless, a1$leaf_expected_g, a1$p_leafless), NA, "results/22_allocation_log_ratios/predicted_mass.csv")
  rec(4, "b", "Predicted dry mass, averaged over watering (geometric-mean scale)", s, x$n[1], "Root: 0% → 100% Fremont", paste(fmt(a0$root_g, a0$root_lo, a0$root_hi), "→", fmt(a1$root_g, a1$root_lo, a1$root_hi)), NA, "results/22_allocation_log_ratios/predicted_mass.csv") }

hdr("4d. Why the dry:fresh ratio comes from 116 plants")
harv <- ds %>% filter(alive_01_10_2021 == 1, !is.na(dryshoot_new)) %>% mutate(scored = !is.na(amf) | !is.na(emf) | !is.na(dse),
  has_dry = !is.na(root_dryweight), has_fresh = !is.na(wet_root) & wet_root > 0)
rc <- harv %>% count(scored, has_dry, has_fresh) %>% arrange(scored, desc(has_dry), desc(has_fresh)); print(as.data.frame(rc), row.names = FALSE)
readr::write_excel_csv(rc, file.path(OUT, "root_mass_counts.csv"))
ud <- harv %>% filter(!scored, has_dry)
cat("\nUnscored plants with measured dry root mass:", nrow(ud), "| with fresh root mass too (ratio training set):", sum(ud$has_fresh), "| without:", sum(!ud$has_fresh), "\n")
cat("In allocation sample (b) (unscored, all model covariates):", nrow(samp[[2]]), "\n")
for (v in c("treatment", "pop", "myco_cat", "new_class")) { t <- table(ud[[v]], ud$has_fresh); cat(sprintf("\nFresh root mass recorded, by %s (FALSE / TRUE):\n", v)); print(t)
  cat("  Fisher exact P =", signif(fisher.test(t, simulate.p.value = TRUE, B = 1e5)$p.value, 3), "\n") }
fm <- glm(has_fresh ~ per_f + treatment, data = ud, family = binomial); cat("\nLogistic: fresh mass recorded ~ ancestry + drought:\n"); print(round(coef(summary(fm)), 4))
rec(4, "d", "Dry:fresh root ratio training set", "unscored harvested survivors", nrow(ud), "Plants with dry root mass / with fresh mass too",
    sprintf("%d unscored plants have measured dry root mass (%d in allocation sample b after covariates); only %d of them also have a whole-root fresh mass, which the ratio needs; %d have no fresh mass recorded",
            nrow(ud), nrow(samp[[2]]), sum(ud$has_fresh), sum(!ud$has_fresh)), NA, "results/22_allocation_log_ratios/root_mass_counts.csv")

hdr("4e. Total new biomass with one pooled dry:fresh ratio")
unsc <- ud %>% filter(has_fresh) %>% mutate(ratio = root_dryweight / wet_root); pooled <- mean(unsc$ratio)
cat("Pooled mean ratio:", round(pooled, 4), "(n =", nrow(unsc), ")\n")
ap <- al %>% mutate(root_pool = ifelse(scored, wet_root * pooled, root_est), total_pool = root_pool + leaf + new_stem, stem_root_pool = root_pool + new_stem, stem_root = root_est + new_stem)
pp <- list()
for (y in c("log(total_new)", "log(total_pool)", "log(stem_root)", "log(stem_root_pool)")) {
  f <- fit_red(y, ap); lab <- c("log(total_new)" = "Total new biomass, predicted ratio (reported)", "log(total_pool)" = "Total new biomass, pooled ratio",
                                "log(stem_root)" = "New stem + root, predicted ratio", "log(stem_root_pool)" = "New stem + root, pooled ratio")[[y]]
  for (tm in c("pf50", "I(pf50^2)", "treatmentDrought")) { r <- row_t(f$m, tm, c(pf50 = 50, `I(pf50^2)` = 2500, treatmentDrought = 1)[[tm]])
    term <- c(pf50 = "Linear ancestry (per 1%)", `I(pf50^2)` = "Quadratic ancestry", treatmentDrought = "Drought")[[tm]]
    pp[[length(pp) + 1]] <- tibble(response = lab, n = nobs(f$m), term = term, string = r$s, P = r$p, back = ifelse(tm == "treatmentDrought", sprintf("%s%%", pc(100 * (exp(r$b) - 1))), ""), lrt = lrt_s(f$lrt))
    rec(4, "e", lab, "(a) all 414 plants", nobs(f$m), term, paste0(r$s, ifelse(tm == "treatmentDrought", sprintf(" (%s%%)", pc(100 * (exp(r$b) - 1))), "")), r$p, "results/22_allocation_log_ratios/pooled_ratio.csv") } }
pp <- bind_rows(pp); print(as.data.frame(pp), row.names = FALSE, right = FALSE); readr::write_excel_csv(pp, file.path(OUT, "pooled_ratio.csv"))

write_rec("strings.csv"); cat("\nDone.\n"); sink()
