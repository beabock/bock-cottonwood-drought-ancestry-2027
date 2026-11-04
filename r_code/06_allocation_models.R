# ==============================================================================
# ALLOCATION MODELS. Two data conventions:
#   1. Missing leaf mass on a harvested survivor = no leaves at harvest -> 0 g.
#   2. Plants scored for colonization had a root subsample removed after fresh weighing
#      and not included in dried root mass -> dried root mass is underestimated. Their root
#      dry mass is estimated as whole-root fresh mass x the dry:fresh ratio predicted from
#      UNSCORED plants.
# Does not change any manuscript result, figure or results string. Supersedes
# an earlier exploratory version (not in this repository).
# ==============================================================================
# Run from the project root (R 4.6.1):  Rscript r_code/06_allocation_models.R
# Writes to results/06_allocation_models/: output.txt, allocation_strings.md/.csv,
#   allocation_root_ratio.csv, allocation_data.csv, allocation_leaf_zero_counts.csv,
#   allocation_diagnostics.pdf
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(lme4); library(lmerTest); library(glmmTMB); library(DHARMa); library(tidyverse)
out <- "results/06_allocation_models"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("ALLOCATION (exploratory) run:", "|", R.version.string,
    "| lme4", as.character(packageVersion("lme4")), "| glmmTMB", as.character(packageVersion("glmmTMB")), "\n")
set.seed(1998)
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, pf50 = per_fc / 50, geno_id = paste(pop, type, geno, sep = "_"))
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# --- AJB strings; ancestry fitted per 50 points (pf50) for stability, reported per 1% ---
EN <- "–"; BETA <- "β"; TIMES <- "×"; CHI <- "χ²"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
coefs <- function(m) {        # coefficient table in per-1% ancestry units, uniform column names
  cs <- coef(summary(m)); if (is.list(cs) && "cond" %in% names(cs)) cs <- cs$cond
  cs <- as.data.frame(cs); rn <- rownames(cs); div <- ifelse(grepl("pf50\\^2", rn), 2500, ifelse(grepl("pf50", rn), 50, 1))
  out <- data.frame(term = gsub("pf50", "per_fc", rn), est = cs[, "Estimate"] / div,
                    stat = cs[, grep("value", colnames(cs))[1]], df = if ("df" %in% colnames(cs)) cs[, "df"] else NA,
                    p = cs[, grep("^Pr", colnames(cs))[1]], stat_name = ifelse("df" %in% colnames(cs), "t", "z"))
  rownames(out) <- out$term; out }
sstr <- function(cf, t) { r <- cf[t, ]
  if (r$stat_name == "t") sprintf("%s = %s, t = %s, df = %s, %s", BETA, num(r$est), dec(r$stat), dff(r$df), P(r$p))
  else sprintf("%s = %s, z = %s, %s", BETA, num(r$est), dec(r$stat), P(r$p)) }
lab <- c(per_fc = "Linear ancestry", "I(per_fc^2)" = "Quadratic ancestry", treatmentDrought = "Drought",
         "per_fc:treatmentDrought" = paste("Ancestry", TIMES, "drought"), "treatmentDrought:I(per_fc^2)" = paste("Quadratic ancestry", TIMES, "drought"))

# ------------------------------------------------------------------------------
hdr("1. Leaf recode and root dry mass estimation")
surv <- ds %>% filter(alive_01_10_2021 == 1)
harv <- surv %>% filter(!is.na(dryshoot_new))     # harvested survivors (new stem weighed)
cat("Survivors:", nrow(surv), "| harvested survivors with new stem mass:", nrow(harv), "\n")
cat("Leaf mass missing among harvested survivors, recoded to 0 g:", sum(is.na(harv$dryleaf_weight)),
    "| already 0 g:", sum(harv$dryleaf_weight == 0, na.rm = TRUE), "\n")
harv <- harv %>% mutate(leaf_recoded = is.na(dryleaf_weight), leaf = coalesce(dryleaf_weight, 0),
                        scored = !is.na(amf) | !is.na(emf) | !is.na(dse))
unsc <- harv %>% filter(!scored, !is.na(root_dryweight), !is.na(wet_root), wet_root > 0) %>% mutate(ratio = root_dryweight / wet_root)
cat("Unscored plants with dry and fresh root mass (ratio training set):", nrow(unsc), "\n")
m_cls <- lm(ratio ~ new_class * treatment, data = unsc)
m_anc <- lm(ratio ~ (per_fc + I(per_fc^2)) * treatment, data = unsc)
aics <- c(cross_type_x_drought = AIC(m_cls), ancestry_x_drought = AIC(m_anc))
m_0 <- lm(ratio ~ 1, data = unsc)
cat("Reference: intercept-only (single mean ratio", round(coef(m_0), 4), ") AIC =", round(AIC(m_0), 1),
    "| F test ancestry x drought vs intercept-only P =", signif(anova(m_0, m_anc)$`Pr(>F)`[2], 3),
    "| cross type x drought vs intercept-only P =", signif(anova(m_0, m_cls)$`Pr(>F)`[2], 3), "\n")
cat("Dry:fresh ratio models on the ratio scale (unscored plants). AIC:", paste(names(aics), round(aics, 1), sep = " = ", collapse = "; "), "\n")
best <- if (aics["ancestry_x_drought"] <= aics["cross_type_x_drought"]) m_anc else m_cls
best_name <- names(which.min(aics)); cat("Chosen (lower AIC):", best_name, "| adj. R2 =", round(summary(best)$adj.r.squared, 3), "\n")
print(round(coef(summary(best)), 5))
grid <- expand_grid(new_class = factor(levels(ds$new_class), levels = levels(ds$new_class)), treatment = factor(levels(ds$treatment), levels = levels(ds$treatment))) %>%
  left_join(ds %>% group_by(new_class) %>% summarise(per_fc = mean(per_fc), .groups = "drop"), by = "new_class")
grid$predicted_ratio <- predict(best, newdata = grid)
grid <- grid %>% left_join(unsc %>% count(new_class, treatment, name = "n_unscored"), by = c("new_class", "treatment"))
cat("\nPredicted dry:fresh ratio by cross type x drought (at each cross type's mean ancestry):\n"); print(as.data.frame(grid %>% mutate(across(c(per_fc, predicted_ratio), ~ round(., 4)))), row.names = FALSE)
write.csv(grid, file.path(out, "allocation_root_ratio.csv"), row.names = FALSE)
harv <- harv %>% mutate(ratio_pred = predict(best, newdata = .),
                        root_est = ifelse(scored, wet_root * ratio_pred, root_dryweight),
                        root_estimated = scored & !is.na(root_est))
sc <- harv %>% filter(scored)
cat("\nScored harvested plants:", nrow(sc), "| with fresh root mass (estimable):", sum(!is.na(sc$wet_root)),
    "| without fresh root mass (root mass missing):", sum(is.na(sc$wet_root)), "\n")
cmp <- sc %>% filter(!is.na(root_est), !is.na(root_dryweight))
cat("Estimated vs recorded dry root mass in scored plants: median ratio", round(median(cmp$root_est / cmp$root_dryweight), 3),
    "(implies ~", round(100 * (1 - median(cmp$root_dryweight / cmp$root_est)), 1), "% of root mass removed as subsample);",
    "estimate < recorded in", sum(cmp$root_est < cmp$root_dryweight), "of", nrow(cmp), "plants\n")

al <- harv %>% mutate(total_new = root_est + leaf + dryshoot_new, log_total_new = log(total_new),
                      LMF = leaf / total_new, RMF = root_est / total_new, SMF = dryshoot_new / total_new,
                      log_root_shoot = log(root_est / (leaf + dryshoot_new)), log_new_stem = log(dryshoot_new)) %>%
  filter(!is.na(root_est)) %>% filter(if_all(c(pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), ~ !is.na(.x)))
write.csv(al %>% select(tree_idold, new_class, per_f, treatment, scored, root_estimated, leaf_recoded, leaf, root_recorded = root_dryweight, root_est,
                        new_stem = dryshoot_new, total_new, LMF, RMF, SMF), file.path(out, "allocation_data.csv"), row.names = FALSE)
samp <- list(`(a) estimated root mass for scored plants` = al, `(b) scored plants excluded` = filter(al, !scored))
for (s in names(samp)) cat(sprintf("Sample %s: n = %d (zero-leaf plants %d; estimated root mass %d)\n", s, nrow(samp[[s]]), sum(samp[[s]]$leaf == 0), sum(samp[[s]]$root_estimated)))

# ------------------------------------------------------------------------------
cov <- "og_drc_mm + myco_cat + soil"; re <- "(1 | pop) + (1 | pop:geno_id)"
full_rhs <- paste("pf50 * treatment + I(pf50^2) * treatment +", cov, "+", re); red_rhs <- paste("pf50 + treatment + I(pf50^2) +", cov, "+", re)
diag_lmm <- function(m) { r <- resid(m, scaled = TRUE); sprintf("max |scaled resid| = %.2f; %d beyond |3|; Shapiro-Wilk P = %s; singular = %s",
                                                               max(abs(r)), sum(abs(r) > 3), formatC(shapiro.test(r)$p.value, digits = 2, format = "g"), isSingular(m)) }
diag_glmm <- function(m) { sr <- simulateResiduals(m, n = 500, plot = FALSE)
  sprintf("DHARMa uniformity P = %s; dispersion P = %s; outliers P = %s", formatC(testUniformity(sr, plot = FALSE)$p.value, digits = 2, format = "g"),
          formatC(testDispersion(sr, plot = FALSE)$p.value, digits = 2, format = "g"), formatC(testOutliers(sr, plot = FALSE, type = "bootstrap")$p.value, digits = 2, format = "g")) }
tmb_ok <- function(m) !is.null(m) && isTRUE(m$sdr$pdHess) && !anyNA(coef(summary(m))$cond)   # converged with usable SEs
FALLBACK_HOW <- "LMM on logit((y(n-1)+0.5)/n) (lmerTest, Satterthwaite df)"
fit_one <- function(y, d, family, force_fallback = FALSE) {
  if (family == "lmm") {
    mf <- lmerTest::lmer(as.formula(paste(y, "~", full_rhs)), data = d); mr <- lmerTest::lmer(as.formula(paste(y, "~", red_rhs)), data = d)
    list(full = mf, red = mr, lrt = anova(mr, mf, refit = TRUE), how = "LMM (lmerTest, Satterthwaite df)")
  } else {
    if (!force_fallback) {
      fit <- function(rhs, ctrl) glmmTMB(as.formula(paste(y, "~", rhs)), data = d, family = ordbeta(), control = ctrl)
      for (ctrl in list(glmmTMBControl(), glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))) {   # default, then BFGS
        mf <- tryCatch(fit(full_rhs, ctrl), error = function(e) NULL); mr <- tryCatch(fit(red_rhs, ctrl), error = function(e) NULL)
        if (tmb_ok(mf) && tmb_ok(mr)) return(list(full = mf, red = mr, lrt = anova(mr, mf), how = "ordered beta GLMM (glmmTMB ordbeta, logit link; Wald z)"))
      }
    }
    n <- nrow(d); d[[paste0(y, "_lg")]] <- qlogis((d[[y]] * (n - 1) + 0.5) / n)
    mf <- lmerTest::lmer(as.formula(paste0(y, "_lg ~ ", full_rhs)), data = d); mr <- lmerTest::lmer(as.formula(paste0(y, "_lg ~ ", red_rhs)), data = d)
    list(full = mf, red = mr, lrt = anova(mr, mf, refit = TRUE), how = FALLBACK_HOW)
  }
}
resp <- list(log_total_new = list(lab = "log total new biomass", fam = "lmm"), LMF = list(lab = "leaf mass fraction", fam = "frac"),
             RMF = list(lab = "root mass fraction", fam = "frac"), SMF = list(lab = "stem mass fraction", fam = "frac"),
             log_root_shoot = list(lab = "log root:shoot", fam = "lmm"))
rows <- list(); diag_store <- list()
hdr("2. Ancestry x drought models: (a) estimated root mass for scored plants; (b) scored plants excluded")
for (y in names(resp)) {
 fits <- lapply(samp, function(d) fit_one(y, d, resp[[y]]$fam))
 if (length(unique(sapply(fits, `[[`, "how"))) > 1) {        # keep (a) and (b) comparable: same method in both samples
   cat(sprintf("\nNOTE: %s - ordered beta did not converge in every sample; both samples refit with the logit LMM so (a) and (b) are comparable.\n", resp[[y]]$lab))
   fits <- lapply(samp, function(d) fit_one(y, d, resp[[y]]$fam, force_fallback = TRUE))
 }
 for (s in names(samp)) {
  f <- fits[[s]]
  chi <- f$lrt$Chisq[2]; dfc <- if ("Chi Df" %in% names(f$lrt)) f$lrt$`Chi Df`[2] else f$lrt$Df[2]; pl <- f$lrt$`Pr(>Chisq)`[2]
  use <- if (pl < 0.05) f$full else f$red; cf <- coefs(use)
  dg <- if (inherits(use, "glmmTMB")) diag_glmm(use) else diag_lmm(use); diag_store[[paste(resp[[y]]$lab, s)]] <- use
  terms <- if (pl < 0.05) names(lab) else c("per_fc", "I(per_fc^2)", "treatmentDrought")
  cat(sprintf("\n--- %s | %s | n = %d | %s | LRT ancestry x drought %s(%d) = %.2f, %s -> %s ---\n", resp[[y]]$lab, s, nobs(use), f$how, CHI, dfc, chi, P(pl), ifelse(pl < 0.05, "FULL", "REDUCED")))
  for (t in terms) cat(sprintf("  %-28s %s\n", lab[t], sstr(cf, t)))
  cat("  diagnostics:", dg, "\n")
  rows[[length(rows) + 1]] <- tibble(section = "2", response = resp[[y]]$lab, sample = s, n = nobs(use), method = f$how, model = ifelse(pl < 0.05, "full", "reduced"),
                                     term = c(paste("LRT ancestry", TIMES, "drought"), lab[terms]),
                                     string = c(sprintf("%s(%d) = %s, %s", CHI, dfc, dec(chi), P(pl)), sapply(terms, function(t) sstr(cf, t))),
                                     est = c(NA, cf[terms, "est"]), p = c(pl, cf[terms, "p"]), diagnostics = dg)
 }
}

# ------------------------------------------------------------------------------
hdr("3. Same-sample check: reduced new stem model on the step-2 samples")
for (s in names(samp)) {
  m <- lmerTest::lmer(as.formula(paste("log_new_stem ~", red_rhs)), data = samp[[s]]); cf <- coefs(m)
  cat(sprintf("\n%s (n = %d): linear %s | quadratic %s | drought %s\n  diagnostics: %s\n", s, nobs(m), sstr(cf, "per_fc"), sstr(cf, "I(per_fc^2)"), sstr(cf, "treatmentDrought"), diag_lmm(m)))
  rows[[length(rows) + 1]] <- tibble(section = "3", response = "log new stem mass (reduced)", sample = s, n = nobs(m), method = "LMM (lmerTest, Satterthwaite df)", model = "reduced",
                                     term = lab[c("per_fc", "I(per_fc^2)", "treatmentDrought")], string = sapply(c("per_fc", "I(per_fc^2)", "treatmentDrought"), function(t) sstr(cf, t)),
                                     est = cf[c("per_fc", "I(per_fc^2)", "treatmentDrought"), "est"], p = cf[c("per_fc", "I(per_fc^2)", "treatmentDrought"), "p"], diagnostics = diag_lmm(m))
}
cat("\nReported (manuscript) reduced new stem model, all 429 plants: linear beta = -0.00528, quadratic beta = -0.000145, drought beta = -0.112 (formatted_results_main_models_df.md)\n")

# ------------------------------------------------------------------------------
hdr("4. Leaf drop: probability of zero leaf mass ~ ancestry x drought (binomial GLMM)")
lz <- harv %>% mutate(zero_leaf = as.integer(leaf == 0)) %>% filter(if_all(c(pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), ~ !is.na(.x)))
cnt <- lz %>% group_by(cross_type = new_class, treatment) %>% summarise(n = n(), zero_leaf = sum(zero_leaf), .groups = "drop") %>% mutate(pct = round(100 * zero_leaf / n, 1)) %>% relocate(n, .after = zero_leaf)
cat("Harvested survivors with covariates:", nrow(lz), "| zero leaf mass:", sum(lz$zero_leaf), "\n\nCounts per cross type x drought:\n"); print(as.data.frame(cnt), row.names = FALSE)
write.csv(cnt, file.path(out, "allocation_leaf_zero_counts.csv"), row.names = FALSE)
zf <- glmer(as.formula(paste("zero_leaf ~", full_rhs)), data = lz, family = binomial, control = ctl)
zr <- glmer(as.formula(paste("zero_leaf ~", red_rhs)), data = lz, family = binomial, control = ctl)
zl <- anova(zr, zf); pz <- zl$`Pr(>Chisq)`[2]; zu <- if (pz < 0.05) zf else zr; cz <- coefs(zu)
cat(sprintf("\nLRT ancestry x drought %s(%d) = %.2f, %s -> %s model; n = %d; singular = %s; genotype SD = %.2f\n", CHI, zl$Df[2], zl$Chisq[2], P(pz),
            ifelse(pz < 0.05, "FULL", "REDUCED"), nobs(zu), isSingular(zu), attr(VarCorr(zu)$`pop:geno_id`, "stddev")))
tz <- if (pz < 0.05) names(lab) else c("per_fc", "I(per_fc^2)", "treatmentDrought")
for (t in tz) cat(sprintf("  %-28s %s (odds ratio %s per unit)\n", lab[t], sstr(cz, t), num(exp(cz[t, "est"]))))
dz <- diag_glmm(zu); cat("  diagnostics:", dz, "\n")
cat("  Odds ratio for drought:", num(exp(cz["treatmentDrought", "est"])), "| odds ratio per 10 points of Fremont ancestry (linear term, at 50%):", num(exp(10 * cz["per_fc", "est"])), "\n")
rows[[length(rows) + 1]] <- tibble(section = "4", response = "P(zero leaf mass)", sample = "harvested survivors", n = nobs(zu), method = "binomial GLMM (lme4 glmer; Wald z)",
                                   model = ifelse(pz < 0.05, "full", "reduced"), term = c(paste("LRT ancestry", TIMES, "drought"), lab[tz]),
                                   string = c(sprintf("%s(%d) = %s, %s", CHI, zl$Df[2], dec(zl$Chisq[2]), P(pz)), sapply(tz, function(t) sstr(cz, t))),
                                   est = c(NA, cz[tz, "est"]), p = c(pz, cz[tz, "p"]), diagnostics = dz)

# ------------------------------------------------------------------------------
hdr("5. Colonized subset (APPROXIMATE: every plant's root mass is estimated)")
col <- al %>% filter(!is.na(amf), !is.na(emf), !is.na(dse))
rep77 <- ds %>% filter(!is.na(amf), !is.na(emf), !is.na(dse), !is.na(dryshoot_new), !is.na(og_drc_mm), !is.na(myco_cat), !is.na(soil))
cat("Colonized plants here:", nrow(col), "| reported colonization-biomass sample: 77 | missing here (no fresh root mass):",
    paste(setdiff(rep77$tree_idold, col$tree_idold), collapse = ", "), "\n")
crhs <- paste("scale(amf) * treatment + scale(emf) * treatment + scale(dse) * treatment + pf50 + I(pf50^2) +", cov, "+", re)
tb <- c(AMF = "scale(amf):treatmentDrought", EMF = "treatmentDrought:scale(emf)", DSE = "treatmentDrought:scale(dse)")
for (y in c("log_total_new", "RMF", "LMF")) {
  isfrac <- y != "log_total_new"
  if (isfrac) {
    for (ctrl in list(glmmTMBControl(), glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))) {   # default, then BFGS
      m <- tryCatch(glmmTMB(as.formula(paste(y, "~", crhs)), data = col, family = ordbeta(), control = ctrl), error = function(e) NULL)
      if (tmb_ok(m)) break
    }
  } else m <- lmerTest::lmer(as.formula(paste(y, "~", crhs)), data = col)
  how <- if (isfrac) "ordered beta GLMM (glmmTMB ordbeta; Wald z)" else "LMM (lmerTest, Satterthwaite df)"
  if (isfrac && !tmb_ok(m)) { n <- nrow(col); col[[paste0(y, "_lg")]] <- qlogis((col[[y]] * (n - 1) + 0.5) / n)
    m <- lmerTest::lmer(as.formula(paste0(y, "_lg ~ ", crhs)), data = col); how <- "FALLBACK: LMM on squeezed logit" }
  cf <- coefs(m); dg <- if (inherits(m, "glmmTMB")) diag_glmm(m) else diag_lmm(m); diag_store[[paste(y, "colonized")]] <- m
  cat(sprintf("\n--- %s (APPROXIMATE) | n = %d | %s ---\n", c(log_total_new = "log total new biomass", RMF = "root mass fraction", LMF = "leaf mass fraction")[y], nobs(m), how))
  for (k in names(tb)) cat(sprintf("  %s %s drought: %s\n", k, TIMES, sstr(cf, tb[k])))
  cat("  diagnostics:", dg, "\n")
  rows[[length(rows) + 1]] <- tibble(section = "5 (approximate)", response = c(log_total_new = "log total new biomass", RMF = "root mass fraction", LMF = "leaf mass fraction")[y],
                                     sample = "colonized, root mass estimated", n = nobs(m), method = how, model = "colonization x drought",
                                     term = paste(names(tb), TIMES, "drought"), string = sapply(tb, function(t) sstr(cf, t)), est = cf[tb, "est"], p = cf[tb, "p"], diagnostics = dg)
}

# ------------------------------------------------------------------------------
hdr("Agreement between samples (a) and (b), step 2")
st <- bind_rows(rows)
agree <- st %>% filter(section == "2") %>% select(response, term, sample, est, p) %>%
  mutate(sample = ifelse(grepl("^\\(a\\)", sample), "a", "b")) %>%
  pivot_wider(names_from = sample, values_from = c(est, p)) %>%
  mutate(agree = ifelse(is.na(p_a) | is.na(p_b), "term only in one model", ifelse((p_a < 0.05) == (p_b < 0.05) & (is.na(est_a) | sign(est_a) == sign(est_b)), "yes", "NO")))
print(as.data.frame(agree %>% mutate(across(c(est_a, est_b, p_a, p_b), ~ signif(., 3)))), right = FALSE)
write.csv(st, file.path(out, "allocation_strings.csv"), row.names = FALSE)
md <- c("# Allocation (exploratory): AJB-style strings", "",
        "Generated by `r_code/06_allocation_models.R` (R 4.6.1). Not part of the manuscript results. Missing leaf mass on harvested survivors is recoded to 0 g; root dry mass of plants scored for colonization is estimated from fresh mass × the dry:fresh ratio predicted from unscored plants. Ancestry coefficients are per 1% Fremont ancestry. LMM rows give t with Satterthwaite df; GLMM rows give Wald z.", "",
        "| Section | Response | Sample | n | Method | Model | Term | String | Diagnostics |", "|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(st))) md <- c(md, sprintf("| %s | %s | %s | %d | %s | %s | %s | %s | %s |", st$section[i], st$response[i], st$sample[i], st$n[i], st$method[i], st$model[i], st$term[i], st$string[i], gsub("\\|", "/", st$diagnostics[i])))
con <- file(file.path(out, "allocation_strings.md"), "w", encoding = "UTF-8"); writeLines(enc2utf8(md), con); close(con)

pdf(file.path(out, "allocation_diagnostics.pdf"), width = 8, height = 4)
for (nm in names(diag_store)) { m <- diag_store[[nm]]
  if (inherits(m, "glmmTMB")) plot(simulateResiduals(m, n = 500), title = nm) else { par(mfrow = c(1, 2)); r <- resid(m, scaled = TRUE)
    qqnorm(r, main = nm, cex.main = 0.7); qqline(r); plot(fitted(m), r, xlab = "Fitted", ylab = "Scaled residual", main = nm, cex.main = 0.7); abline(h = 0, lty = 2) } }
invisible(dev.off())

cat("\nDone.\n")
sink()
