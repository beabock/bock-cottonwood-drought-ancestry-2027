# ==============================================================================
# Do the root, stem and leaf fraction models agree, and does the drought effect on
# root allocation depend on the imputed root mass?
# Data: results/06_allocation_models/allocation_data.csv (written by r_code/06_allocation_models.R; leaf = 0 g when leafless;
#       root_est = recorded dry root mass, or fresh mass x predicted dry:fresh ratio for colonization-scored plants).
# Reported fraction models: reduced ordered beta GLMMs, frac ~ pf50 + treatment + I(pf50^2) + og_drc_mm + myco_cat + soil
#       + (1 | pop) + (1 | pop:geno_id)  (06_allocation_models.R). Reported models are not changed.
# Run from the project root:  Rscript r_code/14_allocation_fraction_checks.R
# ==============================================================================
suppressMessages({ source("r_code/helpers/utilities_data_prep.R"); library(tidyverse); library(lmerTest); library(glmmTMB) })
out <- "results/14_allocation_fraction_checks"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("ALLOCATION FRACTION CHECKS run:", "|", R.version.string, "| glmmTMB", as.character(packageVersion("glmmTMB")), "\n")
EN <- "–"; dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
set.seed(1998)

ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, pf50 = per_fc / 50, geno_id = paste(pop, type, geno, sep = "_"))
al <- read_csv("results/06_allocation_models/allocation_data.csv", show_col_types = FALSE) %>%
  select(tree_idold, scored, root_est, new_stem, leaf, total_new, LMF, RMF, SMF) %>%
  inner_join(ds %>% select(tree_idold, per_f, per_fc, pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id, wet_root, root_dryweight), by = "tree_idold")
stopifnot(nrow(al) == 414)
samp <- list(`(a) estimated root mass` = al, `(b) measured root mass only` = filter(al, !scored))
cov <- "og_drc_mm + myco_cat + soil"; re <- "(1 | pop) + (1 | pop:geno_id)"
red <- paste("pf50 + treatment + I(pf50^2) +", cov, "+", re)
ctlist <- list(glmmTMBControl(), glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
fit_ob <- function(y, d) { for (ctrl in ctlist) { m <- tryCatch(glmmTMB(as.formula(paste(y, "~", red)), data = d, family = ordbeta(), control = ctrl), error = function(e) NULL)
  if (!is.null(m) && isTRUE(m$sdr$pdHess)) return(m) }; stop("ordbeta did not converge: ", y) }
psi_of <- function(m) { p <- m$fit$par; unname(p[names(p) == "psi"]) }
ob_mean <- function(eta, psi) { p0 <- 1 - plogis(eta - psi[1]); p1 <- plogis(eta - psi[2]); p1 + (1 - p0 - p1) * plogis(eta) }   # E[y], incl. exact 0 and 1
zstr <- function(m, t, scale = 1) { cs <- coef(summary(m))$cond; sprintf("β = %s, z = %s, %s", num(cs[t, 1] / scale), dec(cs[t, 3]), P(cs[t, 4])) }
tstr_m <- function(m, t, scale = 1) { cs <- coef(summary(m)); sprintf("β = %s, t = %s, df = %s, %s", num(cs[t, 1] / scale), dec(cs[t, "t value"]), dff(cs[t, "df"]), P(cs[t, "Pr(>|t|)"])) }

# ---------------------------------------------------------------- (a)
hdr("3a. Predicted fractions by watering (fixed effects; each plant's own ancestry and covariates, averaged), and their sum")
pa <- list()
for (s in names(samp)) {
  d <- samp[[s]]
  for (y in c("RMF", "SMF", "LMF")) {
    m <- fit_ob(y, d); psi <- psi_of(m); b <- fixef(m)$cond
    for (trt in levels(d$treatment)) {
      nd <- d %>% mutate(treatment = factor(trt, levels = levels(d$treatment)))
      X <- model.matrix(delete.response(terms(as.formula(paste("y ~ pf50 + treatment + I(pf50^2) +", cov)))), nd); eta <- as.vector(X %*% b[colnames(X)])
      pa[[length(pa) + 1]] <- tibble(sample = s, fraction = y, treatment = trt, expected = mean(ob_mean(eta, psi)), beta_part_only = mean(plogis(eta)),
                                     p_zero = mean(1 - plogis(eta - psi[1])), p_one = mean(plogis(eta - psi[2])), observed_mean = mean(d[[y]][d$treatment == trt]))
    } } }
pa <- bind_rows(pa)
print(as.data.frame(pa %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)
sums <- pa %>% group_by(sample, treatment) %>% summarise(sum_expected = sum(expected), sum_beta_part = sum(beta_part_only), sum_observed = sum(observed_mean), .groups = "drop")
cat("\nSum of the three predicted fractions:\n"); print(as.data.frame(sums %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)
dlt <- pa %>% select(sample, fraction, treatment, expected) %>% pivot_wider(names_from = treatment, values_from = expected) %>% mutate(drought_minus_ww = Drought - `Well-watered`)
cat("\nDrought minus well-watered, expected fraction (percentage points):\n"); print(as.data.frame(dlt %>% mutate(across(where(is.numeric), ~ round(100 * .x, 2)))), row.names = FALSE)
cat("Sum of the three drought - well-watered differences (should be 0 if the models were consistent):\n")
print(dlt %>% group_by(sample) %>% summarise(sum_of_differences_pp = round(100 * sum(drought_minus_ww), 2)))
write.csv(pa, file.path(out, "predicted_fractions.csv"), row.names = FALSE)

# ---------------------------------------------------------------- (b)
hdr("3b. Compositional analysis: additive log-ratios (reference = new stem) among plants with leaves, + leafless hurdle")
alr <- list()
for (s in names(samp)) {
  d <- samp[[s]] %>% filter(leaf > 0, root_est > 0, new_stem > 0) %>% mutate(alr_root = log(root_est / new_stem), alr_leaf = log(leaf / new_stem))
  for (y in c("alr_root", "alr_leaf")) {
    m <- lmerTest::lmer(as.formula(paste(y, "~", red)), data = d)
    full <- lmerTest::lmer(as.formula(paste(y, "~ pf50 * treatment + I(pf50^2) * treatment +", cov, "+", re)), data = d)
    lrt <- anova(m, full, refit = TRUE)
    cs <- coef(summary(m))
    for (tm in c("pf50", "treatmentDrought")) {
      sc <- if (tm == "pf50") 50 else 1
      alr[[length(alr) + 1]] <- tibble(sample = s, n = nobs(m), response = y, term = ifelse(tm == "pf50", "Linear ancestry (per 1%)", "Drought"),
        estimate = cs[tm, 1] / sc, t_value = cs[tm, "t value"], df = cs[tm, "df"], P = cs[tm, "Pr(>|t|)"], string = tstr_m(m, tm, sc),
        ratio_change = ifelse(tm == "treatmentDrought", sprintf("drought × %.3f", exp(cs[tm, 1])), sprintf("%s%% per 10 points", dec(100 * (exp(10 * cs[tm, 1] / sc) - 1), 1))),
        lrt_ancestry_x_drought = sprintf("χ²(%d) = %s, %s", lrt$Df[2], dec(lrt$Chisq[2]), P(lrt$`Pr(>Chisq)`[2])))
    } } }
alr <- bind_rows(alr)
print(as.data.frame(alr %>% select(sample, n, response, term, string, ratio_change, lrt_ancestry_x_drought)), right = FALSE, row.names = FALSE)
s6 <- read_csv("results/07_allocation_appendix_S7/appendix_S7_source.csv", show_col_types = FALSE) %>% filter(response == "P(no leaves at harvest)", term == "Drought")
cat("\nHurdle part (leafless at harvest, binomial GLMM, Appendix S7 source): (a) n = 425:", s6$string_estimated_root, "| (b) n = 253:", s6$string_measured_root_only, "\n")
write.csv(alr, file.path(out, "alr_models.csv"), row.names = FALSE)

# ---------------------------------------------------------------- (c)
hdr("3c. Raw medians by watering treatment")
med <- map_dfr(names(samp), function(s) samp[[s]] %>% group_by(treatment) %>%
  summarise(n = n(), leafless = sum(leaf == 0), RMF = median(RMF), SMF = median(SMF), LMF = median(LMF), root_g = median(root_est), new_stem_g = median(new_stem),
            leaf_g = median(leaf), leaf_g_leafy_only = median(leaf[leaf > 0]), .groups = "drop") %>% mutate(sample = s, .before = 1))
print(as.data.frame(med %>% mutate(across(where(is.double), ~ round(.x, 4)))), row.names = FALSE)
write.csv(med, file.path(out, "raw_medians.csv"), row.names = FALSE)

# ---------------------------------------------------------------- (e)
hdr("3e. Dry:fresh root ratio model and a pooled-ratio refit")
harv <- ds %>% filter(alive_01_10_2021 == 1, !is.na(dryshoot_new)) %>% mutate(scored = !is.na(amf) | !is.na(emf) | !is.na(dse))
unsc <- harv %>% filter(!scored, !is.na(root_dryweight), !is.na(wet_root), wet_root > 0) %>% mutate(ratio = root_dryweight / wet_root)
m_anc <- lm(ratio ~ (per_fc + I(per_fc^2)) * treatment, data = unsc)       # 06_allocation_models.R (the model used)
m_nod <- lm(ratio ~ per_fc + I(per_fc^2), data = unsc)
cat("Training set: n =", nrow(unsc), "unscored plants | pooled mean ratio =", round(mean(unsc$ratio), 4), "\n")
print(round(coef(summary(m_anc)), 5))
a <- anova(m_nod, m_anc); cat(sprintf("Drought terms jointly (drought + drought x linear + drought x quadratic ancestry): F = %.2f, df = %d, %d, P = %.3f\n", a$F[2], a$Df[2], a$Res.Df[2], a$`Pr(>F)`[2]))
cat(sprintf("Raw mean ratio: well-watered %.4f (n = %d), drought %.4f (n = %d); Welch t P = %.3f\n", mean(unsc$ratio[unsc$treatment == "Well-watered"]), sum(unsc$treatment == "Well-watered"),
            mean(unsc$ratio[unsc$treatment == "Drought"]), sum(unsc$treatment == "Drought"), t.test(ratio ~ treatment, data = unsc)$p.value))
pooled <- mean(unsc$ratio)
ap <- al %>% mutate(root_pool = ifelse(scored, wet_root * pooled, root_est), total_pool = root_pool + leaf + new_stem,
                    RMF_pool = root_pool / total_pool, log_rs_pool = log(root_pool / (leaf + new_stem)), log_rs = log(root_est / (leaf + new_stem)))
cat("Scored plants whose root mass is re-estimated with the pooled ratio:", sum(ap$scored), "| median ratio of pooled to predicted estimate:",
    round(median(ap$root_pool[ap$scored] / ap$root_est[ap$scored]), 3), "\n")
cmp <- list()
for (v in list(c("RMF", "RMF_pool", "Root mass fraction (ordered beta)"), c("log_rs", "log_rs_pool", "Root:shoot ratio (log, LMM)"))) {
  for (k in 1:2) {
    y <- v[k]
    m <- if (grepl("RMF", y)) fit_ob(y, ap) else lmerTest::lmer(as.formula(paste(y, "~", red)), data = ap)
    st <- function(t, sc) if (grepl("RMF", y)) zstr(m, t, sc) else tstr_m(m, t, sc)
    pv <- function(t) if (grepl("RMF", y)) coef(summary(m))$cond[t, 4] else coef(summary(m))[t, "Pr(>|t|)"]
    cmp[[length(cmp) + 1]] <- tibble(response = v[3], root_mass = ifelse(k == 1, "predicted ratio (reported)", "pooled ratio"), n = nobs(m),
                                     drought = st("treatmentDrought", 1), drought_P = pv("treatmentDrought"), linear_ancestry = st("pf50", 50), ancestry_P = pv("pf50"))
  } }
cmp <- bind_rows(cmp); print(as.data.frame(cmp %>% select(-drought_P, -ancestry_P)), right = FALSE, row.names = FALSE)
write.csv(cmp, file.path(out, "pooled_ratio.csv"), row.names = FALSE)
sink()
