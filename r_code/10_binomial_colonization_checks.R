# ==============================================================================
# BINOMIAL COLONIZATION AND SURVIVAL CHECKS
#   A. Colonization LMM refit as a binomial GLMM on raw grid-line intersect counts
#      (data/raw_data/gh_amscore_raw.csv, gh_emscore_raw.csv). Does the EMF result hold?
#   B. Survival GLMM with population as a random intercept (as in the other models)
#      instead of a fixed effect.
# Changes no reported result.
# ==============================================================================
# Run from the project root (R 4.6.1):  Rscript r_code/10_binomial_colonization_checks.R
# Writes: results/10_binomial_colonization_checks/output.txt, colonization_counts_glmm.csv, survival_population.csv
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(lme4); library(lmerTest); library(glmmTMB); library(emmeans); library(tidyverse)
out <- "results/10_binomial_colonization_checks"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("BINOMIAL COLONIZATION CHECKS run:", "|", R.version.string, "| lme4", as.character(packageVersion("lme4")),
    "| glmmTMB", as.character(packageVersion("glmmTMB")), "\n")
EN <- "–"; BETA <- "β"; TIMES <- "×"; CHI <- "χ²"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
zstr <- function(b, z, p) sprintf("%s = %s, z = %s, %s", BETA, num(b), dec(z), P(p))
tstr <- function(b, t, df, p) sprintf("%s = %s, t = %s, df = %s, %s", BETA, num(b), dec(t), dff(df), P(p))
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, pf50 = per_fc / 50, geno_id = paste(pop, type, geno, sep = "_"))

# ------------------------------------------------------------------------------
hdr("A. Colonization as binomial counts (grid-line intersects)")
# ------------------------------------------------------------------------------
# raw files have one duplicated plant ID, coalesced to its first non-missing row as in 01_merge_datasets.R
first_nonmissing <- function(d) d %>% group_by(tree_idold) %>% summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop")
am <- read.csv("data/raw_data/gh_amscore_raw.csv", fileEncoding = "UTF-8-BOM") %>% first_nonmissing() %>%
  transmute(tree_idold, amf_k = arb + ves + hyph, dse_k = dse, am_n = arb + ves + hyph + dse + root)
em <- read.csv("data/raw_data/gh_emscore_raw.csv", fileEncoding = "UTF-8-BOM") %>% first_nonmissing() %>%
  transmute(tree_idold, emf_k = brown_tip + dark_finger + dark_thumb + white_tip + brown_bunch + dead_myco,
            em_n = root_alive + root_dead + dark_finger + brown_bunch + dark_thumb + white_tip + brown_tip + dead_myco)
base <- ds %>% filter(!is.na(amf), !is.na(emf), !is.na(dse)) %>%
  select(tree_idold, tree_id, pop, geno_id, per_f, pf50, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>% na.omit() %>%
  mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm))) %>% left_join(am, by = "tree_idold") %>% left_join(em, by = "tree_idold")
stopifnot(nrow(base) == 78, !anyNA(base$am_n), !anyNA(base$em_n))
dev <- with(base, c(amf = max(abs(100 * amf_k / am_n - amf)), dse = max(abs(100 * dse_k / am_n - dse)), emf = max(abs(100 * emf_k / em_n - emf))))
cat("Plants: 78 | intersects per plant: AMF/DSE slides median", median(base$am_n), "(", paste(range(base$am_n), collapse = "-"), "); EMF median",
    median(base$em_n), "(", paste(range(base$em_n), collapse = "-"), ")\n")
cat("Max |100 x count / total - analysed percentage|:", paste(names(dev), signif(dev, 3), sep = " = ", collapse = "; "), "(0 = counts reproduce the percentages)\n")
cat("Note: AMF and DSE are scored as mutually exclusive categories at each intersect on the same slides, so their counts share one denominator.\n")
lc <- bind_rows(base %>% transmute(tree_id, pop, geno_id, per_f, pf50, treatment, og_drc_mm_z, soil, myco_cat, fungal_type = "amf", k = amf_k, n = am_n),
                base %>% transmute(tree_id, pop, geno_id, per_f, pf50, treatment, og_drc_mm_z, soil, myco_cat, fungal_type = "dse", k = dse_k, n = am_n),
                base %>% transmute(tree_id, pop, geno_id, per_f, pf50, treatment, og_drc_mm_z, soil, myco_cat, fungal_type = "emf", k = emf_k, n = em_n)) %>%
  mutate(fungal_type = factor(fungal_type, levels = c("amf", "dse", "emf")), myco_cat = factor(myco_cat), obs = factor(seq_len(n())))
f <- cbind(k, n - k) ~ pf50 * fungal_type * treatment + og_drc_mm_z * fungal_type + soil
g_bin <- glmer(update(f, . ~ . + (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id) + (1 | obs)), data = lc, family = binomial, control = ctl)
g_bb <- glmmTMB(update(f, . ~ . + (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id)), data = lc, family = betabinomial())
g_plain <- glmer(update(f, . ~ . + (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id)), data = lc, family = binomial, control = ctl)
od <- sum(residuals(g_plain, type = "pearson")^2) / df.residual(g_plain)
cat(sprintf("\nOverdispersion without an observation-level random effect: Pearson chi2/df = %.1f, so the reported binomial GLMM includes one (1 | obs).\n", od))
cat("Binomial GLMM (logit link) with observation-level random effect: singular =", isSingular(g_bin), "| convergence messages:",
    paste(g_bin@optinfo$conv$lme4$messages, collapse = "; "), "\n"); print(VarCorr(g_bin))
cat("Beta-binomial GLMM (glmmTMB): converged =", isTRUE(g_bb$sdr$pdHess), "\n")
rows <- list()
for (nm in c("Binomial GLMM + observation-level RE (lme4)", "Beta-binomial GLMM (glmmTMB)")) {
  m <- if (grepl("^Binomial", nm)) g_bin else g_bb
  cf <- if (inherits(m, "glmmTMB")) coef(summary(m))$cond else coef(summary(m))
  it <- cf["pf50:fungal_typeemf", ]
  et <- as.data.frame(summary(emtrends(m, ~ fungal_type, var = "pf50"), infer = c(TRUE, TRUE)))
  et <- et %>% mutate(slope = pf50.trend / 50, lo = asymp.LCL / 50, hi = asymp.UCL / 50, z = z.ratio, p = p.value,
                      or10 = exp(10 * slope), or10_lo = exp(10 * lo), or10_hi = exp(10 * hi))
  cat(sprintf("\n%s, n = %d plants (%d rows)\n  EMF vs. AMF ancestry slope, well-watered (interaction, per 1%% Fremont, logit scale): %s\n", nm, 78, nrow(lc),
              zstr(it[1] / 50, it[3], it[4])))
  for (i in seq_len(nrow(et))) cat(sprintf("  %s ancestry slope (averaged over watering treatments): %s; odds ratio per 10 points of Fremont ancestry %.3f (95%% CI %.3f-%.3f)\n",
                                           toupper(et$fungal_type[i]), zstr(et$slope[i], et$z[i], et$p[i]), et$or10[i], et$or10_lo[i], et$or10_hi[i]))
  rows[[nm]] <- bind_rows(tibble(model = nm, quantity = "EMF vs. AMF ancestry slope, well-watered (logit, per 1%)", string = zstr(it[1] / 50, it[3], it[4]), p = it[4]),
                          tibble(model = nm, quantity = paste(toupper(et$fungal_type), "ancestry slope, averaged over treatments (logit, per 1%)"),
                                 string = sprintf("%s; OR per 10 points %.3f (95%% CI %.3f–%.3f)", mapply(zstr, et$slope, et$z, et$p), et$or10, et$or10_lo, et$or10_hi), p = et$p)) }
readr::write_excel_csv(bind_rows(rows), file.path(out, "colonization_counts_glmm.csv"))
cat("\nReported Gaussian LMM for comparison (Figure 3): EMF slope –0.124 percentage points per 1% (t = –4.80, P < 0.001); ancestry × EMF t = –3.51, P < 0.001; AMF and DSE slopes P = 0.68 and 0.62.\n")

# ------------------------------------------------------------------------------
hdr("B. Survival GLMM: population fixed (reported) vs. random")
# ------------------------------------------------------------------------------
fx <- glmer(alive_01_10_2021 ~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat + (1 | geno_id), data = ds, family = binomial, control = ctl)
fxf <- glmer(alive_01_10_2021 ~ pf50 * treatment + I(pf50^2) * treatment + pop + og_drc_mm + soil + myco_cat + (1 | geno_id), data = ds, family = binomial, control = ctl)
rn <- glmer(alive_01_10_2021 ~ pf50 + treatment + I(pf50^2) + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id), data = ds, family = binomial, control = ctl)
rnf <- glmer(alive_01_10_2021 ~ pf50 * treatment + I(pf50^2) * treatment + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id), data = ds, family = binomial, control = ctl)
div <- c(pf50 = 50, "I(pf50^2)" = 2500, treatmentDrought = 1)
lab <- c(pf50 = "Linear ancestry (per 1%)", "I(pf50^2)" = "Quadratic ancestry", treatmentDrought = "Drought")
sv <- list()
for (nm in c("Population fixed (reported)", "Population random")) {
  m <- if (grepl("fixed", nm)) fx else rn; mf <- if (grepl("fixed", nm)) fxf else rnf
  l <- anova(m, mf); cf <- coef(summary(m))
  vc <- as.data.frame(VarCorr(m))
  cat(sprintf("\n%s: n = %d | AIC = %.1f | singular = %s | random-effect SDs: %s\n", nm, nobs(m), AIC(m), isSingular(m),
              paste(vc$grp, signif(vc$sdcor, 3), sep = " = ", collapse = "; ")))
  cat(sprintf("  ancestry x drought (linear + quadratic) LRT %s(%d) = %.2f, %s\n", CHI, l$Df[2], l$Chisq[2], P(l$`Pr(>Chisq)`[2])))
  for (t in names(lab)) cat(sprintf("  %-26s %s\n", lab[t], zstr(cf[t, 1] / div[t], cf[t, 3], cf[t, 4])))
  sv[[nm]] <- tibble(model = nm, n = nobs(m), AIC = AIC(m), term = c("Ancestry × drought (LRT, 2 df)", unname(lab)),
                     string = c(sprintf("%s(%d) = %s, %s", CHI, l$Df[2], dec(l$Chisq[2]), P(l$`Pr(>Chisq)`[2])), sapply(names(lab), function(t) zstr(cf[t, 1] / div[t], cf[t, 3], cf[t, 4]))),
                     p = c(l$`Pr(>Chisq)`[2], cf[names(lab), 4])) }
readr::write_excel_csv(bind_rows(sv), file.path(out, "survival_population.csv"))
pv <- as.data.frame(VarCorr(rn)) %>% filter(grp == "pop")
cat(sprintf("\nPopulation SD with population random: %.3f (6 populations). Population fixed-effect LRT in the reported model: ", pv$sdcor))
fx0 <- update(fx, . ~ . - pop); lp <- anova(fx0, fx); cat(sprintf("%s(%d) = %.2f, %s\n", CHI, lp$Df[2], lp$Chisq[2], P(lp$`Pr(>Chisq)`[2])))
sink()
