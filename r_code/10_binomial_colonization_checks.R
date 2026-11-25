# ==============================================================================
# SURVIVAL GLMM WITH POPULATION RANDOM: the reported survival GLMM (population fixed, genotype random)
# refit with population as a random intercept (as in the other models) instead of a fixed effect.
# Source of the Methods sentence "Treating population as a random intercept instead gave the same conclusions".
# Changes no reported result. (The file name is historical: a binomial-count refit of the colonization
# model that used to be part A of this script was not reported in the manuscript and was removed.)
# ==============================================================================
# Run from the project root (R 4.6.1):  Rscript r_code/10_binomial_colonization_checks.R
# Writes: results/10_binomial_colonization_checks/output.txt, survival_population.csv
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(lme4); library(lmerTest); library(tidyverse)
out <- "results/10_binomial_colonization_checks"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("SURVIVAL POPULATION CHECK run:", "|", R.version.string, "| lme4", as.character(packageVersion("lme4")), "\n")
EN <- "–"; BETA <- "β"; TIMES <- "×"; CHI <- "χ²"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
zstr <- function(b, z, p) sprintf("%s = %s, z = %s, %s", BETA, num(b), dec(z), P(p))
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, pf50 = per_fc / 50, geno_id = paste(pop, type, geno, sep = "_"))

# ------------------------------------------------------------------------------
hdr("Survival GLMM: population fixed (reported) vs. random")
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
