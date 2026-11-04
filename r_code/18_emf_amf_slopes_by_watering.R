# ==============================================================================
# Does the EMF - AMF ancestry-slope difference depend on watering?
# The test is the ancestry x fungal type (EMF vs AMF) x drought coefficient, per_f:fungal_typeemf:treatmentDrought
# (= drought-minus-well-watered change in the EMF - AMF slope difference). Fits as in r_code/13_emf_ancestry_contrasts.R:
# reported lmer (78 plants), nlme sensitivity (02a_main_models.R), no-control lmer (05_control_inoculum_sensitivity.R),
# and genotype means computed within genotype x watering (OLS; one point per genotype x watering x guild).
# Run from the project root:  Rscript r_code/18_emf_amf_slopes_by_watering.R
# ==============================================================================
suppressMessages({ source("r_code/helpers/utilities_data_prep.R"); library(tidyverse); library(lmerTest); library(emmeans); library(nlme) })
out <- "results/18_emf_amf_slopes_by_watering"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
cat("EMF AND AMF SLOPES BY WATERING run:", "|", R.version.string, "\n\n")
EN <- "–"; dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
tstr <- function(b, t, df, p) sprintf("β = %s, t = %s, df = %s, %s", num(b), dec(t), dff(df), P(p))
ds <- load_and_prepare_data() %>% mutate(geno_id = paste(pop, type, geno, sep = "_"))
mk_dl <- function(d) d %>% filter(!is.na(amf), !is.na(emf), !is.na(dse)) %>%
  select(tree_id, pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>% na.omit() %>%
  mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm)), myco_cat = droplevels(factor(myco_cat))) %>%
  pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>%
  mutate(fungal_type = factor(fungal_type, levels = c("amf", "dse", "emf")), ft = as.integer(fungal_type))
fit_lmer <- function(dl) lmerTest::lmer(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil +
                                          (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id), data = dl)
rows <- list()
add <- function(fit, n, contrast, est, t, df, p) rows[[length(rows) + 1]] <<- tibble(fit, n_plants = n, contrast, estimate = est, t, df, P = p, string = tstr(est, t, df, p))
L3 <- "per_f:fungal_typeemf:treatmentDrought"; L2 <- "per_f:fungal_typeemf"
emf_amf <- function(m, df_method = NULL, data = NULL) {
  a <- list(m, specs = ~ fungal_type, var = "per_f"); if (!is.null(df_method)) a$lmer.df <- df_method; if (!is.null(data)) a$data <- data
  avg <- as.data.frame(summary(contrast(do.call(emtrends, a), list(`emf - amf` = c(-1, 0, 1)))))
  a$specs <- ~ fungal_type | treatment
  byw <- as.data.frame(summary(contrast(do.call(emtrends, a), list(`emf - amf` = c(-1, 0, 1)))))
  list(avg = avg, byw = byw) }
lmer_block <- function(m, lab) {
  n <- n_distinct(model.frame(m)$tree_id); cs <- coef(summary(m)); c3 <- emf_amf(m, "satterthwaite")
  add(lab, n, "EMF − AMF slope, averaged over watering", c3$avg$estimate, c3$avg$t.ratio, c3$avg$df, c3$avg$p.value)
  for (w in c("Well-watered", "Drought")) { r <- filter(c3$byw, treatment == w); add(lab, n, paste0("EMF − AMF slope, ", tolower(w)), r$estimate, r$t.ratio, r$df, r$p.value) }
  add(lab, n, "Difference (drought − well-watered): ancestry × EMF-vs-AMF × drought", cs[L3, 1], cs[L3, "t value"], cs[L3, "df"], cs[L3, "Pr(>|t|)"]) }

dl <- mk_dl(ds); col <- fit_lmer(dl); lmer_block(col, "Reported LMM (Satterthwaite)")
nl <- nlme::lme(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil + myco_cat, random = ~ 1 | pop/geno_id,
                correlation = nlme::corSymm(form = ~ ft | pop/geno_id/tree_id), weights = nlme::varIdent(form = ~ 1 | fungal_type), data = dl, method = "REML",
                control = nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
tt <- summary(nl)$tTable; cn <- emf_amf(nl, data = dl); n78 <- n_distinct(dl$tree_id)
add("nlme (within-plant correlation; containment df)", n78, "EMF − AMF slope, averaged over watering", cn$avg$estimate, cn$avg$t.ratio, cn$avg$df, cn$avg$p.value)
add("nlme (within-plant correlation; containment df)", n78, "EMF − AMF slope, well-watered", tt[L2, "Value"], tt[L2, "t-value"], tt[L2, "DF"], tt[L2, "p-value"])
r <- filter(cn$byw, treatment == "Drought"); add("nlme (within-plant correlation; containment df)", n78, "EMF − AMF slope, drought", r$estimate, r$t.ratio, r$df, r$p.value)
add("nlme (within-plant correlation; containment df)", n78, "Difference (drought − well-watered): ancestry × EMF-vs-AMF × drought", tt[L3, "Value"], tt[L3, "t-value"], tt[L3, "DF"], tt[L3, "p-value"])
coln <- fit_lmer(mk_dl(filter(ds, myco_cat != "None"))); lmer_block(coln, "Reported LMM without control plants (Satterthwaite)")
gm <- dl %>% group_by(geno_id, pop, treatment, fungal_type) %>% summarise(per_f = mean(per_f), colonization = mean(colonization), n_plants = n(), .groups = "drop")
og <- lm(colonization ~ per_f * fungal_type * treatment, data = gm); cg <- coef(summary(og)); dfg <- og$df.residual
ng <- nrow(distinct(gm, geno_id, treatment))
gtr <- function(L, lab) { b <- coef(og); V <- vcov(og); L <- L[names(b)]; L[is.na(L)] <- 0; e <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% V %*% L))
  add("Genotype × watering means, OLS (78 plants; one point per genotype × watering × guild)", n78, lab, e, e / se, dfg, 2 * pt(-abs(e / se), dfg)) }
nm <- names(coef(og)); z <- setNames(rep(0, length(nm)), nm)
L_ww <- z; L_ww["per_f:fungal_typeemf"] <- 1
L_d <- z; L_d["per_f:fungal_typeemf"] <- 1; L_d["per_f:fungal_typeemf:treatmentDrought"] <- 1
L_avg <- z; L_avg["per_f:fungal_typeemf"] <- 1; L_avg["per_f:fungal_typeemf:treatmentDrought"] <- 0.5
gtr(L_avg, "EMF − AMF slope, averaged over watering"); gtr(L_ww, "EMF − AMF slope, well-watered"); gtr(L_d, "EMF − AMF slope, drought")
gtr(setNames(as.numeric(nm == "per_f:fungal_typeemf:treatmentDrought"), nm), "Difference (drought − well-watered): ancestry × EMF-vs-AMF × drought")
cat("Genotype x watering cells:", ng, "| rows:", nrow(gm), "| residual df:", dfg, "\n\n")
tab <- bind_rows(rows)
print(as.data.frame(tab %>% select(fit, n_plants, contrast, string)), right = FALSE, row.names = FALSE)
write.csv(tab, file.path(out, "emf_amf_by_watering.csv"), row.names = FALSE)
sink()
