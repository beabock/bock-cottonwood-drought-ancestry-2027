# ==============================================================================
# EMF ANCESTRY CONTRASTS: exact estimates for every version of the EMF ancestry result.
# Refits the reported colonization models with their published code (08_sample_flow_extra_tests.R,
# 02a_main_models.R, 05_control_inoculum_sensitivity.R) and the genotype-mean check
# (08_sample_flow_extra_tests.R). No model is changed; only contrasts are extracted.
# Run from the project root:  Rscript r_code/13_emf_ancestry_contrasts.R
# ==============================================================================
suppressMessages({ source("r_code/helpers/utilities_data_prep.R"); library(tidyverse); library(lmerTest); library(emmeans); library(nlme) })
out <- "results/13_emf_ancestry_contrasts"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
cat("EMF ANCESTRY CONTRASTS run:", "|", R.version.string, "| lme4", as.character(packageVersion("lme4")),
    "| lmerTest", as.character(packageVersion("lmerTest")), "| emmeans", as.character(packageVersion("emmeans")), "| nlme", as.character(packageVersion("nlme")), "\n\n")
emm_options(lmerTest.limit = 5000, pbkrtest.limit = 5000)
ds <- load_and_prepare_data() %>% mutate(geno_id = paste(pop, type, geno, sep = "_"))
mk_dl <- function(d) d %>% filter(!is.na(amf), !is.na(emf), !is.na(dse)) %>%
  select(tree_id, pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>% na.omit() %>%
  mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm)), myco_cat = droplevels(factor(myco_cat))) %>%
  pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>%
  mutate(fungal_type = factor(fungal_type, levels = c("amf", "dse", "emf")), ft = as.integer(fungal_type))
fit_lmer <- function(dl) lmerTest::lmer(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil +
                                          (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id), data = dl)
rows <- list()
add <- function(contrast, model, n, est, t, df, p, where) rows[[length(rows) + 1]] <<- tibble(contrast, model, n_plants = n, estimate = est, t, df, P = p, manuscript_location = where)
# contrasts from emtrends: EMF - AMF slope, averaged over watering and within each watering treatment
emf_amf <- function(m, df_method = NULL, data = NULL) {
  args <- list(m, specs = ~ fungal_type, var = "per_f"); if (!is.null(df_method)) args$lmer.df <- df_method; if (!is.null(data)) args$data <- data
  avg <- as.data.frame(summary(contrast(do.call(emtrends, args), list(`emf - amf` = c(-1, 0, 1)))))
  args$specs <- ~ fungal_type | treatment
  byw <- as.data.frame(summary(contrast(do.call(emtrends, args), list(`emf - amf` = c(-1, 0, 1)))))
  list(avg = avg, byw = byw) }

# ---- reported lmer, n = 78 ----
dl <- mk_dl(ds); col <- fit_lmer(dl); n78 <- n_distinct(dl$tree_id); cs <- coef(summary(col))
cat("Reported lmer: plants =", n78, "| coefficient per_f:fungal_typeemf (EMF - AMF slope, well-watered):\n"); print(cs["per_f:fungal_typeemf", ], digits = 7)
c1 <- emf_amf(col, "satterthwaite"); print(c1)
add("EMF − AMF ancestry slope, averaged over watering", "Reported LMM (Satterthwaite)", n78, c1$avg$estimate, c1$avg$t.ratio, c1$avg$df, c1$avg$p.value, "Not currently reported")
w <- filter(c1$byw, treatment == "Well-watered")
add("EMF − AMF ancestry slope, well-watered only", "Reported LMM (Satterthwaite)", n78, cs["per_f:fungal_typeemf", "Estimate"], cs["per_f:fungal_typeemf", "t value"],
    cs["per_f:fungal_typeemf", "df"], cs["per_f:fungal_typeemf", "Pr(>|t|)"], "Results colonization paragraph and Figure 3 caption ('Fremont ancestry × EMF interaction'); Appendix S4 table, all plants")
stopifnot(abs(w$estimate - cs["per_f:fungal_typeemf", "Estimate"]) < 1e-8)
d <- filter(c1$byw, treatment == "Drought")
add("EMF − AMF ancestry slope, drought only", "Reported LMM (Satterthwaite)", n78, d$estimate, d$t.ratio, d$df, d$p.value, "Not currently reported")
sl <- as.data.frame(summary(emtrends(col, ~ fungal_type, var = "per_f", lmer.df = "satterthwaite")))
slk <- as.data.frame(summary(emtrends(col, ~ fungal_type, var = "per_f", lmer.df = "kenward-roger")))
cat("\nSlopes, Satterthwaite:\n"); print(sl); cat("Slopes, Kenward-Roger:\n"); print(slk)
e <- filter(sl, fungal_type == "emf"); add("EMF ancestry slope (per 1% Fremont), averaged over watering", "Reported LMM (Satterthwaite)", n78, e$per_f.trend, e$per_f.trend / e$SE, e$df, 2 * pt(-abs(e$per_f.trend / e$SE), e$df), "Figure 3 caption; Results colonization paragraph ('–0.124 ... about 12 percentage points')")
e <- filter(slk, fungal_type == "emf"); add("EMF ancestry slope (per 1% Fremont), averaged over watering", "Reported LMM (Kenward–Roger)", n78, e$per_f.trend, e$per_f.trend / e$SE, e$df, 2 * pt(-abs(e$per_f.trend / e$SE), e$df), "Results colonization paragraph ('t = –4.71, df = 142')")
ck <- emf_amf(col, "kenward-roger"); cat("\nEMF - AMF, Kenward-Roger:\n"); print(ck)

# ---- genotype means (08_sample_flow_extra_tests.R) ----
gm <- dl %>% group_by(geno_id, pop, fungal_type) %>% summarise(per_f = mean(per_f), colonization = mean(colonization), .groups = "drop")
g1 <- coef(summary(lm(colonization ~ per_f, data = filter(gm, fungal_type == "emf"))))["per_f", ]
cat("\nGenotype means, EMF (OLS):\n"); print(g1, digits = 7)
add("EMF ancestry slope (per 1% Fremont)", "Genotype means, OLS (one point per genotype)", n78, g1["Estimate"], g1["t value"], nrow(filter(gm, fungal_type == "emf")) - 2, g1["Pr(>|t|)"], "Results colonization paragraph ('β = –0.155, t = –4.60, df = 35')")
# genotype-mean EMF - AMF difference in slopes (same plants; one point per genotype x guild)
g2 <- coef(summary(lm(colonization ~ per_f * fungal_type, data = filter(gm, fungal_type != "dse") %>% droplevels())))
add("EMF − AMF ancestry slope", "Genotype means, OLS (per_f × guild, AMF and EMF only)", n78, g2["per_f:fungal_typeemf", "Estimate"], g2["per_f:fungal_typeemf", "t value"],
    2 * n_distinct(gm$geno_id) - 4, g2["per_f:fungal_typeemf", "Pr(>|t|)"], "Not currently reported")

# ---- nlme sensitivity (02a_main_models.R) ----
nl <- nlme::lme(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil + myco_cat,
                random = ~ 1 | pop/geno_id, correlation = nlme::corSymm(form = ~ ft | pop/geno_id/tree_id),
                weights = nlme::varIdent(form = ~ 1 | fungal_type), data = dl, method = "REML",
                control = nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
tt <- summary(nl)$tTable; cat("\nnlme coefficient per_f:fungal_typeemf:\n"); print(tt["per_f:fungal_typeemf", ], digits = 7)
add("EMF − AMF ancestry slope, well-watered only", "nlme: within-plant correlation + guild variances (containment df)", n78, tt["per_f:fungal_typeemf", "Value"], tt["per_f:fungal_typeemf", "t-value"],
    tt["per_f:fungal_typeemf", "DF"], tt["per_f:fungal_typeemf", "p-value"], "Results colonization paragraph ('with within-plant correlation ... t = –2.98, df = 180')")
cn <- emf_amf(nl, data = dl); cat("\nnlme EMF - AMF via emtrends (emmeans default df for lme):\n"); print(cn)
add("EMF − AMF ancestry slope, averaged over watering", "nlme (emmeans containment df)", n78, cn$avg$estimate, cn$avg$t.ratio, cn$avg$df, cn$avg$p.value, "Not currently reported")
sn <- as.data.frame(summary(emtrends(nl, ~ fungal_type, var = "per_f", data = dl))); e <- filter(sn, fungal_type == "emf")
add("EMF ancestry slope (per 1% Fremont), averaged over watering", "nlme (emmeans containment df)", n78, e$per_f.trend, e$per_f.trend / e$SE, e$df, 2 * pt(-abs(e$per_f.trend / e$SE), e$df), "Not currently reported")

# ---- without control plants (05_control_inoculum_sensitivity.R) ----
dln <- mk_dl(filter(ds, myco_cat != "None")); coln <- fit_lmer(dln); nn <- n_distinct(dln$tree_id); csn <- coef(summary(coln))
cat("\nNo-control lmer: plants =", nn, "\n"); print(csn["per_f:fungal_typeemf", ], digits = 7)
add("EMF − AMF ancestry slope, well-watered only", "Reported LMM without control plants (Satterthwaite)", nn, csn["per_f:fungal_typeemf", "Estimate"], csn["per_f:fungal_typeemf", "t value"],
    csn["per_f:fungal_typeemf", "df"], csn["per_f:fungal_typeemf", "Pr(>|t|)"], "Results colonization paragraph ('n = 57; t = –2.94, df = 155'); Appendix S4 table, inoculated only")
c2 <- emf_amf(coln, "satterthwaite"); print(c2)
add("EMF − AMF ancestry slope, averaged over watering", "Reported LMM without control plants (Satterthwaite)", nn, c2$avg$estimate, c2$avg$t.ratio, c2$avg$df, c2$avg$p.value, "Not currently reported")
s2 <- as.data.frame(summary(emtrends(coln, ~ fungal_type, var = "per_f", lmer.df = "satterthwaite"))); e <- filter(s2, fungal_type == "emf")
add("EMF ancestry slope (per 1% Fremont), averaged over watering", "Reported LMM without control plants (Satterthwaite)", nn, e$per_f.trend, e$per_f.trend / e$SE, e$df, 2 * pt(-abs(e$per_f.trend / e$SE), e$df), "Not currently reported")

tab <- bind_rows(rows)
cat("\n\nSummary table:\n"); print(as.data.frame(tab %>% mutate(across(c(estimate, t, df, P), ~ signif(.x, 6)))), right = FALSE)
readr::write_excel_csv(tab, file.path(out, "emf_contrasts_raw.csv"))
sink()
