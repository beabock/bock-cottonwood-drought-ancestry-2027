# ==============================================================================
# MAIN MODELS: every manuscript model refit with mother-tree (genotype) random
# intercepts, population in the PERMANOVA, population random intercept in the
# aggregated survival model, and the stomatal model as described in the Methods.
# ==============================================================================
# Run from the project root:  Rscript r_code/02a_main_models.R
# Writes: results/02_main_models/main_models_output.txt, results/02_main_models/genotype_map.csv, results/02_main_models/convergence.csv
#
# Genotype ID = pop_type_geno (inferred mother tree; per_f is identical within
# each ID). Because the ID already contains pop, (1 | pop:geno_id) is genotype
# nested in population. Survival keeps population as a fixed effect (as in the
# earlier model) and adds (1 | geno_id); a pop-random version is also fit.
# Initial diameter is og_drc_mm everywhere (not used in the PERMANOVA).
# ==============================================================================

stopifnot(file.exists("data/edited_data/mega_ds.csv"))
source("r_code/helpers/utilities_data_prep.R")
library(vegan); library(pairwiseAdonis)
library(lme4); library(lmerTest); library(glmmTMB); library(tidyverse); library(emmeans)

dir.create("results", showWarnings = FALSE)
sink("results/02_main_models/main_models_output.txt", split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("MAIN MODELS run:", "\n")

ds <- load_and_prepare_data() %>%
  mutate(per_fc = per_f - 50, pf50 = per_fc / 50,
         geno_id = paste(pop, type, geno, sep = "_"))
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# convergence / singularity log -------------------------------------------------
conv_log <- list()
reg <- function(name, m) {
  msg <- tryCatch(if (inherits(m, "merMod")) paste(m@optinfo$conv$lme4$messages, collapse = "; ") else "", error = function(e) "")
  sing <- tryCatch(if (inherits(m, "merMod")) isSingular(m) else NA, error = function(e) NA)
  vc <- tryCatch(as.data.frame(VarCorr(m)), error = function(e) NULL)
  gv <- if (!is.null(vc) && "grp" %in% names(vc)) vc$sdcor[vc$grp == "pop:geno_id" | vc$grp == "geno_id"] else NA
  pv <- if (!is.null(vc) && "grp" %in% names(vc)) vc$sdcor[vc$grp == "pop"] else NA
  conv_log[[length(conv_log) + 1]] <<- tibble(model = name, converged = (msg == "" | grepl("singular", msg)),
                                              singular = sing, messages = msg,
                                              genotype_sd = if (length(gv)) gv[1] else NA, pop_sd = if (length(pv)) pv[1] else NA)
}

# ------------------------------------------------------------------------------
# Genotype map
# ------------------------------------------------------------------------------
ds <- ds %>% mutate(has_amf = !is.na(amf), has_emf = !is.na(emf), has_dse = !is.na(dse),
                    all3 = has_amf & has_emf & has_dse)
gmap <- ds %>% group_by(geno_id, pop, type, geno) %>%
  summarise(n_cuttings = n(), per_f = round(mean(per_f), 3), n_distinct_per_f = n_distinct(round(per_f, 4)),
            new_class = paste(sort(unique(as.character(new_class))), collapse = "/"),
            n_alive = sum(alive_01_10_2021 == 1, na.rm = TRUE), n_with_biomass = sum(!is.na(dryshoot_new)),
            n_all3_colonization = sum(all3), n_with_porometer = sum(!is.na(por_avg)),
            tree_idold_examples = paste(head(tree_idold, 3), collapse = ";"), .groups = "drop") %>%
  arrange(pop, type, geno)
write.csv(gmap, "results/02_main_models/genotype_map.csv", row.names = FALSE)
hdr("Genotype map (saved to results/02_main_models/genotype_map.csv)")
cat("genotypes:", nrow(gmap), "| cuttings per genotype:", paste(range(gmap$n_cuttings), collapse = "-"),
    "| genotypes with >1 distinct per_f:", sum(gmap$n_distinct_per_f > 1),
    "| genotypes spanning >1 crosstype label:", sum(grepl("/", gmap$new_class)), "\n")
print(table(gmap$new_class))

# ------------------------------------------------------------------------------
# PERMANOVA (population as a term; no diameter). Runs first so the RNG stream is
# reproducible (seed 1998 set in utilities_data_prep.R).
# ------------------------------------------------------------------------------
hdr("PERMANOVA: crosstype and population in the same model; strata = inoculum category")
fc <- ds %>% filter(new_class %in% c("Narrowleaf", "F1 Hybrid", "Fremont")) %>%
  select(new_class, pop, treatment, amf, emf, dse, soil, myco_cat) %>% na.omit()
fmat <- fc %>% select(amf, emf, dse); fenv <- fc %>% select(new_class, pop, treatment, soil, myco_cat)
permA <- adonis2(fmat ~ new_class + pop + treatment + soil, strata = fenv$myco_cat, data = fenv,
                 method = "euclidean", permutations = 999, by = "terms")
fcb <- ds %>% select(per_f, pop, treatment, amf, emf, dse, soil, myco_cat) %>% na.omit()
fmatb <- fcb %>% select(amf, emf, dse); fenvb <- fcb %>% select(per_f, pop, treatment, soil, myco_cat)
permB <- adonis2(fmatb ~ per_f + pop + treatment + soil, strata = fenvb$myco_cat, data = fenvb,
                 method = "euclidean", permutations = 999, by = "terms")
pw <- pairwise.adonis(fmat, fenv$new_class, sim.method = "euclidean", p.adjust.m = "bonferroni")
cat("\nModel A4 (categorical, 3 focal crosstypes): fungal_matrix ~ new_class + pop + treatment + soil; n =", nrow(fc),
    "; sequential terms; permutations restricted within inoculum category\n"); print(permA)
cat("\nModel B4 (continuous ancestry, all crosstypes): fungal_matrix ~ per_f + pop + treatment + soil; n =", nrow(fcb), "\n"); print(permB)
cat("\nPairwise crosstype (crosstype only, no strata/covariates):\n"); print(pw)
permA_m <- adonis2(fmat ~ new_class + pop + treatment + soil, strata = fenv$myco_cat, data = fenv,
                   method = "euclidean", permutations = 999, by = "margin")
permB_m <- adonis2(fmatb ~ per_f + pop + treatment + soil, strata = fenvb$myco_cat, data = fenvb,
                   method = "euclidean", permutations = 999, by = "margin")
cat("\nOrder-independent (marginal) tests, Model A4:\n"); print(permA_m)
cat("\nOrder-independent (marginal) tests, Model B4:\n"); print(permB_m)
cat("\nNote: genotype cannot be a PERMANOVA term (it determines ancestry, so it is collinear with crosstype/per_f).\n",
    "Permutations here still treat trees as exchangeable within inoculum category.\n", sep = "")

# ------------------------------------------------------------------------------
# Colonization LMM
# ------------------------------------------------------------------------------
hdr("COLONIZATION LMM: + genotype (pop:geno_id) and tree_id")
dl <- ds %>% filter(all3) %>%
  select(tree_id, pop, geno_id, per_f, new_class, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>% na.omit() %>%
  mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm))) %>%
  pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>%
  mutate(fungal_type = factor(fungal_type, levels = c("amf", "dse", "emf")), ft = as.integer(fungal_type),
         myco_cat = factor(myco_cat))
cat("trees:", n_distinct(dl$tree_id), "| genotypes:", n_distinct(dl$geno_id), "| rows:", nrow(dl), "\n")
fx <- colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil
col4 <- lmerTest::lmer(update(fx, . ~ . + (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id)), data = dl)
reg("colonization lmer (+pop, myco_cat, pop:geno_id, tree_id)", col4)
print(summary(col4))
col4n <- tryCatch(nlme::lme(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil + myco_cat,
                            random = ~ 1 | pop/geno_id,
                            correlation = nlme::corSymm(form = ~ ft | pop/geno_id/tree_id),
                            weights = nlme::varIdent(form = ~ 1 | fungal_type), data = dl, method = "REML",
                            control = nlme::lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)),
                  error = function(e) { cat("nlme (pop/geno_id + corSymm) DID NOT CONVERGE:", conditionMessage(e), "\n"); NULL })
if (!is.null(col4n)) {
  cat("\nnlme sensitivity (random ~1|pop/geno_id; corSymm within tree; varIdent by fungal type):\n")
  print(summary(col4n)$tTable); print(col4n$modelStruct$corStruct); print(col4n$modelStruct$varStruct)
  print(nlme::VarCorr(col4n))
}
cat("\nAncestry slope by fungal type (lmer):\n")
print(summary(emtrends(col4, ~ fungal_type, var = "per_f", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE)))
if (!is.null(col4n)) { cat("\nAncestry slope by fungal type (nlme):\n")
  print(summary(emtrends(col4n, ~ fungal_type, var = "per_f", data = dl), infer = c(TRUE, TRUE))) }

# ------------------------------------------------------------------------------
# Biomass (full, reduced) + vertex
# ------------------------------------------------------------------------------
hdr("BIOMASS: log(dryshoot_new), (1|pop) + (1|pop:geno_id)")
bio_full4 <- lmerTest::lmer(log(dryshoot_new) ~ per_fc * treatment + I(per_fc^2) * treatment + og_drc_mm + myco_cat + soil +
                              (1 | pop) + (1 | pop:geno_id), data = ds)
bio_red4 <- lmerTest::lmer(log(dryshoot_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil +
                             (1 | pop) + (1 | pop:geno_id), data = ds)
reg("biomass full", bio_full4); reg("biomass reduced", bio_red4)
cat("\n--- FULL (ancestry x drought) ---\n"); print(summary(bio_full4))
cat("\n--- REDUCED (no interactions) ---\n"); print(summary(bio_red4))
bio_lrt4 <- anova(bio_red4, bio_full4, refit = TRUE)
cat("\n--- LRT ancestry (linear + quadratic) x drought ---\n"); print(bio_lrt4)

hdr("BIOMASS PEAK: vertex = 50 - b/(2c) from the reduced model")
vertex_fun <- function(m) { b <- fixef(m); unname(50 - b["per_fc"] / (2 * b["I(per_fc^2)"])) }
bio_red_l <- lme4::lmer(log(dryshoot_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil +
                          (1 | pop) + (1 | pop:geno_id), data = ds)
v_hat <- vertex_fun(bio_red_l)
bb <- fixef(bio_red_l)["per_fc"]; cc <- fixef(bio_red_l)["I(per_fc^2)"]
V <- as.matrix(vcov(bio_red_l))[c("per_fc", "I(per_fc^2)"), c("per_fc", "I(per_fc^2)")]
gr <- c(-1 / (2 * cc), bb / (2 * cc^2)); se_d <- sqrt(as.numeric(t(gr) %*% V %*% gr))
cat("b =", signif(bb, 4), " c =", signif(cc, 4), "\nVertex =", round(v_hat, 2), "% Fremont (Q =", round(1 - v_hat / 100, 3), ")\n")
cat("Delta-method 95% CI:", round(v_hat - 1.96 * se_d, 1), "to", round(v_hat + 1.96 * se_d, 1), "\n")
vc_fun <- function(m) { b <- fixef(m); c(vertex = unname(50 - b["per_fc"] / (2 * b["I(per_fc^2)"])), c = unname(b["I(per_fc^2)"])) }
bs <- lme4::bootMer(bio_red_l, vc_fun, nsim = 2000, type = "parametric", use.u = FALSE, seed = 1998)
vb <- bs$t[, 1]; cb <- bs$t[, 2]
cat("Parametric bootstrap (2000 sims, bootMer, use.u = FALSE): median", round(median(vb), 1), "| 95% percentile CI (all sims)",
    paste(round(quantile(vb, c(0.025, 0.975)), 1), collapse = " to "), "| outside 0-100:", sum(vb < 0 | vb > 100), "of 2000\n")
cat("Refits where the quadratic is NOT concave (c >= 0, so there is no peak):", sum(cb >= 0), "of 2000 (",
    round(100 * mean(cb >= 0), 1), "%)\n")
vb_c <- vb[cb < 0]
cat("95% percentile CI among concave refits only:", paste(round(quantile(vb_c, c(0.025, 0.975)), 1), collapse = " to "),
    "| share of concave refits with vertex inside 0-100:", round(100 * mean(vb_c >= 0 & vb_c <= 100), 1), "%\n")
cat("Original estimate: c =", signif(cc, 3), "(SE", signif(sqrt(V[2, 2]), 3), "), t =", round(cc / sqrt(V[2, 2]), 2), "\n")

# ------------------------------------------------------------------------------
# Survival (full, reduced) -- glmer on ancestry rescaled by 50 (glmmTMB does not
# start with a genotype term); estimates converted back to per-1% units.
# ------------------------------------------------------------------------------
hdr("SURVIVAL: binomial GLMM, population fixed + (1 | geno_id)")
tidy_surv <- function(m) {
  cs <- coef(summary(m)); rn <- rownames(cs)
  div <- ifelse(grepl("pf50\\^2", rn), 2500, ifelse(grepl("pf50", rn), 50, 1))
  cs[, "Estimate"] <- cs[, "Estimate"] / div; cs[, "Std. Error"] <- cs[, "Std. Error"] / div
  rownames(cs) <- gsub("pf50", "per_fc", rn, fixed = TRUE); cs
}
surv_full4 <- glmer(alive_01_10_2021 ~ pf50 * treatment + I(pf50^2) * treatment + pop + og_drc_mm + soil + myco_cat + (1 | geno_id),
                    data = ds, family = binomial, control = ctl)
surv_red4 <- glmer(alive_01_10_2021 ~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat + (1 | geno_id),
                   data = ds, family = binomial, control = ctl)
reg("survival full (pop fixed, geno_id)", surv_full4); reg("survival reduced (pop fixed, geno_id)", surv_red4)
cat("n =", nobs(surv_full4), "\n")
cat("\n--- FULL: raw glmer summary (ancestry in units of 50 points) ---\n"); print(summary(surv_full4))
cat("\n--- FULL: coefficients converted to per-1% ancestry units ---\n"); print(round(tidy_surv(surv_full4), 5))
cat("\n--- REDUCED: raw glmer summary ---\n"); print(summary(surv_red4))
cat("\n--- REDUCED: coefficients converted to per-1% ancestry units ---\n"); print(round(tidy_surv(surv_red4), 5))
sf <- tidy_surv(surv_full4); sr <- tidy_surv(surv_red4)
surv_red4p <- glmer(alive_01_10_2021 ~ pf50 + treatment + I(pf50^2) + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id),
                    data = ds, family = binomial, control = ctl)
surv_full4p <- glmer(alive_01_10_2021 ~ pf50 * treatment + I(pf50^2) * treatment + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id),
                     data = ds, family = binomial, control = ctl)
reg("survival full (pop random + pop:geno_id)", surv_full4p); reg("survival reduced (pop random + pop:geno_id)", surv_red4p)
cat("\n--- Alternative: population RANDOM + genotype nested (reduced), per-1% units ---\n"); print(round(tidy_surv(surv_red4p), 5))
print(VarCorr(surv_red4p))
cat("--- Alternative (full) interaction terms:\n"); print(round(tidy_surv(surv_full4p)[c("per_fc:treatmentDrought", "treatmentDrought:I(per_fc^2)"), ], 5))

# ------------------------------------------------------------------------------
# Colonization x biomass (n = 77) with genotype
# ------------------------------------------------------------------------------
hdr("COLONIZATION x BIOMASS: + (1 | pop:geno_id)")
fb_pf4 <- lmerTest::lmer(log(dryshoot_new) ~ scale(amf) * treatment + scale(emf) * treatment + scale(dse) * treatment +
                           per_fc + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
fb_nv4 <- lmerTest::lmer(log(dryshoot_new) ~ scale(amf) * treatment + scale(emf) * treatment + scale(dse) * treatment +
                           og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
reg("colonization-biomass ancestry-controlled", fb_pf4); reg("colonization-biomass naive", fb_nv4)
cat("--- ancestry-controlled ---\n"); print(summary(fb_pf4)); cat("\n--- naive ---\n"); print(summary(fb_nv4))
zt <- ds %>% filter(!is.na(emf), !is.na(amf), !is.na(dse), !is.na(dryshoot_new), !is.na(og_drc_mm), !is.na(myco_cat)) %>%
  mutate(z = as.numeric(scale(log(dryshoot_new)))) %>% arrange(desc(abs(z)))
cat("\nTrees with |z| > 3 for log(dryshoot_new) among these", nrow(zt), ":", sum(abs(zt$z) > 3), "(max |z| =", round(max(abs(zt$z)), 2), ", tree", zt$tree_idold[1], ")\n")
cat("Max |scaled residual|, ancestry-controlled model:", round(max(abs(resid(fb_pf4, scaled = TRUE))), 2), "\n")

# ------------------------------------------------------------------------------
# Aggregated survival with population random intercept
# ------------------------------------------------------------------------------
hdr("AGGREGATED SURVIVAL: weighted binomial, population random intercept")
gv <- c("pop", "new_class", "treatment", "myco_cat", "soil")
psurv <- ds %>% group_by(across(all_of(gv))) %>%
  summarise(n_total = n(), survival_rate = sum(alive_01_10_2021, na.rm = TRUE) / n(), .groups = "drop") %>% filter(n_total >= 3)
pcol <- ds %>% filter(alive_01_10_2021 == 1) %>% group_by(across(all_of(gv))) %>%
  summarise(n_surv = n(), n_emf = sum(!is.na(emf)), n_amf = sum(!is.na(amf)), n_dse = sum(!is.na(dse)),
            mean_emf = mean(emf, na.rm = TRUE), mean_amf = mean(amf, na.rm = TRUE), mean_dse = mean(dse, na.rm = TRUE),
            mean_per_f = mean(per_f, na.rm = TRUE), mean_drc_mm = mean(og_drc_mm, na.rm = TRUE), .groups = "drop")
sdf <- inner_join(psurv, pcol, by = gv) %>% mutate(mean_per_fc = mean_per_f - 50)
use <- sdf %>% filter(is.finite(mean_emf), is.finite(mean_dse), is.finite(mean_amf))
cat("cells >=3 trees:", nrow(psurv), "| with a survivor:", nrow(sdf), "| usable (all three fungal means):", nrow(use), "\n")
cat("Group size (trees):\n"); print(table(use$n_total))
cat("\nGroups whose colonization mean rests on ONE scored tree (of", nrow(use), "groups):\n")
print(tibble(EMF = sum(use$n_emf == 1), AMF = sum(use$n_amf == 1), DSE = sum(use$n_dse == 1),
             any_fungus_one_tree = sum(use$n_emf == 1 | use$n_amf == 1 | use$n_dse == 1),
             all_three_one_tree = sum(use$n_emf == 1 & use$n_amf == 1 & use$n_dse == 1),
             `AMF/DSE one tree and EMF one tree` = sum(use$n_amf == 1 & use$n_emf == 1)))
cat("Scored trees per group: EMF median", median(use$n_emf), "(range", paste(range(use$n_emf), collapse = "-"),
    "); AMF median", median(use$n_amf), "(range", paste(range(use$n_amf), collapse = "-"), ")\n")
agg_pf4 <- glmer(survival_rate ~ (scale(mean_emf) + scale(mean_dse) + scale(mean_amf)) * treatment +
                   mean_per_fc + I(mean_per_fc^2) + mean_drc_mm + soil + myco_cat + (1 | pop),
                 data = sdf, weights = n_total, family = binomial, control = ctl)
agg_nv4 <- glmer(survival_rate ~ (scale(mean_emf) + scale(mean_dse) + scale(mean_amf)) * treatment +
                   mean_drc_mm + soil + myco_cat + (1 | pop), data = sdf, weights = n_total, family = binomial, control = ctl)
reg("aggregated survival ancestry-controlled (1|pop)", agg_pf4); reg("aggregated survival naive (1|pop)", agg_nv4)
cat("\nn groups in model:", nobs(agg_pf4), "\n")
cat("Population random intercept: converged =", length(agg_pf4@optinfo$conv$lme4$messages) == 0 || all(grepl("singular", agg_pf4@optinfo$conv$lme4$messages)),
    "| singular (pop variance at 0) =", isSingular(agg_pf4), "\n")
cat("\n--- ancestry-controlled ---\n"); print(summary(agg_pf4)); cat("\n--- naive ---\n"); print(summary(agg_nv4))
agg_glm_pf <- glmmTMB(survival_rate ~ (scale(mean_emf) + scale(mean_dse) + scale(mean_amf)) * treatment +
                        mean_per_fc + I(mean_per_fc^2) + mean_drc_mm + soil + myco_cat, data = sdf, family = binomial(), weights = n_total)
agg_glm_stats <- glm(formula(agg_glm_pf), data = sdf, weights = n_total, family = binomial)
ct <- coef(summary(agg_glm_pf))$cond["scale(mean_amf):treatmentDrought", ]
cat("\nSame AMF x drought term from three engines (identical coefficient, population variance is 0):\n")
print(data.frame(engine = c("glmmTMB, no RE", "stats::glm, no RE", "glmer + (1|pop)"),
                 estimate = c(ct["Estimate"], coef(agg_glm_stats)["scale(mean_amf):treatmentDrought"], fixef(agg_pf4)["scale(mean_amf):treatmentDrought"]),
                 se = c(ct["Std. Error"], coef(summary(agg_glm_stats))["scale(mean_amf):treatmentDrought", 2], coef(summary(agg_pf4))["scale(mean_amf):treatmentDrought", 2]),
                 p = c(ct["Pr(>|z|)"], coef(summary(agg_glm_stats))["scale(mean_amf):treatmentDrought", 4], coef(summary(agg_pf4))["scale(mean_amf):treatmentDrought", 4])), digits = 4)

# ------------------------------------------------------------------------------
# Stomatal conductance
# ------------------------------------------------------------------------------
hdr("STOMATAL CONDUCTANCE: model as described in the Methods + alternatives")
re_p <- "(1 | pop) + (1 | pop:geno_id)"
fp <- function(rhs) as.formula(paste("por_avg ~", rhs, "+ og_drc_mm + soil + myco_cat +", re_p))
poro_M <- lmerTest::lmer(fp("per_fc * treatment + I(per_fc^2) * treatment"), data = ds)   # Methods model
poro_A <- lmerTest::lmer(fp("per_fc * treatment + I(per_fc^2)"), data = ds)              # without the inoculum covariate
poro_R <- lmerTest::lmer(fp("per_fc + treatment + I(per_fc^2)"), data = ds)              # additive (reduced)
poro_L <- lmerTest::lmer(fp("per_fc * treatment"), data = ds)                            # linear ancestry x drought only
poro_LR <- lmerTest::lmer(fp("per_fc + treatment"), data = ds)                           # linear, additive
for (nm in c("M", "A", "R", "L", "LR")) reg(paste("stomatal", nm), get(paste0("poro_", nm)))
cat("n =", nobs(poro_M), "\n")
cat("\n--- M: Methods model (linear + quadratic ancestry x drought) ---\n"); print(summary(poro_M))
cat("\n--- R: additive (no interactions) ---\n"); print(summary(poro_R))
cat("\n--- L: linear ancestry x drought, no quadratic ---\n"); print(summary(poro_L))
cat("\n--- Model comparison (ML refits) ---\n")
print(anova(poro_LR, poro_L, poro_R, poro_A, poro_M, refit = TRUE))
cat("\nScaled residual extremes: M", round(max(abs(resid(poro_M, scaled = TRUE))), 2), "| L", round(max(abs(resid(poro_L, scaled = TRUE))), 2), "\n")
cat("\n--- Key terms: A (common quadratic, linear ancestry x drought), R (additive), L (linear, x drought) ---\n")
for (nm in c("A", "R", "L")) {
  cf <- coef(summary(get(paste0("poro_", nm)))); keep <- intersect(c("per_fc", "I(per_fc^2)", "treatmentDrought", "per_fc:treatmentDrought"), rownames(cf))
  cat("Model", nm, "\n"); print(round(cf[keep, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], 5))
}
poro_Alog <- lmerTest::lmer(update(formula(poro_A), log(por_avg) ~ .), data = ds)
cat("Model A on log(por_avg):\n"); print(round(coef(summary(poro_Alog))[c("per_fc", "I(per_fc^2)", "treatmentDrought", "per_fc:treatmentDrought"), c("Estimate", "t value", "Pr(>|t|)")], 5))
cat("Slope of conductance on ancestry (per 1% Fremont) by treatment, model A (at the midpoint, per_fc = 0):\n")
cfA <- fixef(poro_A); cat("  well-watered:", signif(cfA["per_fc"], 3), "| drought:", signif(cfA["per_fc"] + cfA["per_fc:treatmentDrought"], 3), "\n")
poro_Mlog <- lmerTest::lmer(update(formula(poro_M), log(por_avg) ~ .), data = ds)
poro_Llog <- lmerTest::lmer(update(formula(poro_L), log(por_avg) ~ .), data = ds)
cat("\n--- log(por_avg) diagnostics: Methods model terms ---\n"); print(round(coef(summary(poro_Mlog))[c("per_fc", "I(per_fc^2)", "treatmentDrought", "per_fc:treatmentDrought", "treatmentDrought:I(per_fc^2)"), c("Estimate", "t value", "Pr(>|t|)")], 5))
cat("log(por_avg), linear model terms:\n"); print(round(coef(summary(poro_Llog))[c("per_fc", "treatmentDrought", "per_fc:treatmentDrought"), c("Estimate", "t value", "Pr(>|t|)")], 5))
cat("Distinct genotypes with porometer data:", n_distinct(ds$geno_id[!is.na(ds$por_avg)]), "\n")

# ------------------------------------------------------------------------------
# Supp Fig 2: mixed-model replacement for the long-format ANOVA
# ------------------------------------------------------------------------------
hdr("APPENDIX S2: inoculum and soil effects on colonization (mixed model, Type II Satterthwaite F)")
dall <- ds %>% select(tree_id, pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>%
  pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>%
  filter(!is.na(colonization)) %>% na.omit() %>% mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm)))
cat("rows:", nrow(dall), "| trees:", n_distinct(dall$tree_id), "| genotypes:", n_distinct(dall$geno_id), "\n")
sf2 <- lmerTest::lmer(colonization ~ myco_cat * fungal_type + soil + treatment + per_f + og_drc_mm_z +
                        (1 | pop) + (1 | pop:geno_id) + (1 | tree_id), data = dall)
reg("Supp Fig 2 mixed model", sf2)
print(VarCorr(sf2)); sf2_a <- anova(sf2, type = 2); print(sf2_a)

# ------------------------------------------------------------------------------
# Allocation (exploratory; not reported in the manuscript)
# ------------------------------------------------------------------------------
hdr("ALLOCATION (exploratory): leaf mass fraction with genotype")
ds <- ds %>% mutate(lmf_A = dryleaf_weight / (dryleaf_weight + dry_stem), lmf_B = ifelse(dry_stem > 0 & dryleaf_weight <= dry_stem, dryleaf_weight / dry_stem, NA))
for (y in c("lmf_A", "lmf_B")) {
  d <- ds %>% filter(!is.na(.data[[y]]), !is.na(og_drc_mm))
  mf <- lmerTest::lmer(as.formula(paste(y, "~ per_fc * treatment + I(per_fc^2) * treatment + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id)")), data = d)
  mr <- lmerTest::lmer(as.formula(paste(y, "~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id)")), data = d)
  reg(paste("allocation", y, "full"), mf); l <- anova(mr, mf, refit = TRUE)
  cat("\n---", y, "n =", nobs(mf), "| LRT ancestry x drought p =", round(l$`Pr(>Chisq)`[2], 3), "\n")
  print(round(coef(summary(mr))[c("per_fc", "I(per_fc^2)", "treatmentDrought"), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], 5))
}

# ------------------------------------------------------------------------------
# Convergence summary
# ------------------------------------------------------------------------------
conv <- bind_rows(conv_log)
write.csv(conv, "results/02_main_models/convergence.csv", row.names = FALSE)
hdr("CONVERGENCE / SINGULARITY SUMMARY (saved to results/02_main_models/convergence.csv)")
print(as.data.frame(conv %>% mutate(messages = substr(messages, 1, 70), genotype_sd = round(genotype_sd, 3), pop_sd = round(pop_sd, 3))), right = FALSE)

cat("\n--- Session ---\n"); print(sessionInfo()$R.version$version.string)
sink()
cat("Done. See results/02_main_models/main_models_output.txt, results/02_main_models/genotype_map.csv, results/02_main_models/convergence.csv\n")
