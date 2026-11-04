# ==============================================================================
# MAIN MODELS, FOLLOW-UP: genotype-level PERMANOVA, log-scale stomatal model,
# and paste-ready result strings. Sources 02a_main_models.R first (about 90 s) and
# reuses its fitted objects.
# ==============================================================================
# Run from the project root:  Rscript r_code/02b_main_models_followup.R
# Writes: results/02_main_models/main_models_followup_output.txt, results/02_main_models/formatted_results_main_models.md
# Does not modify the manuscript or any existing script.
#
# Why a custom permutation engine: permute::how(plots = genotype) requires every
# genotype to have the same number of trees, and adonis2 stops with "Design must
# be balanced if permuting 'strata'". Genotypes here have 1-4 trees. The engine
# below implements the same idea directly:
#   * genotype-level terms (ancestry, population): whole genotypes are permuted
#     (ancestry among genotypes WITHIN population; population labels among all
#     genotypes), trees stay attached to their genotype;
#   * tree-level terms (drought, soil): trees are permuted WITHIN genotype;
#   * statistic: order-independent (marginal) pseudo-F and R2 from the
#     multivariate linear model on the raw AMF/EMF/DSE matrix, which equals the
#     Euclidean PERMANOVA. Inoculum strata are not used (genotypes contain trees
#     of several inoculum levels).
# ==============================================================================

source("r_code/02a_main_models.R")      # objects: ds, permA, permA_m, permB_m, col4, bio_*, surv_*, agg_*, fb_*, poro_*, sf2_a, ...
library(permute)
sink("results/02_main_models/main_models_followup_output.txt", split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("MAIN MODELS FOLLOW-UP run:", "\n")
BETA <- "β"; TIMES <- "×"; CHI2 <- "χ²"; R2 <- "R²"
set.seed(1998)

# ------------------------------------------------------------------------------
# 1. Colonization subset: genotypes per crosstype
# ------------------------------------------------------------------------------
hdr("1. Colonization subset (the 78 trees): trees and genotypes per crosstype")
d78 <- ds %>% filter(!is.na(per_f), !is.na(pop), !is.na(treatment), !is.na(amf), !is.na(emf), !is.na(dse),
                     !is.na(soil), !is.na(myco_cat)) %>%
  mutate(pop = factor(pop), treatment = factor(treatment), soil = factor(soil))
stopifnot(nrow(d78) == 78)
gcount <- d78 %>% group_by(new_class) %>% summarise(trees = n(), genotypes = n_distinct(geno_id), .groups = "drop")
print(as.data.frame(gcount)); cat("Total: ", nrow(d78), "trees,", n_distinct(d78$geno_id), "genotypes. Trees per genotype:\n"); print(table(table(d78$geno_id)))
foc <- d78 %>% filter(new_class %in% c("Narrowleaf", "F1 Hybrid", "Fremont")) %>% mutate(new_class = droplevels(new_class))
cat("Focal three crosstypes:", nrow(foc), "trees,", n_distinct(foc$geno_id), "genotypes\n")

# ------------------------------------------------------------------------------
# 2. (a) whole-genotype permutation test
# ------------------------------------------------------------------------------
hdr("2a. Genotype-level PERMANOVA (whole genotypes permuted; trees within genotype fixed together)")
cat("permute::how(plots = genotype) check:\n")
chk_how <- tryCatch(adonis2(foc[, c("amf", "emf", "dse")] ~ new_class, data = foc, method = "euclidean",
                            permutations = how(plots = Plots(strata = factor(foc$geno_id), type = "free"),
                                               within = Within(type = "none"), nperm = 99)),
                    error = function(e) conditionMessage(e))
cat("  adonis2 with plots = genotype ->", if (is.character(chk_how)) chk_how else "ran", "\n")
cat("  (genotypes have 1-4 trees each; the custom engine below implements the same permutation)\n")

perm_engine <- function(dat, anc, B = 4999) {
  dat <- droplevels(dat)
  Y <- as.matrix(dat[, c("amf", "emf", "dse")]); Yc <- scale(Y, scale = FALSE); tss <- sum(Yc^2)
  fm <- as.formula(paste("~", anc, "+ pop + treatment + soil"))
  X0 <- model.matrix(fm, dat); asg <- attr(X0, "assign"); terms <- c(anc, "pop", "treatment", "soil")
  statk <- function(X, k) {
    rf <- sum(qr.resid(qr(X), Yc)^2); rr <- sum(qr.resid(qr(X[, asg != k, drop = FALSE]), Yc)^2)
    c(F = ((rr - rf) / sum(asg == k)) / (rf / (nrow(X) - ncol(X))), R2 = (rr - rf) / tss)
  }
  gt <- unique(dat[, c("geno_id", "pop", anc)])
  mapg <- function(d, g, v) { d[[v]] <- g[[v]][match(d$geno_id, g$geno_id)]; d }
  perm_geno <- function(d, k) {
    g <- gt
    if (k == 1) { for (p in unique(g$pop)) { i <- which(g$pop == p); g[[anc]][i] <- g[[anc]][i[sample.int(length(i))]] } }   # within population
    else g$pop <- g$pop[sample.int(nrow(g))]                                                                                   # among all genotypes
    mapg(d, g, terms[k])
  }
  perm_within <- function(d, v) { lv <- levels(d[[v]]); d[[v]] <- factor(ave(as.character(d[[v]]), d$geno_id, FUN = function(x) x[sample.int(length(x))]), levels = lv); d }
  perm_free <- function(d, v) { d[[v]] <- d[[v]][sample.int(nrow(d))]; d }
  out <- lapply(seq_along(terms), function(k) {
    obs <- statk(X0, k); v <- terms[k]
    fg <- replicate(B, { d <- if (k <= 2) perm_geno(dat, k) else perm_within(dat, v); statk(model.matrix(fm, d), k)["F"] })
    ft <- replicate(B, statk(model.matrix(fm, perm_free(dat, v)), k)["F"])
    data.frame(term = v, df = sum(asg == k), R2 = unname(obs["R2"]), F = unname(obs["F"]),
               p_genotype_level = (1 + sum(fg >= obs["F"] - 1e-9)) / (B + 1),
               p_trees_free = (1 + sum(ft >= obs["F"] - 1e-9)) / (B + 1))
  })
  bind_rows(out)
}
cat("\nSanity check: engine marginal R2 / F vs adonis2 by = 'margin'\n")
eng_foc <- perm_engine(foc, "new_class"); eng_all <- perm_engine(d78, "per_f")
cmp <- rbind(data.frame(model = "cross type", term = eng_foc$term, engine_R2 = eng_foc$R2, engine_F = eng_foc$F,
                        adonis_R2 = as.data.frame(permA_m)[eng_foc$term, "R2"], adonis_F = as.data.frame(permA_m)[eng_foc$term, "F"]),
             data.frame(model = "continuous", term = eng_all$term, engine_R2 = eng_all$R2, engine_F = eng_all$F,
                        adonis_R2 = as.data.frame(permB_m)[eng_all$term, "R2"], adonis_F = as.data.frame(permB_m)[eng_all$term, "F"]))
print(cmp, digits = 4, row.names = FALSE)
cat("\n--- Cross type version: marginal terms, n =", nrow(foc), "trees,", n_distinct(foc$geno_id), "genotypes (Narrowleaf, F1 Hybrid, Fremont) ---\n")
print(eng_foc, digits = 4, row.names = FALSE)
cat("\n--- Continuous ancestry version: marginal terms, n =", nrow(d78), "trees,", n_distinct(d78$geno_id), "genotypes ---\n")
print(eng_all, digits = 4, row.names = FALSE)
cat("p_genotype_level: whole genotypes permuted (ancestry within population; population among genotypes); drought and soil permuted within genotype.\n",
    "p_trees_free: same statistic, trees permuted freely (no genotype structure, no inoculum strata) for comparison.\n", sep = "")

# ------------------------------------------------------------------------------
# 2b. PERMANOVA on genotype means
# ------------------------------------------------------------------------------
hdr("2b. PERMANOVA on genotype means of AMF, EMF, DSE (genotypes exchangeable; 9999 permutations)")
gm <- d78 %>% group_by(geno_id, pop, new_class, per_f) %>%
  summarise(n_trees = n(), amf = mean(amf), emf = mean(emf), dse = mean(dse), .groups = "drop")
gmf <- gm %>% filter(new_class %in% c("Narrowleaf", "F1 Hybrid", "Fremont")) %>% mutate(new_class = droplevels(new_class))
cat("genotype means:", nrow(gm), "genotypes (all), ", nrow(gmf), "genotypes (focal); trees per genotype 1-", max(gm$n_trees), "\n")
print(as.data.frame(gmf %>% group_by(new_class) %>% summarise(genotypes = n(), trees = sum(n_trees), .groups = "drop")))
set.seed(1998)
gm_run <- function(dat, rhs, label) {
  cat("\n---", label, "---\n")
  m <- as.matrix(dat[, c("amf", "emf", "dse")])
  f <- as.formula(paste("m ~", rhs))
  s <- adonis2(f, data = dat, method = "euclidean", permutations = 9999, by = "terms")
  g <- adonis2(f, data = dat, method = "euclidean", permutations = 9999, by = "margin")
  cat("sequential:\n"); print(s); cat("marginal:\n"); print(g)
  list(seq = as.data.frame(s), mar = as.data.frame(g))
}
gb_cls1 <- gm_run(gmf, "new_class", "genotype means, cross type only")
gb_cls2 <- gm_run(gmf, "new_class + pop", "genotype means, cross type + population")
gb_pf1 <- gm_run(gm, "per_f", "genotype means, continuous ancestry only")
gb_pf2 <- gm_run(gm, "per_f + pop", "genotype means, continuous ancestry + population")
mm <- as.matrix(gmf[, c("amf", "emf", "dse")])
gm_pw <- pairwise.adonis(mm, gmf$new_class, sim.method = "euclidean", p.adjust.m = "bonferroni")
cat("\nPairwise crosstype on genotype means (Bonferroni):\n"); print(gm_pw)

# ------------------------------------------------------------------------------
# 3. Stomatal conductance, log scale
# ------------------------------------------------------------------------------
hdr("3. log(por_avg) ~ linear + quadratic ancestry + drought + linear ancestry x drought + covariates")
poro_log <- lmerTest::lmer(log(por_avg) ~ per_fc * treatment + I(per_fc^2) + og_drc_mm + soil + myco_cat +
                             (1 | pop) + (1 | pop:geno_id), data = ds)
print(summary(poro_log))
cat("n =", nobs(poro_log), "| genotypes:", ngrps(poro_log)["pop:geno_id"], "| singular:", isSingular(poro_log), "\n")
bad <- which(ds$por_avg <= 0)
cat("por_avg non-NA:", sum(!is.na(ds$por_avg)), "| non-positive values (cannot be logged):", length(bad),
    "| tree(s):", paste(ds$tree_idold[bad], collapse = ","), "| value(s):", paste(round(ds$por_avg[bad], 4), collapse = ","), "\n")
b <- fixef(poro_log); V <- as.matrix(vcov(poro_log))
s_ww <- unname(b["per_fc"]); se_ww <- unname(sqrt(V["per_fc", "per_fc"]))
s_d <- unname(b["per_fc"] + b["per_fc:treatmentDrought"])
se_d <- unname(sqrt(V["per_fc", "per_fc"] + V["per_fc:treatmentDrought", "per_fc:treatmentDrought"] + 2 * V["per_fc", "per_fc:treatmentDrought"]))
bt <- function(s, se) c(pct_per_10 = 100 * (exp(10 * s) - 1), lo = 100 * (exp(10 * (s - 1.96 * se)) - 1), hi = 100 * (exp(10 * (s + 1.96 * se)) - 1),
                        ratio_100_vs_0 = exp(100 * s), r_lo = exp(100 * (s - 1.96 * se)), r_hi = exp(100 * (s + 1.96 * se)))
slope_tab <- rbind(`well-watered` = c(slope_log = s_ww, se = se_ww, bt(s_ww, se_ww)), drought = c(slope_log = s_d, se = se_d, bt(s_d, se_d)))
cat("\nAncestry slope on the log scale per 1% Fremont ancestry, back-transformed. Because the quadratic is symmetric about 50%,\n",
    "the linear coefficient is also the average slope from 0% to 100% (Fremont end vs narrowleaf end ratio = exp(100 x slope)).\n", sep = "")
print(round(slope_tab, 5))
em <- summary(emmeans(poro_log, ~ per_fc | treatment, at = list(per_fc = c(-50, 0, 50)), type = "response", lmer.df = "asymptotic"))
cat("\nPredicted gsw (back-transformed, covariates at means / equal-weighted levels):\n"); print(em)

# ------------------------------------------------------------------------------
# 4-5 helper values
# ------------------------------------------------------------------------------
hdr("4. Survival: LRT for the two ancestry x drought terms (glmer)")
surv_lrt <- anova(surv_red4, surv_full4)
print(surv_lrt)
hdr("5. Aggregated survival with stats::glm (ancestry-controlled and naive)")
agg_glm_naive <- glm(survival_rate ~ (scale(mean_emf) + scale(mean_dse) + scale(mean_amf)) * treatment + mean_drc_mm + soil + myco_cat,
                     data = sdf, weights = n_total, family = binomial)
print(summary(agg_glm_stats)); print(summary(agg_glm_naive))
cat("\nGroups (of", nrow(use), ") whose colonization mean rests on a single scored tree: EMF", sum(use$n_emf == 1), ", AMF", sum(use$n_amf == 1),
    ", DSE", sum(use$n_dse == 1), "\n")
hdr("Fig 5 caption check: refit without the tree with |z| > 3")
fb_pf4_x <- lmerTest::lmer(formula(fb_pf4), data = ds %>% filter(tree_idold != zt$tree_idold[1]))
cat("excluded tree:", zt$tree_idold[1], "(z =", round(zt$z[1], 2), "); n =", nobs(fb_pf4_x), "\n")
print(round(coef(summary(fb_pf4_x))[c("treatmentDrought:scale(emf)", "scale(amf):treatmentDrought", "treatmentDrought:scale(dse)"), c("Estimate", "t value", "Pr(>|t|)")], 4))

# ------------------------------------------------------------------------------
# 6. Paste-ready strings
# ------------------------------------------------------------------------------
hdr("6. Writing results/02_main_models/formatted_results_main_models.md")
fp <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
fb <- function(x, d = 3) formatC(signif(x, d), format = "fg")
fz <- function(x) sprintf("%.2f", x)
fr <- function(x) sprintf("%.1f", x)
cmx <- function(m, t) { cs <- coef(summary(m)); if (is.list(cs)) cs <- cs$cond; cs[t, ] }
bet <- function(m, t, tc = "t") { r <- cmx(m, t); sprintf("%s = %s, %s = %s, %s", BETA, fb(r["Estimate"]), tc, fz(r[grep("value", names(r))[1]]), fp(r[grep("^Pr", names(r))[1]])) }
pv <- function(m, t) unname(cmx(m, t)[grep("^Pr", colnames(coef(summary(m))))[1]])
ssr <- tidy_surv(surv_red4); ssf <- tidy_surv(surv_full4)
rows <- list()
add <- function(section, loc, stat, string, ms_p = NA, p = NA, note = "") {
  flag <- if (!is.na(ms_p) && !is.na(p) && ((ms_p < 0.05) != (p < 0.05))) paste0("**SIGNIFICANCE CHANGED** (manuscript p ", ifelse(ms_p < 0.001, "< 0.001", ms_p), " -> ", fp(p), ")") else ""
  rows[[length(rows) + 1]] <<- data.frame(section = section, location = loc, statistic = stat, string = string, flag = flag, note = note)
}
S1 <- "Results: Biomass and Survival Vary with Host Ancestry"
add(S1, "P1 / Fig 1 caption", "Survival ancestry x drought, full model",
    sprintf("linear and quadratic ancestry %s drought terms, p = %.2f and p = %.2f", TIMES, ssf["per_fc:treatmentDrought", "Pr(>|z|)"], ssf["treatmentDrought:I(per_fc^2)", "Pr(>|z|)"]),
    0.89, ssf["per_fc:treatmentDrought", "Pr(>|z|)"], "glmer, genotype random intercept; both terms non-significant (manuscript 0.89 and 0.47)")
add(S1, "P1 / Fig 2 caption", "Biomass LRT ancestry x drought",
    sprintf("likelihood ratio test [LRT] for ancestry %s drought treatment, %s = %s, df = %d, p = %.2f", TIMES, CHI2, fz(bio_lrt4$Chisq[2]), bio_lrt4$Df[2], bio_lrt4$`Pr(>Chisq)`[2]),
    0.58, bio_lrt4$`Pr(>Chisq)`[2], "full vs reduced model")
add(S1, "P1 / Fig 1 caption", "Survival joint LRT for both interaction terms (new)",
    sprintf("%s = %s, df = %d, p = %.2f", CHI2, fz(surv_lrt$Chisq[2]), surv_lrt$Df[2], surv_lrt$`Pr(>Chisq)`[2]), NA, NA, "optional; mirrors the biomass LRT")
add(S1, "P1 / Fig 1 caption", "Survival linear ancestry (reduced)", sprintf("%s = %s, z = %s, %s", BETA, fb(ssr["per_fc", "Estimate"]), fz(ssr["per_fc", "z value"]), fp(ssr["per_fc", "Pr(>|z|)"])),
    0.028, ssr["per_fc", "Pr(>|z|)"], "manuscript's -2.21 / 0.028 were the full model's well-watered term; reduced-model values shown")
add(S1, "P1", "Survival quadratic ancestry (reduced)", sprintf("%s = %s, z = %s, %s", BETA, fb(ssr["I(per_fc^2)", "Estimate"], 2), fz(ssr["I(per_fc^2)", "z value"]), fp(ssr["I(per_fc^2)", "Pr(>|z|)"])),
    0.40, ssr["I(per_fc^2)", "Pr(>|z|)"])
add(S1, "P1 / Fig 1 caption", "Survival drought (reduced)", sprintf("%s = %s, z = %s, %s", BETA, fb(ssr["treatmentDrought", "Estimate"], 4), fz(ssr["treatmentDrought", "z value"]), fp(ssr["treatmentDrought", "Pr(>|z|)"])),
    0.0005, ssr["treatmentDrought", "Pr(>|z|)"])
add(S1, "P1 / Fig 2 caption", "Biomass linear ancestry (reduced)", paste0("linear term: ", bet(bio_red4, "per_fc")), 0.0005, pv(bio_red4, "per_fc"))
add(S1, "P1 / Fig 2 caption", "Biomass quadratic ancestry (reduced)", paste0("quadratic term: ", bet(bio_red4, "I(per_fc^2)")), 0.0005, pv(bio_red4, "I(per_fc^2)"),
    "still significant, but p is no longer < 0.001")
add(S1, "P1 / Fig 2 caption", "Biomass drought (reduced)", bet(bio_red4, "treatmentDrought"), 0.045, pv(bio_red4, "treatmentDrought"),
    "wording such as 'drought reduced biomass gain overall' (Abstract, Results, Fig 2 caption, Discussion, Conclusions) is no longer supported at p < 0.05")
add(S1, "Abstract / P1 / Fig 2 / Conclusions", "Biomass peak (31%)", "Not reported (CI extends outside 0-100%). Vertex = 31.8%; bootstrap 95% CI -32.2 to 42.4 (results/02_main_models/main_models_output.txt)", NA, NA, "the CI includes values outside 0-100%")

S2 <- "Results: Host Ancestry, Not Drought, Structures Fungal Colonization (PERMANOVA)"
Gc <- function(t) eng_foc[eng_foc$term == t, ]; Ga <- function(t) eng_all[eng_all$term == t, ]
Tc <- function(t) as.data.frame(permA)[t, ]; Ta <- function(t) as.data.frame(permB)[t, ]
add(S2, "P2", "(G) Cross type, whole-genotype permutations",
    sprintf("Host cross type explained %s%% of variance in fungal composition and was the only significant predictor (PERMANOVA with whole-genotype permutations, %d genotypes: F = %s, %s)",
            fr(100 * Gc("new_class")$R2), n_distinct(foc$geno_id), fz(Gc("new_class")$F), fp(Gc("new_class")$p_genotype_level)),
    0.003, Gc("new_class")$p_genotype_level, "marginal R2 / F; population, drought and soil in the same model")
add(S2, "P2", "(T) Cross type, tree-level permutations (model A)",
    sprintf("%s = %s%%, F = %s, %s", R2, fr(100 * Tc("new_class")$R2), fz(Tc("new_class")$F), fp(Tc("new_class")$`Pr(>F)`)), 0.003, Tc("new_class")$`Pr(>F)`,
    "permutations within inoculum category, trees exchangeable")
add(S2, "P2", "(G) Continuous ancestry, whole-genotype permutations",
    sprintf("A parallel model using continuous percent Fremont ancestry across the full sample, including backcrosses (n = %d trees, %d genotypes), gave a consistent result (%s = %s%%, F = %s, %s)",
            nrow(d78), n_distinct(d78$geno_id), R2, fr(100 * Ga("per_f")$R2), fz(Ga("per_f")$F), fp(Ga("per_f")$p_genotype_level)),
    0.001, Ga("per_f")$p_genotype_level)
add(S2, "P2", "(T) Continuous ancestry, tree-level (model B)", sprintf("%s = %s%%, F = %s, %s", R2, fr(100 * Ta("per_f")$R2), fz(Ta("per_f")$F), fp(Ta("per_f")$`Pr(>F)`)), 0.001, Ta("per_f")$`Pr(>F)`)
add(S2, "P2", "(G) Population, same model as cross type",
    sprintf("Population explained %s%% of variance (%s)", fr(100 * Gc("pop")$R2), fp(Gc("pop")$p_genotype_level)), 0.352, Gc("pop")$p_genotype_level, "population is now tested in the same model as cross type")
add(S2, "P2", "(T) Population, tree-level", sprintf("%s = %s%%, %s", R2, fr(100 * Tc("pop")$R2), fp(Tc("pop")$`Pr(>F)`)), 0.352, Tc("pop")$`Pr(>F)`)
add(S2, "P2 / Supp Fig 3", "(G) Drought, permuted within genotype",
    sprintf("%s = %s%%, F = %s, %s", R2, fr(100 * Gc("treatment")$R2), fz(Gc("treatment")$F), fp(Gc("treatment")$p_genotype_level)), 0.79, Gc("treatment")$p_genotype_level)
add(S2, "P2", "(G) Soil legacy, permuted within genotype",
    sprintf("%s = %s%%, F = %s, %s", R2, fr(100 * Gc("soil")$R2), fz(Gc("soil")$F), fp(Gc("soil")$p_genotype_level)), 0.27, Gc("soil")$p_genotype_level)
add(S2, "P2", "(G) Continuous version: soil legacy, permuted within genotype",
    sprintf("%s = %s%%, F = %s, %s (tree-level permutation: %s)", R2, fr(100 * Ga("soil")$R2), fz(Ga("soil")$F), fp(Ga("soil")$p_genotype_level), fp(Ga("soil")$p_trees_free)),
    0.27, Ga("soil")$p_genotype_level, "NOT ROBUST: P < 0.05 only in some specifications (see results/11_soil_legacy_permanova/output.txt and results/16_soil_legacy_specifications/soil_specifications.csv); treat as exploratory")
add(S2, "P2", "(G) Continuous version: drought, permuted within genotype",
    sprintf("%s = %s%%, F = %s, %s", R2, fr(100 * Ga("treatment")$R2), fz(Ga("treatment")$F), fp(Ga("treatment")$p_genotype_level)), 0.79, Ga("treatment")$p_genotype_level)
add(S2, "P2 / Supp Fig 3 caption", "(T) Drought and soil, tree-level (model A)",
    sprintf("drought: %s = %s%%, F = %s, %s; soil: %s = %s%%, F = %s, %s", R2, fr(100 * Tc("treatment")$R2), fz(Tc("treatment")$F), fp(Tc("treatment")$`Pr(>F)`),
            R2, fr(100 * Tc("soil")$R2), fz(Tc("soil")$F), fp(Tc("soil")$`Pr(>F)`)), 0.79, Tc("treatment")$`Pr(>F)`, "Supp Fig 3 caption still says p = 0.784")
gmr <- function(res, term) res$seq[term, ]
add(S2, "P2 (alternative)", "(M) Genotype means, cross type only",
    sprintf("PERMANOVA on genotype means (%d genotypes): %s = %s%%, F = %s, %s", nrow(gmf), R2, fr(100 * gmr(gb_cls1, "new_class")$R2), fz(gmr(gb_cls1, "new_class")$F), fp(gmr(gb_cls1, "new_class")$`Pr(>F)`)),
    0.003, gmr(gb_cls1, "new_class")$`Pr(>F)`)
add(S2, "P2 (alternative)", "(M) Genotype means, cross type + population",
    sprintf("cross type: %s = %s%%, F = %s, %s; population: %s = %s%%, %s", R2, fr(100 * gmr(gb_cls2, "new_class")$R2), fz(gmr(gb_cls2, "new_class")$F), fp(gmr(gb_cls2, "new_class")$`Pr(>F)`),
            R2, fr(100 * gmr(gb_cls2, "pop")$R2), fp(gmr(gb_cls2, "pop")$`Pr(>F)`)), 0.003, gmr(gb_cls2, "new_class")$`Pr(>F)`, "sequential; marginal tests in results/02_main_models/main_models_followup_output.txt")
add(S2, "P2 (alternative)", "(M) Genotype means, continuous ancestry",
    sprintf("continuous ancestry (%d genotypes): %s = %s%%, F = %s, %s", nrow(gm), R2, fr(100 * gmr(gb_pf1, "per_f")$R2), fz(gmr(gb_pf1, "per_f")$F), fp(gmr(gb_pf1, "per_f")$`Pr(>F)`)),
    0.001, gmr(gb_pf1, "per_f")$`Pr(>F)`)
pwp <- gm_pw$p.adjusted
add(S2, "P3 pairwise", "(M) Pairwise on genotype means (Bonferroni)",
    sprintf("Fremont vs narrowleaf p = %.3f; Fremont vs F1 hybrids p = %.3f; narrowleaf vs F1 hybrids p = %.2f", pwp[1], pwp[2], pwp[3]), 0.051, pwp[2],
    "genotype counts: 15 narrowleaf, 4 F1, 15 Fremont")
add(S2, "P3 pairwise", "(T) Pairwise on trees (unchanged model)",
    sprintf("Fremont vs narrowleaf p = %.3f; Fremont vs F1 hybrids p = %.3f; narrowleaf vs F1 hybrids p = %.2f", pw$p.adjusted[1], pw$p.adjusted[2], pw$p.adjusted[3]), 0.051, pw$p.adjusted[2],
    "permutation noise; run in a different RNG state than the manuscript's 0.051")
ct4 <- coef(summary(col4)); nl <- summary(col4n)$tTable
add(S2, "P3 / Fig 3 caption", "EMF x ancestry (lmer + tree_id + genotype)",
    sprintf("Fremont ancestry %s EMF interaction, t = %s, %s", TIMES, fz(ct4["per_f:fungal_typeemf", "t value"]), fp(ct4["per_f:fungal_typeemf", "Pr(>|t|)"])), 0.0005, ct4["per_f:fungal_typeemf", "Pr(>|t|)"],
    "tree and genotype variances are 0, so identical to the model without them")
add(S2, "P3 / Methods", "EMF x ancestry (nlme sensitivity)",
    sprintf("t = %s, %s (nlme, unstructured correlation among fungal types within trees, separate variance per fungal type, genotype nested in population)", fz(nl["per_f:fungal_typeemf", "t-value"]), fp(nl["per_f:fungal_typeemf", "p-value"])),
    0.0005, nl["per_f:fungal_typeemf", "p-value"])
sl <- summary(emtrends(col4, ~ fungal_type, var = "per_f", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE))
add(S2, "Fig 3 caption", "Ancestry slope by fungal type (lmer)",
    sprintf("EMF slope = %s per 1%% Fremont ancestry (%s); AMF = %s (%s); DSE = %s (%s)", fb(sl$per_f.trend[3]), fp(sl$p.value[3]), fb(sl$per_f.trend[1]), fp(sl$p.value[1]), fb(sl$per_f.trend[2]), fp(sl$p.value[2])),
    NA, NA, "supports 'AMF and DSEs showed no significant relationship with ancestry'")

S3 <- "Results: Survival associations (aggregated groups; stats::glm, ancestry-controlled)"
gl <- coef(summary(agg_glm_stats)); gn <- coef(summary(agg_glm_naive))
gs <- function(tb, t, tc = "z") sprintf("%s = %s, z = %s, %s", BETA, fb(tb[t, "Estimate"], 3), fz(tb[t, "z value"]), fp(tb[t, "Pr(>|z|)"]))
add(S3, "Survival associations / Fig 4 caption", "AMF x drought (ancestry-controlled)", paste0("The AMF ", TIMES, " drought interaction was significant and negative, controlling for host ancestry (linear and quadratic) among survivors in each group (", gs(gl, "scale(mean_amf):treatmentDrought"), ")"),
    0.018, gl["scale(mean_amf):treatmentDrought", "Pr(>|z|)"])
add(S3, "Survival associations / Fig 4 caption", "AMF x drought (naive)", gs(gn, "scale(mean_amf):treatmentDrought"), 0.010, gn["scale(mean_amf):treatmentDrought", "Pr(>|z|)"])
add(S3, "Survival associations", "EMF and DSE x drought", sprintf("EMF %s drought and DSE %s drought interactions were not significant (%s and %s, respectively)", TIMES, TIMES,
    fp(gl["scale(mean_emf):treatmentDrought", "Pr(>|z|)"]), fp(gl["scale(mean_dse):treatmentDrought", "Pr(>|z|)"])), 0.179, gl["scale(mean_emf):treatmentDrought", "Pr(>|z|)"])
add(S3, "Survival associations", "EMF x drought detail", gs(gl, "scale(mean_emf):treatmentDrought"), 0.179, gl["scale(mean_emf):treatmentDrought", "Pr(>|z|)"])
add(S3, "Survival associations", "DSE x drought detail", gs(gl, "scale(mean_dse):treatmentDrought"), 0.325, gl["scale(mean_dse):treatmentDrought", "Pr(>|z|)"])
add(S3, "Survival associations / Limitations", "Ancestry terms", sprintf("Percent Fremont ancestry (linear or quadratic) did not significantly predict survival in this model when fungal colonization was included (%s and %s, respectively)",
    fp(gl["mean_per_fc", "Pr(>|z|)"]), fp(gl["I(mean_per_fc^2)", "Pr(>|z|)"])), 0.74, gl["mean_per_fc", "Pr(>|z|)"], "Limitations paragraph quotes the linear p")
add(S3, "Methods / Survival associations", "Single-tree groups",
    sprintf("Colonization means rested on a single scored tree in %d (EMF), %d (AMF) and %d (DSE) of %d groups", sum(use$n_emf == 1), sum(use$n_amf == 1), sum(use$n_dse == 1), nrow(use)), NA, NA, "describe as a limitation; groups have 3-6 trees but 1-5 scored trees")
add(S3, "Methods", "Model description", "weighted binomial generalized linear models (glm; a population random intercept was estimated at zero variance and dropped)", NA, NA, "Methods currently say GLMMs")

S4 <- "Results: Biomass associations (n = 77, with genotype)"
fbs <- function(m, t) bet(m, t)
add(S4, "Biomass associations / Fig 5 caption", "EMF x drought", paste0("(LMM, n = ", nobs(fb_pf4), "; EMF ", TIMES, " drought: ", fbs(fb_pf4, "treatmentDrought:scale(emf)"), ")"), 0.004, pv(fb_pf4, "treatmentDrought:scale(emf)"))
add(S4, "Biomass associations / Fig 5 caption", "EMF x drought (naive)", fbs(fb_nv4, "treatmentDrought:scale(emf)"), 0.003, pv(fb_nv4, "treatmentDrought:scale(emf)"))
add(S4, "Biomass associations", "AMF x drought", fbs(fb_pf4, "scale(amf):treatmentDrought"), 0.617, pv(fb_pf4, "scale(amf):treatmentDrought"))
add(S4, "Biomass associations", "DSE x drought", fbs(fb_pf4, "treatmentDrought:scale(dse)"), 0.984, pv(fb_pf4, "treatmentDrought:scale(dse)"))
add(S4, "Fig 5 caption", "Outlier statement", sprintf("One tree (%s; standardized log biomass z = %s) exceeded |z| = 3; excluding it gives EMF %s drought %s", zt$tree_idold[1], fz(zt$z[1]), TIMES, bet(fb_pf4_x, "treatmentDrought:scale(emf)")),
    NA, NA, "current caption says none met the threshold")

S5 <- "Results: Stomatal conductance (log scale; Methods model without quadratic x drought)"
nm <- function(t) bet(poro_log, t)
add(S5, "Methods", "Model", "log-transformed stomatal conductance with linear and quadratic ancestry, drought, and linear ancestry x drought, plus initial cutting diameter, soil legacy and inoculation category as covariates and population and genotype (nested in population) as random intercepts", NA, NA, "quadratic x drought removed (LRT p = 0.96 on the raw scale)")
add(S5, "Stomatal conductance", "n", sprintf("n = %d trees from %d genotypes", nobs(poro_log), ngrps(poro_log)["pop:geno_id"]), NA, NA,
    sprintf("was n = 204; one tree (%s) has a non-positive conductance value (%s) that cannot be log-transformed", paste(ds$tree_idold[bad], collapse = ","), paste(round(ds$por_avg[bad], 4), collapse = ",")))
add(S5, "Stomatal conductance", "Linear ancestry", paste0("linear ancestry term: ", nm("per_fc")), 0.0005, pv(poro_log, "per_fc"), "simple effect in well-watered trees")
add(S5, "Stomatal conductance", "Quadratic ancestry", paste0("quadratic ancestry term: ", nm("I(per_fc^2)")), 0.043, pv(poro_log, "I(per_fc^2)"))
add(S5, "Stomatal conductance", "Drought", paste0("drought (well-watered vs drought at 50% ancestry): ", nm("treatmentDrought")), 0.072, pv(poro_log, "treatmentDrought"))
add(S5, "Stomatal conductance", "Ancestry x drought", paste0("linear ancestry ", TIMES, " drought: ", nm("per_fc:treatmentDrought")), 0.073, pv(poro_log, "per_fc:treatmentDrought"),
    "ancestry slope is shallower under drought; the text currently says the interaction was not significant")
inoc <- cmx(poro_log, grep("myco_cat", rownames(coef(summary(poro_log))), value = TRUE))
add(S5, "Stomatal conductance", "Inoculation category", sprintf("mycorrhizal inoculation category (all p > %.2f)", floor(min(inoc[, "Pr(>|t|)"]) * 100) / 100), NA, NA)
add(S5, "Stomatal conductance", "Covariates", sprintf("initial cutting diameter: %s; soil legacy: %s", bet(poro_log, "og_drc_mm"), bet(poro_log, "soilTam Legacy")), NA, NA)
add(S5, "Stomatal conductance", "Slope, well-watered", sprintf("In well-watered trees, stomatal conductance declined by %s%% (95%% CI %s to %s%%) per 10 percentage points of Fremont ancestry; Fremont-end trees had %s times the conductance of narrowleaf-end trees (95%% CI %s to %s)",
    fr(-slope_tab[1, "pct_per_10"]), fr(-slope_tab[1, "hi"]), fr(-slope_tab[1, "lo"]), sprintf("%.2f", slope_tab[1, "ratio_100_vs_0"]), sprintf("%.2f", slope_tab[1, "r_lo"]), sprintf("%.2f", slope_tab[1, "r_hi"])), NA, NA,
    "back-transformed from the log scale; average slope over 0-100% (quadratic cancels)")
add(S5, "Stomatal conductance", "Slope, drought", sprintf("Under drought the decline was %s%% (95%% CI %s to %s%%) per 10 percentage points; ratio Fremont-end vs narrowleaf-end = %s (95%% CI %s to %s)",
    fr(-slope_tab[2, "pct_per_10"]), fr(-slope_tab[2, "hi"]), fr(-slope_tab[2, "lo"]), sprintf("%.2f", slope_tab[2, "ratio_100_vs_0"]), sprintf("%.2f", slope_tab[2, "r_lo"]), sprintf("%.2f", slope_tab[2, "r_hi"])), NA, NA)

S6 <- "Supplemental figure captions"
fa <- as.data.frame(sf2_a)
add(S6, "Supp Fig 2 caption", "Inoculum, soil", sprintf("mixed model (tree, genotype and population random intercepts): inoculum F(3, %.0f) = %s, %s; soil F(1, %.0f) = %s (p = %.2f)", fa["myco_cat", "DenDF"], fz(fa["myco_cat", "F value"]), fp(fa["myco_cat", "Pr(>F)"]),
    fa["soil", "DenDF"], fz(fa["soil", "F value"]), fa["soil", "Pr(>F)"]), 0.03, fa["myco_cat", "Pr(>F)"],
    "replaces the ANOVA F(3, 328); denominator df are Satterthwaite; soil was p = 0.06, now p = 0.54 (both non-significant)")
add(S6, "Supp Fig 2 caption", "Inoculum x fungal type", sprintf("inoculum %s fungal type F(6, %.0f) = %s, %s", TIMES, fa["myco_cat:fungal_type", "DenDF"], fz(fa["myco_cat:fungal_type", "F value"]), fp(fa["myco_cat:fungal_type", "Pr(>F)"])), NA, NA, "optional")
add(S6, "Supp Fig 3 caption", "Drought, PERMANOVA", sprintf("PERMANOVA: F = %s, %s (whole-genotype permutations; tree-level p = %.2f)", fz(Gc("treatment")$F), fp(Gc("treatment")$p_genotype_level), Tc("treatment")$`Pr(>F)`), 0.784, Gc("treatment")$p_genotype_level)

S7 <- "Sample sizes"
add(S7, "Methods / Results", "n", sprintf("colonization: %d trees from %d genotypes (%d survivors scored for all three fungi); survival models n = %d; biomass n = %d; colonization-biomass n = %d; aggregated survival: %d groups; stomatal n = %d; %d genotypes overall",
    nrow(d78), n_distinct(d78$geno_id), nrow(d78), nobs(surv_red4), nobs(bio_red4), nobs(fb_pf4), nobs(agg_glm_stats), nobs(poro_log), nrow(gmap)), NA, NA)
add(S7, "Methods", "Genotypes per cross type (78 trees)", paste(sprintf("%s: %d %s (%d trees)", gcount$new_class, gcount$genotypes, ifelse(gcount$genotypes == 1, "genotype", "genotypes"), gcount$trees), collapse = "; "), NA, NA)

tab <- bind_rows(rows)
lines <- c("# Formatted results: main models", "",
           "Generated by `r_code/02b_main_models_followup.R`. The manuscript is not changed. Rows marked **SIGNIFICANCE CHANGED** cross p = 0.05 relative to the current manuscript.",
           "Labels: (G) genotype-level permutation test, (T) tree-level permutation test, (M) PERMANOVA on genotype means. Ancestry coefficients are per 1% Fremont ancestry.",
           "Aggregated survival strings use `stats::glm`. Mixed-model numbers use a genotype random intercept nested in population (survival: population fixed, genotype random).", "")
for (sct in unique(tab$section)) {
  lines <- c(lines, paste0("## ", sct), "", "| Location | Statistic | Paste-ready string | Flag / note |", "|---|---|---|---|")
  for (i in which(tab$section == sct)) {
    fl <- paste(c(if (nzchar(tab$flag[i])) tab$flag[i], if (nzchar(tab$note[i])) tab$note[i]), collapse = "; ")
    lines <- c(lines, sprintf("| %s | %s | %s | %s |", tab$location[i], tab$statistic[i], gsub("|", "\\|", tab$string[i], fixed = TRUE), fl))
  }
  lines <- c(lines, "")
}
con <- file("results/02_main_models/formatted_results_main_models.md", "w", encoding = "UTF-8"); writeLines(enc2utf8(lines), con); close(con)
cat("Wrote", nrow(tab), "strings;", sum(nzchar(tab$flag)), "flagged as significance changes\n")
print(tab %>% filter(nzchar(flag)) %>% select(location, statistic, flag), right = FALSE)
sink()
cat("Done. See results/02_main_models/main_models_followup_output.txt and results/02_main_models/formatted_results_main_models.md\n")
