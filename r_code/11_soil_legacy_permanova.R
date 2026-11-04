# ==============================================================================
# Robustness of the genotype-level soil-legacy result in the continuous-ancestry
# PERMANOVA (marginal R2 = 2.5%, within-genotype permutation p = 0.031; tree-level
# p = 0.109). Standalone: does not source 02a_main_models.R.
# Run from the project root:  Rscript r_code/11_soil_legacy_permanova.R
# Writes results/11_soil_legacy_permanova/output.txt
# ==============================================================================
stopifnot(file.exists("data/edited_data/mega_ds.csv"))
source("r_code/helpers/utilities_data_prep.R")
library(vegan); library(permute); library(lme4); library(lmerTest); library(tidyverse)
sink("results/11_soil_legacy_permanova/output.txt", split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("SOIL CHECK run:", "\n")

ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, geno_id = paste(pop, type, geno, sep = "_"))
d78 <- ds %>% filter(!is.na(per_f), !is.na(pop), !is.na(treatment), !is.na(amf), !is.na(emf), !is.na(dse), !is.na(soil), !is.na(myco_cat),
                     !is.na(og_drc_mm)) %>%
  mutate(pop = factor(pop), treatment = factor(treatment), soil = factor(soil), myco_cat = factor(myco_cat))
stopifnot(nrow(d78) == 78)
foc <- d78 %>% filter(new_class %in% c("Narrowleaf", "F1 Hybrid", "Fremont")) %>% mutate(new_class = droplevels(new_class))

# marginal F / R2 for one term, with a within-group permutation of soil -----------
soil_test <- function(dat, rhs = "per_f + pop + treatment + soil", term = "soil", B = 9999, group = "geno_id", seed = 1) {
  set.seed(seed)
  dat <- droplevels(dat); Y <- scale(as.matrix(dat[, c("amf", "emf", "dse")]), scale = FALSE); tss <- sum(Y^2)
  fm <- as.formula(paste("~", rhs)); X0 <- model.matrix(fm, dat); asg <- attr(X0, "assign")
  k <- match(term, attr(terms(fm), "term.labels"))
  stat <- function(X) { rf <- sum(qr.resid(qr(X), Y)^2); rr <- sum(qr.resid(qr(X[, asg != k, drop = FALSE]), Y)^2)
    c(F = ((rr - rf) / sum(asg == k)) / (rf / (nrow(X) - ncol(X))), R2 = (rr - rf) / tss) }
  obs <- stat(X0); key <- if (is.null(group)) rep(1, nrow(dat)) else dat[[group]]
  lv <- levels(dat[[term]])
  fp <- replicate(B, { d <- dat; d[[term]] <- factor(ave(as.character(d[[term]]), key, FUN = function(x) x[sample.int(length(x))]), levels = lv); stat(model.matrix(fm, d))["F"] })
  c(R2 = unname(obs["R2"]), F = unname(obs["F"]), p = (1 + sum(fp >= obs["F"] - 1e-9)) / (B + 1), mc_se = sqrt(0.05 * 0.95 / B))
}

# ------------------------------------------------------------------------------
hdr("0. Was inoculum used as strata in the genotype-level scheme?")
cat("No. In r_code/02b_main_models_followup.R the genotype-level permutations (whole genotypes for ancestry and population; trees within genotype\n",
    "for drought and soil) use no inoculum strata, and neither does the genotype-means PERMANOVA. Inoculum strata (strata = myco_cat) were used\n",
    "only in the tree-level (T) rows from 02a_main_models.R. The 'p_trees_free' column in the follow-up output is also unstratified.\n", sep = "")
cat("Inoculum composition within genotypes of the 78 trees (genotypes with trees in >1 inoculum category):",
    sum(tapply(d78$myco_cat, d78$geno_id, function(x) n_distinct(x) > 1)), "of", n_distinct(d78$geno_id), "\n")

# ------------------------------------------------------------------------------
hdr("1. How much within-genotype soil contrast exists?")
sg <- d78 %>% group_by(geno_id) %>% summarise(n = n(), n_ww_soil = sum(soil == "No Tam Legacy"), n_tam = sum(soil == "Tam Legacy"), .groups = "drop")
inf <- sg %>% filter(n_ww_soil > 0, n_tam > 0)
cat("genotypes:", nrow(sg), "| singletons (no within-genotype information):", sum(sg$n == 1),
    "| genotypes with BOTH soils (informative):", nrow(inf), "covering", sum(inf$n), "of 78 trees\n")
cat("Number of distinct within-genotype soil arrangements (permutation space):", format(prod(choose(sg$n, sg$n_tam)), scientific = TRUE, digits = 3), "\n")
cat("Soil split overall:", paste(names(table(d78$soil)), table(d78$soil), collapse = " / "), "\n")

# ------------------------------------------------------------------------------
hdr("2. Replication: Monte Carlo error and seeds (continuous ancestry, within-genotype permutation of soil)")
reps <- sapply(1:5, function(s) soil_test(d78, B = 19999, seed = s))
colnames(reps) <- paste("seed", 1:5); print(round(reps, 4))
cat("p range across seeds:", paste(round(range(reps["p", ]), 4), collapse = " to "), "| Monte Carlo SE at B = 19999 is about", round(sqrt(0.03 * 0.97 / 19999), 4), "\n")

# ------------------------------------------------------------------------------
hdr("3. Same term under other specifications (within-genotype permutation, B = 9999 unless noted)")
spec <- list(
  `continuous ancestry, 78 trees (reported)` = function() soil_test(d78),
  `tree-level free permutation, same model` = function() soil_test(d78, group = NULL),
  `cross type (3 focal; reported)` = function() soil_test(foc, rhs = "new_class + pop + treatment + soil"),
  `continuous ancestry, 3 focal cross types only` = function() soil_test(foc, rhs = "per_f + pop + treatment + soil"),
  `5-level cross type, 78 trees` = function() soil_test(d78, rhs = "new_class + pop + treatment + soil"),
  `continuous + inoculum category as covariate` = function() soil_test(d78, rhs = "per_f + pop + treatment + myco_cat + soil"),
  `continuous + initial diameter as covariate` = function() soil_test(d78, rhs = "per_f + pop + treatment + og_drc_mm + soil"),
  `continuous, soil permuted within genotype x inoculum` = function() { d <- d78; d$gi <- paste(d$geno_id, d$myco_cat); soil_test(d, group = "gi") },
  `continuous, no population term` = function() soil_test(d78, rhs = "per_f + treatment + soil"))
tab <- t(sapply(spec, function(f) f())); print(round(tab, 4))

# ------------------------------------------------------------------------------
hdr("4. Genotype as a fixed blocking term (test of soil within genotype only; permutations within genotype)")
Ym <- as.matrix(d78[, c("amf", "emf", "dse")])
for (f in list(Ym ~ geno_id + soil, Ym ~ geno_id + treatment + soil, Ym ~ geno_id + treatment + myco_cat + soil)) {
  set.seed(1)
  a <- adonis2(f, data = d78, method = "euclidean", by = "margin",
               permutations = how(blocks = factor(d78$geno_id), nperm = 9999))
  cat("\n", deparse(f), "\n"); print(a[rownames(a) == "soil", ])
}

# ------------------------------------------------------------------------------
hdr("5. Leave-one-genotype-out: does p depend on a single genotype? (B = 1999 each)")
lo <- map_dfr(unique(d78$geno_id), function(g) {
  r <- soil_test(filter(d78, geno_id != g), B = 1999, seed = 11)
  tibble(dropped = g, n_trees_dropped = sum(d78$geno_id == g), R2 = r["R2"], F = r["F"], p = r["p"])
})
print(as.data.frame(lo %>% arrange(desc(p)) %>% head(8) %>% mutate(across(c(R2, F, p), ~ round(., 4)))), row.names = FALSE)
cat("p range over", nrow(lo), "leave-one-genotype-out fits:", paste(round(range(lo$p), 3), collapse = " to "),
    "| fits with p >= 0.05:", sum(lo$p >= 0.05), "| median p:", round(median(lo$p), 3), "\n")

# ------------------------------------------------------------------------------
hdr("6. Which fungus drives it, and is it visible in univariate models?")
cat("Raw means by soil (all 78 trees):\n"); print(round(d78 %>% group_by(soil) %>% summarise(n = n(), AMF = mean(amf), EMF = mean(emf), DSE = mean(dse)) %>% as.data.frame() %>% column_to_rownames("soil"), 2))
cat("\nWithin-genotype soil differences (Tam Legacy minus No Tam Legacy; genotypes with both soils; mean of genotype-level differences):\n")
wd <- d78 %>% group_by(geno_id, soil) %>% summarise(across(c(amf, emf, dse), mean), .groups = "drop") %>%
  pivot_wider(names_from = soil, values_from = c(amf, emf, dse)) %>% drop_na()
for (v in c("amf", "emf", "dse")) { dd <- wd[[paste0(v, "_Tam Legacy")]] - wd[[paste0(v, "_No Tam Legacy")]]
  cat(sprintf("  %s: mean diff %.2f (SE %.2f), paired t p = %.3f, n genotypes = %d\n", toupper(v), mean(dd), sd(dd) / sqrt(length(dd)), t.test(dd)$p.value, length(dd))) }
cat("\nUnivariate mixed models, y ~ per_f + treatment + soil + (1|pop) + (1|pop:geno_id):\n")
for (v in c("amf", "emf", "dse")) {
  m <- lmerTest::lmer(as.formula(paste(v, "~ per_f + treatment + soil + (1 | pop) + (1 | pop:geno_id)")), data = d78)
  cf <- coef(summary(m))["soilTam Legacy", ]
  cat(sprintf("  %s: soil beta = %.2f, t = %.2f, p = %.3f\n", toupper(v), cf["Estimate"], cf["t value"], cf["Pr(>|t|)"]))
}
dl <- d78 %>% select(tree_id, pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>%
  pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>% mutate(fungal_type = factor(fungal_type))
l0 <- lmerTest::lmer(colonization ~ per_f * fungal_type + treatment + soil + (1 | pop) + (1 | pop:geno_id) + (1 | tree_id), data = dl)
l1 <- lmerTest::lmer(colonization ~ per_f * fungal_type + treatment + soil * fungal_type + (1 | pop) + (1 | pop:geno_id) + (1 | tree_id), data = dl)
cat("\nLong-format LMM: soil main effect (averaged over fungi) p =", round(coef(summary(l0))["soilTam Legacy", "Pr(>|t|)"], 3), "\n")
cat("Soil x fungal type interaction (does soil shift the AMF/EMF/DSE mix?):\n"); print(anova(l0, l1, refit = TRUE))

# ------------------------------------------------------------------------------
hdr("7. Multiplicity")
cat("The genotype-level tables contain 8 tests (4 terms x 2 versions). Soil p = 0.031 would be 0.25 after Bonferroni over all 8 and 0.063 over the 2 soil tests.\n",
    "Soil is not a hypothesis stated in the Introduction (it is a covariate/nuisance factor), so it is exploratory.\n", sep = "")
sink()
cat("Done. See results/11_soil_legacy_permanova/output.txt\n")
