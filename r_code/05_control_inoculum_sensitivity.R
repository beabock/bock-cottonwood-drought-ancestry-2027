# ==============================================================================
# CONTROL-INOCULUM SENSITIVITY: control-inoculum sensitivity (contamination question), plants whose
# new biomass exceeds leaf mass, and treatment / inoculum / cross type splits of the
# stomatal conductance and colonization subsets. Models are those of 02a.
# ==============================================================================
# Run from the project root:  Rscript r_code/05_control_inoculum_sensitivity.R   (about 3 min)
# Writes: results/05_control_inoculum_sensitivity/output.txt, results/05_control_inoculum_sensitivity/control_sensitivity.csv, results/05_control_inoculum_sensitivity/newbiomass_gt_leaf.csv,
#         (Appendix S4 is built from control_sensitivity.csv by r_code/30_appendix_tables.R)
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(lme4); library(lmerTest); library(tidyverse); library(vegan)
sink("results/05_control_inoculum_sensitivity/output.txt", split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("CONTROL INOCULUM SENSITIVITY run:", "\n")
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, geno_id = paste(pop, type, geno, sep = "_"))

EN <- "–"; BETA <- "β"; TIMES <- "×"; R2 <- "R²"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
tstr <- function(cs, term) sprintf("%s = %s, t = %s, df = %s, %s", BETA, num(cs[term, "Estimate"]), dec(cs[term, "t value"]), dff(cs[term, "df"]), P(cs[term, "Pr(>|t|)"]))
Fstr <- function(r2, F, df1, df2, p) sprintf("%s = %s%%, F = %s, df = %s, %s, %s", R2, dec(100 * r2, 1), dec(F), dff(df1), dff(df2), P(p))

res <- list()
add <- function(analysis, data, n, term, string, p, sign_val) res[[length(res) + 1]] <<- tibble(analysis, data, n, term, string, p, sign = sign(sign_val))
sets <- list(`All plants (reported)` = ds, `Without control-inoculum plants` = filter(ds, myco_cat != "None"))

# ------------------------------------------------------------------------------
# 1a. Colonization LMM: Fremont ancestry x EMF
# ------------------------------------------------------------------------------
hdr("1a. Colonization LMM: Fremont ancestry x EMF interaction")
for (nm in names(sets)) {
  dl <- sets[[nm]] %>% filter(!is.na(amf), !is.na(emf), !is.na(dse)) %>%
    select(tree_id, pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>% na.omit() %>%
    mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm)), myco_cat = droplevels(factor(myco_cat))) %>%
    pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>%
    mutate(fungal_type = factor(fungal_type, levels = c("amf", "dse", "emf")))
  m <- lmerTest::lmer(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil +
                        (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id), data = dl)
  cs <- coef(summary(m)); n <- n_distinct(dl$tree_id)
  cat(sprintf("%s: n = %d trees (%d genotypes)\n  %s\n", nm, n, n_distinct(dl$geno_id), tstr(cs, "per_f:fungal_typeemf")))
  add("Colonization LMM", nm, n, paste0("Fremont ancestry ", TIMES, " EMF"), tstr(cs, "per_f:fungal_typeemf"), cs["per_f:fungal_typeemf", "Pr(>|t|)"], cs["per_f:fungal_typeemf", "Estimate"])
}

# ------------------------------------------------------------------------------
# 1b. Colonization-biomass LMM: AMF, EMF, DSE x drought
# ------------------------------------------------------------------------------
hdr("1b. Colonization-biomass LMM: fungus x drought")
terms_b <- c(EMF = "treatmentDrought:scale(emf)", AMF = "scale(amf):treatmentDrought", DSE = "treatmentDrought:scale(dse)")
for (nm in names(sets)) {
  d <- sets[[nm]] %>% mutate(myco_cat = droplevels(factor(myco_cat)))
  m <- lmerTest::lmer(log(dryshoot_new) ~ scale(amf) * treatment + scale(emf) * treatment + scale(dse) * treatment +
                        per_fc + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = d)
  cs <- coef(summary(m)); n <- nobs(m)
  cat(sprintf("%s: n = %d\n", nm, n))
  for (k in names(terms_b)) {
    cat(sprintf("  %s %s drought: %s\n", k, TIMES, tstr(cs, terms_b[k])))
    add("Colonization-biomass LMM", nm, n, paste(k, TIMES, "drought"), tstr(cs, terms_b[k]), cs[terms_b[k], "Pr(>|t|)"], cs[terms_b[k], "Estimate"])
  }
}

# ------------------------------------------------------------------------------
# 1c. PERMANOVA, genotype-level permutations (engine from 02b_main_models_followup.R)
# ------------------------------------------------------------------------------
hdr("1c. PERMANOVA with whole-genotype permutations (marginal tests; 4999 permutations)")
perm_engine <- function(dat, anc, B = 4999) {
  dat <- droplevels(dat)
  Yc <- scale(as.matrix(dat[, c("amf", "emf", "dse")]), scale = FALSE); tss <- sum(Yc^2)
  fm <- as.formula(paste("~", anc, "+ pop + treatment + soil"))
  X0 <- model.matrix(fm, dat); asg <- attr(X0, "assign"); k <- 1
  statk <- function(X) { rf <- sum(qr.resid(qr(X), Yc)^2); rr <- sum(qr.resid(qr(X[, asg != k, drop = FALSE]), Yc)^2)
    c(F = ((rr - rf) / sum(asg == k)) / (rf / (nrow(X) - ncol(X))), R2 = (rr - rf) / tss) }
  gt <- unique(dat[, c("geno_id", "pop", anc)])
  obs <- statk(X0)
  fg <- replicate(B, { g <- gt
    for (p in unique(g$pop)) { i <- which(g$pop == p); g[[anc]][i] <- g[[anc]][i[sample.int(length(i))]] }
    d <- dat; d[[anc]] <- g[[anc]][match(d$geno_id, g$geno_id)]; statk(model.matrix(fm, d))["F"] })
  tibble(R2 = obs["R2"], F = obs["F"], df1 = sum(asg == 1), df2 = nrow(X0) - ncol(X0), p = (1 + sum(fg >= obs["F"] - 1e-9)) / (B + 1),
         n_trees = nrow(dat), n_genotypes = n_distinct(dat$geno_id))
}
for (nm in names(sets)) {
  d78 <- sets[[nm]] %>% filter(!is.na(per_f), !is.na(pop), !is.na(treatment), !is.na(amf), !is.na(emf), !is.na(dse), !is.na(soil), !is.na(myco_cat)) %>%
    mutate(pop = factor(pop), treatment = factor(treatment), soil = factor(soil))
  foc <- d78 %>% filter(new_class %in% c("Narrowleaf", "F1 Hybrid", "Fremont")) %>% mutate(new_class = droplevels(new_class))
  set.seed(1998)
  for (v in list(list(d = foc, anc = "new_class", lab = "Cross type (3 focal)"), list(d = d78, anc = "per_f", lab = "Continuous ancestry"))) {
    r <- perm_engine(v$d, v$anc); s <- Fstr(r$R2, r$F, r$df1, r$df2, r$p)
    cat(sprintf("%s, %s: %d trees, %d genotypes\n  %s\n", nm, v$lab, r$n_trees, r$n_genotypes, s))
    add("PERMANOVA (genotype-level)", nm, r$n_trees, v$lab, s, r$p, 1)
  }
}
cat("\nNote: P values for 'All plants' are recomputed here with a fresh permutation stream, so they can differ in the third decimal from formatted_results_main_models_df.md.\n")

hdr("1. Summary: does any conclusion change without control plants?")
tab <- bind_rows(res) %>% group_by(analysis, term) %>%
  mutate(sig = p < 0.05, conclusion_changes = (sig != sig[data == "All plants (reported)"]) | (sign != sign[data == "All plants (reported)"])) %>% ungroup()
print(as.data.frame(tab %>% select(analysis, term, data, n, string, conclusion_changes)), right = FALSE)
write.csv(tab, "results/05_control_inoculum_sensitivity/control_sensitivity.csv", row.names = FALSE)

# The Appendix S4 table is built from control_sensitivity.csv by r_code/30_appendix_tables.R.

hdr("1. Colonization range by guild (min / median / max, %), all scored plants")
rng <- function(d, lab) d %>% pivot_longer(c(amf, emf, dse), names_to = "guild", values_to = "col") %>% filter(!is.na(col)) %>%
  group_by(guild = toupper(guild)) %>% summarise(n = n(), min = min(col), median = median(col), max = max(col), zero = sum(col == 0), .groups = "drop") %>% mutate(set = lab)
rt <- bind_rows(rng(ds, "With controls (all inoculum treatments)"), rng(filter(ds, myco_cat != "None"), "Without controls"),
                rng(filter(ds, myco_cat == "None"), "Controls only")) %>% mutate(guild = factor(guild, levels = c("AMF", "EMF", "DSE"))) %>% arrange(set, guild)
print(as.data.frame(rt %>% mutate(across(c(min, median, max), ~ round(., 1)))), right = FALSE)

# ------------------------------------------------------------------------------
# 2. New biomass > leaf mass
# ------------------------------------------------------------------------------
hdr("2. Plants whose new biomass (dryshoot_new) exceeds dry leaf mass")
b <- ds %>% filter(!is.na(dryleaf_weight), !is.na(dryshoot_new))
big <- b %>% filter(dryshoot_new > dryleaf_weight) %>%
  transmute(tree_idold, cross_type = new_class, pct_fremont = round(per_f, 1), treatment, pop, genotype = geno_id, inoculum = myco_cat, soil,
            new_biomass_g = dryshoot_new, leaf_g = dryleaf_weight, stem_g = dry_stem, ratio_new_to_leaf = round(dryshoot_new / dryleaf_weight, 2)) %>%
  arrange(cross_type, pop, tree_idold)
cat("Plants with both values:", nrow(b), "| new biomass > leaf mass:", nrow(big), "\n\n")
print(as.data.frame(big), right = FALSE)
write.csv(big, "results/05_control_inoculum_sensitivity/newbiomass_gt_leaf.csv", row.names = FALSE)
cmp <- function(v) {
  t <- table(b[[v]], b$dryshoot_new > b$dryleaf_weight)
  cat(sprintf("\nBy %s (n with new > leaf / n total, %%):\n", v))
  print(data.frame(level = rownames(t), flagged = t[, "TRUE"], total = rowSums(t), pct = round(100 * t[, "TRUE"] / rowSums(t), 1)), row.names = FALSE)
  ft <- tryCatch(fisher.test(t, simulate.p.value = TRUE, B = 1e5)$p.value, error = function(e) NA)
  cat("  Fisher exact test (simulated) P =", signif(ft, 3), "\n")
}
for (v in c("new_class", "pop", "treatment", "myco_cat", "soil")) cmp(v)
cat("\nGenotypes contributing flagged plants:\n"); print(sort(table(big$genotype), decreasing = TRUE))
cat("Median ratio new/leaf among flagged plants:", median(big$ratio_new_to_leaf), "\n")

# ------------------------------------------------------------------------------
# 3. Splits of the stomatal conductance and colonization subsets
# ------------------------------------------------------------------------------
hdr("3. How the stomatal conductance (204) and colonization (78) plants split")
gsw <- ds %>% filter(!is.na(por_avg), por_avg > 0)
col <- ds %>% filter(!is.na(amf), !is.na(emf), !is.na(dse), !is.na(per_f), !is.na(soil), !is.na(treatment), !is.na(myco_cat))
cat("Stomatal conductance plants (por_avg > 0):", nrow(gsw), "| colonization plants (all three guilds scored):", nrow(col), "\n")
for (v in c("treatment", "myco_cat", "new_class")) {
  cat(sprintf("\n-- %s --\n", v))
  print(data.frame(level = levels(factor(ds[[v]])),
                   conductance = as.vector(table(factor(gsw[[v]], levels = levels(factor(ds[[v]]))))),
                   colonization = as.vector(table(factor(col[[v]], levels = levels(factor(ds[[v]]))))),
                   all_plants = as.vector(table(factor(ds[[v]], levels = levels(factor(ds[[v]])))))), row.names = FALSE)
}
cat("\nConductance plants, cross type x treatment:\n"); print(table(gsw$new_class, gsw$treatment))
cat("\nColonization plants, cross type x treatment:\n"); print(table(col$new_class, col$treatment))
cat("\nColonization plants, inoculum x treatment:\n"); print(table(col$myco_cat, col$treatment))
cat("\nOverlap: plants in both subsets:", length(intersect(gsw$tree_idold, col$tree_idold)), "\n")
sink()
