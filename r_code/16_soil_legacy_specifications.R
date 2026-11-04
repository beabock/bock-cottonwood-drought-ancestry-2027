# ==============================================================================
# Soil legacy on colonization across every specification, under the cutoff labels.
# R2 / F / P come from the outputs of earlier steps: results/11_soil_legacy_permanova/output.txt (r_code/11_soil_legacy_permanova.R),
# results/02_main_models/main_models_followup_output.txt (genotype-level engine), results/02_main_models/main_models_output.txt (tree-level adonis2), results/02_main_models/formatted_results_main_models_df.md (Appendix S2 LMM).
# This script only rebuilds each design to report n and df (residual df = n - columns of the model matrix) and
# assembles the table; it does not rerun any permutation test.
# Run from the project root:  Rscript r_code/16_soil_legacy_specifications.R
# ==============================================================================
suppressMessages({ source("r_code/helpers/utilities_data_prep.R"); library(tidyverse) })
out <- "results/16_soil_legacy_specifications"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, geno_id = paste(pop, type, geno, sep = "_"))
d78 <- ds %>% filter(!is.na(per_f), !is.na(pop), !is.na(treatment), !is.na(amf), !is.na(emf), !is.na(dse), !is.na(soil), !is.na(myco_cat), !is.na(og_drc_mm)) %>%
  mutate(pop = factor(pop), treatment = factor(treatment), soil = factor(soil), myco_cat = factor(myco_cat))
foc <- d78 %>% filter(new_class %in% c("Narrowleaf", "F1 Hybrid", "Fremont")) %>% mutate(new_class = droplevels(new_class))
rdf <- function(d, rhs) nrow(d) - ncol(model.matrix(as.formula(paste("~", rhs)), droplevels(d)))
sc <- readLines("results/11_soil_legacy_permanova/output.txt", encoding = "UTF-8")
grab <- function(lab) { l <- sc[startsWith(sc, lab)][1]; v <- as.numeric(tail(strsplit(trimws(sub(lab, "", l, fixed = TRUE)), "\\s+")[[1]], 3)); setNames(v, c("R2", "F", "P")) }
blk <- function(form) { i <- which(trimws(sc) == form)[1]; l <- sc[i + which(startsWith(sc[(i + 1):(i + 12)], "soil"))[1]]; v <- strsplit(trimws(l), "\\s+")[[1]]; c(R2 = as.numeric(v[4]), F = as.numeric(v[5]), P = as.numeric(v[6])) }
fu <- readLines("results/02_main_models/main_models_followup_output.txt", encoding = "UTF-8"); v4 <- readLines("results/02_main_models/main_models_output.txt", encoding = "UTF-8")
eng <- function(after) { i <- grep(after, fu, fixed = TRUE)[1]; l <- fu[i + which(grepl("^\\s+soil", fu[(i + 1):(i + 8)]))[1]]
  v <- strsplit(trimws(l), "\\s+")[[1]]; c(R2 = as.numeric(v[3]), F = as.numeric(v[4]), P = as.numeric(v[5]), P_free = as.numeric(v[6])) }
adon <- function(call_pat, nth) { i <- grep(call_pat, v4, fixed = TRUE)[nth]; l <- v4[i + which(grepl("^soil", v4[(i + 1):(i + 10)]))[1]]
  v <- strsplit(trimws(l), "\\s+")[[1]]; c(R2 = as.numeric(v[4]), F = as.numeric(v[5]), P = as.numeric(v[6])) }

rows <- list()
add <- function(spec, scheme, perms, n, df2, v, src, reported = "", lab = "") rows[[length(rows) + 1]] <<- tibble(specification = spec, permutation_scheme = scheme, permutations = perms,
  n = n, df = sprintf("1, %d", df2), R2_pct = 100 * v[["R2"]], F = v[["F"]], P = v[["P"]], source = src, reported = reported)
g_ct <- eng("--- Cross type version"); g_pf <- eng("--- Continuous ancestry version")
add("Cross type (3 focal), population + drought + soil (marginal)", "Genotype-level engine: soil permuted among cuttings within genotype", 4999, nrow(foc), rdf(foc, "new_class + pop + treatment + soil"),
    g_ct, "results/02_main_models/main_models_followup_output.txt; r_code/02b_main_models_followup.R", "REPORTED (Results, colonization paragraph)")
add("Same model, cuttings permuted freely (no genotype structure)", "Tree-level, free", 4999, nrow(foc), rdf(foc, "new_class + pop + treatment + soil"),
    c(R2 = g_ct[["R2"]], F = g_ct[["F"]], P = g_ct[["P_free"]]), "results/02_main_models/main_models_followup_output.txt (p_trees_free)")
add("Continuous ancestry, population + drought + soil (marginal)", "Genotype-level engine: soil permuted within genotype", 4999, nrow(d78), rdf(d78, "per_f + pop + treatment + soil"),
    g_pf, "results/02_main_models/main_models_followup_output.txt; r_code/02b_main_models_followup.R", "formatted_results_main_models_df.md")
add("Same model, cuttings permuted freely", "Tree-level, free", 4999, nrow(d78), rdf(d78, "per_f + pop + treatment + soil"),
    c(R2 = g_pf[["R2"]], F = g_pf[["F"]], P = g_pf[["P_free"]]), "results/02_main_models/main_models_followup_output.txt (p_trees_free)")
add("Cross type (3 focal), adonis2 sequential", "Tree-level, strata = inoculum", 999, nrow(foc), rdf(foc, "new_class + pop + treatment + soil"),
    adon("adonis2(formula = fmat ~ new_class", 1), "results/02_main_models/main_models_output.txt; r_code/02a_main_models.R")
add("Cross type (3 focal), adonis2 marginal", "Tree-level, strata = inoculum", 999, nrow(foc), rdf(foc, "new_class + pop + treatment + soil"),
    adon("adonis2(formula = fmat ~ new_class", 2), "results/02_main_models/main_models_output.txt; r_code/02a_main_models.R")
add("Continuous ancestry, adonis2 sequential", "Tree-level, strata = inoculum", 999, nrow(d78), rdf(d78, "per_f + pop + treatment + soil"),
    adon("adonis2(formula = fmatb ~ per_f", 1), "results/02_main_models/main_models_output.txt; r_code/02a_main_models.R")
add("Continuous ancestry, adonis2 marginal", "Tree-level, strata = inoculum", 999, nrow(d78), rdf(d78, "per_f + pop + treatment + soil"),
    adon("adonis2(formula = fmatb ~ per_f", 2), "results/02_main_models/main_models_output.txt; r_code/02a_main_models.R")
specs <- list(
  list("continuous ancestry, 78 trees (reported)", "Continuous ancestry, 78 plants (soil-check engine)", d78, "per_f + pop + treatment + soil", "Within genotype"),
  list("tree-level free permutation, same model", "Continuous ancestry, cuttings permuted freely", d78, "per_f + pop + treatment + soil", "Tree-level, free"),
  list("cross type (3 focal; reported)", "Cross type (3 focal)", foc, "new_class + pop + treatment + soil", "Within genotype"),
  list("continuous ancestry, 3 focal cross types only", "Continuous ancestry, 3 focal cross types only", foc, "per_f + pop + treatment + soil", "Within genotype"),
  list("5-level cross type, 78 trees", "5-level cross type, 78 plants", d78, "new_class + pop + treatment + soil", "Within genotype"),
  list("continuous + inoculum category as covariate", "Continuous + inoculum covariate", d78, "per_f + pop + treatment + myco_cat + soil", "Within genotype"),
  list("continuous + initial diameter as covariate", "Continuous + initial diameter covariate", d78, "per_f + pop + treatment + og_drc_mm + soil", "Within genotype"),
  list("continuous, soil permuted within genotype x inoculum", "Continuous, soil permuted within genotype × inoculum", d78, "per_f + pop + treatment + soil", "Within genotype × inoculum"),
  list("continuous, no population term", "Continuous, no population term", d78, "per_f + treatment + soil", "Within genotype"))
for (s in specs) add(s[[2]], s[[5]], 9999, nrow(s[[3]]), rdf(s[[3]], s[[4]]), grab(s[[1]]), "results/11_soil_legacy_permanova/output.txt section 3; r_code/11_soil_legacy_permanova.R")
for (f in c("Ym ~ geno_id + soil", "Ym ~ geno_id + treatment + soil", "Ym ~ geno_id + treatment + myco_cat + soil")) {
  rhs <- sub("Ym ~ ", "", f, fixed = TRUE)
  add(paste("Genotype as fixed blocking term:", rhs), "adonis2, permutations within genotype blocks", 9999, nrow(d78), rdf(d78, rhs), blk(f), "results/11_soil_legacy_permanova/output.txt section 4; r_code/11_soil_legacy_permanova.R") }
rs <- readLines("results/02_main_models/formatted_results_main_models_df.md", encoding = "UTF-8"); l <- rs[grepl("Inoculum / soil / interaction", rs)]
m <- regmatches(l, regexpr("soil F\\(1, [0-9.]+\\) = [0-9.]+, P = [0-9.]+", l))
tab <- bind_rows(rows) %>% mutate(P_below_0.05 = P < 0.05)
cat("Soil legacy specifications (cutoff labels)\n\n"); print(as.data.frame(tab %>% mutate(R2_pct = round(R2_pct, 2), F = round(F, 3))), right = FALSE, row.names = FALSE)
cat("\nUnivariate colonization LMM (Appendix S2, all scored plants, Type II):", m, "(formatted_results_main_models_df.md; not a PERMANOVA)\n")
lo <- readLines("results/11_soil_legacy_permanova/output.txt"); cat(lo[grepl("p range over", lo)], "\n")
write.csv(tab, file.path(out, "soil_specifications.csv"), row.names = FALSE)
