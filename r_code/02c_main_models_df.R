# ==============================================================================
# MAIN MODELS, DEGREES OF FREEDOM: degrees of freedom for every reported statistic (AJB format),
# aggregated-survival grouping sensitivity, and a check of dry_shoottotal.
# No model is changed: all reported models are those of 02a/02b, obtained by
# sourcing 02b_main_models_followup.R (which sources 02a_main_models.R).
# ==============================================================================
# Run from the project root:  Rscript r_code/02c_main_models_df.R   (about 8 min)
# Writes: results/02_main_models/formatted_results_main_models_df.md, results/02_main_models/main_models_df_output.txt
# Note: sourcing 02a and 02b rewrites their results/ outputs; they are
# deterministic, so only the "run:" timestamp lines change.
# ==============================================================================

source("r_code/02b_main_models_followup.R")
sink("results/02_main_models/main_models_df_output.txt", split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("MAIN MODELS WITH DF run:", "\n")

# ------------------------------------------------------------------------------
# AJB formatting helpers: capital P, en dash for minus, df in parentheses
# ------------------------------------------------------------------------------
EN <- "–"; BETA <- "β"; TIMES <- "×"; CHI <- "χ²"; R2 <- "R²"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))   # d significant digits, trailing zeros kept
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
tstr <- function(t, df) sprintf("t(%s) = %s", dff(df), dec(t))
zstr <- function(z) sprintf("z = %s", dec(z))
Fstr <- function(F, df1, df2) sprintf("F(%s, %s) = %s", dff(df1), dff(df2), dec(F))
chistr <- function(x, df) sprintf("%s(%s) = %s", CHI, dff(df), dec(x))
lmt <- function(m, term) { cs <- coef(summary(m)); sprintf("%s = %s, %s, %s", BETA, num(cs[term, "Estimate"]), tstr(cs[term, "t value"], cs[term, "df"]), P(cs[term, "Pr(>|t|)"])) }
glz <- function(tab, term, d = 3) sprintf("%s = %s, %s, %s", BETA, num(tab[term, "Estimate"], d), zstr(tab[term, "z value"]), P(tab[term, "Pr(>|z|)"]))
resid_df <- function(dat, rhs) nrow(dat) - ncol(model.matrix(as.formula(paste("~", rhs)), droplevels(dat)))

rows <- list()
add <- function(section, loc, stat, string, ms_p = NA, p = NA, note = "") {
  flag <- if (!is.na(ms_p) && !is.na(p) && ((ms_p < 0.05) != (p < 0.05))) paste0("**SIGNIFICANCE CHANGED** (manuscript P ", ifelse(ms_p < 0.001, "< 0.001", ms_p), " → ", P(p), ")") else ""
  rows[[length(rows) + 1]] <<- data.frame(section = section, location = loc, statistic = stat, string = string, flag = flag, note = note)
}
pcol <- function(m, term) coef(summary(m))[term, "Pr(>|t|)"]

# ------------------------------------------------------------------------------
# 1. Strings with df
# ------------------------------------------------------------------------------
S1 <- "Results: Biomass and Survival Vary with Host Ancestry"
add(S1, "P1 / Fig 1 caption", "Survival ancestry x drought (full model)",
    sprintf("linear and quadratic ancestry %s drought terms, %s and %s (joint LRT %s, %s)", TIMES,
            P(ssf["per_fc:treatmentDrought", "Pr(>|z|)"]), P(ssf["treatmentDrought:I(per_fc^2)", "Pr(>|z|)"]),
            chistr(surv_lrt$Chisq[2], surv_lrt$Df[2]), P(surv_lrt$`Pr(>Chisq)`[2])), 0.89, ssf["per_fc:treatmentDrought", "Pr(>|z|)"], "GLMM: Wald z tests have no df")
add(S1, "P1 / Fig 2 caption", "Biomass LRT ancestry x drought",
    sprintf("LRT for ancestry %s drought treatment, %s, %s", TIMES, chistr(bio_lrt4$Chisq[2], bio_lrt4$Df[2]), P(bio_lrt4$`Pr(>Chisq)`[2])), 0.58, bio_lrt4$`Pr(>Chisq)`[2])
add(S1, "P1 / Fig 1 caption", "Survival linear ancestry (reduced)", glz(ssr, "per_fc"), 0.028, ssr["per_fc", "Pr(>|z|)"])
add(S1, "P1", "Survival quadratic ancestry (reduced)", glz(ssr, "I(per_fc^2)", 2), 0.40, ssr["I(per_fc^2)", "Pr(>|z|)"])
add(S1, "P1 / Fig 1 caption", "Survival drought (reduced)", glz(ssr, "treatmentDrought", 4), 0.0005, ssr["treatmentDrought", "Pr(>|z|)"])
add(S1, "P1 / Fig 2 caption", "Biomass linear ancestry (reduced)", lmt(bio_red4, "per_fc"), 0.0005, pcol(bio_red4, "per_fc"))
add(S1, "P1 / Fig 2 caption", "Biomass quadratic ancestry (reduced)", lmt(bio_red4, "I(per_fc^2)"), 0.0005, pcol(bio_red4, "I(per_fc^2)"), "still significant, but no longer P < 0.001")
add(S1, "P1 / Fig 2 caption", "Biomass drought (reduced)", lmt(bio_red4, "treatmentDrought"), 0.045, pcol(bio_red4, "treatmentDrought"))

S2 <- "Results: Host Ancestry, Not Drought, Structures Fungal Colonization (PERMANOVA)"
rdf_foc <- resid_df(foc, "new_class + pop + treatment + soil"); rdf_all <- resid_df(d78, "per_f + pop + treatment + soil")
G <- function(eng, t) eng[eng$term == t, ]
gstr <- function(eng, t, rdf) sprintf("%s = %s%%, %s, %s", R2, dec(100 * G(eng, t)$R2, 1), Fstr(G(eng, t)$F, G(eng, t)$df, rdf), P(G(eng, t)$p_genotype_level))
add(S2, "P2", sprintf("(G) Cross type, whole-genotype permutations (%d genotypes, %d trees)", n_distinct(foc$geno_id), nrow(foc)), gstr(eng_foc, "new_class", rdf_foc), 0.003, G(eng_foc, "new_class")$p_genotype_level, "marginal tests; population, drought and soil in the same model")
add(S2, "P2", "(G) Population (same model)", gstr(eng_foc, "pop", rdf_foc), 0.352, G(eng_foc, "pop")$p_genotype_level)
add(S2, "P2 / Supp Fig 3", "(G) Drought, permuted within genotype", gstr(eng_foc, "treatment", rdf_foc), 0.79, G(eng_foc, "treatment")$p_genotype_level)
add(S2, "P2", "(G) Soil legacy, permuted within genotype", gstr(eng_foc, "soil", rdf_foc), 0.27, G(eng_foc, "soil")$p_genotype_level)
add(S2, "P2", "(G) Continuous ancestry (37 genotypes, 78 trees)", gstr(eng_all, "per_f", rdf_all), 0.001, G(eng_all, "per_f")$p_genotype_level)
add(S2, "P2", "(G) Continuous version: population / drought", sprintf("population: %s; drought: %s", gstr(eng_all, "pop", rdf_all), gstr(eng_all, "treatment", rdf_all)), NA, NA)
add(S2, "P2", "(G) Continuous version: soil legacy", gstr(eng_all, "soil", rdf_all), 0.27, G(eng_all, "soil")$p_genotype_level,
    "NOT ROBUST (P varies across specifications; results/11_soil_legacy_permanova/output.txt, results/16_soil_legacy_specifications/soil_specifications.csv); treat as exploratory")
TA <- as.data.frame(permA); TB <- as.data.frame(permB)
tl <- function(T, t) sprintf("%s = %s%%, %s, %s", R2, dec(100 * T[t, "R2"], 1), Fstr(T[t, "F"], T[t, "Df"], T["Residual", "Df"]), P(T[t, "Pr(>F)"]))
add(S2, "P2", "(T) Cross type, tree-level (sequential; strata = inoculum)", tl(TA, "new_class"), 0.003, TA["new_class", "Pr(>F)"])
add(S2, "P2", "(T) Population / drought / soil, tree-level", sprintf("population: %s; drought: %s; soil: %s", tl(TA, "pop"), tl(TA, "treatment"), tl(TA, "soil")), 0.79, TA["treatment", "Pr(>F)"])
add(S2, "P2", "(T) Continuous ancestry, tree-level", tl(TB, "per_f"), 0.001, TB["per_f", "Pr(>F)"])
gmt <- function(res, t) sprintf("%s = %s%%, %s, %s", R2, dec(100 * res$seq[t, "R2"], 1), Fstr(res$seq[t, "F"], res$seq[t, "Df"], res$seq["Residual", "Df"]), P(res$seq[t, "Pr(>F)"]))
add(S2, "P2 (alternative)", sprintf("(M) Genotype means, cross type only (%d genotypes)", nrow(gmf)), gmt(gb_cls1, "new_class"), 0.003, gb_cls1$seq["new_class", "Pr(>F)"])
add(S2, "P2 (alternative)", "(M) Genotype means, cross type + population", sprintf("cross type: %s; population: %s", gmt(gb_cls2, "new_class"), gmt(gb_cls2, "pop")), 0.003, gb_cls2$seq["new_class", "Pr(>F)"])
add(S2, "P2 (alternative)", "(M) Genotype means, continuous ancestry (37 genotypes)", gmt(gb_pf1, "per_f"), 0.001, gb_pf1$seq["per_f", "Pr(>F)"])
pw_str <- function(pwt, counts) {
  sp <- strsplit(as.character(pwt$pairs), " vs "); out <- character(nrow(pwt))
  for (i in seq_len(nrow(pwt))) out[i] <- sprintf("%s vs. %s: %s, %s", sp[[i]][1], sp[[i]][2],
      Fstr(pwt$F.Model[i], pwt$Df[i], counts[sp[[i]][1]] + counts[sp[[i]][2]] - 2), P(pwt$p.adjusted[i]))
  paste(out, collapse = "; ") }
add(S2, "P3 pairwise", "(M) Pairwise on genotype means (Bonferroni)", pw_str(gm_pw, table(as.character(gmf$new_class))), 0.051, gm_pw$p.adjusted[2])
add(S2, "P3 pairwise", "(T) Pairwise on trees (Bonferroni)", pw_str(pw, table(as.character(foc$new_class))), 0.051, pw$p.adjusted[2])
ct <- coef(summary(col4)); nl <- summary(col4n)$tTable
add(S2, "P3 / Fig 3 caption", "EMF x ancestry (lmer)", sprintf("Fremont ancestry %s EMF interaction, %s, %s", TIMES, tstr(ct["per_f:fungal_typeemf", "t value"], ct["per_f:fungal_typeemf", "df"]), P(ct["per_f:fungal_typeemf", "Pr(>|t|)"])), 0.0005, ct["per_f:fungal_typeemf", "Pr(>|t|)"])
add(S2, "P3 / Methods", "EMF x ancestry (nlme sensitivity)", sprintf("%s, %s", tstr(nl["per_f:fungal_typeemf", "t-value"], nl["per_f:fungal_typeemf", "DF"]), P(nl["per_f:fungal_typeemf", "p-value"])), 0.0005, nl["per_f:fungal_typeemf", "p-value"],
    "nlme reports containment df, not Satterthwaite (lmerTest does not support nlme models)")
add(S2, "Fig 3 caption", "Ancestry slope by fungal type", sprintf("EMF slope = %s per 1%% (%s, %s); AMF = %s (%s, %s); DSE = %s (%s, %s)",
    num(sl$per_f.trend[3]), tstr(sl$t.ratio[3], sl$df[3]), P(sl$p.value[3]), num(sl$per_f.trend[1]), tstr(sl$t.ratio[1], sl$df[1]), P(sl$p.value[1]),
    num(sl$per_f.trend[2]), tstr(sl$t.ratio[2], sl$df[2]), P(sl$p.value[2])), NA, NA, "Satterthwaite df via emmeans")

S3 <- "Results: Survival associations (aggregated groups; stats::glm)"
gl <- coef(summary(agg_glm_stats)); gn <- coef(summary(agg_glm_naive))
add(S3, "Survival associations / Fig 4", "AMF x drought (ancestry-controlled)", glz(gl, "scale(mean_amf):treatmentDrought"), 0.018, gl["scale(mean_amf):treatmentDrought", "Pr(>|z|)"], "GLM: Wald z tests have no df")
add(S3, "Survival associations / Fig 4", "AMF x drought (naive)", glz(gn, "scale(mean_amf):treatmentDrought"), 0.010, gn["scale(mean_amf):treatmentDrought", "Pr(>|z|)"])
add(S3, "Survival associations", "EMF x drought", glz(gl, "scale(mean_emf):treatmentDrought"), 0.179, gl["scale(mean_emf):treatmentDrought", "Pr(>|z|)"])
add(S3, "Survival associations", "DSE x drought", glz(gl, "scale(mean_dse):treatmentDrought"), 0.325, gl["scale(mean_dse):treatmentDrought", "Pr(>|z|)"])
add(S3, "Survival associations / Limitations", "Ancestry terms", sprintf("linear %s; quadratic %s", P(gl["mean_per_fc", "Pr(>|z|)"]), P(gl["I(mean_per_fc^2)", "Pr(>|z|)"])), 0.74, gl["mean_per_fc", "Pr(>|z|)"])

S4 <- "Results: Biomass associations (n = 77)"
add(S4, "Biomass associations / Fig 5", "EMF x drought", lmt(fb_pf4, "treatmentDrought:scale(emf)"), 0.004, pcol(fb_pf4, "treatmentDrought:scale(emf)"))
add(S4, "Biomass associations / Fig 5", "EMF x drought (naive)", lmt(fb_nv4, "treatmentDrought:scale(emf)"), 0.003, pcol(fb_nv4, "treatmentDrought:scale(emf)"))
add(S4, "Biomass associations", "AMF x drought", lmt(fb_pf4, "scale(amf):treatmentDrought"), 0.617, pcol(fb_pf4, "scale(amf):treatmentDrought"))
add(S4, "Biomass associations", "DSE x drought", lmt(fb_pf4, "treatmentDrought:scale(dse)"), 0.984, pcol(fb_pf4, "treatmentDrought:scale(dse)"))
add(S4, "Fig 5 caption", "Without SFN117", lmt(fb_pf4_x, "treatmentDrought:scale(emf)"), NA, NA, "EMF x drought refit without the |z| > 3 tree")

S5 <- "Results: Stomatal conductance (log scale; n = 203)"
add(S5, "Stomatal conductance", "Linear ancestry", lmt(poro_log, "per_fc"), 0.0005, pcol(poro_log, "per_fc"))
add(S5, "Stomatal conductance", "Quadratic ancestry", lmt(poro_log, "I(per_fc^2)"), 0.043, pcol(poro_log, "I(per_fc^2)"))
add(S5, "Stomatal conductance", "Drought (at 50% ancestry)", lmt(poro_log, "treatmentDrought"), 0.072, pcol(poro_log, "treatmentDrought"))
add(S5, "Stomatal conductance", "Ancestry x drought", lmt(poro_log, "per_fc:treatmentDrought"), 0.073, pcol(poro_log, "per_fc:treatmentDrought"))
inoc <- coef(summary(poro_log))[grep("myco_cat", rownames(coef(summary(poro_log)))), ]
add(S5, "Stomatal conductance", "Inoculation category", sprintf("all P > %.2f (t df %s%s%s)", floor(min(inoc[, "Pr(>|t|)"]) * 100) / 100, dff(min(inoc[, "df"])), EN, dff(max(inoc[, "df"]))), NA, NA)
add(S5, "Stomatal conductance", "Slopes (back-transformed)", sprintf("well-watered: %s%% per 10 points (95%% CI %s%s%s%%); drought: %s%% (%s%s%s%%)",
    dec(slope_tab[1, "pct_per_10"], 1), dec(slope_tab[1, "lo"], 1), " to ", dec(slope_tab[1, "hi"], 1),
    dec(slope_tab[2, "pct_per_10"], 1), dec(slope_tab[2, "lo"], 1), " to ", dec(slope_tab[2, "hi"], 1)), NA, NA, "negative values = decline")

S6 <- "Supplemental figure captions"
fa <- as.data.frame(sf2_a)
add(S6, "Supp Fig 2", "Inoculum / soil / interaction (mixed model, Type II, Satterthwaite)", sprintf("inoculum %s, %s; soil %s, %s; inoculum %s fungal type %s, %s",
    Fstr(fa["myco_cat", "F value"], fa["myco_cat", "NumDF"], fa["myco_cat", "DenDF"]), P(fa["myco_cat", "Pr(>F)"]),
    Fstr(fa["soil", "F value"], fa["soil", "NumDF"], fa["soil", "DenDF"]), P(fa["soil", "Pr(>F)"]), TIMES,
    Fstr(fa["myco_cat:fungal_type", "F value"], fa["myco_cat:fungal_type", "NumDF"], fa["myco_cat:fungal_type", "DenDF"]), P(fa["myco_cat:fungal_type", "Pr(>F)"])), 0.03, fa["myco_cat", "Pr(>F)"])
add(S6, "Supp Fig 3", "Drought PERMANOVA (whole-genotype permutations)", sprintf("PERMANOVA: %s, %s (tree-level %s)", Fstr(G(eng_foc, "treatment")$F, 1, rdf_foc), P(G(eng_foc, "treatment")$p_genotype_level), P(TA["treatment", "Pr(>F)"])), 0.784, G(eng_foc, "treatment")$p_genotype_level)

tab <- bind_rows(rows)
lines <- c("# Formatted results: main models with degrees of freedom (AJB format)", "",
           "Generated by `r_code/02c_main_models_df.R`. The manuscript is not changed.",
           paste0("Format: capital P; en dash (", EN, ") for minus signs; t(df), F(df1, df2) and ", CHI, "(df) written inline so subscripts can be applied in Word. ",
                  "t df are Satterthwaite (lmerTest) unless noted. z tests (binomial GLMM/GLM) have no df. ",
                  "PERMANOVA F uses the term df and the residual df of the full model."),
           "Labels: (G) genotype-level permutation test, (T) tree-level permutation test, (M) PERMANOVA on genotype means. Rows marked **SIGNIFICANCE CHANGED** cross P = 0.05 relative to the current manuscript.", "")
for (s in unique(tab$section)) {
  lines <- c(lines, paste0("## ", s), "", "| Location | Statistic | Paste-ready string | Flag / note |", "|---|---|---|---|")
  for (i in which(tab$section == s)) lines <- c(lines, sprintf("| %s | %s | %s | %s |", tab$location[i], tab$statistic[i],
      gsub("|", "\\|", tab$string[i], fixed = TRUE), paste(c(if (nzchar(tab$flag[i])) tab$flag[i], if (nzchar(tab$note[i])) tab$note[i]), collapse = "; ")))
  lines <- c(lines, "")
}
con <- file("results/02_main_models/formatted_results_main_models_df.md", "w", encoding = "UTF-8"); writeLines(enc2utf8(lines), con); close(con)
hdr("1. formatted_results_main_models_df.md written")
cat(nrow(tab), "strings;", sum(nzchar(tab$flag)), "flagged\n")

# ------------------------------------------------------------------------------
# 2. Aggregated survival: grouping sensitivity (stats::glm)
# ------------------------------------------------------------------------------
hdr("2. Aggregated survival sensitivity (AMF, EMF, DSE x drought)")
scored <- function(d) !is.na(d$amf) & !is.na(d$emf) & !is.na(d$dse)
build <- function(gvars) {
  cells <- ds %>% group_by(across(all_of(gvars))) %>%
    summarise(n_total = n(), survival_rate = sum(alive_01_10_2021, na.rm = TRUE) / n(),
              prop_tam = mean(soil == "Tam Legacy"), prop_AM = mean(myco_cat == "AM"), prop_EM = mean(myco_cat == "EM"),
              prop_Both = mean(myco_cat == "Both"), .groups = "drop") %>% filter(n_total >= 3)
  colz <- ds %>% filter(alive_01_10_2021 == 1) %>% mutate(sc = scored(.)) %>% group_by(across(all_of(gvars))) %>%
    summarise(n_scored = sum(sc), mean_emf = mean(emf, na.rm = TRUE), mean_amf = mean(amf, na.rm = TRUE), mean_dse = mean(dse, na.rm = TRUE),
              mean_per_f = mean(per_f, na.rm = TRUE), mean_drc_mm = mean(og_drc_mm, na.rm = TRUE), .groups = "drop")
  inner_join(cells, colz, by = gvars) %>% mutate(mean_per_fc = mean_per_f - 50) %>%
    filter(is.finite(mean_emf), is.finite(mean_amf), is.finite(mean_dse))
}
fit_agg <- function(d, covs, label) {
  base <- tibble(analysis = label, n_groups = nrow(d), n_trees = sum(d$n_total), scored_median = median(d$n_scored),
                 scored_range = if (nrow(d)) paste(range(d$n_scored), collapse = EN) else "", groups_by_treatment = paste(table(d$treatment), collapse = "/"))
  f <- as.formula(paste("survival_rate ~ (scale(mean_emf) + scale(mean_dse) + scale(mean_amf)) * treatment + mean_per_fc + I(mean_per_fc^2) + mean_drc_mm +", covs))
  npar <- NA
  ok <- nrow(d) > 0 && n_distinct(d$treatment) == 2
  if (ok) { npar <- ncol(model.matrix(f, d)); ok <- nrow(d) >= npar + 5 }
  if (!ok) return(bind_cols(base, tibble(note = sprintf("not estimable: %d groups (model has %s parameters; both treatments needed)",
                                                         nrow(d), ifelse(is.na(npar), "15+", npar)))))
  m <- glm(f, data = d, weights = n_total, family = binomial); cs <- coef(summary(m))
  mq <- glm(f, data = d, weights = n_total, family = quasibinomial); cq <- coef(summary(mq))   # SEs scaled by estimated dispersion
  g <- function(t) sprintf("%s = %s, %s, %s", BETA, num(cs[t, "Estimate"]), zstr(cs[t, "z value"]), P(cs[t, "Pr(>|z|)"]))
  q <- function(t) sprintf("%s, %s", tstr(cq[t, "t value"], df.residual(mq)), P(cq[t, "Pr(>|t|)"]))
  bind_cols(base, tibble(AMF_x_drought = g("scale(mean_amf):treatmentDrought"), EMF_x_drought = g("scale(mean_emf):treatmentDrought"),
                         DSE_x_drought = g("scale(mean_dse):treatmentDrought"),
                         AMF_quasi = q("scale(mean_amf):treatmentDrought"), EMF_quasi = q("scale(mean_emf):treatmentDrought"),
                         DSE_quasi = q("scale(mean_dse):treatmentDrought"), AMF_P = cs["scale(mean_amf):treatmentDrought", "Pr(>|z|)"],
                         dispersion = round(sum(residuals(m, type = "pearson")^2) / df.residual(m), 2), n_par = length(coef(m)), note = ""))
}
cur <- build(c("pop", "new_class", "treatment", "myco_cat", "soil"))
coa <- build(c("pop", "new_class", "treatment"))
sens <- bind_rows(
  fit_agg(cur, "soil + myco_cat", "current grouping (pop x cross type x drought x inoculum x soil), all groups"),
  fit_agg(filter(cur, n_scored >= 2), "soil + myco_cat", "current grouping, >= 2 scored plants"),
  fit_agg(filter(cur, n_scored >= 3), "soil + myco_cat", "current grouping, >= 3 scored plants"),
  fit_agg(coa, "prop_tam + prop_AM + prop_EM + prop_Both", "(a) pop x cross type x drought; inoculum and soil as group proportions"),
  fit_agg(filter(coa, n_scored >= 2), "prop_tam + prop_AM + prop_EM + prop_Both", "(a) with >= 2 scored plants"),
  fit_agg(filter(coa, n_scored >= 3), "prop_tam + prop_AM + prop_EM + prop_Both", "(a) with >= 3 scored plants"))
cat("'scored plants' = survivors with AMF, EMF and DSE all scored. Colonization predictors are z-scored over the groups in each model.\n",
    "The 'current grouping, all groups' row therefore differs slightly from the reported model, which z-scored over all 83 candidate groups.\n\n", sep = "")
print(as.data.frame(sens %>% select(-AMF_P)), right = FALSE)
write.csv(sens, "results/02_main_models/survival_sensitivity.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. dry_shoottotal: total aboveground or stem only?
# ------------------------------------------------------------------------------
hdr("4. dry_shoottotal (dry_stem in mega_ds.csv): total aboveground or stem only?")
raw <- read_csv("data/raw_data/gh_drybiomass_aboveground.csv", show_col_types = FALSE)
cat("Raw columns:", paste(names(raw), collapse = ", "), "\n")
cat("Metadata definition (metadata_ghcottieproject_summer2021.xlsx, sheet gh_drybiomass_aboveground):\n",
    "  dry_shoottotal = 'Mass (g) of dry stem biomass as defined as no leaves or roots'\n",
    "  wet_stem       = 'Mass (g) of wet stem as defined as no leaves or roots'\n",
    "  dryleaf_weight = 'Mass (g) of dry leaf biomass'; dryshoot_new = 'Mass (g) of NEW biomass as measured from bud scars'\n", sep = "")
b <- ds %>% filter(!is.na(dry_stem)) %>% select(tree_idold, wet_stem, dry_stem, dryleaf_weight, dryshoot_new)
bb <- b %>% filter(!is.na(wet_stem), wet_stem > 0)
cat("\nRows with dry_shoottotal:", nrow(b), "| with wet_stem too:", nrow(bb), "\n")
cat("dry_shoottotal / wet_stem: median", round(median(bb$dry_stem / bb$wet_stem), 3), "IQR",
    paste(round(quantile(bb$dry_stem / bb$wet_stem, c(0.25, 0.75)), 3), collapse = "-"),
    "| rows where dry_shoottotal > wet_stem:", sum(bb$dry_stem > bb$wet_stem), "\n")
bl <- bb %>% filter(!is.na(dryleaf_weight))
cat("(dry_shoottotal - dryleaf) / wet_stem (the stem dry-matter fraction IF leaves were included): median",
    round(median((bl$dry_stem - bl$dryleaf_weight) / bl$wet_stem), 3), "\n")
m1 <- lm(dry_stem ~ wet_stem + dryleaf_weight, data = bl)
cat("lm(dry_shoottotal ~ wet_stem + dryleaf_weight): leaf coefficient =", round(coef(m1)["dryleaf_weight"], 3),
    "(SE", round(coef(summary(m1))["dryleaf_weight", 2], 3), "); about 1 expected if leaves were included, about 0 if not;",
    "wet_stem coefficient =", round(coef(m1)["wet_stem"], 3), "\n")
bc <- b %>% filter(!is.na(dryleaf_weight))
cat("Rows where dry_shoottotal < dryleaf_weight:", sum(bc$dry_stem < bc$dryleaf_weight), "of", nrow(bc), "\n")
cat("Rows where dryshoot_new > dry_shoottotal:", sum(b$dryshoot_new > b$dry_stem, na.rm = TRUE), "| dryshoot_new > dryleaf_weight:",
    sum(bc$dryshoot_new > bc$dryleaf_weight, na.rm = TRUE), "of", sum(!is.na(bc$dryshoot_new)), "\n")
cat("Median dry leaf / (leaf + dry_shoottotal):", round(median(bc$dryleaf_weight / (bc$dryleaf_weight + bc$dry_stem)), 3), "\n")
sink()
cat("Done. See results/02_main_models/main_models_df_output.txt, results/02_main_models/formatted_results_main_models_df.md, results/02_main_models/survival_sensitivity.csv\n")
