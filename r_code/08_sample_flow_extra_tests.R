# ==============================================================================
# SAMPLE FLOW AND EXTRA TESTS: sample flow (plants per analysis), survival by group, inoculum effects, EMF simple slopes.
# Refits the reported models (02a_main_models.R structure), checks each reproduces the
# reported numbers, then runs the extra tests. Changes no manuscript result.
# ==============================================================================
# Run from the project root (R 4.6.1):  Rscript r_code/08_sample_flow_extra_tests.R
# Writes to results/08_sample_flow_extra_tests/: output.txt, flow_table.csv,
#   survival_by_group.csv, inoculum_effects.csv, emf_simple_slopes.csv
# The parametric bootstrap for item 5 is NOT rerun here; it is in 02a_main_models.R
# (output in results/02_main_models/main_models_output.txt).
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(lme4); library(lmerTest); library(emmeans); library(tidyverse)
out <- "results/08_sample_flow_extra_tests"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("SAMPLE FLOW AND EXTRA TESTS run:", "|", R.version.string, "| lme4", as.character(packageVersion("lme4")),
    "| lmerTest", as.character(packageVersion("lmerTest")), "| emmeans", as.character(packageVersion("emmeans")), "\n")

EN <- "–"; BETA <- "β"; TIMES <- "×"; CHI <- "χ²"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
tstr <- function(b, t, df, p) sprintf("%s = %s, t = %s, df = %s, %s", BETA, num(b), dec(t), dff(df), P(p))
zstr <- function(b, z, p) sprintf("%s = %s, z = %s, %s", BETA, num(b), dec(z), P(p))
check <- function(label, got, want, tol = 0.0005) {
  ok <- abs(got - want) <= tol * max(1, abs(want))
  cat(sprintf("  CHECK %-55s got %s | reported %s | %s\n", label, signif(got, 4), signif(want, 4), ifelse(ok, "OK", "MISMATCH")))
  invisible(ok) }

ds <- load_and_prepare_data() %>%
  mutate(per_fc = per_f - 50, pf50 = per_fc / 50, geno_id = paste(pop, type, geno, sep = "_"),
         has_amf = !is.na(amf), has_emf = !is.na(emf), has_dse = !is.na(dse), all3 = has_amf & has_emf & has_dse,
         scored_any = has_amf | has_emf | has_dse,
         scoring = case_when(all3 ~ "all three", has_emf & !has_amf & !has_dse ~ "EMF only",
                             has_amf & has_dse & !has_emf ~ "AMF and DSE only", scored_any ~ "other", TRUE ~ "not scored"),
         alive = alive_01_10_2021)
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
cat("myco_cat levels (first = reference):", paste(levels(factor(ds$myco_cat)), collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# Reported models (as 02a_main_models.R / 02b_main_models_followup.R / 04_conductance_surfaces.R)
# ------------------------------------------------------------------------------
surv_red <- glmer(alive ~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat + (1 | geno_id),
                  data = ds, family = binomial, control = ctl)
bio_red <- lmerTest::lmer(log(dryshoot_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil +
                            (1 | pop) + (1 | pop:geno_id), data = ds)
fb <- lmerTest::lmer(log(dryshoot_new) ~ scale(amf) * treatment + scale(emf) * treatment + scale(dse) * treatment +
                       per_fc + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
dl <- ds %>% filter(all3) %>%
  select(tree_id, pop, geno_id, per_f, new_class, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>% na.omit() %>%
  mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm))) %>%
  pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>%
  mutate(fungal_type = factor(fungal_type, levels = c("amf", "dse", "emf")), myco_cat = factor(myco_cat))
col <- lmerTest::lmer(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil +
                        (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id), data = dl)
poro <- lmerTest::lmer(log(por_avg) ~ per_fc * treatment + I(per_fc^2) + og_drc_mm + soil + myco_cat +
                         (1 | pop) + (1 | pop:geno_id), data = ds %>% filter(por_avg > 0))
side_mean <- function(d, side) { cc <- as.vector(outer(paste0("l", 1:3, "_", side), c("wet", "dry"), paste, sep = "_"))
  m <- rowMeans(as.matrix(sapply(d[, cc], as.numeric)), na.rm = TRUE); m[is.nan(m)] <- NA; m }
por <- ds %>% filter(!is.na(por_avg)) %>% mutate(adaxial = side_mean(., "t"), abaxial = side_mean(., "b"))
fit_side <- function(y) lmerTest::lmer(as.formula(paste0("log(", y, ") ~ per_fc * treatment + I(per_fc^2) + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id)")),
                                       data = por %>% filter(!is.na(.data[[y]]), .data[[y]] > 0))
ad <- fit_side("adaxial"); ab <- fit_side("abaxial")
cs <- function(m) coef(summary(m))
cat("\nReproduction checks against formatted_results_main_models_df.md:\n")
check("survival linear ancestry (per 1%)", cs(surv_red)["pf50", "Estimate"] / 50, -0.0112, 0.01)
check("survival drought", cs(surv_red)["treatmentDrought", "Estimate"], -1.205, 0.001)
check("biomass linear ancestry", cs(bio_red)["per_fc", "Estimate"], -0.00528, 0.001)
check("biomass quadratic ancestry", cs(bio_red)["I(per_fc^2)", "Estimate"], -0.000145, 0.005)
check("EMF x drought (colonization-biomass)", cs(fb)["treatmentDrought:scale(emf)", "Estimate"], -0.445, 0.002)
check("Fremont ancestry x EMF (colonization)", cs(col)["per_f:fungal_typeemf", "Estimate"], -0.171, 0.003)
check("conductance linear ancestry", cs(poro)["per_fc", "Estimate"], -0.00851, 0.002)
cat("  n: survival", nobs(surv_red), "| biomass", nobs(bio_red), "| colonization plants", n_distinct(dl$tree_id),
    "| colonization-biomass", nobs(fb), "| conductance", nobs(poro), "| adaxial", nobs(ad), "| abaxial", nobs(ab), "\n")

# ------------------------------------------------------------------------------
hdr("1. Flow of plants through the experiment and the models")
# ------------------------------------------------------------------------------
cat("Cuttings in the data set (planted):", nrow(ds), "| unique plant IDs:", n_distinct(ds$tree_idold), "\n")
cat("Survival status (alive_01_10_2021): alive", sum(ds$alive == 1, na.rm = TRUE), "| dead", sum(ds$alive == 0, na.rm = TRUE),
    "| missing", sum(is.na(ds$alive)), "\n")
cat("With new stem mass (dryshoot_new): survivors", sum(ds$alive == 1 & !is.na(ds$dryshoot_new), na.rm = TRUE),
    "| recorded dead", sum(ds$alive == 0 & !is.na(ds$dryshoot_new), na.rm = TRUE), "| status missing", sum(is.na(ds$alive) & !is.na(ds$dryshoot_new)), "\n")
dead_bio <- ds %>% filter(alive == 0, !is.na(dryshoot_new))
if (nrow(dead_bio)) { cat("Plants recorded dead on 1 Oct 2021 but with new stem mass (in the biomass model if covariates complete):\n")
  print(as.data.frame(dead_bio %>% select(tree_idold, new_class, treatment, dryshoot_new, dryleaf_weight, root_dryweight, notes)), row.names = FALSE) }
cat("\nColonization scoring pattern, by survival status:\n"); print(table(scoring = ds$scoring, alive = ds$alive, useNA = "ifany"))
cat("\nPlants with each guild scored (any status): AMF", sum(ds$has_amf), "| EMF", sum(ds$has_emf), "| DSE", sum(ds$has_dse), "\n")
s2 <- ds %>% filter(!is.na(myco_cat)) %>% summarise(AMF = sum(has_amf), EMF = sum(has_emf), DSE = sum(has_dse))
cat("Appendix S2 plot data (guild scored and inoculum known): AMF", s2$AMF, "| EMF", s2$EMF, "| DSE", s2$DSE, "\n")

# first failing criterion = reason for exclusion
reasons <- function(d, crit) {
  why <- rep(NA_character_, nrow(d))
  for (nm in names(crit)) { fail <- is.na(why) & !crit[[nm]](d); why[fail] <- nm }
  why }
cc <- function(v) function(d) rowSums(is.na(as.matrix(d[, v]))) == 0
alloc <- read.csv("results/06_allocation_models/allocation_data.csv")
flow <- list()
add_flow <- function(model, n_rep, pool_desc, d, crit, n_model) {
  w <- reasons(d, crit); tab <- table(factor(w, levels = names(crit)))
  stopifnot(sum(is.na(w)) == n_model)
  flow[[length(flow) + 1]] <<- tibble(model = model, reported_n = n_rep, model_n = n_model, pool = pool_desc, pool_n = nrow(d),
                                      dropped = sum(!is.na(w)),
                                      drop_reasons = paste(sprintf("%s: %d", names(tab)[tab > 0], tab[tab > 0]), collapse = "; "))
  cat(sprintf("\n%s: reported n = %d | refit n = %d | pool = %s (%d)\n", model, n_rep, n_model, pool_desc, nrow(d)))
  print(tab[tab > 0])
  w }
w <- add_flow("Survival (binomial GLMM)", 563, "all cuttings", ds,
              list(`survival status missing` = cc("alive"), `initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), nobs(surv_red))
w <- add_flow("New stem biomass (LMM)", 429, "all cuttings", ds,
              list(`died or survival status missing` = function(d) d$alive %in% 1, `survivor without new stem mass` = cc("dryshoot_new"), `initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), nobs(bio_red))
harv <- ds %>% filter(alive == 1, !is.na(dryshoot_new))
w <- add_flow("Allocation, (a) estimated root mass", 414, "harvested survivors (alive, new stem mass)", harv,
              list(`root mass missing (scored: no fresh mass; unscored: no dry mass)` = function(d) ifelse(d$scored_any, !is.na(d$wet_root), !is.na(d$root_dryweight)),
                   `initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), nrow(alloc))
w <- add_flow("Allocation, (b) measured root mass only", 242, "harvested survivors (alive, new stem mass)", harv,
              list(`colonization-scored (root subsample removed)` = function(d) !d$scored_any,
                   `root dry mass missing` = cc("root_dryweight"), `initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), sum(!alloc$scored))
w <- add_flow("Leaf drop (binomial GLMM)", 425, "harvested survivors (alive, new stem mass)", harv,
              list(`initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), 425)
w <- add_flow("Colonization (LMM)", 78, "all cuttings", ds,
              list(`not scored for all three guilds` = function(d) d$all3, `initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), n_distinct(dl$tree_id))
w <- add_flow("Colonization-biomass (LMM)", 77, "plants with all three guilds scored", ds %>% filter(all3),
              list(`no new stem mass` = cc("dryshoot_new"), `initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), nobs(fb))
if (nobs(fb) < 78) { cat("Scored plant(s) not in the colonization-biomass model:\n")
  print(as.data.frame(ds %>% filter(all3, is.na(dryshoot_new) | is.na(og_drc_mm)) %>% select(tree_idold, alive, dryshoot_new, og_drc_mm, notes)), row.names = FALSE) }
w <- add_flow("Conductance, whole leaf (LMM)", 203, "plants with a porometer mean", por,
              list(`non-positive mean (cannot be logged)` = function(d) d$por_avg > 0, `initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), nobs(poro))
w <- add_flow("Conductance, adaxial (LMM)", 203, "plants with a porometer mean", por,
              list(`no adaxial reading` = function(d) !is.na(d$adaxial), `non-positive adaxial mean` = function(d) d$adaxial > 0,
                   `initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), nobs(ad))
w <- add_flow("Conductance, abaxial (LMM)", 200, "plants with a porometer mean", por,
              list(`no abaxial reading` = function(d) !is.na(d$abaxial), `non-positive abaxial mean` = function(d) d$abaxial > 0,
                   `initial cutting diameter missing` = cc("og_drc_mm"),
                   `other covariate missing` = cc(c("per_f", "treatment", "pop", "soil", "myco_cat", "geno_id"))), nobs(ab))
flow <- bind_rows(flow); readr::write_excel_csv(flow, file.path(out, "flow_table.csv"))
bm <- ds %>% filter(!is.na(dryshoot_new), !is.na(og_drc_mm))
cat("
Biomass model rows by survival status (alive_01_10_2021):
"); print(table(status = bm$alive, useNA = "ifany"))
print(as.data.frame(bm %>% filter(is.na(alive) | alive == 0) %>% select(tree_idold, alive, new_class, treatment, dryshoot_new, dryleaf_weight, root_dryweight, scoring)), row.names = FALSE)
cat("\nPorometer plants: por_avg present", nrow(por), "| non-positive por_avg:", paste(por$tree_idold[por$por_avg <= 0], round(por$por_avg[por$por_avg <= 0], 4), collapse = ", "), "\n")

# ------------------------------------------------------------------------------
hdr("2. Survival (%) by drought x cross type and by inoculum x drought (survival model sample)")
# ------------------------------------------------------------------------------
sd_ <- ds %>% filter(if_all(c(alive, pf50, treatment, pop, og_drc_mm, soil, myco_cat, geno_id), ~ !is.na(.x)))
stopifnot(nrow(sd_) == nobs(surv_red))
stab <- function(d, g) d %>% group_by(across(all_of(g))) %>% summarise(planted = n(), survived = sum(alive), .groups = "drop") %>%
  mutate(survival_pct = round(100 * survived / planted, 1))
s_ct <- stab(sd_, c("new_class", "treatment")); s_in <- stab(sd_, c("myco_cat", "treatment")); s_tr <- stab(sd_, "treatment")
cat("By treatment:\n"); print(as.data.frame(s_tr), row.names = FALSE)
cat("\nBy cross type x drought:\n"); print(as.data.frame(s_ct), row.names = FALSE)
cat("\nBy inoculum x drought:\n"); print(as.data.frame(s_in), row.names = FALSE)
readr::write_excel_csv(bind_rows(s_tr %>% mutate(table = "treatment"), s_ct %>% mutate(table = "cross type x drought"), s_in %>% mutate(table = "inoculum x drought")) %>%
            relocate(table), file.path(out, "survival_by_group.csv"))

# ------------------------------------------------------------------------------
hdr("3. Inoculation effects on performance (reduced models; reference = control)")
# ------------------------------------------------------------------------------
inoc_rows <- list()
inoc_lab <- c(myco_catAM = "AMF vs. control", myco_catEM = "EMF vs. control", myco_catBoth = "Dual vs. control")
# survival (GLMM, Wald z; LRTs)
s_noinoc <- update(surv_red, . ~ . - myco_cat)
s_int <- update(surv_red, . ~ . + myco_cat:treatment)
l_main <- anova(s_noinoc, surv_red); l_int <- anova(surv_red, s_int)
cs_s <- cs(surv_red)
cat("Survival (n =", nobs(surv_red), ")\n"); print(round(cs_s[names(inoc_lab), ], 4))
cat(sprintf("  inoculum main effect LRT %s(%d) = %.2f, %s | inoculum %s drought LRT %s(%d) = %.2f, %s\n",
            CHI, l_main$Df[2], l_main$Chisq[2], P(l_main$`Pr(>Chisq)`[2]), TIMES, CHI, l_int$Df[2], l_int$Chisq[2], P(l_int$`Pr(>Chisq)`[2])))
inoc_rows[["surv"]] <- tibble(response = "Survival", n = nobs(surv_red), term = c(inoc_lab, "Inoculum (main effect)", paste("Inoculum", TIMES, "drought")),
  string = c(sapply(names(inoc_lab), function(t) zstr(cs_s[t, 1], cs_s[t, 3], cs_s[t, 4])),
             sprintf("LRT %s(%d) = %s, %s", CHI, l_main$Df[2], dec(l_main$Chisq[2]), P(l_main$`Pr(>Chisq)`[2])),
             sprintf("LRT %s(%d) = %s, %s", CHI, l_int$Df[2], dec(l_int$Chisq[2]), P(l_int$`Pr(>Chisq)`[2]))),
  p = c(cs_s[names(inoc_lab), 4], l_main$`Pr(>Chisq)`[2], l_int$`Pr(>Chisq)`[2]))
# LMMs: coefficient t tests; Satterthwaite F tests for the inoculum term and inoculum x drought (added to the reduced model)
lmm_inoc <- function(m, label) {
  m_int <- update(m, . ~ . + myco_cat:treatment)
  a <- anova(m, type = 3)["myco_cat", ]; ai <- anova(m_int, type = 3)["treatment:myco_cat", ]
  c1 <- cs(m)
  cat(sprintf("\n%s (n = %d)\n", label, nobs(m))); print(round(c1[names(inoc_lab), ], 4))
  cat(sprintf("  inoculum F(%d, %.1f) = %.2f, %s | inoculum %s drought F(%d, %.1f) = %.2f, %s\n", a$NumDF, a$DenDF, a$`F value`, P(a$`Pr(>F)`),
              TIMES, ai$NumDF, ai$DenDF, ai$`F value`, P(ai$`Pr(>F)`)))
  ctr <- as.data.frame(summary(contrast(emmeans(m_int, ~ myco_cat | treatment, lmer.df = "satterthwaite"), "trt.vs.ctrl"), infer = c(TRUE, TRUE)))
  cat("  Inoculum vs. control within each watering treatment (interaction model; ratio = exp(log-scale difference); CIs and P Dunnett-adjusted):\n")
  print(ctr %>% transmute(treatment, contrast, ratio = round(exp(estimate), 3), lower = round(exp(lower.CL), 3), upper = round(exp(upper.CL), 3),
                          t = round(t.ratio, 2), df = round(df, 1), p_dunnett = round(p.value, 4)), row.names = FALSE)
  tibble(response = label, n = nobs(m), term = c(inoc_lab, "Inoculum (main effect)", paste("Inoculum", TIMES, "drought")),
         string = c(sapply(names(inoc_lab), function(t) tstr(c1[t, "Estimate"], c1[t, "t value"], c1[t, "df"], c1[t, "Pr(>|t|)"])),
                    sprintf("F = %s, df = %d, %s, %s", dec(a$`F value`), a$NumDF, dff(a$DenDF), P(a$`Pr(>F)`)),
                    sprintf("F = %s, df = %d, %s, %s", dec(ai$`F value`), ai$NumDF, dff(ai$DenDF), P(ai$`Pr(>F)`))),
         p = c(c1[names(inoc_lab), "Pr(>|t|)"], a$`Pr(>F)`, ai$`Pr(>F)`)) }
inoc_rows[["bio"]] <- lmm_inoc(bio_red, "New stem biomass (log)")
ta <- ds %>% inner_join(alloc %>% select(tree_idold, total_new, alloc_scored = scored), by = "tree_idold") %>% mutate(log_total_new = log(total_new))
stopifnot(nrow(ta) == 414)
tot_a <- lmerTest::lmer(log_total_new ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ta)
tot_b <- update(tot_a, data = filter(ta, !alloc_scored))
check("total new biomass (a) drought", cs(tot_a)["treatmentDrought", "Estimate"], -0.220, 0.002)
inoc_rows[["tot_a"]] <- lmm_inoc(tot_a, "Total new biomass (log), (a) estimated root mass")
inoc_rows[["tot_b"]] <- lmm_inoc(tot_b, "Total new biomass (log), (b) measured root mass only")
inoc <- bind_rows(inoc_rows); readr::write_excel_csv(inoc, file.path(out, "inoculum_effects.csv"))
cat("\nAJB strings:\n"); print(as.data.frame(inoc %>% select(response, term, string)), right = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
hdr("4. EMF x drought robustness: ancestry x drought added; EMF simple slopes")
# ------------------------------------------------------------------------------
# scale() inside the reported formula is evaluated on every row of ds with the guild scored (before
# rows with other missing values are dropped), so 1 SD = SD over all plants with that guild scored.
zsd <- c(amf = sd(ds$amf, na.rm = TRUE), emf = sd(ds$emf, na.rm = TRUE), dse = sd(ds$dse, na.rm = TRUE))
cat("SD used for the reported 'per SD' scaling (all plants with the guild scored): AMF", round(zsd["amf"], 2), "(n =", sum(ds$has_amf), ") | EMF",
    round(zsd["emf"], 2), "(n =", sum(ds$has_emf), ") | DSE", round(zsd["dse"], 2), "(n =", sum(ds$has_dse), ")\n")
fbd <- ds %>% mutate(amf_z = (amf - mean(amf, na.rm = TRUE)) / zsd["amf"], dse_z = (dse - mean(dse, na.rm = TRUE)) / zsd["dse"],
                     emf_z = (emf - mean(emf, na.rm = TRUE)) / zsd["emf"], emf10 = emf / 10)
cat("SD of EMF among the 77 model plants:", round(sd(fbd$emf[!is.na(fbd$dryshoot_new) & fbd$all3 & !is.na(fbd$og_drc_mm)]), 2), "\n")
f_rep <- log(dryshoot_new) ~ amf_z * treatment + emf_z * treatment + dse_z * treatment + per_fc + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id)
f_rob <- update(f_rep, . ~ . + per_fc:treatment + I(per_fc^2):treatment)
m_rep <- lmerTest::lmer(f_rep, data = fbd); m_rob <- lmerTest::lmer(f_rob, data = fbd)
check("manual z-scoring reproduces EMF x drought", cs(m_rep)["treatmentDrought:emf_z", "Estimate"], -0.445, 0.002)
emf_rows <- list()
for (nm in c("reported", "with ancestry x drought")) {
  m <- if (nm == "reported") m_rep else m_rob
  r <- cs(m)["treatmentDrought:emf_z", ]
  cat(sprintf("\n%s model (n = %d): EMF %s drought %s\n", nm, nobs(m), TIMES, tstr(r["Estimate"], r["t value"], r["df"], r["Pr(>|t|)"])))
  if (nm != "reported") { for (t in grep("per_fc.*Drought|Drought.*per_fc", rownames(cs(m)), value = TRUE)) { r2 <- cs(m)[t, ]
    cat(sprintf("  %s: %s\n", t, tstr(r2["Estimate"], r2["t value"], r2["df"], r2["Pr(>|t|)"]))) }
    lr <- anova(m_rep, m_rob); cat(sprintf("  LRT for the two ancestry x drought terms: %s(%d) = %.2f, %s\n", CHI, lr$Df[2], lr$Chisq[2], P(lr$`Pr(>Chisq)`[2]))) }
  # simple slopes per 10 percentage points of EMF (emf_z = emf/SD, so slope per 10 pp = slope per SD x 10/SD)
  m10 <- lmerTest::lmer(update(formula(m), . ~ . - emf_z - emf_z:treatment + emf10 * treatment), data = fbd)
  stopifnot(abs(cs(m10)["treatmentDrought:emf10", "t value"] - cs(m)["treatmentDrought:emf_z", "t value"]) < 1e-4)
  et <- summary(emtrends(m10, ~ treatment, var = "emf10", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE))
  et <- as.data.frame(et) %>% transmute(model = nm, treatment, slope_log_per_10pp = emf10.trend, SE, df, lower = lower.CL, upper = upper.CL, t = t.ratio, p = p.value,
                                       pct_change = 100 * (exp(slope_log_per_10pp) - 1), pct_lo = 100 * (exp(lower) - 1), pct_hi = 100 * (exp(upper) - 1))
  print(et %>% mutate(across(where(is.numeric), ~ round(., 4))))
  emf_rows[[nm]] <- et }
pc <- function(x) dash(sprintf("%.1f", x))
emf_tab <- bind_rows(emf_rows) %>% mutate(string = sprintf("%s%% per 10 percentage points of EMF (95%% CI %s to %s%%), t = %s, df = %s, %s",
                                                           pc(pct_change), pc(pct_lo), pc(pct_hi), dec(t), sapply(df, dff), sapply(p, P)))
readr::write_excel_csv(emf_tab, file.path(out, "emf_simple_slopes.csv"))
cat("\n"); print(as.data.frame(emf_tab %>% select(model, treatment, string)), right = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
hdr("5. Biomass curvature: fitted peak of the quadratic ancestry curve")
# ------------------------------------------------------------------------------
bb <- fixef(bio_red)["per_fc"]; cq <- fixef(bio_red)["I(per_fc^2)"]
vtx <- unname(50 - bb / (2 * cq))
V <- as.matrix(vcov(bio_red))[c("per_fc", "I(per_fc^2)"), c("per_fc", "I(per_fc^2)")]
gr <- c(-1 / (2 * cq), bb / (2 * cq^2)); se_v <- sqrt(as.numeric(t(gr) %*% V %*% gr))
cat(sprintf("Vertex = 50 - b/(2c) = %.1f%% Fremont ancestry (b = %s, c = %s); delta-method 95%% CI %.1f to %.1f%%\n",
            vtx, signif(bb, 4), signif(cq, 4), vtx - 1.96 * se_v, vtx + 1.96 * se_v))
cat("Parametric bootstrap (r_code/02a_main_models.R; result in results/02_main_models/main_models_output.txt, rerun under R 4.6.1):\n")
bl <- readLines("results/02_main_models/main_models_output.txt"); i <- grep("^Parametric bootstrap", bl); cat(paste0("  ", bl[i:(i + 3)]), sep = "\n")
pr <- emmeans(bio_red, ~ per_fc, at = list(per_fc = c(-50, vtx - 50, 50)), type = "response")
cat("\nPredicted new stem biomass (g; back-transformed, averaged over treatments and covariates) at 0%, the vertex, and 100% Fremont:\n"); print(summary(pr))

# ------------------------------------------------------------------------------
hdr("6. Colonization LMM: why the ancestry slopes have df = 215.6")
# ------------------------------------------------------------------------------
cat("Rows:", nrow(dl), "(", n_distinct(dl$tree_id), "plants x 3 guilds ) | genotypes:", n_distinct(dl$geno_id), "| populations:", n_distinct(dl$pop),
    "| fixed-effect parameters:", length(fixef(col)), "| residual df = rows - fixed =", nrow(dl) - length(fixef(col)), "\n")
cat("\nVariance components:\n"); vc <- as.data.frame(VarCorr(col)); print(vc %>% mutate(share = round(vcov / sum(vcov), 3)) %>% select(grp, vcov, sdcor, share))
cat("Singular fit:", isSingular(col), "\n")
sl_s <- summary(emtrends(col, ~ fungal_type, var = "per_f", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE))
sl_k <- summary(emtrends(col, ~ fungal_type, var = "per_f", lmer.df = "kenward-roger"), infer = c(TRUE, TRUE))
cat("\nAncestry slope by guild, Satterthwaite (as reported):\n"); print(sl_s)
cat("\nAncestry slope by guild, Kenward-Roger:\n"); print(sl_k)
col_g <- update(col, . ~ . + (1 | pop:geno_id:fungal_type))
cat("\nSensitivity: genotype-by-guild random intercept added, (1 | pop:geno_id:fungal_type):\n")
print(as.data.frame(VarCorr(col_g)) %>% select(grp, vcov, sdcor))
sl_g <- summary(emtrends(col_g, ~ fungal_type, var = "per_f", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE)); print(sl_g)
r <- cs(col_g)["per_f:fungal_typeemf", ]; cat("  ancestry x EMF interaction:", tstr(r["Estimate"], r["t value"], r["df"], r["Pr(>|t|)"]), "\n")
cat("\nGenotype-mean check: colonization averaged within genotype, regressed on ancestry (one point per genotype):\n")
gm <- dl %>% group_by(geno_id, pop, fungal_type) %>% summarise(per_f = mean(per_f), colonization = mean(colonization), n = n(), .groups = "drop")
gm_rows <- list()
for (g in c("amf", "dse", "emf")) { d <- filter(gm, fungal_type == g)
  m1 <- lm(colonization ~ per_f, data = d); m2 <- lmerTest::lmer(colonization ~ per_f + (1 | pop), data = d)
  r1 <- coef(summary(m1))["per_f", ]; r2 <- coef(summary(m2))["per_f", ]
  cat(sprintf("  %s: genotypes = %d | OLS %s | with population random intercept %s\n", toupper(g), nrow(d),
              tstr(r1[1], r1[3], df.residual(m1), r1[4]), tstr(r2["Estimate"], r2["t value"], r2["df"], r2["Pr(>|t|)"]))) }

# ------------------------------------------------------------------------------
hdr("7. Reconciling two pairs of numbers")
# ------------------------------------------------------------------------------
cat("7a. AMF x drought in the group-level survival GLM: -1.91 (text) vs -1.83 (Appendix S5)\n")
gv <- c("pop", "new_class", "treatment", "myco_cat", "soil")
psurv <- ds %>% group_by(across(all_of(gv))) %>% summarise(n_total = n(), survival_rate = sum(alive, na.rm = TRUE) / n(), .groups = "drop") %>% filter(n_total >= 3)
pcol <- ds %>% filter(alive == 1) %>% group_by(across(all_of(gv))) %>%
  summarise(mean_emf = mean(emf, na.rm = TRUE), mean_amf = mean(amf, na.rm = TRUE), mean_dse = mean(dse, na.rm = TRUE),
            mean_per_f = mean(per_f, na.rm = TRUE), mean_drc_mm = mean(og_drc_mm, na.rm = TRUE), .groups = "drop")
sdf <- inner_join(psurv, pcol, by = gv) %>% mutate(mean_per_fc = mean_per_f - 50)
fa <- survival_rate ~ (scale(mean_emf) + scale(mean_dse) + scale(mean_amf)) * treatment + mean_per_fc + I(mean_per_fc^2) + mean_drc_mm + soil + myco_cat
g_text <- glm(fa, data = sdf, weights = n_total, family = binomial)
used <- sdf[rownames(model.frame(g_text)) %>% as.integer(), ]
g_s5 <- glm(fa, data = used, weights = n_total, family = binomial)
rt <- coef(summary(g_text))["scale(mean_amf):treatmentDrought", ]; r5 <- coef(summary(g_s5))["scale(mean_amf):treatmentDrought", ]
cat(sprintf("  Text (z-scored over all %d candidate groups with an AMF mean): %s\n", sum(is.finite(sdf$mean_amf)), zstr(rt[1], rt[3], rt[4])))
cat(sprintf("  Appendix S5 (z-scored over the %d groups in the model):      %s\n", nobs(g_s5), zstr(r5[1], r5[3], r5[4])))
cat(sprintf("  SD of group mean AMF: all candidate groups %.3f | model groups %.3f | ratio %.3f; coefficient ratio %.3f\n",
            sd(sdf$mean_amf[is.finite(sdf$mean_amf)]), sd(used$mean_amf), sd(sdf$mean_amf[is.finite(sdf$mean_amf)]) / sd(used$mean_amf), rt[1] / r5[1]))
cat("  Per 1 percentage point of mean AMF (scale-free):", signif(rt[1] / sd(sdf$mean_amf[is.finite(sdf$mean_amf)]), 4), "vs", signif(r5[1] / sd(used$mean_amf), 4), "\n")

cat("\n7b. Fremont ancestry x EMF: -0.171 (Appendix S6) vs EMF slope -0.124 (Figure 3 legend)\n")
b <- fixef(col)
sl <- c(AMF_ww = b["per_f"], AMF_dr = b["per_f"] + b["per_f:treatmentDrought"],
        EMF_ww = b["per_f"] + b["per_f:fungal_typeemf"], EMF_dr = b["per_f"] + b["per_f:fungal_typeemf"] + b["per_f:treatmentDrought"] + b["per_f:fungal_typeemf:treatmentDrought"])
names(sl) <- c("AMF_ww", "AMF_dr", "EMF_ww", "EMF_dr")
cat(sprintf("  Reference guild = AMF, reference treatment = well-watered. Coefficients: per_f = %.4f; per_f:EMF = %.4f; per_f:drought = %.4f; per_f:EMF:drought = %.4f\n",
            b["per_f"], b["per_f:fungal_typeemf"], b["per_f:treatmentDrought"], b["per_f:fungal_typeemf:treatmentDrought"]))
cat(sprintf("  Ancestry slopes: AMF well-watered %.4f, drought %.4f | EMF well-watered %.4f, drought %.4f\n", sl["AMF_ww"], sl["AMF_dr"], sl["EMF_ww"], sl["EMF_dr"]))
cat(sprintf("  EMF slope averaged over watering treatments = (%.4f + %.4f)/2 = %.4f  <- Figure 3 legend (-0.124)\n", sl["EMF_ww"], sl["EMF_dr"], mean(sl[c("EMF_ww", "EMF_dr")])))
cat(sprintf("  EMF minus AMF slope, well-watered = %.4f  <- Appendix S6 'Fremont ancestry x EMF' (-0.171)\n", sl["EMF_ww"] - sl["AMF_ww"]))
ct <- summary(pairs(emtrends(col, ~ fungal_type, var = "per_f", lmer.df = "satterthwaite"), adjust = "none"))
cat("  EMF minus AMF slope averaged over treatments (emtrends contrast):\n"); print(ct)
cat("\nDone.\n")
sink()
