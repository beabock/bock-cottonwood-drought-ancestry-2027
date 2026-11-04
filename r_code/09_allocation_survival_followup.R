# ==============================================================================
# ALLOCATION, PREDICTED SURVIVAL AND FOLLOW-UP TESTS.
# Analyses for the reframed manuscript; everything is written to results/09_allocation_survival_followup/.
#   1. Allocation: root:new-stem core LMM in samples (a) and (b), total new biomass, and stem, leaf and root mass
#      fractions (ordered beta GLMMs). (Section 2, Figure 1, is drawn by r_code/28_figures.R.)
#   3. Conductance vs. EMF colonization and vs. root mass fraction (plants with both).
#   4. Counts and tables.
#   5. Inoculum as a fixed effect in the colonization LMM (and colonization-biomass model).
#   6. EMF x drought with other responses and within narrowleaf plants.
# Uses the alive-only growth data (SFF209 removed) and the allocation data of 06_allocation_models.R.
# ==============================================================================
# Run from the project root (R 4.6.1):  Rscript r_code/09_allocation_survival_followup.R
# Writes results/09_allocation_survival_followup/: output.txt, strings.csv/.md, table_*.csv
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(lme4); library(lmerTest); library(glmmTMB); library(emmeans); library(tidyverse)
out <- "results/09_allocation_survival_followup"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
cat("ALLOCATION AND SURVIVAL FOLLOW-UP run:", "|", R.version.string, "| lme4", as.character(packageVersion("lme4")),
    "| lmerTest", as.character(packageVersion("lmerTest")), "| glmmTMB", as.character(packageVersion("glmmTMB")), "\n")
set.seed(1998)

# --- AJB strings ---
EN <- "–"; BETA <- "β"; TIMES <- "×"; CHI <- "χ²"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
tstr <- function(b, t, df, p) sprintf("%s = %s, t = %s, df = %s, %s", BETA, num(b), dec(t), dff(df), P(p))
zstr <- function(b, z, p) sprintf("%s = %s, z = %s, %s", BETA, num(b), dec(z), P(p))
pc <- function(x) dash(sprintf("%.1f", x))
strings <- list()
add <- function(section, analysis, sample, n, term, string, p = NA) {
  strings[[length(strings) + 1]] <<- tibble(section, analysis, sample, n, term, string, p)
  cat(sprintf("  [%s] %s | %s | n = %s | %s: %s\n", section, analysis, sample, n, term, string)) }
row_t <- function(m, t) { r <- coef(summary(m))[t, ]; list(s = tstr(r["Estimate"], r["t value"], r["df"], r["Pr(>|t|)"]), p = unname(r["Pr(>|t|)"]), b = unname(r["Estimate"])) }


ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, pf50 = per_fc / 50, geno_id = paste(pop, type, geno, sep = "_"))
treatment_levels <- levels(ds$treatment)
cov <- "og_drc_mm + myco_cat + soil"; re <- "(1 | pop) + (1 | pop:geno_id)"

# core LMM (ancestry linear + quadratic x drought; LRT; reduced if P >= 0.05)
core_lmm <- function(y, d, label, sample, section = "1") {
  full <- lmerTest::lmer(as.formula(paste(y, "~ per_fc * treatment + I(per_fc^2) * treatment +", cov, "+", re)), data = d)
  red <- lmerTest::lmer(as.formula(paste(y, "~ per_fc + treatment + I(per_fc^2) +", cov, "+", re)), data = d)
  l <- anova(red, full, refit = TRUE); pl <- l$`Pr(>Chisq)`[2]; m <- if (pl < 0.05) full else red
  add(section, label, sample, nobs(m), paste("Ancestry", TIMES, "drought (LRT, 2 df)"), sprintf("%s(%d) = %s, %s", CHI, l$Df[2], dec(l$Chisq[2]), P(pl)), pl)
  for (t in c(per_fc = "Linear ancestry", "I(per_fc^2)" = "Quadratic ancestry", treatmentDrought = "Drought") %>% imap_chr(~ .y) %>% names()) {
    lab <- c(per_fc = "Linear ancestry", "I(per_fc^2)" = "Quadratic ancestry", treatmentDrought = "Drought")[t]
    r <- row_t(m, t); add(section, label, sample, nobs(m), lab, r$s, r$p) }
  r <- resid(m, scaled = TRUE)
  cat(sprintf("  diagnostics: max |scaled resid| = %.2f; %d beyond |3|; Shapiro-Wilk P = %s; singular = %s; model = %s\n",
              max(abs(r)), sum(abs(r) > 3), formatC(shapiro.test(r)$p.value, digits = 2, format = "g"), isSingular(m), ifelse(pl < 0.05, "full", "reduced")))
  m }

# allocation data (06_allocation_models.R: leaf 0 g when leafless; root mass estimated for colonization-scored plants)
al <- read_csv("results/06_allocation_models/allocation_data.csv", show_col_types = FALSE) %>%
  select(tree_idold, scored, root_est, new_stem, leaf, total_new, LMF, RMF, SMF) %>%
  inner_join(ds %>% select(tree_idold, per_f, per_fc, pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id, new_class), by = "tree_idold") %>%
  mutate(root_stem = root_est / new_stem, log_root_stem = log(root_stem), log_total_new = log(total_new))
stopifnot(nrow(al) == 414)
samp <- list(`(a) estimated root mass` = al, `(b) measured root mass only` = filter(al, !scored))

# ------------------------------------------------------------------------------
hdr("1. Allocation: root:new-stem ratio (core LMM) and new Figure 2")
# ------------------------------------------------------------------------------
cat("Root:new-stem ratio = root dry mass / new stem dry mass (log scale). Sample (a): n = 414; (b): n = 242.\n")
rs_mod <- imap(samp, ~ core_lmm("log(root_stem)", .x, "log root:new-stem ratio", .y))
tot_a <- core_lmm("log(total_new)", samp[[1]], "log total new biomass (Figure 2A model)", names(samp)[1])
emm_rs <- map(rs_mod, ~ summary(emmeans(.x, ~ treatment, type = "response")))
cat("\nRoot:new-stem ratio, back-transformed estimated marginal means by treatment:\n"); print(emm_rs)
for (s in names(samp)) { e <- emm_rs[[s]]
  add("1", "Root:new-stem ratio, model mean", s, nobs(rs_mod[[s]]), "Well-watered / drought (95% CI)",
      sprintf("%.2f (%.2f–%.2f) / %.2f (%.2f–%.2f); drought %s %.2f", e$response[1], e$lower.CL[1], e$upper.CL[1], e$response[2], e$lower.CL[2], e$upper.CL[2], TIMES,
              e$response[2] / e$response[1])) }

# mass fractions: ordered beta GLMMs (reduced; as 06_allocation_models.R), sample (a)
fit_ob <- function(y, d) { f <- as.formula(paste(y, "~ pf50 + treatment + I(pf50^2) +", cov, "+", re))
  for (ctrl in list(glmmTMBControl(), glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))) {
    m <- tryCatch(glmmTMB(f, data = d, family = ordbeta(), control = ctrl), error = function(e) NULL)
    if (!is.null(m) && isTRUE(m$sdr$pdHess)) return(m) }
  stop("ordbeta did not converge for ", y) }
fr <- c(SMF = "Stem", LMF = "Leaf", RMF = "Root")
ob <- map(names(fr), ~ fit_ob(.x, samp[[1]])) %>% set_names(names(fr))
ob_mean <- function(eta, psi) { p0 <- 1 - plogis(eta - psi[1]); p1 <- plogis(eta - psi[2]); p1 + (1 - p0 - p1) * plogis(eta) }
psi_of <- function(m) { p <- m$fit$par; unname(p[names(p) == "psi"]) }
# glmmTMB's predict(type = "response") for ordbeta returns only the beta-component mean plogis(eta); the Figure 2 curves
# (r_code/28_figures.R) use the full expected value, which adds the probabilities of exact 0 and 1 from the cutpoints (psi). Check both parts.
for (nm in names(ob)) { m <- ob[[nm]]; X <- getME(m, "X"); eta <- as.vector(X %*% fixef(m)$cond); psi <- psi_of(m); y <- samp[[1]][[nm]]
  dev <- max(abs(plogis(eta) - predict(m, re.form = NA, type = "response")))
  cat(sprintf("%s: predict() = plogis(eta) (max |difference| %.1g) | predicted share of exact zeros %.3f vs observed %.3f | mean expected fraction %.3f vs observed mean %.3f\n",
              nm, dev, mean(1 - plogis(eta - psi[1])), mean(y == 0), mean(ob_mean(eta, psi)), mean(y)))
  stopifnot(dev < 1e-6) }
stopifnot(all(map_int(ob, nobs) == nrow(samp[[1]])))

# ------------------------------------------------------------------------------
hdr("3. Above-belowground: whole-leaf conductance vs. EMF colonization and root mass fraction")
# ------------------------------------------------------------------------------
ab <- ds %>% filter(!is.na(amf), !is.na(emf), !is.na(dse), por_avg > 0) %>% left_join(al %>% select(tree_idold, RMF), by = "tree_idold") %>%
  filter(if_all(c(per_fc, treatment, og_drc_mm, soil, myco_cat, pop, geno_id), ~ !is.na(.x))) %>% mutate(emf10 = emf / 10, rmf10 = RMF * 10)
cat("Plants with whole-leaf conductance and all three guilds scored:", nrow(ab), "| with root mass fraction:", sum(!is.na(ab$RMF)),
    "(root mass is estimated for all of them, so RMF results are approximate) | genotypes:", n_distinct(ab$geno_id), "\n")
# core conductance structure (as reported): linear ancestry x drought + quadratic ancestry + covariates
crhs <- "per_fc * treatment + I(per_fc^2) + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id)"
for (pr in c(emf10 = "EMF colonization (per 10 percentage points)", rmf10 = "Root mass fraction (per 0.1)")) {}
for (v in c("emf10", "rmf10")) {
  lab <- c(emf10 = "EMF colonization (per 10 percentage points)", rmf10 = "Root mass fraction (per 0.1; approximate)")[v]
  d <- filter(ab, !is.na(.data[[v]]))
  m <- lmerTest::lmer(as.formula(paste("log(por_avg) ~", v, "+", crhs)), data = d)
  mi <- lmerTest::lmer(as.formula(paste("log(por_avg) ~", v, "* treatment +", crhs)), data = d)
  r <- row_t(m, v); ri <- row_t(mi, paste0(v, ":treatmentDrought"))
  add("3", "log whole-leaf conductance", "plants with conductance and colonization", nobs(m), lab,
      sprintf("%s; %s%% change in conductance (95%% CI %s to %s%%)", r$s, pc(100 * (exp(r$b) - 1)),
              pc(100 * (exp(confint(m, parm = v, method = "Wald")[1]) - 1)), pc(100 * (exp(confint(m, parm = v, method = "Wald")[2]) - 1))), r$p)
  add("3", "log whole-leaf conductance", "plants with conductance and colonization", nobs(mi), paste(lab, TIMES, "drought"), ri$s, ri$p)
  cat("  singular:", isSingular(m), "\n") }

# ------------------------------------------------------------------------------
hdr("4. Counts and tables")
# ------------------------------------------------------------------------------
raw <- read.csv("data/edited_data/mega_ds.csv")
uh <- raw %>% filter(tree_idold %in% ds$tree_idold[ds$alive_01_10_2021 %in% 1 & is.na(ds$dryshoot_new)]) %>%
  transmute(plant = tree_idold, new_stem_mass = dryshoot_new, stem_mass = dry_stem, leaf_mass = dryleaf_weight, root_mass = root_dryweight,
            fresh_root_mass = wet_root, notes = coalesce(notes, og_notes)) %>%
  mutate(status = case_when(!is.na(stem_mass) | !is.na(leaf_mass) | !is.na(root_mass) ~ "harvested; new stem mass not recorded",
                            TRUE ~ "no harvest masses recorded")) %>% arrange(status, plant)
readr::write_excel_csv(uh, file.path(out, "table_unharvested_survivors.csv"))
cat("Survivors (alive 1 Oct 2021):", sum(ds$alive_01_10_2021 %in% 1), "| with new stem mass:", sum(ds$alive_01_10_2021 %in% 1 & !is.na(ds$dryshoot_new)),
    "| without:", nrow(uh), "\n"); print(table(uh$status)); print(as.data.frame(uh), row.names = FALSE)

mt <- ds %>% distinct(pop, geno_id, new_class) %>% count(pop, new_class) %>% pivot_wider(names_from = new_class, values_from = n, values_fill = 0) %>%
  select(pop, any_of(levels(ds$new_class))) %>%
  left_join(ds %>% group_by(pop) %>% summarise(mother_trees = n_distinct(geno_id), cuttings = n(), survivors = sum(alive_01_10_2021 %in% 1), .groups = "drop"), by = "pop")
gl <- ds %>% distinct(geno_id, new_class) %>% count(geno_id) %>% filter(n > 1)
cat("\nMother trees by population and cross type (genotypes spanning more than one cross-type label:", nrow(gl), ")\n"); print(as.data.frame(mt), row.names = FALSE)
cat("Totals: mother trees", n_distinct(ds$geno_id), "| narrowleaf", n_distinct(ds$geno_id[ds$new_class == "Narrowleaf"]), "| Fremont",
    n_distinct(ds$geno_id[ds$new_class == "Fremont"]), "| hybrid or backcross", n_distinct(ds$geno_id[!ds$new_class %in% c("Narrowleaf", "Fremont")]), "\n")
cat("Population elevations are not in the data (mega_ds.csv has no elevation column). GPS coordinates by population and cross type are in\n",
    "data/raw_data/WC Garden Collections_Coords.xlsx (sheet 'Sources for GH Expt') for Indian Creek, Little Colorado, Blue River and San Francisco River only.\n", sep = "")
readr::write_excel_csv(mt, file.path(out, "table_mother_trees.csv"))

colz <- ds %>% filter(!is.na(amf), !is.na(emf), !is.na(dse))
ct1 <- colz %>% count(myco_cat, treatment) %>% pivot_wider(names_from = treatment, values_from = n, values_fill = 0) %>% mutate(Total = rowSums(across(where(is.numeric))))
ct2 <- colz %>% count(new_class, treatment) %>% pivot_wider(names_from = treatment, values_from = n, values_fill = 0) %>% mutate(Total = rowSums(across(where(is.numeric))))
cat("\nColonized plants (all three guilds, n =", nrow(colz), "), inoculum x drought:\n"); print(as.data.frame(ct1), row.names = FALSE)
cat("Colonized plants by cross type x drought:\n"); print(as.data.frame(ct2), row.names = FALSE)
readr::write_excel_csv(bind_rows(ct1 %>% rename(group = myco_cat) %>% mutate(table = "inoculum x drought"), ct2 %>% rename(group = new_class) %>% mutate(table = "cross type x drought")) %>%
                         relocate(table) %>% mutate(group = as.character(group)), file.path(out, "table_colonized_plants.csv"))

# reported survival GLMM (population fixed, genotype random), predicted at 0% and 100% Fremont ancestry
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
marg_x <- function(fe_terms, nd, cols) colMeans(model.matrix(fe_terms, data = nd)[, cols, drop = FALSE])
curve_by_treatment <- function(beta, V, fe_terms, base, xgrid, xvar, xfun = identity) {
  map_dfr(xgrid, function(g) map_dfr(treatment_levels, function(trt) {
    nd <- base; nd[[xvar]] <- xfun(g); nd$treatment <- factor(trt, levels = treatment_levels)
    xm <- marg_x(fe_terms, nd, names(beta)); tibble(x = g, treatment = trt, fit = sum(xm * beta), se = sqrt(as.numeric(t(xm) %*% V %*% xm))) })) %>%
    mutate(treatment = factor(treatment, levels = treatment_levels)) }
surv_red <- glmer(alive_01_10_2021 ~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat + (1 | geno_id), data = ds, family = binomial, control = ctl)
sp <- ds %>% filter(if_all(c(alive_01_10_2021, pf50, treatment, pop, og_drc_mm, soil, myco_cat, geno_id), ~ !is.na(.x)))
stopifnot(nrow(sp) == nobs(surv_red))
c1 <- curve_by_treatment(fixef(surv_red), as.matrix(vcov(surv_red)), delete.response(terms(~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat)),
                         sp %>% select(pop, og_drc_mm, soil, myco_cat, treatment), c(0, 100), "pf50", function(g) (g - 50) / 50) %>%
  mutate(y = plogis(fit), lower = plogis(fit - 1.96 * se), upper = plogis(fit + 1.96 * se))
ps <- c1 %>% filter(x %in% c(0, 100)) %>% transmute(fremont_ancestry = paste0(x, "%"), treatment, survival_pct = 100 * y, lower = 100 * lower, upper = 100 * upper)
cat("\nPredicted survival (reported survival GLMM, averaged over population, inoculum, soil and diameter; 95% CI):\n"); print(as.data.frame(ps %>% mutate(across(where(is.numeric), ~ round(., 1)))), row.names = FALSE)
for (i in seq_len(nrow(ps))) add("4", "Predicted survival", "survival model, n = 563", nobs(surv_red), paste(ps$fremont_ancestry[i], "Fremont,", ps$treatment[i]),
                                 sprintf("%.1f%% (95%% CI %.1f–%.1f%%)", ps$survival_pct[i], ps$lower[i], ps$upper[i]))
readr::write_excel_csv(ps, file.path(out, "table_predicted_survival.csv"))

# ------------------------------------------------------------------------------
hdr("5. Inoculum as a fixed effect")
# ------------------------------------------------------------------------------
dl <- colz %>% select(tree_id, pop, geno_id, per_f, new_class, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>% na.omit() %>%
  mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm))) %>% pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>%
  mutate(fungal_type = factor(fungal_type, levels = c("amf", "dse", "emf")), myco_cat = factor(myco_cat))
col_rep <- lmerTest::lmer(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil + (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id), data = dl)
col_fix <- lmerTest::lmer(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil + myco_cat * fungal_type + (1 | pop) + (1 | pop:geno_id) + (1 | tree_id), data = dl)
for (nm in c("reported (inoculum random)", "inoculum fixed + inoculum x fungal type")) {
  m <- if (grepl("^reported", nm)) col_rep else col_fix
  r <- row_t(m, "per_f:fungal_typeemf"); add("5", "Colonization LMM", nm, n_distinct(dl$tree_id), "EMF vs. AMF ancestry slope, well-watered", r$s, r$p)
  et <- as.data.frame(summary(emtrends(m, ~ fungal_type, var = "per_f", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE)))
  for (g in c("emf", "amf", "dse")) { e <- et[et$fungal_type == g, ]
    add("5", "Colonization LMM", nm, n_distinct(dl$tree_id), paste(toupper(g), "ancestry slope (averaged over treatments)"), tstr(e$per_f.trend, e$t.ratio, e$df, e$p.value), e$p.value) } }
a5 <- anova(col_fix, type = 3)
for (t in c("myco_cat", "fungal_type:myco_cat")) add("5", "Colonization LMM", "inoculum fixed + inoculum x fungal type", n_distinct(dl$tree_id), sub("fungal_type:myco_cat", paste("Inoculum", TIMES, "fungal type"), sub("^myco_cat$", "Inoculum", t)),
                                                     sprintf("F = %s, df = %d, %s, %s", dec(a5[t, "F value"]), a5[t, "NumDF"], dff(a5[t, "DenDF"]), P(a5[t, "Pr(>F)"])), a5[t, "Pr(>F)"])
cat("Colonization-biomass model: inoculum category is already a fixed effect in the reported model (log new stem ~ ... + og_drc_mm + myco_cat + soil); no refit needed.\n")

# ------------------------------------------------------------------------------
hdr("6. EMF x drought robustness: other responses; narrowleaf plants only")
# ------------------------------------------------------------------------------
zsd <- c(amf = sd(ds$amf, na.rm = TRUE), emf = sd(ds$emf, na.rm = TRUE), dse = sd(ds$dse, na.rm = TRUE))   # as scale() in the reported model
zmu <- c(amf = mean(ds$amf, na.rm = TRUE), emf = mean(ds$emf, na.rm = TRUE), dse = mean(ds$dse, na.rm = TRUE))
fb <- ds %>% filter(!is.na(amf), !is.na(emf), !is.na(dse)) %>% left_join(al %>% select(tree_idold, log_total_new, log_root_stem), by = "tree_idold") %>%
  mutate(amf_z = (amf - zmu["amf"]) / zsd["amf"], emf_z = (emf - zmu["emf"]) / zsd["emf"], dse_z = (dse - zmu["dse"]) / zsd["dse"])
frhs <- function(anc) paste("amf_z * treatment + emf_z * treatment + dse_z * treatment +", anc, "+ og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id)")
runs <- list(list("log new stem biomass (reported)", "log(dryshoot_new)", fb, "per_fc + I(per_fc^2)"),
             list("log total new biomass (root mass estimated; approximate)", "log_total_new", fb, "per_fc + I(per_fc^2)"),
             list("log root:new-stem ratio (root mass estimated; approximate)", "log_root_stem", fb, "per_fc + I(per_fc^2)"),
             list("log new stem biomass, narrowleaf plants only (linear ancestry only)", "log(dryshoot_new)", filter(fb, new_class == "Narrowleaf"), "per_fc"))
for (r in runs) {
  m <- lmerTest::lmer(as.formula(paste(r[[2]], "~", frhs(r[[4]]))), data = r[[3]])
  x <- row_t(m, "treatmentDrought:emf_z")
  add("6", "Colonization-biomass model", r[[1]], nobs(m), paste("EMF", TIMES, "drought (per SD of EMF, 13.8 percentage points)"), x$s, x$p)
  cat(sprintf("  genotypes = %d; singular = %s; max |scaled resid| = %.2f\n", ngrps(m)["pop:geno_id"], isSingular(m), max(abs(resid(m, scaled = TRUE))))) }

# ------------------------------------------------------------------------------
hdr("7. Conductance drought effect; Figure 4 predictions")
# ------------------------------------------------------------------------------
# (b) Reported whole-leaf conductance model: per_fc = Fremont ancestry - 50, so the drought coefficient is the drought
#     effect at 50% Fremont ancestry; at other ancestry values it is drought + (ancestry - 50) x (ancestry x drought).
pm <- lmerTest::lmer(log(por_avg) ~ per_fc * treatment + I(per_fc^2) + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id), data = filter(ds, por_avg > 0))
stopifnot(nobs(pm) == 203)
b <- fixef(pm); V <- as.matrix(vcov(pm)); dfm <- coef(summary(pm))["treatmentDrought", "df"]
for (a in c(0, 50, 100)) {
  L <- setNames(rep(0, length(b)), names(b)); L["treatmentDrought"] <- 1; L["per_fc:treatmentDrought"] <- a - 50
  est <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% V %*% L)); tv <- est / se
  add("7", "Conductance: drought effect by ancestry", "reported whole-leaf model", nobs(pm), sprintf("Drought at %d%% Fremont ancestry", a),
      sprintf("%s; %s%% conductance under drought", tstr(est, tv, dfm, 2 * pt(-abs(tv), dfm)), pc(100 * (exp(est) - 1))), 2 * pt(-abs(tv), dfm)) }
cat("  (df for the 0% and 100% contrasts use the Satterthwaite df of the drought coefficient; approximate.)\n")
# (c) Figure 4 curves (r_code/28_figures.R): for each grid value of the focal guild, every one of the 77 plants keeps its
#     observed values of the other guilds and covariates, the model-matrix rows are averaged, and exp() is applied. The
#     linear predictor is linear in the other guilds, so this equals holding them at their means over the 77 plants.
zz <- fb %>% filter(!is.na(dryshoot_new), !is.na(og_drc_mm))
f4 <- log(dryshoot_new) ~ (amf_z + emf_z + dse_z) * treatment + per_fc + I(per_fc^2) + og_drc_mm + myco_cat + soil
m4 <- lmer(update(f4, . ~ . + (1 | pop) + (1 | pop:geno_id)), data = zz); stopifnot(nobs(m4) == 77)
fe4 <- delete.response(terms(f4)); g <- 0.5; nd <- zz %>% mutate(emf_z = g, treatment = factor("Drought", levels = treatment_levels))
avg_rows <- colMeans(model.matrix(fe4, nd)); at_means <- model.matrix(fe4, nd %>% mutate(amf_z = mean(zz$amf_z), dse_z = mean(zz$dse_z)))
cat(sprintf("Figure 4 check (EMF panel, drought, EMF z = 0.5): averaged design rows give %.6f on the log scale; other guilds at their means give %.6f\n",
            sum(avg_rows * fixef(m4)), mean(as.vector(at_means %*% fixef(m4)))))

st <- bind_rows(strings)
readr::write_excel_csv(st, file.path(out, "strings.csv"))
md <- c("# Allocation and survival follow-up: AJB-style strings", "", "Generated by `r_code/09_allocation_survival_followup.R`. Exploratory; not reported results. Ancestry coefficients per 1% Fremont ancestry; t tests with Satterthwaite df; GLMM rows give Wald z.", "",
        "| Section | Analysis | Sample | n | Term | String |", "|---|---|---|---|---|---|",
        sprintf("| %s | %s | %s | %s | %s | %s |", st$section, st$analysis, st$sample, st$n, st$term, st$string))
con <- file(file.path(out, "strings.md"), "w", encoding = "UTF-8"); writeLines(enc2utf8(md), con); close(con)
cat("\nDone.\n")
sink()
