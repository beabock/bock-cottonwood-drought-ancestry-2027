# ==============================================================================
# Conductance values of record with measurement day accounted for.
# Three fits of log whole-leaf conductance with the reported fixed structure
#   per_fc * treatment + I(per_fc^2) + og_drc_mm + soil + myco_cat, random (1 | pop) + (1 | pop:geno_id):
#   reported     plant means, 203 plants (r_code/02b_main_models_followup.R)
#   reading      log of each plant-day mean, + phase (wet/dry) + day pair, + (1 | plant) (as in an earlier conductance script not included here)
#   pair_17_22   plant means, plants measured on 17/22 Sep, the only day pair with both watering treatments
# Drought contrasts (lmerTest::contest, Satterthwaite df for t and P; 95% CI = estimate ± 1.96 SE, the convention of the reported conductance CIs):
#   averaged over ancestry = drought + mean(per_fc of the fitted plants) x (ancestry x drought), and at 0, 50 and 100% Fremont.
# Also refits the conductance vs. EMF colonization and root mass fraction tests (r_code/09_allocation_survival_followup.R)
# with day pair as a covariate. The ancestry x drought coefficient is reported but not interpreted.
# Also fits the same reading-level model to each leaf surface, and a joint surface model.
# Writes s10_*.csv, read by r_code/30_appendix_tables.R (Appendix S5 tables and text) and r_code/28_figures.R.
# Run from the project root:  Rscript r_code/19_conductance_models_of_record.R
# ==============================================================================
suppressMessages({ source("r_code/helpers/utilities_data_prep.R"); library(tidyverse); library(lmerTest) })
out <- "results/19_conductance_models_of_record"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
sink(file.path(out, "output.txt"), split = TRUE)
cat("CONDUCTANCE MODELS OF RECORD run:", "|", R.version.string, "| lme4", as.character(packageVersion("lme4")), "| lmerTest", as.character(packageVersion("lmerTest")), "\n")
EN <- "–"; dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
tstr <- function(b, t, df, p) sprintf("β = %s, t = %s, df = %s, %s", num(b), dec(t), dff(df), P(p))
norm <- function(x) toupper(sub("^0+", "", as.character(x)))

# ---- data (as in an earlier conductance script not included here) ----
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, geno_id = paste(pop, type, geno, sep = "_"), key = norm(poro_id))
rd <- read_csv("data/raw_data/gh_poromeasurements_allgroups.csv", col_types = cols(poro_id = "c"), show_col_types = FALSE) %>%
  mutate(key = norm(poro_id), date = as.Date(date, format = "%m/%d/%Y"), phase = water) %>%
  rowwise() %>% mutate(day_mean = mean(c_across(l1_t:l3_b), na.rm = TRUE), adaxial = mean(c(l1_t, l2_t, l3_t), na.rm = TRUE),
                       abaxial = mean(c(l1_b, l2_b, l3_b), na.rm = TRUE)) %>% ungroup()
pair_of <- function(d) case_when(d %in% as.Date(c("2021-09-08", "2021-09-09")) ~ "8/9 Sep", d %in% as.Date(c("2021-09-17", "2021-09-22")) ~ "17/22 Sep",
                                 d %in% as.Date(c("2021-09-20", "2021-09-23")) ~ "20/23 Sep", d %in% as.Date(c("2021-09-21", "2021-09-26")) ~ "21/26 Sep")
LV <- c("8/9 Sep", "17/22 Sep", "20/23 Sep", "21/26 Sep")
plants <- ds %>% filter(!is.na(por_avg)) %>% select(tree_idold, key, pop, geno_id, per_f, per_fc, treatment, og_drc_mm, soil, myco_cat, por_avg, amf, emf, dse)
pd <- rd %>% inner_join(plants, by = "key") %>% mutate(day_pair = factor(pair_of(date), levels = LV))
stopifnot(n_distinct(pd$tree_idold) == 205, !anyNA(pd$day_pair))
plants <- plants %>% left_join(pd %>% group_by(tree_idold) %>% summarise(day_pair = first(day_pair), .groups = "drop"), by = "tree_idold")
d_rep <- plants %>% filter(por_avg > 0, !is.na(og_drc_mm)); stopifnot(nrow(d_rep) == 203)

# ---- fits ----
fe <- "per_fc * treatment + I(per_fc^2) + og_drc_mm + soil + myco_cat"; re <- "(1 | pop) + (1 | pop:geno_id)"
fits <- list(
  reported = lmerTest::lmer(as.formula(paste("log(por_avg) ~", fe, "+", re)), data = d_rep),
  reading = lmerTest::lmer(as.formula(paste("log(day_mean) ~", fe, "+ phase + day_pair +", re, "+ (1 | pop:geno_id:tree_idold)")),
                           data = pd %>% filter(tree_idold %in% d_rep$tree_idold, day_mean > 0)),
  pair_17_22 = lmerTest::lmer(as.formula(paste("log(por_avg) ~", fe, "+", re)), data = filter(d_rep, day_pair == "17/22 Sep")))
lab <- c(reported = "Reported model (plant means)", reading = "Reading-level model (plant-day means; phase and day pair as fixed effects)",
         pair_17_22 = "Plants measured on 17/22 Sep only (plant means)")
mf_plants <- function(k) if (k == "reading") unique(model.frame(fits[[k]])$`pop:geno_id:tree_idold`) else NULL
n_pl <- c(reported = 203, reading = n_distinct(pd$tree_idold[pd$tree_idold %in% d_rep$tree_idold & pd$day_mean > 0]), pair_17_22 = nobs(fits$pair_17_22))
mean_fc <- c(reported = mean(d_rep$per_fc), reading = mean(d_rep$per_fc), pair_17_22 = mean(d_rep$per_fc[d_rep$day_pair == "17/22 Sep"]))

rows <- list()
for (k in names(fits)) {
  m <- fits[[k]]; cs <- coef(summary(m)); b <- fixef(m)
  coef_row <- function(tm, label) { r <- cs[tm, ]; rows[[length(rows) + 1]] <<- tibble(model = k, model_label = lab[[k]], n_plants = n_pl[[k]], n_obs = nobs(m), term = label,
    estimate = r[[1]], t = r[["t value"]], df = r[["df"]], P = r[["Pr(>|t|)"]], pct = NA_real_, lo = NA_real_, hi = NA_real_,
    string = tstr(r[[1]], r[["t value"]], r[["df"]], r[["Pr(>|t|)"]])) }
  coef_row("per_fc", "Linear ancestry (per 1% Fremont; well-watered)"); coef_row("I(per_fc^2)", "Quadratic ancestry")
  coef_row("per_fc:treatmentDrought", "Ancestry × drought (reported, not interpreted)")
  slope_row <- function(m, b, w_dr = 0.5) {   # linear ancestry slope averaged over watering (equal weights); % per 10 points
    L <- setNames(rep(0, length(b)), names(b)); L["per_fc"] <- 1; L["per_fc:treatmentDrought"] <- w_dr
    ct <- contest(m, L, joint = FALSE); se <- ct$`Std. Error`; e <- ct$Estimate
    tibble(estimate = e, t = ct$`t value`, df = ct$df, P = ct$`Pr(>|t|)`, pct = 100 * (exp(10 * e) - 1), lo = 100 * (exp(10 * (e - 1.96 * se)) - 1), hi = 100 * (exp(10 * (e + 1.96 * se)) - 1),
           string = sprintf("%s%% per 10 percentage points (95%% CI %s to %s%%), t = %s, df = %s, %s", dec(100 * (exp(10 * e) - 1), 1), dec(100 * (exp(10 * (e - 1.96 * se)) - 1), 1),
                            dec(100 * (exp(10 * (e + 1.96 * se)) - 1), 1), dec(ct$`t value`), dff(ct$df), P(ct$`Pr(>|t|)`))) }
  rows[[length(rows) + 1]] <- bind_cols(tibble(model = k, model_label = lab[[k]], n_plants = n_pl[[k]], n_obs = nobs(m), term = "Ancestry slope, averaged over watering"), slope_row(m, b))
  if (k == "reading") { coef_row("phasewet", "Phase: after watering vs. end of cycle"); for (lv in LV[-1]) coef_row(paste0("day_pair", lv), paste("Day pair", lv, "vs. 8/9 Sep")) }
  for (a in list(list("Drought, averaged over ancestry", mean_fc[[k]]), list("Drought at 0% Fremont", -50), list("Drought at 50% Fremont", 0), list("Drought at 100% Fremont", 50))) {
    L <- setNames(rep(0, length(b)), names(b)); L["treatmentDrought"] <- 1; L["per_fc:treatmentDrought"] <- a[[2]]
    ct <- contest(m, L, joint = FALSE); q <- 1.96   # Wald CI, as the reported conductance CIs
    pct <- 100 * (exp(ct$Estimate) - 1); lo <- 100 * (exp(ct$Estimate - q * ct$`Std. Error`) - 1); hi <- 100 * (exp(ct$Estimate + q * ct$`Std. Error`) - 1)
    rows[[length(rows) + 1]] <- tibble(model = k, model_label = lab[[k]], n_plants = n_pl[[k]], n_obs = nobs(m),
      term = if (a[[1]] == "Drought, averaged over ancestry") sprintf("%s (at mean ancestry, %.1f%% Fremont)", a[[1]], a[[2]] + 50) else a[[1]],
      estimate = ct$Estimate, t = ct$`t value`, df = ct$df, P = ct$`Pr(>|t|)`, pct = pct, lo = lo, hi = hi,
      string = sprintf("%s%% (95%% CI %s to %s%%), t = %s, df = %s, %s", dec(pct, 1), dec(lo, 1), dec(hi, 1), dec(ct$`t value`), dff(ct$df), P(ct$`Pr(>|t|)`)))
  }
}
tab <- bind_rows(rows)
cat("\nConductance fits (singular: ", paste(names(fits), sapply(fits, isSingular), sep = " = ", collapse = "; "), ")\n\n", sep = "")
print(as.data.frame(tab %>% select(model, n_plants, n_obs, term, string)), right = FALSE, row.names = FALSE)
write.csv(tab, file.path(out, "conductance_models.csv"), row.names = FALSE)

# ---- surface-specific reading-level models (as the whole-leaf model of record) and the joint surface model ----
pd_s <- pd %>% filter(tree_idold %in% d_rep$tree_idold)
srows <- list()
for (sf in c("adaxial", "abaxial")) {
  d <- pd_s %>% filter(!is.na(.data[[sf]]), .data[[sf]] > 0)
  m <- lmerTest::lmer(as.formula(paste0("log(", sf, ") ~ ", fe, " + phase + day_pair + ", re, " + (1 | pop:geno_id:tree_idold)")), data = d)
  b <- fixef(m); cs <- coef(summary(m))
  L <- setNames(rep(0, length(b)), names(b)); L["treatmentDrought"] <- 1; L["per_fc:treatmentDrought"] <- mean(d_rep$per_fc)
  ct <- contest(m, L, joint = FALSE)
  srows[[length(srows) + 1]] <- bind_rows(
    bind_cols(tibble(surface = sf, n_plants = n_distinct(d$tree_idold), n_obs = nobs(m), term = "Ancestry slope, averaged over watering"), slope_row(m, b)),
    tibble(surface = sf, n_plants = n_distinct(d$tree_idold), n_obs = nobs(m), term = "Drought, averaged over ancestry (at mean ancestry)", estimate = ct$Estimate, t = ct$`t value`, df = ct$df, P = ct$`Pr(>|t|)`,
           pct = 100 * (exp(ct$Estimate) - 1), lo = 100 * (exp(ct$Estimate - 1.96 * ct$`Std. Error`) - 1), hi = 100 * (exp(ct$Estimate + 1.96 * ct$`Std. Error`) - 1),
           string = sprintf("%s%% (95%% CI %s to %s%%), t = %s, df = %s, %s", dec(100 * (exp(ct$Estimate) - 1), 1), dec(100 * (exp(ct$Estimate - 1.96 * ct$`Std. Error`) - 1), 1),
                            dec(100 * (exp(ct$Estimate + 1.96 * ct$`Std. Error`) - 1), 1), dec(ct$`t value`), dff(ct$df), P(ct$`Pr(>|t|)`))),
    tibble(surface = sf, n_plants = n_distinct(d$tree_idold), n_obs = nobs(m), term = "Ancestry × drought (reported, not interpreted)", estimate = cs["per_fc:treatmentDrought", 1],
           t = cs["per_fc:treatmentDrought", "t value"], df = cs["per_fc:treatmentDrought", "df"], P = cs["per_fc:treatmentDrought", "Pr(>|t|)"],
           string = tstr(cs["per_fc:treatmentDrought", 1], cs["per_fc:treatmentDrought", "t value"], cs["per_fc:treatmentDrought", "df"], cs["per_fc:treatmentDrought", "Pr(>|t|)"])))
}
long <- pd_s %>% select(tree_idold, pop, geno_id, per_fc, treatment, og_drc_mm, soil, myco_cat, phase, day_pair, date, adaxial, abaxial) %>%
  pivot_longer(c(adaxial, abaxial), names_to = "surface", values_to = "gsw") %>% filter(!is.na(gsw), gsw > 0) %>% mutate(surface = factor(surface, levels = c("abaxial", "adaxial")))
jm <- lmerTest::lmer(log(gsw) ~ surface * (per_fc * treatment + I(per_fc^2)) + og_drc_mm + soil + myco_cat + phase + day_pair +
                       (1 | pop) + (1 | pop:geno_id) + (1 | pop:geno_id:tree_idold) + (1 | pop:geno_id:tree_idold:phase), data = long)
bj <- fixef(jm); L <- setNames(rep(0, length(bj)), names(bj)); L["surfaceadaxial:per_fc"] <- 1; L["surfaceadaxial:per_fc:treatmentDrought"] <- 0.5
ct <- contest(jm, L, joint = FALSE)
srows[[length(srows) + 1]] <- tibble(surface = "adaxial − abaxial (joint model)", n_plants = n_distinct(long$tree_idold), n_obs = nobs(jm),
  term = "Difference in ancestry slope between surfaces, averaged over watering", estimate = ct$Estimate, t = ct$`t value`, df = ct$df, P = ct$`Pr(>|t|)`,
  string = tstr(ct$Estimate, ct$`t value`, ct$df, ct$`Pr(>|t|)`))
st <- bind_rows(srows)
cat("
Surface-specific reading-level models (phase + day pair; plant random) and joint surface model (plant and plant-day random):
")
print(as.data.frame(st %>% select(surface, n_plants, n_obs, term, string)), right = FALSE, row.names = FALSE)
write.csv(st, file.path(out, "surface_models.csv"), row.names = FALSE)

# ---- plant-level means for Figures S5.1 and S5.2 (r_code/28_figures.R): mean of all readings per plant and surface ----
pm <- rd %>% inner_join(d_rep %>% select(tree_idold, key, per_f, treatment, day_pair), by = "key") %>% group_by(tree_idold, per_f, treatment, day_pair) %>%
  summarise(whole_leaf = mean(c(l1_t, l2_t, l3_t, l1_b, l2_b, l3_b), na.rm = TRUE), adaxial = mean(c(l1_t, l2_t, l3_t), na.rm = TRUE),
            abaxial = mean(c(l1_b, l2_b, l3_b), na.rm = TRUE), n_days = n_distinct(date), .groups = "drop")
stopifnot(nrow(pm) == 203, max(abs(pm$whole_leaf - d_rep$por_avg[match(pm$tree_idold, d_rep$tree_idold)])) < 1e-9)
write.csv(pm, file.path(out, "plant_means.csv"), row.names = FALSE)
p89 <- pm %>% filter(day_pair == "8/9 Sep")   # Appendix S5: every plant measured on 8/9 Sep was well-watered with > 50% Fremont ancestry
cat(sprintf("Plants measured on 8/9 Sep: %d; all well-watered: %s; lowest Fremont ancestry %.1f%% (all above 50%%: %s)\n",
            nrow(p89), all(p89$treatment == "Well-watered"), min(p89$per_f), all(p89$per_f > 50)))

# ---- design table: plants by day pair x watering x ancestry tercile (203 model plants) ----
cuts <- quantile(d_rep$per_f, c(1/3, 2/3))
dt <- d_rep %>% mutate(tercile = cut(per_f, c(-Inf, cuts, Inf), labels = sprintf(c("≤ %.1f%%", "%.1f–%.1f%%", "> %.1f%%"), c(cuts[1], cuts[1], cuts[2]), c(NA, cuts[2], NA)) %>% sub("NA", "", .))) %>%
  count(day_pair, treatment, tercile) %>% complete(day_pair = factor(LV, levels = LV), treatment, tercile, fill = list(n = 0)) %>%
  pivot_wider(names_from = tercile, values_from = n) %>% mutate(Total = rowSums(across(where(is.numeric))))
cat("\nPlants by day pair x watering x ancestry tercile (Fremont ancestry terciles of the 203 model plants; cut points",
    paste(sprintf("%.1f%%", cuts), collapse = " and "), "):\n"); print(as.data.frame(dt), row.names = FALSE)
write.csv(dt, file.path(out, "daypair_watering_tercile.csv"), row.names = FALSE)

# ---- conductance vs. EMF colonization and root mass fraction, with day pair as a covariate ----
al <- read_csv("results/06_allocation_models/allocation_data.csv", show_col_types = FALSE) %>% select(tree_idold, RMF)
ab <- d_rep %>% filter(!is.na(amf), !is.na(emf), !is.na(dse)) %>% left_join(al, by = "tree_idold") %>% mutate(emf10 = emf / 10, rmf10 = RMF * 10)
cat("\nPlants with conductance and all three guilds:", nrow(ab), "| with root mass fraction:", sum(!is.na(ab$RMF)), "\n")
cat("Day pairs among them:\n"); print(table(ab$day_pair, ab$treatment))
crhs <- paste(fe, "+", re)
rows2 <- list()
for (v in c("emf10", "rmf10")) for (dp in c(FALSE, TRUE)) {
  d <- filter(ab, !is.na(.data[[v]])); extra <- if (dp) " + day_pair" else ""
  m <- lmerTest::lmer(as.formula(paste("log(por_avg) ~", v, extra, "+", crhs)), data = d)
  mi <- lmerTest::lmer(as.formula(paste("log(por_avg) ~", v, "* treatment", extra, "+", crhs)), data = d)
  r <- coef(summary(m))[v, ]; ri <- coef(summary(mi))[paste0(v, ":treatmentDrought"), ]; q <- 1.96
  rows2[[length(rows2) + 1]] <- tibble(predictor = c(emf10 = "EMF colonization (per 10 percentage points)", rmf10 = "Root mass fraction (per 0.1; root mass estimated)")[[v]],
    day_pair_covariate = dp, n = nobs(m), estimate = r["Estimate"], t = r["t value"], df = r["df"], P = r["Pr(>|t|)"],
    pct = 100 * (exp(r["Estimate"]) - 1), lo = 100 * (exp(r["Estimate"] - q * r["Std. Error"]) - 1), hi = 100 * (exp(r["Estimate"] + q * r["Std. Error"]) - 1),
    string = sprintf("β = %s, t = %s, df = %s, %s; %s%% change in conductance (95%% CI %s to %s%%)", num(r["Estimate"]), dec(r["t value"]), dff(r["df"]), P(r["Pr(>|t|)"]),
                     dec(100 * (exp(r["Estimate"]) - 1), 1), dec(100 * (exp(r["Estimate"] - q * r["Std. Error"]) - 1), 1), dec(100 * (exp(r["Estimate"] + q * r["Std. Error"]) - 1), 1)),
    x_drought = tstr(ri["Estimate"], ri["t value"], ri["df"], ri["Pr(>|t|)"]), singular = isSingular(m))
}
t2 <- bind_rows(rows2)
cat("\n"); print(as.data.frame(t2 %>% select(predictor, day_pair_covariate, n, string, x_drought, singular)), right = FALSE, row.names = FALSE)
cat("(Without day pair these reproduce r_code/09_allocation_survival_followup.R section 3, which uses the Wald normal CI; here the CI uses the t distribution.)\n")
write.csv(t2, file.path(out, "conductance_vs_emf_rmf.csv"), row.names = FALSE)
sink()
