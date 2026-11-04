# ==============================================================================
# Colonization on every scored plant, the Q3 test, inoculum contrasts,
# and source-site propagules (EMF vs. ancestry in control plants).
# Single-guild LMMs keep the covariates and random effects of the reported colonization model
#   (r_code/08_sample_flow_extra_tests.R: per_f * treatment, initial diameter (z), soil; random population, inoculum,
#   genotype within population; the plant random intercept is not needed with one observation per plant):
#   colonization ~ per_f * treatment + og_drc_mm_z + soil + (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id)
# Slopes per 1% Fremont ancestry from emtrends (Satterthwaite df), against zero, averaged over watering and in each
# watering level. Because ancestry is a genotype-level predictor, each slope is also tested with genotype-level df
# (number of genotypes in that sample - 2), and as an OLS slope on genotype means (one point per genotype, or per
# genotype x watering).
# Run from the project root:  Rscript r_code/21_colonization_all_scored_plants.R
# Writes results/21_colonization_all_scored_plants/: output.txt, strings.csv, guild_slopes.csv, appendix_S2_samples.csv, inoculum_medians.csv
# ==============================================================================
source("r_code/helpers/model_helpers.R")
OUT <- "results/21_colonization_all_scored_plants"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "output.txt"), split = TRUE); run_stamp("Colonization, all scored plants")
src <- "results/21_colonization_all_scored_plants/strings.csv"
ds <- ds %>% mutate(all3 = !is.na(amf) & !is.na(emf) & !is.na(dse), inoc = factor(ifelse(myco_cat == "None", "Control", "Inoculated"), levels = c("Control", "Inoculated")))

hdr("2a. Which plants and models Appendix S2 uses")
# Figure S2 (r_code/28_figures.R): every plant with a score for that guild and an inoculum category
fs2 <- ds %>% pivot_longer(c(amf, emf, dse), names_to = "guild", values_to = "col") %>% filter(!is.na(col), !is.na(myco_cat))
# Appendix S2 LMM (r_code/02a_main_models.R): long data, rows with complete covariates
dall <- ds %>% select(tree_id, alive_01_10_2021, pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse, all3) %>%
  pivot_longer(c(amf, emf, dse), names_to = "guild", values_to = "colonization") %>% filter(!is.na(colonization)) %>%
  filter(if_all(c(pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat), ~ !is.na(.x)))
a2 <- bind_rows(fs2 %>% group_by(guild) %>% summarise(plants = n_distinct(tree_id), all_three = n_distinct(tree_id[all3]), partial = plants - all_three,
                                                      alive_at_harvest = n_distinct(tree_id[alive_01_10_2021 %in% 1]), genotypes = n_distinct(geno_id), .groups = "drop") %>% mutate(use = "Figure S2 boxes and points"),
                dall %>% group_by(guild) %>% summarise(plants = n_distinct(tree_id), all_three = n_distinct(tree_id[all3]), partial = plants - all_three,
                                                       alive_at_harvest = n_distinct(tree_id[alive_01_10_2021 %in% 1]), genotypes = n_distinct(geno_id), .groups = "drop") %>% mutate(use = "Appendix S2 LMM (inoculum, soil, inoculum × guild)")) %>%
  mutate(guild = toupper(guild)) %>% relocate(use)
print(as.data.frame(a2), row.names = FALSE)
cat("Appendix S2 LMM: rows =", nrow(dall), "| plants =", n_distinct(dall$tree_id), "| genotypes =", n_distinct(dall$geno_id),
    "| model: colonization ~ myco_cat * fungal_type + soil + treatment + per_f + og_drc_mm_z + (1 | pop) + (1 | pop:geno_id) + (1 | tree_id); Type II F, Satterthwaite\n")
readr::write_excel_csv(a2, file.path(OUT, "appendix_S2_samples.csv"))
rec(2, "a", "Appendix S2 samples", "Figure S2 and S2 LMM", n_distinct(dall$tree_id), "Plants per guild",
    sprintf("Figure S2 and the S2 LMM both use every scored plant: AMF %d, DSE %d, EMF %d plants (EMF: %d with all three guilds + %d EMF only; AMF/DSE: %d + %d AMF and DSE only); LMM %d rows from %d plants, %d genotypes",
            a2$plants[a2$use != "Figure S2 boxes and points" & a2$guild == "AMF"], a2$plants[a2$use != "Figure S2 boxes and points" & a2$guild == "DSE"],
            a2$plants[a2$use != "Figure S2 boxes and points" & a2$guild == "EMF"], a2$all_three[a2$use != "Figure S2 boxes and points" & a2$guild == "EMF"],
            a2$partial[a2$use != "Figure S2 boxes and points" & a2$guild == "EMF"], a2$all_three[a2$use != "Figure S2 boxes and points" & a2$guild == "AMF"],
            a2$partial[a2$use != "Figure S2 boxes and points" & a2$guild == "AMF"], nrow(dall), n_distinct(dall$tree_id), n_distinct(dall$geno_id)), NA, "results/21_colonization_all_scored_plants/appendix_S2_samples.csv")

hdr("2b, 2d. Ancestry slope for each guild against zero; ancestry x drought (Q3)")
one_guild <- function(g, d, lab) {
  d <- d %>% filter(!is.na(.data[[g]])) %>% filter(if_all(c(pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat), ~ !is.na(.x))) %>%
    mutate(colonization = .data[[g]], og_drc_mm_z = as.numeric(scale(og_drc_mm)))
  m <- lmerTest::lmer(colonization ~ per_f * treatment + og_drc_mm_z + soil + (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id), data = d)
  G <- n_distinct(d$geno_id); Gw <- d %>% group_by(treatment) %>% summarise(G = n_distinct(geno_id)) %>% deframe()
  avg <- as.data.frame(summary(emtrends(m, ~ 1, var = "per_f", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE)))
  byw <- as.data.frame(summary(emtrends(m, ~ treatment, var = "per_f", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE)))
  gm <- d %>% group_by(geno_id) %>% summarise(per_f = mean(per_f), y = mean(colonization), .groups = "drop")
  gw <- d %>% group_by(geno_id, treatment) %>% summarise(per_f = mean(per_f), y = mean(colonization), .groups = "drop")
  ols <- function(dd) { cs <- coef(summary(lm(y ~ per_f, data = dd)))["per_f", ]; tibble(b = cs[1], t = cs[3], df = nrow(dd) - 2, p = cs[4]) }
  iw <- row_t(m, "per_f:treatmentDrought")
  rows <- bind_rows(
    tibble(scope = "averaged over watering", b = avg$per_f.trend, se = avg$SE, t = avg$t.ratio, df = avg$df, p = avg$p.value, G = G),
    tibble(scope = as.character(byw$treatment), b = byw$per_f.trend, se = byw$SE, t = byw$t.ratio, df = byw$df, p = byw$p.value, G = Gw[as.character(byw$treatment)])) %>%
    mutate(p_geno_df = 2 * pt(-abs(t), G - 2),
           ols = list(ols(gm), ols(filter(gw, treatment == "Well-watered")), ols(filter(gw, treatment == "Drought"))),
           guild = toupper(g), sample = lab, n = nrow(d), n_ww = sum(d$treatment == "Well-watered"), n_dr = sum(d$treatment == "Drought"), singular = isSingular(m),
           geno_var = as.data.frame(VarCorr(m)) %>% filter(grp == "pop:geno_id") %>% pull(vcov)) %>% unnest_wider(ols, names_sep = "_")
  list(rows = rows, inter = iw, m = m, d = d) }
samples <- list(`78 plants with all three guilds` = filter(ds, all3), `all plants scored for the guild` = ds)
gs <- list(); q3 <- list()
for (s in names(samples)) for (g in c("emf", "amf", "dse")) {
  r <- one_guild(g, samples[[s]], s); gs[[length(gs) + 1]] <- r$rows
  for (i in seq_len(nrow(r$rows))) { x <- r$rows[i, ]
    rec(2, "b", sprintf("%s ancestry slope (per 1%% Fremont), single-guild LMM", toupper(g)), sprintf("%s (well-watered %d, drought %d; %d genotypes)", s, x$n_ww, x$n_dr, x$G), x$n,
        paste("Slope,", x$scope), sprintf("%s; genotype-level df (%d genotypes − 2): t = %s, df = %d, %s; genotype-mean OLS: %s", tstr(x$b, x$t, x$df, x$p), x$G, dec(x$t), x$G - 2, P(x$p_geno_df),
                                         tstr(x$ols_b, x$ols_t, x$ols_df, x$ols_p)), x$p_geno_df, "results/21_colonization_all_scored_plants/guild_slopes.csv") }
  rec(2, "d", sprintf("%s colonization: ancestry × drought (Q3)", toupper(g)), s, r$rows$n[1], paste("Ancestry", TIMES, "drought (drought − well-watered slope)"), r$inter$s, r$inter$p, src)
  cat(sprintf("    %s / %s: genotype variance %.3g; singular %s\n", s, g, r$rows$geno_var[1], r$rows$singular[1])) }
gs <- bind_rows(gs); readr::write_excel_csv(gs %>% select(-geno_var), file.path(OUT, "guild_slopes.csv"))

# Q3 in the reported three-guild model: EMF slope, drought minus well-watered (emtrends contrast)
dl <- ds %>% filter(all3) %>% select(tree_id, pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat, amf, emf, dse) %>% na.omit() %>%
  mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm)), myco_cat = factor(myco_cat)) %>% pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>%
  mutate(fungal_type = factor(fungal_type, levels = c("amf", "dse", "emf")))
col <- lmerTest::lmer(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil + (1 | pop) + (1 | myco_cat) + (1 | pop:geno_id) + (1 | tree_id), data = dl)
et <- emtrends(col, ~ treatment | fungal_type, var = "per_f", lmer.df = "satterthwaite")
cq <- as.data.frame(summary(contrast(et, list(`drought - well-watered` = c(-1, 1))))) %>% filter(fungal_type == "emf")
sq <- as.data.frame(summary(et, infer = c(TRUE, TRUE))) %>% filter(fungal_type == "emf")
rec(2, "d", "EMF colonization: ancestry × drought (Q3), reported three-guild model", "78 plants", 78, "EMF slope, drought − well-watered", tstr(cq$estimate, cq$t.ratio, cq$df, cq$p.value), cq$p.value, src)
for (i in 1:2) rec(2, "d", "EMF ancestry slope by watering, reported three-guild model", "78 plants", 78, paste("EMF slope,", sq$treatment[i]), tstr(sq$per_f.trend[i], sq$t.ratio[i], sq$df[i], sq$p.value[i]), sq$p.value[i], src)

hdr("2e. Colonization in control vs. inoculated plants (all scored plants)")
iq <- function(x) sprintf("median %.1f%% (IQR %.1f–%.1f%%, n = %d)", median(x), quantile(x, 0.25), quantile(x, 0.75), length(x))
med <- list()
for (k in list(c("emf", "EM", "EMF inoculum"), c("emf", "EM|Both", "EMF or dual inoculum"), c("amf", "AM", "AMF inoculum"), c("amf", "AM|Both", "AMF or dual inoculum"))) {
  d <- ds %>% filter(!is.na(.data[[k[1]]]), myco_cat == "None" | grepl(paste0("^(", k[2], ")$"), myco_cat)) %>% mutate(grp = ifelse(myco_cat == "None", "Control", k[3]), y = .data[[k[1]]])
  w <- wilcox.test(y ~ grp, data = d, exact = FALSE)
  dd <- d %>% filter(if_all(c(pop, geno_id, treatment, og_drc_mm, soil, per_f), ~ !is.na(.x))) %>% mutate(grp = factor(grp, levels = c("Control", k[3])))
  lm1 <- lmerTest::lmer(y ~ grp + per_f + treatment + soil + (1 | pop) + (1 | pop:geno_id), data = dd); rl <- row_t(lm1, paste0("grp", k[3]))
  W <- unname(w$statistic); ctl_x <- d$y[d$grp == "Control"]; trt_x <- d$y[d$grp != "Control"]
  med[[length(med) + 1]] <- tibble(guild = toupper(k[1]), comparison = paste("Control vs.", k[3]), control = iq(ctl_x), inoculated = iq(trt_x), wilcoxon_W = W, wilcoxon_P = w$p.value,
                                   lmm = rl$s, lmm_n = nobs(lm1))
  rec(2, "e", sprintf("%s colonization, control vs. %s", toupper(k[1]), k[3]), "all plants scored for the guild", nrow(d), "Medians and test",
      sprintf("control %s; %s %s; Wilcoxon rank-sum W = %s, %s; LMM (inoculum + ancestry + drought + soil, random population and genotype): difference %s", iq(ctl_x), k[3], iq(trt_x), format(W, big.mark = ""), P(w$p.value), rl$s),
      w$p.value, "results/21_colonization_all_scored_plants/inoculum_medians.csv") }
med <- bind_rows(med); print(as.data.frame(med), row.names = FALSE); readr::write_excel_csv(med, file.path(OUT, "inoculum_medians.csv"))

hdr("3a. Control plants only: EMF colonization vs. Fremont ancestry")
for (s in names(samples)) {
  d <- samples[[s]] %>% filter(myco_cat == "None", !is.na(emf)) %>% filter(if_all(c(pop, geno_id, per_f, treatment, og_drc_mm, soil), ~ !is.na(.x))) %>%
    mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm)))
  m <- lmerTest::lmer(emf ~ per_f * treatment + og_drc_mm_z + soil + (1 | pop) + (1 | pop:geno_id), data = d)
  a <- as.data.frame(summary(emtrends(m, ~ 1, var = "per_f", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE))); G <- n_distinct(d$geno_id)
  gm <- d %>% group_by(geno_id) %>% summarise(per_f = mean(per_f), y = mean(emf), .groups = "drop"); o <- coef(summary(lm(y ~ per_f, data = gm)))["per_f", ]
  rec(3, "a", "EMF ancestry slope (per 1% Fremont), control plants only", sprintf("controls, %s (well-watered %d, drought %d; %d genotypes)", s, sum(d$treatment == "Well-watered"), sum(d$treatment == "Drought"), G), nrow(d),
      "Slope, averaged over watering", sprintf("%s; genotype-level df: t = %s, df = %d, %s; genotype-mean OLS: %s", tstr(a$per_f.trend, a$t.ratio, a$df, a$p.value), dec(a$t.ratio), G - 2,
                                                P(2 * pt(-abs(a$t.ratio), G - 2)), tstr(o[1], o[3], nrow(gm) - 2, o[4])), 2 * pt(-abs(a$t.ratio), G - 2), src)
  cat("    singular:", isSingular(m), "\n") }

hdr("3b. Does the ancestry effect on EMF differ between control and inoculated plants?")
for (s in names(samples)) {
  d <- samples[[s]] %>% filter(!is.na(emf)) %>% filter(if_all(c(pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat), ~ !is.na(.x))) %>% mutate(og_drc_mm_z = as.numeric(scale(og_drc_mm)))
  m <- lmerTest::lmer(emf ~ per_f * treatment + per_f * inoc + og_drc_mm_z + soil + (1 | pop) + (1 | pop:geno_id), data = d)
  r <- row_t(m, "per_f:inocInoculated"); sl <- as.data.frame(summary(emtrends(m, ~ inoc, var = "per_f", lmer.df = "satterthwaite"), infer = c(TRUE, TRUE)))
  rec(3, "b", "EMF: ancestry × inoculation (control vs. any inoculum)", sprintf("%s (control %d, inoculated %d)", s, sum(d$inoc == "Control"), sum(d$inoc == "Inoculated")), nrow(d),
      "Ancestry × inoculated", sprintf("%s; slopes: control %s, inoculated %s", r$s, tstr(sl$per_f.trend[1], sl$t.ratio[1], sl$df[1], sl$p.value[1]), tstr(sl$per_f.trend[2], sl$t.ratio[2], sl$df[2], sl$p.value[2])), r$p, src) }

write_rec("strings.csv"); cat("\nDone.\n"); sink()
