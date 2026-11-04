# ==============================================================================
# APPENDIX TABLES: the tables of Appendices S1 and S4-S7, one CSV per manuscript table, in results/appendices/.
#   S1  Table S1.1 mother trees by population and cross type; Table S1.2 collection points
#         (collection points and elevations from r_code/12_map_appendix_S1.R)
#   S4  main colonization analyses repeated without control plants
#         (r_code/05_control_inoculum_sensitivity.R, r_code/18_emf_amf_slopes_by_watering.R)
#   S5  stomatal conductance: Tables S5.1-S5.4 (r_code/19_conductance_models_of_record.R)
#   S6  Table S6.2 numbers of plants in each analysis (r_code/08_sample_flow_extra_tests.R and the counts below).
#         Table S6.1 (each analysis and how it departs from the core model) is descriptive and is given in the manuscript only.
#   S7  allocation of new growth and leaf drop (r_code/07_allocation_appendix_S7.R, r_code/22_allocation_log_ratios.R)
# Run last, after every other step:  Rscript r_code/30_appendix_tables.R
# CSVs use readr::write_excel_csv (UTF-8 with BOM, so Excel shows β, χ² and en dashes correctly).
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(tidyverse)
out <- "results/appendices"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
EN <- "–"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
rd <- function(f) read_csv(f, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
blank_rep <- function(x) ifelse(!is.na(lag(x)) & x == lag(x), "", x)        # show a grouping label once
yn <- function(x) ifelse(x, "Yes", "No")
p_of <- function(s) { v <- str_match(s, "P ([=<]) ([0-9.]+)"); as.numeric(v[, 3]) }   # P from an AJB string
# Sig. change: P < 0.05 in one analysis and not the other; Sign change: estimates of opposite sign;
# Concl. change: a change in significance, or a significant result in both with opposite signs
changes <- function(p1, p2, e1, e2) {
  sig <- (p1 < 0.05) != (p2 < 0.05); sgn <- sign(e1) != sign(e2)
  tibble(`Sig. change` = yn(sig), `Sign change` = ifelse(is.na(sgn), "–", yn(sgn)),
         `Concl. change` = yn(sig | (!is.na(sgn) & sgn & p1 < 0.05 & p2 < 0.05))) }
write_tab <- function(d, f) { readr::write_excel_csv(d, file.path(out, f)); cat("written:", file.path(out, f), "\n") }
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, geno_id = paste(pop, type, geno, sep = "_"))
states <- c(AZ = "Arizona", UT = "Utah", NM = "New Mexico")

# ------------------------------------------------------------------------------
# Appendix S1: Tables S1.1 and S1.2
# ------------------------------------------------------------------------------
cpts <- rd(file.path(out, "appendix_S1_collection_points.csv"))
stopifnot(setequal(cpts$Code, unique(ds$pop)))
erng <- cpts %>% group_by(Code) %>% summarise(Population = first(Population), State = states[first(State)],
                                               `Elevation (m)` = paste0(min(`Elevation (m)`), EN, max(`Elevation (m)`)), .groups = "drop")
s1 <- ds %>% distinct(pop, geno_id, new_class) %>% count(pop, new_class) %>%
  pivot_wider(names_from = new_class, values_from = n, values_fill = 0) %>% select(pop, any_of(levels(ds$new_class))) %>%
  left_join(ds %>% group_by(pop) %>% summarise(`Mother trees` = n_distinct(geno_id), Cuttings = n(), .groups = "drop"), by = "pop") %>%
  rename(Code = pop, `F₁ hybrid` = `F1 Hybrid`, `Backcross narrowleaf` = `Backcross Narrowleaf`) %>%
  left_join(erng, by = "Code") %>% relocate(Population, State, Code) %>% arrange(Population)
s1 <- bind_rows(s1 %>% mutate(across(where(is.numeric), as.character)),
                s1 %>% summarise(across(where(is.numeric), sum)) %>% mutate(across(everything(), as.character), Population = "Total", State = "", Code = "",
                                                                             `Elevation (m)` = paste0(min(cpts$`Elevation (m)`), EN, max(cpts$`Elevation (m)`))))
stopifnot(s1$`Mother trees`[s1$Population == "Total"] == "43", s1$Cuttings[s1$Population == "Total"] == "574")
write_tab(s1, "appendix_S1_table1_mother_trees.csv")

date_range <- function(x) {   # "03/03/2021 - 03/06/2021" -> "3–6 Mar 2021" (all ranges fall within one month)
  d <- str_split_fixed(x, " - ", 2); a <- as.Date(d[, 1], "%m/%d/%Y"); b <- as.Date(d[, 2], "%m/%d/%Y")
  stopifnot(format(a, "%Y-%m") == format(b, "%Y-%m"))
  sprintf("%d%s%d %s", as.integer(format(a, "%d")), EN, as.integer(format(b, "%d")), format(b, "%b %Y")) }
s1b <- cpts %>% transmute(Population, State = states[State], `Cross type collected`, `Genotypes collected`,
                          `Collection dates` = date_range(`Collection dates`),
                          `Latitude (°N)` = sprintf("%.5f", Latitude), `Longitude (°W)` = sprintf("%.5f", -Longitude), `Elevation (m)`) %>%
  mutate(across(c(Population, State, `Collection dates`), blank_rep))
write_tab(s1b, "appendix_S1_table2_collection_points.csv")

# ------------------------------------------------------------------------------
# Appendix S4: main colonization analyses with and without control plants
# ------------------------------------------------------------------------------
sl <- rd("results/18_emf_amf_slopes_by_watering/emf_amf_by_watering.csv") %>%
  filter(fit %in% c("Reported LMM (Satterthwaite)", "Reported LMM without control plants (Satterthwaite)")) %>%
  mutate(data = ifelse(grepl("without control", fit), "inoc", "all"),
         Test = recode(contrast, "EMF − AMF slope, averaged over watering" = "EMF − AMF ancestry slope, averaged over watering",
                       "EMF − AMF slope, well-watered" = "EMF − AMF ancestry slope, well-watered plants",
                       "EMF − AMF slope, drought" = "EMF − AMF ancestry slope, drought plants",
                       "Difference (drought − well-watered): ancestry × EMF-vs-AMF × drought" = "Change in EMF − AMF slope with drought (ancestry × fungal type × drought)"),
         Analysis = "Colonization LMM") %>%
  select(Analysis, Test, data, n = n_plants, string, p = P, est = estimate)
cs <- rd("results/05_control_inoculum_sensitivity/control_sensitivity.csv") %>% filter(analysis != "Colonization LMM") %>%
  mutate(data = ifelse(data == "All plants (reported)", "all", "inoc"),
         Analysis = recode(analysis, "Colonization-biomass LMM" = "Colonization–biomass LMM", "PERMANOVA (genotype-level)" = "PERMANOVA (genotype-level permutations)"),
         Test = recode(term, "Cross type (3 focal)" = "Cross type (narrowleaf, F₁ hybrid, Fremont)", "Continuous ancestry" = "Continuous Fremont ancestry"),
         est = ifelse(grepl("PERMANOVA", analysis), NA, sign)) %>%
  select(Analysis, Test, data, n, string, p, est)
s4 <- bind_rows(sl, cs) %>% pivot_wider(names_from = data, values_from = c(n, string, p, est))
s4 <- bind_cols(s4 %>% transmute(Analysis, Test, `n, all` = n_all, `n, inoc.` = n_inoc, `All plants` = string_all, `Inoculated plants only` = string_inoc),
                changes(s4$p_all, s4$p_inoc, s4$est_all, s4$est_inoc)) %>% mutate(Analysis = blank_rep(Analysis))
stopifnot(nrow(s4) == 9)
write_tab(s4, "appendix_S4_control_sensitivity.csv")

# ------------------------------------------------------------------------------
# Appendix S5: Tables S5.1-S5.4 (stomatal conductance)
# ------------------------------------------------------------------------------
cm <- rd("results/19_conductance_models_of_record/conductance_models.csv")
nm <- cm %>% distinct(model, n_plants, n_obs)
col_lab <- c(reading = sprintf("Model of record: plant-day means (%d plants, %d plant-days)", nm$n_plants[nm$model == "reading"], nm$n_obs[nm$model == "reading"]),
             reported = sprintf("Plant means (%d plants)", nm$n_plants[nm$model == "reported"]),
             pair_17_22 = sprintf("Plants measured on 17 and 22 September (%d plants)", nm$n_plants[nm$model == "pair_17_22"]))
term_lab <- c("Ancestry slope, averaged over watering" = "Ancestry, averaged over watering (per 10 percentage points)",
              "Drought, averaged over ancestry" = "Drought, averaged over ancestry (at mean ancestry)",
              "Drought at 0% Fremont" = "Drought at 0% Fremont ancestry", "Drought at 50% Fremont" = "Drought at 50% Fremont ancestry",
              "Drought at 100% Fremont" = "Drought at 100% Fremont ancestry",
              "Linear ancestry (per 1% Fremont; well-watered)" = "Linear ancestry (per 1%, well-watered)", "Quadratic ancestry" = "Quadratic ancestry",
              "Ancestry × drought (reported, not interpreted)" = "Ancestry × drought (reported, not interpreted)",
              "Phase: after watering vs. end of cycle" = "Reading after watering vs. end of cycle",
              "Day pair 17/22 Sep vs. 8/9 Sep" = "Dates 17/22 vs. 8/9 September", "Day pair 20/23 Sep vs. 8/9 Sep" = "Dates 20/23 vs. 8/9 September",
              "Day pair 21/26 Sep vs. 8/9 Sep" = "Dates 21/26 vs. 8/9 September")
s5_1 <- cm %>% mutate(Term = term_lab[sub(" \\(at mean ancestry, [0-9.]+% Fremont\\)", "", term)]) %>%
  select(model, Term, string) %>% pivot_wider(names_from = model, values_from = string, values_fill = "–") %>%
  mutate(Term = factor(Term, levels = unname(term_lab))) %>% arrange(Term) %>% mutate(Term = as.character(Term)) %>%
  select(Term, all_of(names(col_lab))) %>% rename(!!!setNames(names(col_lab), col_lab))
stopifnot(!anyNA(s5_1$Term), nrow(s5_1) == length(term_lab))
write_tab(s5_1, "appendix_S5_table1_conductance_models.csv")

sm <- rd("results/19_conductance_models_of_record/surface_models.csv")
s5_2 <- sm %>% mutate(Surface = recode(surface, adaxial = "Adaxial", abaxial = "Abaxial", "adaxial − abaxial (joint model)" = "Adaxial minus abaxial, joint model"),
                      `Plants (observations)` = sprintf("%d (%d)", n_plants, n_obs),
                      col = case_when(startsWith(term, "Ancestry slope") | startsWith(term, "Difference") ~ "anc",
                                      startsWith(term, "Drought") ~ "dr", TRUE ~ "x"),
                      string = ifelse(startsWith(term, "Difference"), paste("Difference in ancestry slope:", string), string)) %>%
  select(Surface, `Plants (observations)`, col, string) %>% pivot_wider(names_from = col, values_from = string, values_fill = "–") %>%
  transmute(Surface, `Plants (observations)`, `Ancestry, averaged over watering (per 10 percentage points)` = anc,
            `Drought, averaged over ancestry` = dr, `Ancestry × drought (reported, not interpreted)` = x)
stopifnot(nrow(s5_2) == 3)
write_tab(s5_2, "appendix_S5_table2_surface_models.csv")

s5_3 <- rd("results/19_conductance_models_of_record/conductance_vs_emf_rmf.csv") %>%
  transmute(Predictor = predictor, `Pair of dates as covariate` = yn(day_pair_covariate), Plants = n,
            `Effect on conductance` = string, `Interaction with drought` = x_drought)
write_tab(s5_3, "appendix_S5_table3_conductance_vs_emf_rmf.csv")

s5_4 <- rd("results/19_conductance_models_of_record/daypair_watering_tercile.csv") %>%
  rename(`Pair of dates` = day_pair, Watering = treatment) %>% rename_with(~ paste("Fremont ancestry", .x), starts_with("≤"))
write_tab(s5_4, "appendix_S5_table4_plants_by_dates.csv")

# ------------------------------------------------------------------------------
# Appendix S6: Table S6.2, numbers of plants in each analysis
# ------------------------------------------------------------------------------
fl <- rd("results/08_sample_flow_extra_tests/flow_table.csv")
s6 <- fl %>% transmute(Analysis = sub(" \\((LMM|binomial GLMM)\\)", "", model) %>% sub("^Allocation, ", "Allocation ", .) %>% sub("Colonization-biomass", "Colonization–biomass", .),
                       `Model` = coalesce(str_match(model, "\\((LMM|binomial GLMM)\\)")[, 2], "See Appendix S7") %>% recode("binomial GLMM" = "Binomial GLMM"),
                       `Eligible plants` = paste0(recode(pool, "all cuttings" = "All cuttings", "harvested survivors (alive, new stem mass)" = "Harvested survivors", "plants with all three guilds scored" = "Plants with all three guilds scored", "plants with a porometer mean" = "Plants with porometer readings"), " (", pool_n, ")"), `Plants analyzed` = model_n, Excluded = dropped,
                       `Reasons for exclusion (plants)` = gsub(": ([0-9]+)", " (\\1)", drop_reasons) %>% str_replace_all("(^|\\); )([a-z])", function(m) toupper(m)))
s6_row <- function(a, m, pool, pool_n, n, why) tibble(Analysis = a, Model = m, `Eligible plants` = sprintf("%s (%d)", pool, pool_n), `Plants analyzed` = n, Excluded = pool_n - n, `Reasons for exclusion (plants)` = why)
d78 <- ds %>% filter(!is.na(per_f), !is.na(pop), !is.na(treatment), !is.na(amf), !is.na(emf), !is.na(dse), !is.na(soil), !is.na(myco_cat)); stopifnot(nrow(d78) == 78)
foc78 <- d78 %>% filter(new_class %in% c("Narrowleaf", "F1 Hybrid", "Fremont"))
bc <- table(droplevels(d78$new_class[!d78$new_class %in% c("Narrowleaf", "F1 Hybrid", "Fremont")]))
pm <- rd("results/19_conductance_models_of_record/plant_means.csv")
dp <- rd("results/19_conductance_models_of_record/daypair_watering_tercile.csv")
ce <- rd("results/19_conductance_models_of_record/conductance_vs_emf_rmf.csv") %>% filter(!day_pair_covariate)
alr <- rd("results/14_allocation_fraction_checks/alr_models.csv") %>% distinct(sample, n)
n_rd <- cm$n_obs[cm$model == "reading"][1]; n_pd <- sum(pm$n_days); n1722 <- sum(dp$Total[dp$`day_pair` == "17/22 Sep"])
n_pos <- sum(ds$por_avg > 0, na.rm = TRUE); n_any <- sum(!is.na(ds$por_avg))
hv <- ds %>% filter(alive_01_10_2021 %in% 1, !is.na(dryshoot_new)) %>% mutate(sc = !is.na(amf) | !is.na(emf) | !is.na(dse))
# EMF models on all EMF-scored plants: same filters as r_code/21_colonization_all_scored_plants.R (all plants; controls)
# and r_code/20_emf_drought_influence.R (EMF-only colonization-biomass model)
e <- ds %>% filter(!is.na(emf))
cov_ok <- complete.cases(e[c("pop", "per_f", "treatment", "og_drc_mm", "soil", "myco_cat")])
fit161 <- e$alive_01_10_2021 %in% 1 & !is.na(e$dryshoot_new) & !is.na(e$og_drc_mm)
stopifnot(sum(cov_ok) == 165, sum(fit161) == 161, sum(cov_ok & e$myco_cat == "None") == 43)
why <- function(x) paste(sprintf("%s (%d)", names(x)[x > 0], x[x > 0]), collapse = "; ")
r165 <- c("Initial cutting diameter missing" = sum(is.na(e$og_drc_mm)), "Other covariate missing" = sum(!cov_ok & !is.na(e$og_drc_mm)))
r161 <- c("Died or survival status missing" = sum(!e$alive_01_10_2021 %in% 1), "Survivor without new stem mass" = sum(e$alive_01_10_2021 %in% 1 & is.na(e$dryshoot_new)),
          "Initial cutting diameter missing" = sum(e$alive_01_10_2021 %in% 1 & !is.na(e$dryshoot_new) & is.na(e$og_drc_mm)))
r43 <- c("Inoculated (AMF, EMF or dual)" = sum(e$myco_cat != "None"), "Covariate missing" = sum(e$myco_cat == "None" & !cov_ok))
stopifnot(sum(r165) == nrow(e) - 165, sum(r161) == nrow(e) - 161, sum(r43) == nrow(e) - 43)
s6 <- s6 %>% mutate(`Reasons for exclusion (plants)` = ifelse(Analysis == "Conductance, whole leaf",
                      sprintf("Non-positive mean (cannot be logged) (%d), leaving %d with a positive mean; Initial cutting diameter missing (1)", n_any - n_pos, n_pos), `Reasons for exclusion (plants)`))
s6 <- bind_rows(s6,
  s6_row("Leaf drop, measured root mass only", "Binomial GLMM", "Harvested survivors", nrow(hv), sum(!hv$sc & !is.na(hv$og_drc_mm)),
         sprintf("Colonization-scored (root subsample removed) (%d); Initial cutting diameter missing (%d)", sum(hv$sc), sum(!hv$sc & is.na(hv$og_drc_mm)))),
  s6_row("Allocation log-ratios, plants with leaves (a)", "LMM", "Allocation (a)", 414, alr$n[grepl("^\\(a\\)", alr$sample)], sprintf("No leaves at harvest (%d)", 414 - alr$n[grepl("^\\(a\\)", alr$sample)])),
  s6_row("Allocation log-ratios, plants with leaves (b)", "LMM", "Allocation (b)", 242, alr$n[grepl("^\\(b\\)", alr$sample)], sprintf("No leaves at harvest (%d)", 242 - alr$n[grepl("^\\(b\\)", alr$sample)])),
  s6_row("PERMANOVA, cross type", "PERMANOVA", "Plants with all three guilds scored", 78, nrow(foc78),
         paste(sprintf("%s (%d)", c("Backcross Narrowleaf" = "Backcross to narrowleaf", "Backcross Fremont" = "Backcross to Fremont")[names(bc)], as.integer(bc)), collapse = "; ")),
  s6_row("PERMANOVA, continuous ancestry", "PERMANOVA", "Plants with all three guilds scored", 78, 78, ""),
  s6_row("Colonization and PERMANOVA, inoculated plants only (Appendix S4)", "LMM; PERMANOVA", "Plants with all three guilds scored", 78, sum(d78$myco_cat != "None"),
         sprintf("Control inoculum (%d); cross type PERMANOVA: %d of the %d focal plants", sum(d78$myco_cat == "None"), sum(foc78$myco_cat != "None"), nrow(foc78))),
  s6_row("Conductance, whole leaf, plant-day model of record", "LMM", "Plants in the whole-leaf model", 203, cm$n_plants[cm$model == "reading"][1],
         sprintf("None; %d of %d plant-day means used (non-positive plant-day mean (%d))", n_rd, n_pd, n_pd - n_rd)),
  s6_row("Conductance, plants measured on 17 and 22 September", "LMM", "Plants in the whole-leaf model", 203, n1722, sprintf("Measured on another pair of dates (%d)", 203 - n1722)),
  s6_row("Conductance vs. EMF colonization", "LMM", "Plants in the whole-leaf model", 203, ce$n[grepl("^EMF", ce$predictor)], sprintf("Not scored for all three guilds (%d)", 203 - ce$n[grepl("^EMF", ce$predictor)])),
  s6_row("Conductance vs. root mass fraction", "LMM", "Plants in the whole-leaf model", 203, ce$n[grepl("^Root", ce$predictor)],
         sprintf("Not scored for all three guilds (%d); No root mass fraction (%d)", 203 - ce$n[grepl("^EMF", ce$predictor)], ce$n[grepl("^EMF", ce$predictor)] - ce$n[grepl("^Root", ce$predictor)])),
  s6_row("EMF colonization vs. ancestry, all EMF-scored plants", "LMM", "Plants scored for EMF", nrow(e), sum(cov_ok), why(r165)),
  s6_row("EMF-only colonization–biomass, EMF-scored survivors", "LMM", "Plants scored for EMF", nrow(e), sum(fit161), why(r161)),
  s6_row("EMF colonization vs. ancestry, control plants only", "LMM", "Plants scored for EMF", nrow(e), sum(cov_ok & e$myco_cat == "None"), why(r43)))
write_tab(s6, "appendix_S6_table2_sample_sizes.csv")

# ------------------------------------------------------------------------------
# Appendix S7: allocation of new growth and leaf drop, estimated root mass vs. measured root mass only
# ------------------------------------------------------------------------------
lrt_lab <- "Ancestry × drought (likelihood ratio test, 2 df)"
a7 <- rd("results/07_allocation_appendix_S7/appendix_S7_source.csv") %>%
  transmute(Response = recode(response, "log total new biomass" = "Total new biomass (log)", "log root:shoot" = "Root:shoot ratio (log)",
                              "log root:new-stem ratio" = "Root:new-stem ratio (log)", "P(no leaves at harvest)" = "No leaves at harvest (probability)"),
            Model = method, Term = ifelse(grepl("^Ancestry × drought", term), lrt_lab, term),
            n_a = n_estimated_root, s_a = string_estimated_root, n_b = n_measured_root_only, s_b = string_measured_root_only,
            p_a = p_estimated_root, p_b = p_measured_root_only, e_a = est_estimated_root, e_b = est_measured_root_only)
# log ratios in plants that kept their leaves (r_code/22_allocation_log_ratios.R)
lr <- rd("results/22_allocation_log_ratios/log_ratios.csv") %>% filter(response %in% c("log(root:new stem)", "log(leaf:new stem)")) %>%
  mutate(s = ifelse(startsWith(sample, "(a)"), "a", "b"),
         Response = recode(response, "log(root:new stem)" = "Root:new-stem ratio (log), plants with leaves", "log(leaf:new stem)" = "Leaf:new-stem ratio (log), plants with leaves"),
         Term = sub(" \\(per 1%\\)", "", term))
lr <- bind_rows(lr %>% distinct(s, Response, n, string = lrt_ancestry_x_drought) %>% mutate(Term = lrt_lab, P = p_of(string), estimate = NA),
                lr %>% select(s, Response, n, Term, string, P, estimate)) %>%
  pivot_wider(id_cols = c(Response, Term), names_from = s, values_from = c(n, string, P, estimate)) %>%
  transmute(Response, Model = "LMM", Term = factor(Term, levels = c(lrt_lab, "Linear ancestry", "Drought")),
            n_a, s_a = string_a, n_b, s_b = string_b, p_a = P_a, p_b = P_b, e_a = estimate_a, e_b = estimate_b) %>%
  arrange(desc(Response), Term) %>% mutate(Term = as.character(Term))
leaf <- a7$Response == "No leaves at harvest (probability)"
s7 <- bind_rows(a7[!leaf, ], lr, a7[leaf, ])
s7 <- bind_cols(s7 %>% transmute(Response, Model, Term, `n (estimated root mass)` = n_a, `Estimated root mass` = s_a,
                                 `n (measured root mass only)` = n_b, `Measured root mass only` = s_b),
                changes(s7$p_a, s7$p_b, s7$e_a, s7$e_b)) %>%
  mutate(Model = ifelse(Response == lag(Response, default = ""), "", Model), Response = blank_rep(Response))
stopifnot(nrow(s7) == 34)
write_tab(s7, "appendix_S7_allocation.csv")
