# ==============================================================================
# Appendix S5 consistency.
#   8b  surface-specific reading-level values (results/19_conductance_models_of_record/surface_models.csv) and the surface difference in the
#       plant-mean joint model (r_code/30_appendix_tables.R, Table S5.2), averaged over watering and well-watered only,
#       so the two published surface-difference values can be compared like for like.
#   8c  counts: plants with readings, non-positive surface means, plant-days, and the excluded whole-leaf plant.
# Data as r_code/19_conductance_models_of_record.R (raw porometer file joined to the plants).
# Run from the project root:  Rscript r_code/26_conductance_counts.R
# Writes results/26_conductance_counts/: output.txt, strings.csv, surface_reading_level.csv, conductance_counts.csv
# ==============================================================================
source("r_code/helpers/model_helpers.R")
OUT <- "results/26_conductance_counts"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "output.txt"), split = TRUE); run_stamp("Conductance counts")
src <- "results/26_conductance_counts/conductance_counts.csv"
norm <- function(x) toupper(sub("^0+", "", as.character(x)))
d0 <- ds %>% mutate(key = norm(poro_id))
rd <- read_csv("data/raw_data/gh_poromeasurements_allgroups.csv", col_types = cols(poro_id = "c"), show_col_types = FALSE) %>%
  mutate(key = norm(poro_id), date = as.Date(date, format = "%m/%d/%Y")) %>% rowwise() %>%
  mutate(day_mean = mean(c_across(l1_t:l3_b), na.rm = TRUE), adaxial = mean(c(l1_t, l2_t, l3_t), na.rm = TRUE), abaxial = mean(c(l1_b, l2_b, l3_b), na.rm = TRUE),
         n_neg = sum(c_across(l1_t:l3_b) <= 0, na.rm = TRUE)) %>% ungroup()
plants <- d0 %>% filter(!is.na(por_avg)) %>% select(tree_idold, key, og_drc_mm, por_avg)
pd <- rd %>% inner_join(plants, by = "key"); stopifnot(n_distinct(pd$tree_idold) == 205)
pm <- pd %>% group_by(tree_idold, og_drc_mm, por_avg) %>% summarise(n_days = n_distinct(date), adaxial = mean(c(l1_t, l2_t, l3_t), na.rm = TRUE),
                                                                    abaxial = mean(c(l1_b, l2_b, l3_b), na.rm = TRUE), n_neg = sum(n_neg), .groups = "drop")
cnt <- list(); k <- function(stage, n, note = "") cnt[[length(cnt) + 1]] <<- tibble(stage, n, note)
k("Plants with porometer readings", nrow(pm))
k("  whole-leaf mean > 0", sum(pm$por_avg > 0), paste(pm$tree_idold[pm$por_avg <= 0], collapse = ", "))
k("  ... and initial diameter recorded (whole-leaf model)", sum(pm$por_avg > 0 & !is.na(pm$og_drc_mm)), paste(pm$tree_idold[pm$por_avg > 0 & is.na(pm$og_drc_mm)], collapse = ", "))
an <- pm %>% filter(por_avg > 0, !is.na(og_drc_mm))
for (s in c("adaxial", "abaxial")) {
  k(sprintf("Non-positive %s mean, among all %d plants with readings", s, nrow(pm)), sum(pm[[s]] <= 0, na.rm = TRUE), paste(pm$tree_idold[which(pm[[s]] <= 0)], collapse = ", "))
  k(sprintf("Non-positive %s mean, among the 203 analyzed plants", s), sum(an[[s]] <= 0, na.rm = TRUE), paste(an$tree_idold[which(an[[s]] <= 0)], collapse = ", ")) }
ex <- pm %>% filter(por_avg <= 0)
k("Excluded whole-leaf plant: negative readings / all readings", sum(ex$n_neg), sprintf("%s: whole-leaf mean %.4f; %d of %d readings <= 0", ex$tree_idold, ex$por_avg, sum(ex$n_neg),
                                                                                      sum(!is.na(as.matrix(pd %>% filter(tree_idold %in% ex$tree_idold) %>% select(l1_t:l3_b))))))
pda <- pd %>% filter(tree_idold %in% an$tree_idold)
k("Analyzed plants measured on both days", sum(an$n_days == 2)); k("Analyzed plants measured on one day", sum(an$n_days == 1))
k("Plants with one day only, among all 205 with readings", sum(pm$n_days == 1))
k("Plant-days among the 203 analyzed plants", nrow(pda), sprintf("= 2 x %d + %d", sum(an$n_days == 2), sum(an$n_days == 1)))
k("  plant-day means <= 0 (excluded from the reading-level model)", sum(pda$day_mean <= 0), paste(sprintf("%s %s (%.4f)", pda$tree_idold[pda$day_mean <= 0], pda$date[pda$day_mean <= 0], pda$day_mean[pda$day_mean <= 0]), collapse = "; "))
k("  plant-days in the reading-level model", sum(pda$day_mean > 0))
k("Surface reading level: adaxial plant-days > 0 / plants", sum(pda$adaxial > 0, na.rm = TRUE), sprintf("%d plants", n_distinct(pda$tree_idold[pda$adaxial > 0])))
k("Surface reading level: abaxial plant-days > 0 / plants", sum(pda$abaxial > 0, na.rm = TRUE), sprintf("%d plants", n_distinct(pda$tree_idold[which(pda$abaxial > 0)])))
cnt <- bind_rows(cnt); print(as.data.frame(cnt), right = FALSE, row.names = FALSE); readr::write_excel_csv(cnt, file.path(OUT, "conductance_counts.csv"))
for (i in seq_len(nrow(cnt))) rec(8, "c", "Conductance counts", "porometer plants", NA, trimws(cnt$stage[i]), paste0(cnt$n[i], ifelse(cnt$note[i] == "", "", paste0(" (", cnt$note[i], ")"))), NA, src)

hdr("8b. Surface-specific values and the surface difference, by model")
s10s <- read_csv("results/19_conductance_models_of_record/surface_models.csv", show_col_types = FALSE)
print(as.data.frame(s10s %>% select(surface, n_plants, n_obs, term, string)), right = FALSE, row.names = FALSE)
# plant-mean joint model as Table S5.2 (r_code/30_appendix_tables.R)
side_mean <- function(d, side) { cc <- as.vector(outer(paste0("l", 1:3, "_", side), c("wet", "dry"), paste, sep = "_"))
  m <- rowMeans(as.matrix(sapply(d[, cc], as.numeric)), na.rm = TRUE); m[is.nan(m)] <- NA; m }
por <- ds %>% filter(!is.na(por_avg)) %>% mutate(adaxial = side_mean(., "t"), abaxial = side_mean(., "b"))
long <- por %>% select(tree_id, pop, geno_id, per_fc, treatment, og_drc_mm, soil, myco_cat, adaxial, abaxial) %>% pivot_longer(c(adaxial, abaxial), names_to = "surface", values_to = "gsw") %>%
  filter(!is.na(gsw), gsw > 0) %>% mutate(surface = factor(surface, levels = c("abaxial", "adaxial")))
joint <- lmerTest::lmer(log(gsw) ~ surface * (per_fc * treatment + I(per_fc^2)) + og_drc_mm + soil + myco_cat + (1 | pop) + (1 | pop:geno_id) + (1 | tree_id), data = long)
ww <- lin(joint, c(`surfaceadaxial:per_fc` = 1)); avg <- lin(joint, c(`surfaceadaxial:per_fc` = 1, `surfaceadaxial:per_fc:treatmentDrought` = 0.5))
stopifnot(abs(ww$est + 0.0237) < 0.00005)
sd8 <- bind_rows(
  s10s %>% transmute(model = "Reading level (model of record): plant-day means, phase and day pair fixed, plant random", surface, n_plants, n_obs, term, string),
  tibble(model = "Plant-mean joint model (Table S5.2)", surface = "adaxial − abaxial", n_plants = n_distinct(model.frame(joint)$tree_id), n_obs = nobs(joint),
         term = c("Difference in ancestry slope, well-watered (Table S5.2 row)", "Difference in ancestry slope, averaged over watering"), string = c(ww$s, avg$s)))
print(as.data.frame(sd8 %>% select(model, surface, term, string)), right = FALSE, row.names = FALSE)
readr::write_excel_csv(sd8, file.path(OUT, "surface_reading_level.csv"))
for (i in seq_len(nrow(sd8))) rec(8, "b", sd8$model[i], sd8$surface[i], sprintf("%s plants, %s obs", sd8$n_plants[i], sd8$n_obs[i]), sd8$term[i], sd8$string[i], NA, "results/26_conductance_counts/surface_reading_level.csv")
write_rec("strings.csv"); cat("\nDone.\n"); sink()
