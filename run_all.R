# ==============================================================================
# Run every analysis behind the manuscript, in order, from the repository root:
#   Rscript run_all.R              # full run (about 20 min; step 12 needs internet)
#   Rscript run_all.R --skip-map   # offline: keeps the committed Appendix S1 map and collection-point files
#   Rscript run_all.R --from=12    # resume at step 12 after a failure (earlier outputs must be in place)
# Each script runs in its own R process (several scripts open sinks, set seeds and attach packages that would
# leak between scripts in one session). Stops at the first failure.
# Writes run_all_log.txt (exit code and seconds per script; not committed).
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
skip_map <- "--skip-map" %in% args
from <- sub("^--from=", "", grep("^--from=", args, value = TRUE))
stopifnot(file.exists("run_all.R"), dir.exists("r_code"), dir.exists("data"))

steps <- c(
  "r_code/01_merge_datasets.R",
  "r_code/02c_main_models_df.R",                  # sources 02b_main_models_followup.R, which sources 02a_main_models.R
  "r_code/03_conductance_leaf_sides.R",
  "r_code/04_conductance_surfaces.R",             # after 02: adds its block to results/02_main_models/formatted_results_main_models_df.md
  "r_code/05_control_inoculum_sensitivity.R",
  "r_code/06_allocation_models.R",
  "r_code/07_allocation_appendix_S7.R",
  "r_code/08_sample_flow_extra_tests.R",
  "r_code/09_allocation_survival_followup.R",
  "r_code/10_binomial_colonization_checks.R",
  "r_code/11_soil_legacy_permanova.R",            # about 6 min (9999-permutation PERMANOVAs)
  "r_code/12_map_appendix_S1.R",                  # downloads elevations and rivers (internet; about 5 min)
  "r_code/13_emf_ancestry_contrasts.R",
  "r_code/14_allocation_fraction_checks.R",
  "r_code/15_total_biomass_without_leaves.R",
  "r_code/16_soil_legacy_specifications.R",
  "r_code/17_inoculation_contrasts.R",
  "r_code/18_emf_amf_slopes_by_watering.R",
  "r_code/19_conductance_models_of_record.R",
  "r_code/20_emf_drought_influence.R",
  "r_code/21_colonization_all_scored_plants.R",
  "r_code/22_allocation_log_ratios.R",
  "r_code/23_new_stem_curvature.R",               # 2000 bootstrap refits
  "r_code/24_soil_missing_stems_survival_weights.R",
  "r_code/25_drought_rack_replication.R",
  "r_code/26_conductance_counts.R",
  "r_code/27_small_counts.R",
  "r_code/28_figures.R",                          # every figure except the Appendix S1 map
  "r_code/29_sensitivity_checks.R",
  "r_code/30_appendix_tables.R")                  # Appendix S1 and S4-S7 tables
if (skip_map) steps <- setdiff(steps, "r_code/12_map_appendix_S1.R")
if (length(from)) steps <- steps[as.integer(substr(basename(steps), 1, 2)) >= as.integer(from)]

# output folders (git does not keep empty folders, and several scripts open a sink without creating one)
out_dirs <- unique(sub("^02[abc]_.*", "02_main_models", sub("\\.R$", "", basename(steps))))
out_dirs <- setdiff(out_dirs, c("01_merge_datasets", "12_map_appendix_S1", "30_appendix_tables"))
for (d in c("data/edited_data", "figures", "results/appendices", file.path("results", out_dirs)))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

rscript <- file.path(R.home("bin"), "Rscript")
log_file <- "run_all_log.txt"
cat(R.version.string, if (length(from)) paste("| resumed at step", from), "\n", file = log_file, append = length(from) > 0)
t_all <- Sys.time()
for (s in steps) {
  t0 <- Sys.time()
  rc <- system2(rscript, shQuote(s), stdout = "", stderr = "")
  sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  line <- sprintf("%-55s exit %d  %5d s", s, rc, sec)
  cat(line, "\n"); cat(line, "\n", file = log_file, append = TRUE)
  if (rc != 0) stop("failed: ", s, call. = FALSE)
}
tot <- sprintf("TOTAL %d s", round(as.numeric(difftime(Sys.time(), t_all, units = "secs"))))
cat(tot, "\n"); cat(tot, "\n", file = log_file, append = TRUE)
