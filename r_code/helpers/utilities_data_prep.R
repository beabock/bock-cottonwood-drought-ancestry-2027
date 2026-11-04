# ==============================================================================
# UTILITIES AND DATA PREPARATION SCRIPT
# ==============================================================================
# This script handles all data loading, cleaning, and preparation tasks
# that are common across analysis scripts. All analysis scripts should
# source this file to ensure consistent data preparation.
# ==============================================================================

# Load required libraries for data preparation
library(dplyr)
library(ggplot2)
library(tidyr)
library(vegan)
library(patchwork)
library(purrr)

# Note: Individual analysis scripts will set up their own output sinks
# This script only handles data preparation

set.seed(1998)  # For reproducibility

# Watering-treatment colors (the figures themselves use r_code/helpers/figure_style.R, which defines the same colors)
drought_palette <- c("Well-watered" = "#3B9AB2", "Drought" = "#F2AD00")

# ==============================================================================
# DATA LOADING AND PREPARATION FUNCTION
# ==============================================================================

load_and_prepare_data <- function() {
  # Load data
  ds <- read.csv("data/edited_data/mega_ds.csv") %>%
    relocate(contains("note"),
             contains("id"),
             "treat",
             .after = last_col()) %>%
    mutate(
      
      new_class = dplyr::recode(new_class,
                         "NAR" = "Narrowleaf",
                         "BCN" = "Backcross Narrowleaf",
                         "F1H" = "F1 Hybrid",
                         "BCF" = "Backcross Fremont",
                         "FRE" = "Fremont",
                         .default = NA_character_),
      fremont_ancestry = case_when(
        per_f >= .70 ~ "Fremont-like",
        per_f < .70 ~ "Non-Fremont",
        TRUE ~ NA_character_
      ),
      fremont_ancestry = factor(fremont_ancestry, levels = c("Non-Fremont", "Fremont-like")),
      new_class = factor(new_class, levels = c("Narrowleaf", "Backcross Narrowleaf", "F1 Hybrid", "Backcross Fremont", "Fremont")),
      treatment = case_when(
        moist == "Well-Watered" ~ "Well-watered",
        moist == "Droughted" ~ "Drought",
        TRUE ~ NA_character_
      ),
      treatment = factor(treatment, levels = c("Well-watered", "Drought")),
      # Mycorrhizal categories
      myco_cat = case_when(
        myco == "None (Background)" ~ "None",
        myco == "AM inoc" ~ "AM",
        myco == "EM inoc" ~ "EM",
        myco == "Both" ~ "Both",
        TRUE ~ NA_character_
      ),
      myco_cat = factor(myco_cat, levels = c("None", "AM", "EM", "Both")),
      # Add size classes for group-level analyses
      height_class = cut(og_h_m, breaks = quantile(og_h_m, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                         labels = c("Low", "Medium", "High"), include.lowest = TRUE),
      drc_class = cut(og_drc_mm, breaks = quantile(og_drc_mm, probs = c(0, 1/4, 1/2, 3/4, 1), na.rm = TRUE),
                      labels = c("Tiny", "Small", "Medium", "Large"), include.lowest = TRUE),
      final_drc_class = cut(drc, breaks = quantile(drc, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                            labels = c("Small", "Medium", "Large"), include.lowest = TRUE)
    ) %>%
    # Use rowwise for calculations that need to be performed on each row
    rowwise() %>%
    mutate(
      # Calculate additional physiological variables for drought adaptation analysis
      # Stomatal density: average across all measured sides and leaves
      # --- DATA CLEANING FIX ---
      # Ensure physiological columns are numeric to prevent warnings in H3 loop
      across(matches("stom_dens|por_wet"), as.numeric),
      across(matches("^l_[0-9]+$"), as.numeric),
      stomatal_density_avg = mean(c_across(matches("^l1_.*_s[12]d$")), na.rm = TRUE),
      # Root architecture metrics: total length and surface area
      root_total_length = length,
      root_total_sa = sa,
      # Root diameter class distributions (ratios for normalization) - use corrected column names
      fine_root_ratio_length = (x0l_0_5l + x0_5l_1l + x1l_1_5l + x1_5_2l) / root_total_length,
      srl = root_total_length / root_dryweight,
      fine_srl = (x0l_0_5l + x0_5l_1l + x1l_1_5l + x1_5_2l)  / root_dryweight,
      coarse_root_ratio_length = (x2l_2_5l + x2_5_3l + x3l_3_5l + x3_5l_4l + x4_4_5l + x4_5l) / root_total_length,
      fine_root_ratio_sa = (x0l_0_5l + x0_5l_1l + x1l_1_5l + x1_5_2l)  / root_total_sa,
      coarse_root_ratio_sa = (x2sa_2_5sa + x2_5_3sa + x3sa_3_5sa + x3_5sa_4sa + x4_4_5sa + x4_5sa) / root_total_sa,
      # Leaf morphology: extract and average leaf areas (assumes columns named l_1, l_2, l_3, etc.)
      leaf_area_avg = mean(c_across(matches("^l_[0-9]+$") & where(is.numeric)), na.rm = TRUE),
      # Porometer leaf health score (average health across leaves, 1=healthy, 2=miscolored)
      across(starts_with("poro_l"), as.numeric),
      poro_health_avg = mean(c_across(starts_with("poro_l")), na.rm = TRUE),
      # Growth metrics
      resprout_count = branch_rs,
      flush_count = flush,
      per_f = per_f *100)%>%
    ungroup() # Return to standard dataframe operations

  # --- data fixes (2026-10-27) ---
  # Colonization can be scored only on survivors: SFF209 was recorded dead but has an EMF score, so it is dropped.
  # Growth (harvest mass) variables are kept only for plants recorded alive at harvest; this removes ICF908 and
  # JCF118 (recorded dead but weighed) and BRN910, ICN916 and BRF413 (no survival status). Survival is unaffected.
  growth_vars <- intersect(c("dryshoot_new", "dryleaf_weight", "root_dryweight", "wet_root", "dry_stem", "wet_stem",
                             "total_shoot", "total_biomass", "root_biomass_allocation"), names(ds))
  ds <- ds %>%
    mutate(across(c(emf, white_tip, brown_tip, dead_myco), ~ if_else(tree_idold == "SFF209", NA_real_, as.numeric(.x))),
           across(all_of(growth_vars), ~ if_else(alive_01_10_2021 %in% 1, as.numeric(.x), NA_real_)))

  return(ds)
}
