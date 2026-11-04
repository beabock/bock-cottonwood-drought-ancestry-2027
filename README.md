# Hybrid ancestry predicts survival, growth allocation, and ectomycorrhizal colonization of cottonwood cuttings regardless of drought

Beatrice M. Bock, Elena A. Schaefer, Maya Scull, Hillary F. Cooper, Arthur Keith, Gerard J. Allan, Kevin R. Hultine,
Christopher E. Doughty, Thomas G. Whitham, and Catherine A. Gehring

Data and R code for the manuscript (American Journal of Botany, in review). This repository holds the raw and merged
data, the scripts that produce every statistic, figure and appendix table reported in the paper, their outputs, and the
final figures.

**Repository:** https://github.com/beabock/bock-cottonwood-drought-ancestry-2027

**Archived version:** Zenodo, DOI: [ZENODO DOI PLACEHOLDER]

## Contents

| Folder or file | What it is |
|---|---|
| `data/raw_data/` | Raw greenhouse data (survival, biomass, colonization scores, porometer readings, root scans, leaf areas, soil CHN) and the collection record |
| `data/edited_data/winrhizaodata_final.csv` | Per-plant WinRHIZO root totals. **Hand-edited:** compiled by hand from the per-scan rows in `data/raw_data/gh_rootscandata_final.csv`; no script builds it. `r_code/01_merge_datasets.R` reads it as is. |
| `data/edited_data/mega_ds.csv` | Merged plant table (574 cuttings), written by `r_code/01_merge_datasets.R`. Cross types follow the Q-score cutoffs of Scull et al. (2025). |
| `data/DATA_DICTIONARY.md` | Every column of every data file |
| `r_code/` | Analysis scripts, numbered in run order; `r_code/helpers/` holds the files the scripts source (data preparation, shared model helpers, figure style) |
| `run_all.R` | Runs every script in order |
| `results/NN_name/` | Output of script `r_code/NN_name.R`: a text log (`output.txt` or `*_output.txt`) and the tables it writes |
| `results/appendices/` | The tables of Appendices S1 and S4–S7, one CSV per table (see below) |
| `results/sessionInfo.txt` | R session information recorded for the manuscript analyses |
| `figures/` | Final figures: `Fig1`–`Fig4`, `AppendixS1_map`, `FigS2`, `FigS3`, `AppendixS5_daypair` (Figure S5.1) and `AppendixS5_surfaces_daypair` (Figure S5.2), each as vector PDF and 600-dpi PNG |

## Running the analyses

From the repository root:

```bash
Rscript run_all.R
```

`run_all.R` runs each script in its own R process (the scripts open sinks, set seeds and attach packages that would leak
between scripts in one session) and stops at the first failure. A full run takes about 20 minutes on a Windows 11
laptop; most of it is the 9999-permutation PERMANOVAs in step 11 (about 6 minutes) and the downloads in step 12 (about
5 minutes). After a failure, `Rscript run_all.R --from=12` resumes at step 12.

**Internet access.** Step 12 (`r_code/12_map_appendix_S1.R`) downloads point elevations (USGS Elevation Point Query
Service), terrain tiles (`elevatr`), Natural Earth rivers and USGS NHD flowlines to build the Appendix S1 map and the
collection-point elevations. The elevation service sometimes times out; the script queries failed points again up to
five times. To run offline, use `Rscript run_all.R --skip-map`; the committed
`results/appendices/appendix_S1_collection_points.csv`, `appendix_S1_elevation_ranges.csv` and
`figures/AppendixS1_map.*` are then used as they are.

### Run order

Script `r_code/NN_name.R` writes its output to `results/NN_name/` (the three step-02 scripts write to
`results/02_main_models/`). The exceptions are the publication files: the appendix tables in `results/appendices/`
(steps 12 and 30) and the figures in `figures/` (steps 12 and 28). Scripts 02a and 02b are not run on their own: 02c
sources 02b, which sources 02a.

| Step | Script | What it does |
|---|---|---|
| 01 | `01_merge_datasets.R` | Builds the merged plant table `data/edited_data/mega_ds.csv` |
| 02 | `02a_main_models.R`, `02b_main_models_followup.R`, `02c_main_models_df.R` | Survival, biomass, colonization and PERMANOVA models; genotype-level permutation tests; degrees of freedom |
| 03 | `03_conductance_leaf_sides.R` | Stomatal conductance by leaf side and measurement day |
| 04 | `04_conductance_surfaces.R` | Surface-specific conductance (plant means); survival grouping sensitivity |
| 05 | `05_control_inoculum_sensitivity.R` | Main colonization results without control plants (Appendix S4) |
| 06 | `06_allocation_models.R` | Allocation of new growth, leaf loss and estimated root mass |
| 07 | `07_allocation_appendix_S7.R` | Allocation models with estimated and with measured root mass (source of Appendix S7) |
| 08 | `08_sample_flow_extra_tests.R` | Sample flow (Table S6.2), survival by group, inoculum effects, EMF simple slopes |
| 09 | `09_allocation_survival_followup.R` | Root:new-stem ratio, mass fractions, predicted survival, counts |
| 10 | `10_binomial_colonization_checks.R` | Colonization as a binomial GLMM on raw counts; survival with population random |
| 11 | `11_soil_legacy_permanova.R` | Robustness of the soil-legacy PERMANOVA result (about 6 min) |
| 12 | `12_map_appendix_S1.R` | Appendix S1 map and collection-point elevations (internet; about 5 min) |
| 13 | `13_emf_ancestry_contrasts.R` | Estimates for each version of the EMF ancestry result |
| 14 | `14_allocation_fraction_checks.R` | Agreement of the stem, leaf and root fraction models |
| 15 | `15_total_biomass_without_leaves.R` | Total new biomass without leaves |
| 16 | `16_soil_legacy_specifications.R` | Soil legacy on colonization across specifications |
| 17 | `17_inoculation_contrasts.R` | Inoculum-versus-control contrasts and their multiplicity |
| 18 | `18_emf_amf_slopes_by_watering.R` | EMF minus AMF ancestry slope by watering treatment |
| 19 | `19_conductance_models_of_record.R` | Conductance models with measurement dates accounted for (Appendix S5) |
| 20 | `20_emf_drought_influence.R` | EMF × drought on new stem biomass: influence of single plants and genotypes |
| 21 | `21_colonization_all_scored_plants.R` | Single-guild colonization models on every scored plant |
| 22 | `22_allocation_log_ratios.R` | Log ratios, predicted masses, pooled dry:fresh ratio |
| 23 | `23_new_stem_curvature.R` | Bootstrap of the new-stem peak; inventory of quadratic terms |
| 24 | `24_soil_missing_stems_survival_weights.R` | Soil legacy, survivors without new stem mass, survivorship weighting |
| 25 | `25_drought_rack_replication.R` | Drought tested against the number of racks |
| 26 | `26_conductance_counts.R` | Appendix S5 counts and consistency |
| 27 | `27_small_counts.R` | Holm family of inoculum contrasts; mother-tree counts |
| 28 | `28_figures.R` | Every figure except the Appendix S1 map |
| 29 | `29_sensitivity_checks.R` | EMF-scored plant counts, predicted EMF colonization, quadratic-term summary, models refit without plant SFF111 |
| 30 | `30_appendix_tables.R` | The tables of Appendices S1 and S4–S7 |

Steps 20–29 share `r_code/helpers/model_helpers.R`; each writes `output.txt` and a `strings.csv` with the formatted
values it computes.

### Where each figure and appendix table comes from

| Item | Script | File |
|---|---|---|
| Figures 1–4 | `28_figures.R` | `figures/Fig1`–`Fig4` |
| Appendix S1 map | `12_map_appendix_S1.R` | `figures/AppendixS1_map` |
| Tables S1.1, S1.2 | `30_appendix_tables.R` | `results/appendices/appendix_S1_table1_mother_trees.csv`, `appendix_S1_table2_collection_points.csv` |
| Appendix S2 and S3 figures | `28_figures.R` | `figures/FigS2`, `FigS3` |
| Appendix S4 table | `05`, `18`, `30_appendix_tables.R` | `results/appendices/appendix_S4_control_sensitivity.csv` |
| Tables S5.1–S5.4 | `19`, `30_appendix_tables.R` | `results/appendices/appendix_S5_table1_conductance_models.csv`, `appendix_S5_table2_surface_models.csv`, `appendix_S5_table3_conductance_vs_emf_rmf.csv`, `appendix_S5_table4_plants_by_dates.csv` |
| Figures S5.1, S5.2 | `28_figures.R` | `figures/AppendixS5_daypair`, `figures/AppendixS5_surfaces_daypair` |
| Table S6.2 | `08`, `30_appendix_tables.R` | `results/appendices/appendix_S6_table2_sample_sizes.csv` |
| Appendix S7 table | `07`, `22`, `30_appendix_tables.R` | `results/appendices/appendix_S7_allocation.csv` |

Table S6.1 (each analysis and how it departs from the core model) describes the models and is given in the manuscript
only.

All figures share `r_code/helpers/figure_style.R`: `theme_bw()`, Arial, AJB sizes, and one color rule (blue =
well-watered, gold = drought, with circle/solid line and triangle/dashed line as redundant cues; everything else is gray
or black). The blue/gold pair stays distinct under simulated deuteranopia, protanopia and tritanopia.

### Reproducibility

Run-time stamps are not printed, so a rerun on the same R and package versions reproduces the committed outputs byte
for byte (PDFs apart from the creation date that cairo writes).

## Requirements

- **R 4.6.1** (the version used for the manuscript; `results/sessionInfo.txt`).
- Packages, with the versions used for the manuscript analyses (`results/sessionInfo.txt`):
  lme4 2.0-6, lmerTest 3.2-1, glmmTMB 1.1.15.2, nlme 3.1-169, emmeans 2.0.4, vegan 2.7-6, permute 0.9-10,
  pairwiseAdonis 0.4.1, Matrix 1.7-5, ggplot2 4.0.3, ragg 1.5.2, tidyverse 2.0.0 (dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0,
  purrr 1.2.2, stringr 1.6.0, forcats 1.0.1, tibble 3.3.1, lubridate 1.9.5), readxl 1.6.0, janitor 2.2.1,
  patchwork 1.3.2, systemfonts 1.3.2, MASS 7.3-65.
- Not in `sessionInfo.txt` (loaded only by individual scripts): DHARMa 0.5.0 (step 06), and for the map (step 12)
  sf 1.1.3, terra 1.9.50, elevatr 0.99.1, ggspatial 1.1.11, rnaturalearth 1.2.0, nhdplusTools 1.5.2, ggrepel 0.9.8,
  maps 3.4.3.
- `pairwiseAdonis` is not on CRAN. Install it from GitHub at the commit used for the manuscript analyses:

  ```r
  remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis@cb190f7668a0c82c0b0853927db239e7b9ec3e83")
  ```
- vegan ≥ 2.7 changed the default of `adonis2()` to an omnibus test; every script passes `by = "terms"` or
  `by = "margin"` explicitly.
- Figures use Arial (`ragg` and `cairo_pdf`); other fonts change the rendering slightly.

## Licenses

- Code (`r_code/`, `run_all.R`): MIT License, see `LICENSE`.
- Data, statistical output and figures (`data/`, `results/`, `figures/`): Creative Commons Attribution 4.0
  International (CC BY 4.0), see `LICENSE-DATA`.

## Citation

Bock, B. M., E. A. Schaefer, M. Scull, H. F. Cooper, A. Keith, G. J. Allan, K. R. Hultine, C. E. Doughty, T. G. Whitham,
and C. A. Gehring. Hybrid ancestry predicts survival, growth allocation, and ectomycorrhizal colonization of cottonwood
cuttings regardless of drought. American Journal of Botany (in review). Data and code:
https://github.com/beabock/bock-cottonwood-drought-ancestry-2027; archived on Zenodo, DOI: [ZENODO DOI PLACEHOLDER].

## Contact

Beatrice M. Bock, Beatrice.bock@usda.gov (ORCID 0000-0003-2240-9360)
