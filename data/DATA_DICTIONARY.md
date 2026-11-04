# Data dictionary

Every column of every file in `data/`. Definitions come from the project metadata workbook
(`raw_data/metadata_ghcottieproject_summer2021.xlsx`), from how `r_code/01_merge_datasets.R` and the analysis
scripts use each column, and from the values themselves.

One detail is not known: the magnification or field area of the stomatal counts in `gh_stomdens.csv` (not used in the manuscript).

Conventions used throughout:

- **Plant IDs.** Each cutting has an *old* ID (`tree_idold`, 6–7 characters, e.g. `BRF203`) and a *new* ID (`tree_id`,
  7 characters, e.g. `BRF0203`). Both are population code + parental type + genotype number + replicate number. The new ID
  zero-pads the genotype number to two digits. Note that in `gh_allinfo.csv` the column named `tree_id` holds the **old** ID
  (the merge script renames it to `tree_idold`).
- **Population codes (`pop`).** BR Blue River (Alpine, AZ); IC Indian Creek (Monticello, UT); JC Jack's Canyon (Winslow, AZ);
  LC Little Colorado River (Springerville, AZ); OC Oak Creek (Sedona, AZ); SF San Francisco River (Reserve, NM).
- **Parental type (`type`).** Field classification of the mother tree from leaf morphology: F Fremont, H hybrid, N narrowleaf.
  The analyses use the Q-score classes in `new_xtype.csv` instead.
- **Mother tree (`old_class`).** Population + type + genotype number, e.g. `BRF2`. This is the genotype (mother-tree) ID.
- **Treatment code (`treat`).** Three letters: watering (D drought, W well-watered) + soil legacy (L tamarisk legacy, N no
  legacy / cottonwood soil) + inoculum (A AMF, E EMF, B both, N none). `?`, `???` or a trailing `?` mark a code that was
  uncertain on the data sheet; the merge script drops these plants.
- Missing values are empty cells or `NA`.

---

## raw_data/

### gh_allinfo.csv
One row per cutting (900 rows): pre-treatment and harvest measurements transcribed from the handwritten harvest sheets.
Read by `01_merge_datasets.R` and by several robustness scripts.

| Column | Meaning |
|---|---|
| `tree_id` | Old plant ID (despite the name; see conventions). One ID (`BRF413`) appears twice with conflicting metadata and is dropped by the merge script. |
| `pop` | Population code. |
| `type` | Parental type of the mother tree (F, H, N). |
| `geno` | Genotype number of the mother tree within population × type. |
| `tree` | Replicate (cutting) number within the mother tree. |
| `og_h(m)` | Height before treatments, May 2021, cm (the "(m)" in the header is a typo). |
| `og_drc(mm)` | Basal (root-collar) stem diameter before treatments, May 2021, mm. This is the "initial cutting diameter" covariate in every model. |
| `og_leaf` | Leaf health score in May 2021: 1 good (green, healthy) to 5 poor (yellow or pink leaves); `?` unknown. |
| `og_notes` | Notes from May 2021 (e.g. resprout, number of leaves). |
| `alive_preinnoc` | Alive in May 2021, before inoculation (1 yes, 0 no; `MIA` = cutting missing, not found). |
| `treat` | Three-letter treatment code. |
| `moist` | Watering treatment: D drought (watered every fourth day), W well-watered (every other day); `?` unknown. |
| `soil` | Soil legacy: L tamarisk legacy, N no legacy (cottonwood soil); `?` unknown. |
| `myco` | Inoculum: A AMF, E EMF, B both (dual), N none (control); `?`/`B?` uncertain. |
| `alive_01_10_2021` | Alive at harvest, 1 October 2021 (1 yes, 0 no). This is the survival response. |
| `drc` | Basal stem diameter at harvest, mm. |
| `h` | Height at harvest, cm. |
| `wet_root` | Fresh mass of the whole root system at harvest, g (gently blotted). Used with the dry:fresh ratio to estimate root dry mass of plants whose roots were subsampled for colonization. |
| `wet_stem` | Fresh mass of the stem (no leaves or roots) at harvest, g. |
| `flush` | Number of growth flushes, counted from bud scars; `?` unknown. |
| `branch_rs` | Number of branches or resprouts from the cutting. |
| `poro_id` | Porometer plant ID; links to `gh_poromeasurements_allgroups.csv`, `gh_stomdens.csv` and the soil workbook. |
| `poro_l1`, `poro_l2`, `poro_l3` | Condition of porometer leaf 1, 2 and 3: 1 green and healthy, 2 discoloured; `?` unknown, `MIA` leaf missing, 3 = replacement leaf (the original labelled leaf was lost and another leaf was measured). |
| `notes` | Harvest notes (e.g. pink stems, leaves only on resprouts). |

### gh_cuttingskey.csv
Key from old to new plant IDs with treatment assignments (946 rows; 100 rows are exact duplicates, removed by the merge script).

| Column | Meaning |
|---|---|
| `tree_idold` | Old plant ID. |
| `old_class` | Mother tree (population + type + genotype), e.g. `BRF2`. Links to `new_xtype.csv`. |
| `pop` | Population code. |
| `type` | Parental type (F, H, N). |
| `geno` | Genotype number. |
| `tree` | Replicate number. |
| `tree_id` | New plant ID (7 characters). Links to `winrhizaodata_final.csv` and `gh_rootscandata_final.csv`. |
| `treat` | Three-letter treatment code. |
| `moist` | Watering treatment (D, W). |
| `soil` | Soil legacy (L, N). |
| `myco` | Inoculum (A, E, B, N). |

### new_xtype.csv
Ancestry of each mother tree from microsatellite Q-scores (Structure, K = 2; Scull et al. 2025). One row per mother tree (183 rows, including trees not used in the experiment).

| Column | Meaning |
|---|---|
| `old_class` | Mother tree ID (population + type + genotype). |
| `per_f` | Proportion of Fremont ancestry, 1 − Q (0–1; the scripts multiply by 100 for percent Fremont ancestry). |
| `per_n` | Proportion of narrowleaf ancestry, Q (0–1). |
| `new_class` | Cross type from the Q cutoffs: `FRE` Fremont (Q < 0.1), `BCF` backcross to Fremont (0.1–0.4), `F1H` F₁ hybrid (0.4–0.6), `BCN` backcross to narrowleaf (0.6–0.9), `NAR` narrowleaf (Q > 0.9). |
| `old_xtype` | Field classification before genotyping: `FRE`, `HYB`, `NAR`. |

### gh_drybiomass_aboveground.csv
Aboveground dry mass at harvest (515 rows).

| Column | Meaning |
|---|---|
| `tree_identification` | Old-style label from the original data sheet (e.g. `BRF2-03`); an artifact kept from the original labelling, not used. The header in the file has a trailing space (`tree_identification `). |
| `tree_idold` | Old plant ID. |
| `pop` | Population code. |
| `type` | Parental type. |
| `geno` | Genotype number. |
| `indiv` | Replicate (cutting) number. |
| `wet_stem` | Fresh stem mass (no leaves or roots), g. |
| `flush` | Number of flushes from bud scars; `?` unknown. |
| `branch_rs` | Number of branches or resprouts. |
| `dryleaf_weight` | Dry mass of all leaves present at harvest, g. Blank for harvested survivors with no leaves; the analyses set leaf mass to 0 g for those plants. |
| `dry_shoottotal` | Dry mass of the whole stem (original cutting plus new growth; no leaves or roots), g. Renamed `dry_stem` in `mega_ds.csv`. |
| `dryshoot_new` | Dry mass of new stem growth produced during the experiment (separated at the bud scars), g. The "new stem biomass" response. |

### gh_rootdryweight_FINAL.csv

| Column | Meaning |
|---|---|
| `tree_idold` | Old plant ID. |
| `root_dryweight` | Dry mass of the root system, g (65 °C to constant mass). For plants scored for colonization this excludes the fresh subsample removed for scoring. |

### gh_amscore_raw.csv
Arbuscular mycorrhizal and DSE scoring by the grid-line intersect method (about 100 intersections per sample, compound microscope, 200×). Counts per sample.

| Column | Meaning |
|---|---|
| `tree_idold` | Old plant ID (one duplicated ID; the merge script keeps the first non-missing value of each column). |
| `arb` | Intersections with arbuscules. |
| `ves` | Intersections with vesicles. |
| `hyph` | Intersections with AMF hyphae only. |
| `dse` | Intersections with dark septate endophytes (melanized septate hyphae, microsclerotia). |
| `root` | Intersections with root but no fungal structure. |

The merge script rescales each row to sum to 100 and defines AMF colonization (%) as `arb + ves + hyph`.

### gh_emscore_raw.csv
Ectomycorrhizal scoring of root tips at grid-line intersections under a dissecting microscope (40×). Counts per sample (row totals 92–131). Morphotypes defined by Elena A. Schaefer for Fremont cottonwood.

| Column | Meaning |
|---|---|
| `tree_idold` | Old plant ID. |
| `root_alive` | Live root tips without an ectomycorrhizal mantle (non-colonized). |
| `root_dead` | Dead root tips (all 0). |
| `dark_finger` | EMF morphotype "dark finger". |
| `brown_bunch` | EMF morphotype "brown bunch". |
| `dark_thumb` | EMF morphotype "dark thumb". |
| `white_tip` | EMF morphotype "white tip". |
| `brown_tip` | EMF morphotype "brown tip". |
| `dead_myco` | Dead ectomycorrhizal tips. |

The merge script rescales each row to sum to 100 and defines EMF colonization (%) as the sum of the six mycorrhizal columns (`dark_finger` through `dead_myco`).

### gh_poromeasurements_allgroups.csv
Stomatal conductance (gsw) from a LI-600 porometer, around midday, on three labelled leaves per plant (422 rows: one per plant per measurement day).

| Column | Meaning |
|---|---|
| `poro_id` | Porometer plant ID (links to `gh_allinfo.csv`). |
| `water` | Point in the watering cycle: `dry` end of a watering cycle (driest day), `wet` after the next watering. |
| `date` | Measurement date (M/D/YYYY), 8–26 September 2021. The four pairs of dates are the measurement sessions in Appendix S5. |
| `l1_t`, `l2_t`, `l3_t` | Conductance of leaf 1, 2, 3, top (adaxial) surface, mol m⁻² s⁻¹. Small negative values are instrument readings near zero. |
| `l1_b`, `l2_b`, `l3_b` | Conductance of leaf 1, 2, 3, bottom (abaxial) surface, mol m⁻² s⁻¹. |

### gh_stomdens.csv
Stomatal counts from leaf impressions on a microscopy subsample (60 rows; 33 with counts). Not used in the manuscript.

| Column | Meaning |
|---|---|
| `tree_idold` | Old plant ID. Covers only part of the sample; the merge script links this file by `poro_id` instead. |
| `poro_id` | Porometer plant ID. |
| `l1_t_s1d` | Leaf 1, top (adaxial) surface, stomatal density, field 1. Probably a count per microscope field of view. Magnification or field area not recorded; values (abaxial median about 150) are in the range of *Populus* densities per mm², so either the counts are already per mm² or the field is about 1 mm². |
| `l1_b_s1d` | Leaf 1, bottom (abaxial) surface, stomatal density, field 1 (unit as for `l1_t_s1d`). |
| `l1_t_s2d` | Leaf 1, top surface, stomatal density, field 2 ("side 2" in the metadata; unit as for `l1_t_s1d`). |
| `l1_b_s2d` | Leaf 1, bottom surface, stomatal density, field 2 (unit as for `l1_t_s1d`). |
| `notes` | Observations (trichomes, stomatal size). |

### gh_sa_raw_1027.csv
Area of each leaf at harvest, measured from scans in ImageJ (429 rows; six plants have two rows and are dropped by the merge script). Not used in the manuscript results.

| Column | Meaning |
|---|---|
| `tree_idold` | Old plant ID. |
| `tree_id` | New plant ID. |
| `poro_id` | Porometer plant ID for the subsample measured with the porometer. This column is independent of `poro_id` in `gh_allinfo.csv`; the merge script blanks any `poro_id` that ends up assigned to more than one plant. |
| `type` | Parental type. |
| `treat` | Three-letter treatment code (`?` marks uncertain codes). |
| `l_1` … `l_57` | One-sided area of leaf 1 … 57 of the plant, cm². A plant has as many filled columns as it had leaves. |

### gh_rootscandata_final.csv
WinRHIZO output, one row per root scan (790 rows; large root systems were scanned in several parts, so a plant can have several rows). Not read by any script; `edited_data/winrhizaodata_final.csv` holds the per-plant totals. Diameter classes are 0.5 mm wide (class 1 = 0–0.5 mm … class 10 = > 4.5 mm).

| Column | Meaning |
|---|---|
| `treeidentification` | Scan label (lower case, e.g. `brf_02_06`). |
| `tree_id` | New plant ID. |
| `AvgDiam(mm)` | Average root diameter of the scan, mm. |
| `RootVolume(cm3)` | Root volume of the scan, cm³. |
| `Tips` | Number of root tips (WinRHIZO count). |
| `Forks` | Number of root forks (branch points, WinRHIZO count). |
| `Crossings` | Number of root crossings (overlaps, WinRHIZO count). |
| `l1` … `l10` | Root length in diameter classes 1–10, cm. |
| `sa1` … `sa10` | Root surface area in diameter classes 1–10, cm². |
| `0<.PA.<=0.5000000` … `.PA.>4.5000000` (10 columns) | Root projected area in each diameter class (0–0.5 mm … > 4.5 mm), cm². |
| `0<.V.<=0.5000000` … `.V.>4.5000000` (10 columns) | Root volume in each diameter class, cm³. |
| `0<.T.<=0.5000000` … `.T.>4.5000000` (10 columns) | Number of tips in each diameter class. |

Length, area and volume are in WinRHIZO's default units (cm, cm², cm³).

### gh_soildata_rawfromASU.xlsx
Carbon, hydrogen and nitrogen of soil samples (Perkin Elmer Series II CHN analyzer, ASU METAL core). Not used in the manuscript; `01_merge_datasets.R` joins the "Final Data" sheet to plants by `poro_id`.

**Sheet "Raw Data"** (instrument export; two header rows):

| Column | Meaning |
|---|---|
| `Run` | Sample or standard label (`K1` acetanilide standard, `BLANK`, `1S`–`16S` soil samples). |
| `Run #` | Instrument run number. |
| `Weight` | Sample mass, mg. |
| `Created on` | Instrument field (all 0). |
| `Mode` | Analysis mode (`CHN`). |
| `Comment` | Instrument field (all 0). |
| `Results: Carbon`, `Hydrogen`, `Nitrogen` | Measured content, % of sample mass (for blanks, raw signal). |
| `Signals: ZR`, `CR`, `HR`, `NR` | Detector signals (zero reading, carbon, hydrogen, nitrogen). |
| `Blank Values: Carbon`, `Hydrogen`, `Nitrogen` | Blank signals subtracted. |
| `Fill Time: Seconds` | Instrument fill time, s. |
| `Messages` | Instrument messages. |

**Sheet "Final Data"** (title rows, then one row per soil sample):

| Column | Meaning |
|---|---|
| `Run` | Sample label (`1S`–`16S`). |
| `sample_id` | Plant whose pot soil was sampled (new ID, lower case). Blank for the six field-soil samples `11S`–`16S`. |
| `treat` | Treatment code of that plant, lower case. Renamed `soil_chem_treat` in `mega_ds.csv`. |
| `poro_id` | Porometer ID of that plant (join key). |
| `note` | Sample type code: first part `nt` no tamarisk / `ta` tamarisk soil. Second part: `st` steam-sterilized, `ns` not sterilized; third part: `ad` air-dried, `od` oven-dried. |
| `leg` | Soil source: `nt` no tamarisk, `ta` tamarisk. |
| `w` | Soil-legacy treatment of the potted plant: `leg` legacy, `non` no legacy. |
| `Weight (mg)` | Sample mass, mg. |
| `%Carbon`, `%Hydrogen`, `%Nitrogen` | Content, % of sample mass. |
| `C:N` | Carbon:nitrogen ratio (%Carbon / %Nitrogen). |
| two unnamed columns | Empty. |

**Sheet "Sheet1"**: pivot table of mean C:N by soil source (`leg` = `nt`, `ta`).
**Sheet "QC"**: replicate runs of the acetanilide standard (`Run`, `Weight`, `%Carbon`, `%Hydrogen`, `%Nitrogen`) with their standard deviation (`Std`), mean (`Avg`) and relative standard deviation (`%RSD`).

### WC Garden Collections_Coords.xlsx
Collection record of the 2021 hybrid-garden cuttings (Hillary F. Cooper). Read by `12_map_appendix_S1.R` (sheet "Sheet1" only).

- **Sheet "Metadata"**: one note (which collections were used in the greenhouse experiment).
- **Sheet "Sheet1"** (collection record; the source of Appendix S1): `Populations` (river and nearest town; filled on the first row of each population), `Date Collected` (date range), `Crosstypes (F,H,N)` (Fremont, Hybrid, Narrowleaf), `# of Genotypes` (genotypes collected of that cross type), `GPS coord.` (degrees and decimal minutes), and a notes column (e.g. "Unknown if this data sheet is correct for this experiment"). Two further columns are empty.
- **Sheet "Sources for GH Expt"**: the same columns for the four populations listed there.
- **Sheet "GPS Coords"**: `Populations`, `GPS coord.`, `Lat` and `Lon` (decimal degrees, converted from the coordinates), and a note ("Weber estimated by BB").
- **Sheet "Polygons_GPS"**: `Populations`, `GPS coord.`, `Lat`, `Lon`, and `xmin`, `ymin`, `xmax`, `ymax` (bounding box of each population, decimal degrees), `Notes`.

Known error, left unchanged in the raw workbook: the `Lat`/`Lon` columns of "GPS Coords" and "Polygons_GPS" read the minutes as decimal degrees (e.g. N38° 04.184' → 38.04184 instead of 38.06973). Do not use them; the map script converts the "Sheet1" coordinates itself.

### metadata_ghcottieproject_summer2021.xlsx
The original metadata workbook. Sheet `readmefirst` defines the two plant IDs; each other sheet has two columns, the variable name and its definition, for the data file of the same name (`gh_cuttingskey`, `gh_allinfo`, `gh_stomdens`, `new_xtype`, `gh_poromeasurements_allgroups`, `gh_sa_raw_1027`, `gh_amscore_raw`, `gh_emscore_raw`, `winrhizaodata_final`, `gh_drybiomass_aboveground`). Not read by any script. Where this dictionary differs from the workbook, the dictionary follows the data and code (noted with `BEA:` where unresolved).

---

## edited_data/

### winrhizaodata_final.csv
Per-plant WinRHIZO totals (229 plants). **Hand-edited**: compiled from the per-scan rows of `raw_data/gh_rootscandata_final.csv`; no script builds it. Read by `01_merge_datasets.R`. Not used in the manuscript results.

| Column | Meaning |
|---|---|
| (unnamed first column) | Row number. |
| `tree_id` | New plant ID. |
| `d` | Average root diameter, mm (the metadata workbook says cm, which is wrong). |
| `vol` | Total root volume, cm³. |
| `tips` | Number of root tips (WinRHIZO count). |
| `forks` | Number of root forks (WinRHIZO count). |
| `crossings` | Number of root crossings (WinRHIZO count). |
| `0L-0.5L`, `0.5L-1L`, `1L-1.5L`, `1.5-2L`, `2L-2.5L`, `2.5-3L`, `3L-3.5L`, `3.5L-4L`, `4-4.5L` | Root length in each 0.5-mm diameter class (0–0.5 … 4–4.5 mm), cm. |
| `4.5l+` | **Unreliable, do not use.** Meant to be root length with diameter > 4.5 mm (cm), but in every row it equals `4.5sa+`: the surface-area values were copied into it during the hand edit. `4.5sa+` matches the scan file for all 229 plants; this column matches the scanned length for only 150. No analysis uses it. |
| `0sa-0.5sa`, `0.5sa-1sa`, `1sa-1.5sa`, `1.5-2sa`, `2sa-2.5sa`, `2.5-3sa`, `3sa-3.5sa`, `3.5sa-4sa`, `4-4.5sa` | Root surface area in each 0.5-mm diameter class, cm². |
| `4.5sa+` | Root surface area with diameter > 4.5 mm, cm². |
| `length` | Total root length, cm. |
| `sa` | Total root surface area, cm². |

(The metadata workbook gives the classes in cm and with some typos, e.g. "0.5L-1L: 0 - 1cm"; the classes are 0.5 mm wide, as in the scan file.)

### mega_ds.csv
The merged plant table (574 cuttings with a watering treatment, a known inoculum and a Q-score class), written by `r_code/01_merge_datasets.R` and read by every analysis script through `r_code/helpers/utilities_data_prep.R`. Column names were cleaned with `janitor::clean_names()` (lower case; `(`, `)`, `%`, `-` and `.` replaced; names starting with a digit get an `x` prefix). Columns that are copies of a raw column are listed with that source; see the raw file for the definition.

| Column | Meaning |
|---|---|
| (unnamed first column) | Row number written by `write.csv`. |
| `tree_idold` | Old plant ID (zero-padding removed where the new ID had been used). |
| `old_class` | Mother tree ID (from `gh_cuttingskey.csv`). |
| `per_f`, `per_n`, `new_class`, `old_xtype` | From `new_xtype.csv` (proportion Fremont, proportion narrowleaf, Q-cutoff cross type, field cross type). |
| `og_h_m`, `og_drc_mm`, `og_leaf`, `og_notes`, `alive_preinnoc`, `alive_01_10_2021`, `drc`, `h`, `wet_root`, `wet_stem`, `flush`, `branch_rs`, `poro_id`, `poro_l1`, `poro_l2`, `poro_l3`, `notes` | From `gh_allinfo.csv` (`og_h(m)` → `og_h_m`, `og_drc(mm)` → `og_drc_mm`). `poro_id` is blanked where it was assigned to more than one plant. |
| `tree_id` | New plant ID (from `gh_cuttingskey.csv`). |
| `pop`, `type`, `geno`, `tree`, `treat` | From `gh_allinfo.csv`, filled from `gh_cuttingskey.csv` where missing. |
| `moist` | Watering treatment, recoded: `Droughted`, `Well-Watered`. |
| `soil` | Soil legacy, recoded: `Tam Legacy`, `No Tam Legacy`. |
| `myco` | Inoculum, recoded: `AM inoc`, `EM inoc`, `Both`, `None (Background)`. |
| `tree_identification`, `indiv`, `dryleaf_weight`, `dryshoot_new` | From `gh_drybiomass_aboveground.csv`. |
| `dry_stem` | `dry_shoottotal` from `gh_drybiomass_aboveground.csv` (whole-stem dry mass, g). |
| `root_dryweight` | From `gh_rootdryweight_FINAL.csv`. |
| `l_1` … `l_57` | Leaf areas from `gh_sa_raw_1027.csv`, cm². |
| `x1` | Row number of `winrhizaodata_final.csv` (its unnamed first column). |
| `d`, `vol`, `tips`, `forks`, `crossings`, `length`, `sa` | From `winrhizaodata_final.csv`. |
| `x0l_0_5l`, `x0_5l_1l`, `x1l_1_5l`, `x1_5_2l`, `x2l_2_5l`, `x2_5_3l`, `x3l_3_5l`, `x3_5l_4l`, `x4_4_5l`, `x4_5l` | Root length by diameter class from `winrhizaodata_final.csv` (`0L-0.5L` … `4.5l+`). |
| `x0sa_0_5sa`, `x0_5sa_1sa`, `x1sa_1_5sa`, `x1_5_2sa`, `x2sa_2_5sa`, `x2_5_3sa`, `x3sa_3_5sa`, `x3_5sa_4sa`, `x4_4_5sa`, `x4_5sa` | Root surface area by diameter class from `winrhizaodata_final.csv` (`0sa-0.5sa` … `4.5sa+`). |
| `l1_t_s1d`, `l1_b_s1d`, `l1_t_s2d`, `l1_b_s2d`, `stomdens_notes` | From `gh_stomdens.csv` (`notes` renamed `stomdens_notes`). |
| `date_wet`, `date_dry` | Porometer measurement date after watering (`wet`) and at the end of a watering cycle (`dry`), from `gh_poromeasurements_allgroups.csv`. |
| `l1_t_wet` … `l3_b_dry` (12 columns) | Conductance of leaf 1–3, top (`t`) or bottom (`b`) surface, after watering (`wet`) or at the end of a cycle (`dry`), mol m⁻² s⁻¹. Only one reading per plant and point in the cycle is kept (the first non-missing); the conductance analyses read all days from the raw porometer file. |
| `run`, `sample_id`, `note`, `leg`, `w` | From the "Final Data" sheet of `gh_soildata_rawfromASU.xlsx`. |
| `soil_chem_treat` | `treat` from the soil workbook (lower-case treatment code of the sampled plant). |
| `weight_mg`, `percent_carbon`, `percent_hydrogen`, `percent_nitrogen`, `c_n` | `Weight (mg)`, `%Carbon`, `%Hydrogen`, `%Nitrogen`, `C:N` from the soil workbook. |
| `general_new_class` | Three-level class: `FRE`, `NAR`, or `HYB` for `BCF`, `BCN` and `F1H`. |
| `total_shoot` | `dryleaf_weight + dry_stem`, g (missing if either is missing). |
| `total_biomass` | `root_dryweight + total_shoot`, g. Includes the original cutting; the manuscript's "total new biomass" is computed in the scripts instead. |
| `top_bottom_stomata_s1`, `top_bottom_stomata_s2` | Ratio of top to bottom stomatal density, field 1 and field 2. |
| `top_bottom_stomata_avg` | Mean of the two ratios. |
| `arb`, `ves`, `hyph`, `dse`, `root` | AMF/DSE intersection counts from `gh_amscore_raw.csv`, rescaled to percent of all intersections. `dse` is DSE colonization (%). |
| `amf` | AMF colonization, %: `arb + ves + hyph`. |
| `root_alive`, `root_dead`, `dark_finger`, `brown_bunch`, `dark_thumb`, `white_tip`, `brown_tip`, `dead_myco` | EMF tip counts from `gh_emscore_raw.csv`, rescaled to percent of all tips scored. |
| `fremont_vs_everything` | `FRE` for `FRE` and `BCF`, `NAR_HYB` otherwise. |
| `por_bottoms` | Mean of all bottom-surface conductance readings in this table, mol m⁻² s⁻¹. |
| `por_tops` | Mean of all top-surface readings. |
| `por_drys` | Mean of all end-of-cycle (`dry`) readings. |
| `por_wets` | Mean of all after-watering (`wet`) readings. |
| `por_wet_dry_ratio` | `por_wets / por_drys`. |
| `por_tops_bottoms_ratio` | `por_tops / por_bottoms`. |
| `por_avg` | Mean of `por_bottoms`, `por_tops`, `por_drys` and `por_wets`. |
| `stom_dens_bottoms` | Mean of `l1_b_s1d` and `l1_b_s2d`. |
| `stom_dens_tops` | Mean of `l1_t_s1d` and `l1_t_s2d`. |
| `stom_dens_avg` | Mean of `stom_dens_bottoms` and `stom_dens_tops`. |
| `stom_dens_top_to_bottom` | `stom_dens_tops / stom_dens_bottoms`. |
| `emf` | EMF colonization, %: `brown_tip + dark_finger + dark_thumb + white_tip + brown_bunch + dead_myco`. |
| `total_leaves_sa` | Sum of `l_1` … `l_57` (total leaf area at harvest, cm²; 0 when no leaf was measured). |
| `avg_leaf_sa` | Mean leaf area, cm². |
| `total_leaves` | Number of leaves with an area. |
| `drc_growth` | `drc − og_drc_mm`, mm. |
| `height_growth` | `h − og_h_m` (same unit as the heights). |
| `number_stomata_tops`, `number_stomata_bottoms` | Stomatal density × `total_leaves_sa` for the top and bottom surfaces (unit depends on the stomatal density unit; see `l1_t_s1d`). |
| `number_stomata_total` | Sum of the two. |
| `root_biomass_allocation` | `root_dryweight / total_biomass`. |
| `fine_root_percentage` | Root length < 0.5 mm diameter / total length (a proportion, despite the name). |
| `spec_root_length` | `length / root_dryweight`, cm g⁻¹. |
| `sla` | `total_leaves_sa / total_biomass`, cm² g⁻¹. Not specific leaf area in the usual sense (divides by total plant mass, not leaf mass). |

Columns that are entirely empty after filtering are dropped by the merge script, so they do not appear here.
