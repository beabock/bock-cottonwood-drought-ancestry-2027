#10/30/24
#Merging all of EAS's datasets.

library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(janitor)

theme_set(theme_bw())
# Run from the project root (the folder containing run_all.R);
# all paths below are relative to it.

gh_all <- read_csv("data/raw_data/gh_allinfo.csv")%>%
  rename("tree_idold" = tree_id)
#why are some treatments missing? Did those tree_ids die?
#ID is "tree_id"
gh_all %>%
  group_by(tree_idold) %>%
  filter(n() > 1)
#tree_idold number BRF413 has two entries. Second one has less info in it. They do have contradicting metadata (treatment and myco). Remove both for now, ask EAS.

gh_all <- gh_all %>%
  filter(tree_idold != "BRF413")

gh_all %>%
  group_by(treat)%>%
  summarize(n=n()) #267 NAs in treat here

# View(gh_all %>%
#        filter(is.na(treat)))

cuttingskey <- read_csv("data/raw_data/gh_cuttingskey.csv")%>%
  distinct()
#Metadata (treat, moist, soil, myco) missing in a bunch of rows here too. Why? 100 are duplicated, too. Need to ask eas about this. removing duplicates for now....
#id is "tree_idold". There is also a tree id column though

duplicates <- cuttingskey[duplicated(cuttingskey) | duplicated(cuttingskey, fromLast = TRUE), ]
#100 exact duplicates

mega <- full_join(gh_all, cuttingskey, by = "tree_idold")

mega <- mega %>%
  mutate(
    pop = coalesce(pop.x, pop.y),
    type = coalesce(type.x, type.y),
    geno = coalesce(geno.x, geno.y),
    tree = coalesce(tree.x, tree.y),
    treat = coalesce(treat.x, treat.y),
    moist = coalesce(moist.x, moist.y),
    soil = coalesce(soil.x, soil.y),
    myco = coalesce(myco.x, myco.y)
    # Add more common columns as needed
  ) %>%
  # Optionally, remove the original columns
  select(-ends_with(".x"), -ends_with(".y"))

mega %>%
  group_by(tree_idold) %>%
  filter(n() > 1)

mega %>%
  group_by(treat) %>%
  summarize(n=n())

# View(mega %>%
#   filter(is.na(treat)))

above_mass <- read_csv("data/raw_data/gh_drybiomass_aboveground.csv")
#What is dry_shoot?
#id is "tree_idold".

above_mass %>%
  group_by(tree_idold) %>%
  filter(n() > 1)

mega <- full_join(mega, above_mass)

mega %>%
  group_by(tree_idold) %>%
  filter(n() > 1)

test <- mega %>%
  filter(tree_idold == "BRF404")

mega <- mega %>%
  group_by(tree_idold) %>%
  summarize(across(everything(), ~ first(na.omit(.))), .groups = "drop")
#Seemed to work. nice.

poro <- read_csv("data/raw_data/gh_poromeasurements_allgroups.csv")

poro <- poro %>% pivot_wider(names_from = water, values_from = colnames(poro)[3:9])

below_mass <- read_csv("data/raw_data/gh_rootdryweight_FINAL.csv")#Some masses missing here.
#id is "tree_idold".
mega <- full_join(mega, below_mass)

mega %>%
  group_by(tree_idold) %>%
  filter(n() > 1)

#scanning <- read_csv("data/raw_data/gh_rootscandata_final.csv")
#ID is "tree_id". doesn't jabe treeidold. Could merge with cuttingskey to get around messing with it too much.
#Scanning has multiple rows per plant. Winrhizo has somehow compiled those rows. Probably for plants that are too big for just one scan. Use winRhizo, not scanning.

sa <- read_csv("data/raw_data/gh_sa_raw_1027.csv")
#Looks similar to another dataset. Come back to this.
#Has new and old tree id
#is each column the surface area of each leaf?
#Some plants have more than one sa row (e.g. JCN1514)

print(sa %>%
  group_by(treat)%>%
  summarize(n=n()), n = 22) #49 NAs in this one

dups <- sa %>%
  group_by(tree_idold) %>%
  filter(n() > 1)
#6 have more than one row in sa. Why?
#For now, removing any that have >1 row. Can come back and do something else.
sa <- sa %>%
  filter(!tree_idold %in% dups$tree_idold)

mega <- full_join(mega, sa) #Removing duplicates in sa removed the issue of many-many in the join.
mega %>%
  group_by(tree_idold) %>%
  filter(n() > 1)

mega <- mega %>%
  group_by(tree_idold) %>%
  summarize(across(everything(), ~ first(na.omit(.))), .groups = "drop")

test <- mega %>%
  filter(tree_idold == "BRF213")

poro <- poro %>% mutate(poro_id = as.character(poro_id))

soil <- read_xlsx("data/raw_data/gh_soildata_rawfromASU.xlsx", skip = 2, sheet = 2)%>% #What is sheet 1 and sheet 2? Is sheet 2 the tranasformed data from sheet 1
  mutate(poro_id = as.character(poro_id))%>%
  filter(!is.na(poro_id))
#Has poro id.
soil <- full_join(poro, soil, by = "poro_id")

soil %>%
  group_by(poro_id) %>%
  filter(n() > 1)

# Fix (2026-07-31): gh_stomdens.csv's own tree_idold only
# covers 33 of the ~228 measured poro_ids -- it was never the right link.
# gh_allinfo.csv already gives mega its own poro_id column (joined further
# below), which covers ~225 of them. Drop gh_stomdens's tree_idold so it
# doesn't collide with mega's, and join everything below on poro_id alone
# instead of the old implicit full_join() (which silently required
# tree_idold AND poro_id AND notes AND treat to all match at once --
# this is what was dropping ~85% of the porometer data). gh_stomdens's
# "notes" (trichome/stomata observations from the microscopy sample) and
# the ASU soil sheet's "treat" are metadata in a different context than
# mega's own "treat"/"notes" columns, so rename rather than drop them --
# keeps the info without splitting mega's real treat/notes into .x/.y.
stom <- read_csv("data/raw_data/gh_stomdens.csv")%>%
  mutate(poro_id = as.character(poro_id))%>%
  filter(!is.na(poro_id)) %>%
  select(-tree_idold)
stom <- full_join(stom, soil, by = "poro_id")

stom %>%
  group_by(poro_id) %>%
  filter(n() > 1)

stom <- stom %>% rename(stomdens_notes = notes, soil_chem_treat = treat)

new_xtype <- read_csv("data/raw_data/new_xtype.csv")
#Good to know.
#Would link to other datasets by the old class data, which is the source tree. not tree id.

mega <- full_join(new_xtype, mega)

mega <- mega %>%
  group_by(tree_idold) %>%
  summarize(across(everything(), ~ first(na.omit(.))), .groups = "drop")

mega %>%
  group_by(tree_idold) %>%
  filter(n() > 1)

# meta <- read_xlsx("data/raw_data/metadata_ghcottieproject_summer2021.xlsx") #Need to go through the sheets on this one manually

winrhizo <- read_csv("data/edited_data/winrhizaodata_final.csv")
#Has tree id (new)

mega <- full_join(mega, winrhizo)

#Step 1: Combine all datasets. Need to add in fung scoring data when I get that from EAS

mega %>%
  group_by(tree_idold) %>%
  filter(n() > 1)

# Fix (2026-07-31): join the porometer/stomatal-density data
# onto mega via poro_id (which mega already carries from gh_allinfo.csv)
# instead of the old implicit full_join(), which required tree_idold AND
# poro_id AND notes AND treat to all match simultaneously and dropped
# ~85% of the physiological data as a result.
#
# A handful of poro_id values are ambiguous -- assigned to more than one
# tree_idold in the raw data. Some of these come from gh_allinfo.csv
# itself; a few more were introduced above when gh_sa_raw_1027.csv's own,
# independent poro_id column got blended in via the sa join. We can't tell
# which tree an ambiguous poro_id's measurements actually belong to, so we
# blank it out for those trees rather than risk attaching one tree's
# porometer/stomatal-density reading to a different tree.
mega <- mega %>% mutate(poro_id = as.character(poro_id))

ambiguous_poro_ids <- mega %>%
  filter(!is.na(poro_id)) %>%
  group_by(poro_id) %>%
  filter(n_distinct(tree_idold) > 1) %>%
  pull(poro_id) %>%
  unique()
cat("n ambiguous poro_id values (assigned to >1 tree_idold), zeroed out:",
    length(ambiguous_poro_ids), "\n")

mega <- mega %>%
  mutate(poro_id = if_else(poro_id %in% ambiguous_poro_ids, NA_character_, poro_id))

mega <- full_join(mega, stom, by = "poro_id")

mega <- mega %>%
  group_by(tree_idold) %>%
  summarize(across(everything(), ~ first(na.omit(.))), .groups = "drop")

#Check if old_id matches any new_ids. I think it does for a couple...

mega <- mega %>%
  mutate(tree_idold_modified = gsub("(^[A-Za-z]{3})0", "\\1", tree_idold))

mega <- mega %>%
  mutate(tree_idold = ifelse(tree_idold != tree_idold_modified, tree_idold_modified, tree_idold)) %>%
  select(-tree_idold_modified)

mega <- mega %>%
  group_by(tree_idold) %>%
  summarize(across(everything(), ~ first(na.omit(.))), .groups = "drop")

mega %>%
  group_by(treat)%>%
  summarize(n=n())


#Why is treat NA for 257 observations? Check earlier in code.

mega <- mega %>%
  select(where(~ !all(is.na(.))))%>%
  filter(treat != "???",
         treat != "WLB?",
         !is.na(tree), #2 NAs here
         !is.na(myco),
         !is.na(new_class))%>% #14 NAs for myco.
  # new_class is kept as in new_xtype.csv, which follows the Scull et al. (2025) Q cutoffs
  # (Fremont Q < 0.1, backcross Fremont 0.1-0.4, ...). An earlier override here,
  # new_class = if_else(per_f > 0.8, "FRE", new_class), relabelled BRF2 (Q = 0.106) and
  # BRH1 (Q = 0.143) as Fremont; removed 2026-10-29.
  mutate(treat = ifelse(treat == "wla", "WLA", treat),
         general_new_class = case_when(
           new_class %in% c("BCF", "BCN", "F1H") ~ "HYB",
           TRUE ~ new_class  # Keep the original value for other cases
         ))

mega %>%
  group_by(new_class)%>%
  summarize(n = n())
#no NAs here

mega %>%
  group_by(general_new_class)%>%
  summarize(n = n())


#Removed NAs. May want to be careful with that, lost 257 obs.

mega %>%
  group_by(soil)%>%
  summarize(n=n())

mega <- mega %>%
   mutate(across(c(h, wet_stem, flush, branch_rs, poro_l1, poro_l2, poro_l3, og_leaf ), as.numeric))%>%
  rename("dry_stem" = "dry_shoottotal")%>%
  mutate(
    total_shoot = if_else(is.na(dryleaf_weight) | is.na(dry_stem), NA_real_, dryleaf_weight + dry_stem),
    total_biomass = if_else(is.na(root_dryweight) | is.na(total_shoot), NA_real_, root_dryweight + total_shoot),
    top_bottom_stomata_s1 = l1_t_s1d/l1_b_s1d,
    top_bottom_stomata_s2 = l1_t_s2d/l1_b_s2d,
    top_bottom_stomata_avg = rowMeans(cbind(top_bottom_stomata_s1, top_bottom_stomata_s2), na.rm = TRUE)
  )
  
em <- read_csv("data/raw_data/gh_emscore_raw.csv")

em %>%
  mutate(row_sum = rowSums(select(., -1))) %>%
  filter(row_sum != 100)
#Not all EM counts are at 100, so we'll need to either rarefy or treat them as ratios.

em_adj <- em %>%
  mutate(row_sum = rowSums(select(., -1))) %>%  # Calculate row sums, excluding the first column
  mutate(across(-1, ~ . / row_sum * 100)) %>%   # Adjust all other columns based on the row sum
  select(-row_sum)


em_adj %>%
  mutate(row_sum = rowSums(select(., -1))) %>%
  filter(row_sum < 99 | row_sum > 101)
#Nice, fixed it.

dups <- em_adj %>%
  group_by(tree_idold) %>%
  summarize(n = n()) %>%
  filter(n > 1) %>%
  pull(tree_idold)  # Get the list of duplicated tree_idold values

# Step 2: Coalesce duplicated rows for each duplicated tree_idold
coalesced_data <- em_adj %>%
  filter(tree_idold %in% dups) %>%  # Filter duplicated rows
  group_by(tree_idold) %>%
  summarize(across(everything(), ~ first(na.omit(.))), .groups = "drop")

# Step 3: Combine coalesced rows with the non-duplicated data
em_adj <- em_adj %>%
  filter(!tree_idold %in% dups) %>%  # Keep non-duplicated rows
  bind_rows(coalesced_data)  

am <- read_csv("data/raw_data/gh_amscore_raw.csv")

am %>%
  mutate(row_sum = rowSums(select(., -1))) %>%
  filter(row_sum != 100)
#this one has same problem.

am_adj <- am %>%
  mutate(row_sum = rowSums(select(., -1))) %>%  # Calculate row sums, excluding the first column
  mutate(across(-1, ~ . / row_sum * 100)) %>%   # Adjust all other columns based on the row sum
  select(-row_sum)

am_adj %>%
  mutate(row_sum = rowSums(select(., -1))) %>%
  filter(row_sum < 99 | row_sum > 101)
#Nice, fixed it.

dups <- am_adj %>%
  group_by(tree_idold) %>%
  summarize(n = n()) %>%
  filter(n > 1) %>%
  pull(tree_idold)  # Get the list of duplicated tree_idold values

# Step 2: Coalesce duplicated rows for each duplicated tree_idold
coalesced_data <- am_adj %>%
  filter(tree_idold %in% dups) %>%  # Filter duplicated rows
  group_by(tree_idold) %>%
  summarize(across(everything(), ~ first(na.omit(.))), .groups = "drop")

# Step 3: Combine coalesced rows with the non-duplicated data
am_adj <- am_adj %>%
  filter(!tree_idold %in% dups) %>%  # Keep non-duplicated rows
  bind_rows(coalesced_data)  %>%
  mutate(amf = arb + ves + hyph)

myco <- full_join(am_adj, em_adj)

myco %>% group_by(tree_idold) %>%
  summarize(n = n()) %>%
  filter(n > 1)
#No duplictaes, good.

mega <- full_join(mega, myco)

mega %>%
  group_by(new_class)%>%
  summarize(n = n())
#12 NAs popped up

test <- mega %>%
  filter(is.na(new_class))%>%
  select_if(~ any(!is.na(.)))
#Seems like these 12 have myco data but no other data. Weird. Not sure how to match these.

mega %>%
  group_by(tree_idold)%>%
  summarize(n=n())%>%
  filter(n>1)


#Sort into fremonts and not fremonts
mega <- mega %>%
  mutate(fremont_vs_everything = case_when(
    new_class %in% c("FRE", "BCF") ~ "FRE",
    new_class %in% c("HYB", "NAR", "F1H", "BCN") ~ "NAR_HYB",
    TRUE ~ general_new_class # keep the original class for any other value
  ))

mega %>%
  group_by(general_new_class)%>%
  summarize(n = n())

mega %>%
  group_by(fremont_vs_everything)%>%
  summarize(n = n())

##porometer avg all drys, all wets, then all everything


mega <- mega %>%
  rowwise() %>%
  mutate(
    por_bottoms = mean(c_across(ends_with("_b_wet") | ends_with("_b_dry") & where(is.numeric)), na.rm = TRUE),  # Bottoms
    por_tops = mean(c_across(ends_with("_t_wet") | ends_with("_t_dry") & where(is.numeric)), na.rm = TRUE),    # Tops
    por_drys = mean(c_across(ends_with("_dry") & where(is.numeric)), na.rm = TRUE),                           # Dry measurements
    por_wets = mean(c_across(ends_with("_wet") & where(is.numeric)), na.rm = TRUE),                           # Wet measurements
    por_wet_dry_ratio = por_wets / por_drys,                                                                  # Wet-to-dry ratio
    por_tops_bottoms_ratio = por_tops / por_bottoms,                                                          # Tops-to-bottoms ratio
    por_avg = mean(c(por_bottoms, por_tops, por_drys, por_wets), na.rm = TRUE)                                # Overall average
  ) %>%
  ungroup()


mega <- mega %>%
  rowwise() %>%
  mutate(
    # FIX (2026-07-31): inside rowwise(), l1_b_s1d/l1_b_s2d are already
    # scalars, so mean(x, y, na.rm=TRUE) was silently passing y as mean()'s
    # `trim` argument instead of averaging the two sub-samples -- it
    # returned x (s1d) only, every time, for both sides and for the
    # tops/bottoms average. Use mean(c(...)) instead to actually average.
    stom_dens_bottoms = mean(c(l1_b_s1d, l1_b_s2d), na.rm = TRUE),  # Bottoms
    stom_dens_tops = mean(c(l1_t_s1d, l1_t_s2d), na.rm = TRUE),
    stom_dens_avg = mean(c(stom_dens_bottoms, stom_dens_tops), na.rm = TRUE),
    stom_dens_top_to_bottom = stom_dens_tops/stom_dens_bottoms) %>%
  ungroup()

mega <- mega %>%
  mutate(
    myco = case_when(
      myco == "A" ~ "AM inoc",
      myco == "E" ~ "EM inoc",
      myco == "B" ~ "Both",
      myco == "N" ~ "None (Background)",
      TRUE ~ myco # Keep any other values unchanged
    ),
    soil = case_when(
      soil == "L" ~ "Tam Legacy",
      soil == "N" ~ "No Tam Legacy",
      TRUE ~ soil
    ),
    moist = case_when(
      moist == "D" ~ "Droughted",
      moist == "W" ~ "Well-Watered",
      TRUE ~ moist
    ),
    amf = ves + hyph + arb,
    emf = brown_tip + dark_finger + dark_thumb + white_tip + brown_bunch + dead_myco,
    total_leaves_sa = rowSums(select(., starts_with("l_")), na.rm = TRUE),
    avg_leaf_sa = rowMeans(select(., starts_with("l_")), na.rm = TRUE),
    total_leaves = rowSums(!is.na(select(., starts_with("l_")))),
    drc_growth = drc - `og_drc(mm)`,
    height_growth = h - `og_h(m)`,
    number_stomata_tops = stom_dens_tops*total_leaves_sa,
    number_stomata_bottoms = stom_dens_bottoms*total_leaves_sa,
    number_stomata_total = number_stomata_tops + number_stomata_bottoms,
    root_biomass_allocation = root_dryweight/total_biomass,
    fine_root_percentage = `0L-0.5L`/length,
    spec_root_length = length/root_dryweight,
    sla = total_leaves_sa/total_biomass
  )%>%
  clean_names()%>%
  filter(!is.na(moist)) #Rows that don't have this arent very useful



write.csv(mega, "data/edited_data/mega_ds.csv")

