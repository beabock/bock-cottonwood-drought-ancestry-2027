# ==============================================================================
# Small counts: the Holm family of inoculum contrasts (9a) and mother-tree counts (9b); 9c is in 20_emf_drought_influence.R.
#   9a  the Holm family of inoculum-vs-control contrasts (results/17_inoculation_contrasts/inoculation_contrasts.csv).
#   9b  mother trees: genotypes collected (Appendix S1 collection points, from the collection workbook Sheet1), genotypes in
#       the cutting key, and genotypes in the analysis data. The repo has no record of how the 43 were chosen; this only
#       lays out where the others drop out.
# Run from the project root:  Rscript r_code/27_small_counts.R
# Writes results/27_small_counts/: output.txt, strings.csv, mother_trees.csv
# ==============================================================================
source("r_code/helpers/model_helpers.R")
OUT <- "results/27_small_counts"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "output.txt"), split = TRUE); run_stamp("Small counts")

hdr("9a. Holm family for the inoculum contrasts")
s6 <- read_csv("results/17_inoculation_contrasts/inoculation_contrasts.csv", show_col_types = FALSE)
fam <- s6 %>% group_by(response, n) %>% summarise(contrasts = n(), inocula = n_distinct(contrast), watering = n_distinct(treatment), .groups = "drop")
print(as.data.frame(fam), row.names = FALSE); stopifnot(nrow(s6) == 18, all(fam$contrasts == 6))
rec(9, "a", "Holm family, inoculum vs. control", "results/17_inoculation_contrasts/inoculation_contrasts.csv", 18, "Contrasts",
    sprintf("18 = 3 responses × 3 inocula × 2 watering levels: %s; smallest Holm P = %.3f", paste(sprintf("%s (n = %d)", fam$response, fam$n), collapse = "; "), min(s6$P_holm_all18)),
    min(s6$P_holm_all18), "results/17_inoculation_contrasts/inoculation_contrasts.csv")

hdr("9b. Mother trees collected, in the cutting key, and analysed")
cp <- read_csv("results/appendices/appendix_S1_collection_points.csv", show_col_types = FALSE) %>% transmute(pop = Code, type = substr(`Cross type collected`, 1, 1), collected = `Genotypes collected`)
key <- read_csv("data/raw_data/gh_cuttingskey.csv", show_col_types = FALSE) %>% mutate(assigned = !is.na(moist) & moist != "")
xt <- read_csv("data/raw_data/new_xtype.csv", show_col_types = FALSE)
mt <- key %>% group_by(pop, type, geno) %>% summarise(cuttings_in_key = n(), with_treatment = sum(assigned), .groups = "drop") %>%
  mutate(old_class = paste0(pop, type, geno), has_q_score = old_class %in% xt$old_class,
         in_analysis = paste(pop, type, geno, sep = "_") %in% ds$geno_id)
out_t <- cp %>% left_join(mt %>% group_by(pop, type) %>% summarise(in_key = n(), analysed = sum(in_analysis), not_analysed = paste(old_class[!in_analysis], collapse = ", "), .groups = "drop"), by = c("pop", "type"))
print(as.data.frame(out_t), row.names = FALSE)
readr::write_excel_csv(mt, file.path(OUT, "mother_trees.csv"))
cat("\nMother trees in the cutting key but not analysed:\n"); print(as.data.frame(mt %>% filter(!in_analysis)), row.names = FALSE)
cat("Cuttings per analysed mother tree in the key:", paste(range(mt$cuttings_in_key[mt$in_analysis]), collapse = "–"), "\n")
rec(9, "b", "Mother trees", "collection workbook Sheet1; data/raw_data/gh_cuttingskey.csv; mega_ds.csv", 43, "Collected / in cutting key / analysed",
    sprintf("%d genotypes collected at 18 collection points; %d mother trees (%d cuttings) in the cutting key; %d analysed. Not analysed: %s. No record of how the cutting-key trees were chosen from those collected",
            sum(cp$collected), nrow(mt), sum(mt$cuttings_in_key), sum(mt$in_analysis),
            paste(sprintf("%s (%d cuttings, %d with a watering treatment, Q score %s)", mt$old_class[!mt$in_analysis], mt$cuttings_in_key[!mt$in_analysis], mt$with_treatment[!mt$in_analysis],
                          ifelse(mt$has_q_score[!mt$in_analysis], "present", "absent")), collapse = "; ")), NA, "results/27_small_counts/mother_trees.csv")
write_rec("strings.csv"); cat("\nDone.\n"); sink()
