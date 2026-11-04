# ==============================================================================
# Sensitivity checks and counts reported in the manuscript.
#   a  EMF-scored plants: 166 scored, 165 in the EMF models (one plant lacks initial cutting diameter)
#   b  predicted EMF colonization at 0% and 100% Fremont ancestry and the relative decline (colonization LMM, 78 plants)
#   c  summary of the quadratic-term inventory (results/23_new_stem_curvature/quadratic_inventory.csv)
#   d  the plant with the smallest new stem mass (SFF111, 0.001 g; largest total new biomass)
#   e  the reported growth and allocation models refit without SFF111
# Run from the project root:  Rscript r_code/29_sensitivity_checks.R
# Writes results/29_sensitivity_checks/: output.txt, checks.csv
# ==============================================================================
source("r_code/helpers/model_helpers.R")
OUT <- "results/29_sensitivity_checks"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "output.txt"), split = TRUE); run_stamp("Sensitivity checks")
set.seed(1998)
tl <- levels(ds$treatment)
red <- paste("pf50 + treatment + I(pf50^2) +", cov, "+", re)

hdr("Models of record refit for the checks")
stem_m <- lmerTest::lmer(log(dryshoot_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
stem_pts <- ds %>% filter(if_all(c(dryshoot_new, per_fc, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), ~ !is.na(.x)))
stopifnot(nrow(stem_pts) == nobs(stem_m), nobs(stem_m) == 425)
a414 <- al %>% mutate(root_stem = root_est / new_stem)
tot_m <- lmerTest::lmer(log(total_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = a414)
rs_m <- lmerTest::lmer(log(root_stem) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = a414)
lv <- al %>% filter(leaf > 0, root_est > 0, new_stem > 0) %>% mutate(leaf_stem = leaf / new_stem)
ls_m <- lmerTest::lmer(as.formula(paste("log(leaf_stem) ~", red)), data = lv)
stopifnot(nobs(tot_m) == 414, nobs(rs_m) == 414, nobs(ls_m) == 348)
cat("Linear ancestry (per 1%), as reported:\n  new stem ", row_t(stem_m, "per_fc")$s, "\n  total new ", row_t(tot_m, "per_fc")$s,
    "\n  root:new stem ", row_t(rs_m, "per_fc")$s, "\n  leaf:new stem (plants with leaves) ", row_t(ls_m, "pf50", 50)$s, "\n")
fungal_types <- c("amf", "emf", "dse")
d3 <- ds %>% select(tree_id, pop, per_f, treatment, og_drc_mm, amf, emf, dse, soil, myco_cat) %>% na.omit() %>%
  pivot_longer(c(amf, emf, dse), names_to = "fungal_type", values_to = "colonization") %>%
  mutate(fungal_type = factor(fungal_type, levels = fungal_types), og_drc_mm_z = (og_drc_mm - mean(og_drc_mm)) / sd(og_drc_mm))
col_lmm <- lmer(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil + (1 | pop) + (1 | myco_cat) + (1 | tree_id), data = d3)
stopifnot(nobs(col_lmm) == 234)
fe3 <- delete.response(terms(colonization ~ per_f * fungal_type * treatment + og_drc_mm_z * fungal_type + soil))

hdr("Checks")
# ------------------------------------------------------------------------------
CHK <- list(); chk <- function(item, what, value) { CHK[[length(CHK) + 1]] <<- tibble(item = item, check = what, value = as.character(value)); cat(sprintf("[%s] %s: %s\n", item, what, value)) }
# (a) EMF-scored plants: Figure S2 uses every scored plant; the models need all covariates
emf_all <- ds %>% filter(!is.na(emf))
emf_mod <- emf_all %>% filter(if_all(c(pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat), ~ !is.na(.x)))
drop <- anti_join(emf_all, emf_mod, by = "tree_idold")
miss <- drop %>% select(tree_idold, pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat) %>%
  mutate(missing = pmap_chr(list(pop, geno_id, per_f, treatment, og_drc_mm, soil, myco_cat), function(...) {
    v <- list(...); paste(c("pop", "geno_id", "per_f", "treatment", "og_drc_mm", "soil", "myco_cat")[map_lgl(v, ~ is.na(.x))], collapse = ", ") }))
chk("a", "EMF-scored plants (all three guilds + EMF only)", sprintf("%d (%d + %d)", nrow(emf_all), sum(!is.na(emf_all$amf)), sum(is.na(emf_all$amf))))
chk("a", "EMF-scored plants with every model covariate", nrow(emf_mod))
chk("a", "Plants dropped from the EMF models", paste(sprintf("%s (missing: %s; %s)", miss$tree_idold, miss$missing,
    ifelse(is.na(drop$amf), "EMF only", "all three guilds")), collapse = "; "))
# (b) predicted EMF colonization at 0% and 100% Fremont (reported 78-plant colonization LMM, averaged over watering and covariates)
b3 <- fixef(col_lmm); V3 <- as.matrix(vcov(col_lmm))
xe <- function(pf) rowMeans(sapply(tl, function(trt) colMeans(model.matrix(fe3, d3 %>% mutate(per_f = pf, fungal_type = factor("emf", levels = fungal_types),
                                                                                                treatment = factor(trt, levels = tl)))[, names(b3)])))
x0 <- xe(0); x100 <- xe(100); e0 <- sum(x0 * b3); e100 <- sum(x100 * b3)
se_d <- sqrt(as.numeric(t(x100 - x0) %*% V3 %*% (x100 - x0)))
sims <- MASS::mvrnorm(10000, b3, V3); rel <- 1 - (sims %*% x100) / (sims %*% x0)
chk("b", "Predicted EMF colonization, 0% Fremont", sprintf("%.1f%% (95%% CI %.1f-%.1f)", e0, e0 - 1.96 * sqrt(t(x0) %*% V3 %*% x0), e0 + 1.96 * sqrt(t(x0) %*% V3 %*% x0)))
chk("b", "Predicted EMF colonization, 100% Fremont", sprintf("%.1f%% (95%% CI %.1f-%.1f)", e100, e100 - 1.96 * sqrt(t(x100) %*% V3 %*% x100), e100 + 1.96 * sqrt(t(x100) %*% V3 %*% x100)))
chk("b", "Difference, 100% minus 0% Fremont (percentage points)", sprintf("%.1f (95%% CI %.1f to %.1f)", e100 - e0, e100 - e0 - 1.96 * se_d, e100 - e0 + 1.96 * se_d))
chk("b", "Relative decline from 0% to 100% Fremont", sprintf("%.0f%% (95%% CI %.0f-%.0f%%; 10000 draws from the fixed-effect sampling distribution)",
    100 * (1 - e100 / e0), 100 * quantile(rel, 0.025), 100 * quantile(rel, 0.975)))
# (c) quadratic-term inventory
inv <- read_csv("results/23_new_stem_curvature/quadratic_inventory.csv", show_col_types = FALSE)
qx <- grepl(TIMES, inv$term, fixed = TRUE)
chk("c", "Quadratic terms in the inventory", sprintf("%d in all: %d quadratic main effects and %d quadratic × drought interactions", nrow(inv), sum(!qx), sum(qx)))
chk("c", "Terms with P < 0.05", sprintf("%d of %d: %s", sum(inv$P_lt_0.05), nrow(inv), paste(sprintf("%s [%s]", inv$model[inv$P_lt_0.05], inv$sample[inv$P_lt_0.05]), collapse = "; ")))
chk("c", "Expected by chance at alpha = 0.05", sprintf("%.1f of %d (terms are not independent; several are the same response in two root-mass samples)", 0.05 * nrow(inv), nrow(inv)))
# (d) the smallest new stem mass
lo <- ds %>% filter(!is.na(dryshoot_new)) %>% arrange(dryshoot_new) %>% slice(1:3)
chk("d", "Smallest new stem masses (g)", paste(sprintf("%s %.3f", lo$tree_idold, lo$dryshoot_new), collapse = "; "))
z <- al %>% filter(new_stem == min(new_stem))
chk("d", "Allocation ratios of the 0.001 g plant", sprintf("%s: root:new stem %.0f, leaf:new stem %.0f, leaf %.3f g, root %.3f g, %s, %.1f%% Fremont",
    z$tree_idold, z$root_est / z$new_stem, z$leaf / z$new_stem, z$leaf, z$root_est, as.character(z$treatment), z$per_f))
chk("d", "Rank of this plant's total new biomass", sprintf("%d of %d (%.3f g; median %.3f g)", rank(-al$total_new)[al$tree_idold == z$tree_idold], nrow(al), z$total_new, median(al$total_new)))
# (e) reported allocation models refit without this plant
nz <- function(d) filter(d, tree_idold != z$tree_idold)
refit <- list(`New stem biomass (log), linear ancestry` = list(stem_m, "per_fc", 1), `Total new biomass (log), linear ancestry` = list(tot_m, "per_fc", 1),
              `Root:new stem (log), linear ancestry` = list(rs_m, "per_fc", 1), `Leaf:new stem (log), plants with leaves, linear ancestry` = list(ls_m, "pf50", 50),
              `New stem biomass (log), quadratic ancestry` = list(stem_m, "I(per_fc^2)", 1))
for (nm in names(refit)) { m <- refit[[nm]][[1]]; tm <- refit[[nm]][[2]]; sc <- refit[[nm]][[3]]
  m2 <- update(m, data = nz(model.frame(m) %>% as_tibble() %>% { if (nrow(.) == nrow(stem_pts)) stem_pts else if (nrow(.) == 414) a414 else lv }))
  chk("e", paste(nm, "| reported vs without", z$tree_idold), sprintf("%s | %s (n = %d)", row_t(m, tm, sc)$s, row_t(m2, tm, sc)$s, nobs(m2))) }
readr::write_excel_csv(bind_rows(CHK), file.path(OUT, "checks.csv"))
cat("\nDone.\n"); sink()
