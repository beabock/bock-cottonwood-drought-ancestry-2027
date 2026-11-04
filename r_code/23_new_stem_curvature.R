# ==============================================================================
# New stem curvature.
#   5a  Reruns the parametric bootstrap of the reported reduced new stem model (r_code/02a_main_models.R: bootMer,
#       2000 simulations, use.u = FALSE, seed 1998), keeping the linear (b) and quadratic (c) coefficients of every refit.
#       With x = Fremont ancestry - 50, the slope is b + 2cx, so the fitted curve declines across the whole gradient
#       (0-100% Fremont) when b - 100c <= 0 and b + 100c <= 0; it has a peak inside the gradient when c < 0 and 0 <= 50 - b/(2c) <= 100.
#   5b  Inventory of quadratic ancestry terms tested in the reported models, from the files of record (and refits where no file holds them).
# Run from the project root:  Rscript r_code/23_new_stem_curvature.R   (about 10 min)
# Writes results/23_new_stem_curvature/: output.txt, strings.csv, bootstrap_draws.csv, quadratic_inventory.csv
# ==============================================================================
source("r_code/helpers/model_helpers.R")
OUT <- "results/23_new_stem_curvature"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "output.txt"), split = TRUE); run_stamp("New stem curvature")
src <- "results/23_new_stem_curvature/strings.csv"

hdr("5a. Bootstrap draws of the new stem curve")
bio_red_l <- lme4::lmer(log(dryshoot_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
bc <- function(m) { b <- fixef(m); c(vertex = unname(50 - b["per_fc"] / (2 * b["I(per_fc^2)"])), c = unname(b["I(per_fc^2)"]), b = unname(b["per_fc"])) }
bs <- lme4::bootMer(bio_red_l, bc, nsim = 2000, type = "parametric", use.u = FALSE, seed = 1998)
dr <- as_tibble(bs$t) %>% mutate(concave = c < 0, peak_inside = concave & vertex >= 0 & vertex <= 100,
                                 decline_whole = (b - 100 * c) <= 0 & (b + 100 * c) <= 0, increase_whole = (b - 100 * c) >= 0 & (b + 100 * c) >= 0)
readr::write_excel_csv(dr, file.path(OUT, "bootstrap_draws.csv"))
q <- quantile(dr$vertex, c(0.025, 0.975)); cat(sprintf("Reproduction: 95%% percentile CI of the vertex %.1f to %.1f (reported –43.1 to 42.5); outside 0–100: %d (reported 152)\n", q[1], q[2], sum(dr$vertex < 0 | dr$vertex > 100)))
stopifnot(abs(q[1] + 43.1) < 0.05, abs(q[2] - 42.5) < 0.05)
cs <- fixef(bio_red_l); cat(sprintf("Fitted: b = %.5f, c = %.7f, vertex %.1f%%; decline over whole gradient for the fit: %s\n", cs["per_fc"], cs["I(per_fc^2)"], bc(bio_red_l)["vertex"],
                                    (cs["per_fc"] - 100 * cs["I(per_fc^2)"]) <= 0 & (cs["per_fc"] + 100 * cs["I(per_fc^2)"]) <= 0))
sh <- c(concave = mean(dr$concave), peak_inside = mean(dr$peak_inside), decline_whole = mean(dr$decline_whole), increase_whole = mean(dr$increase_whole),
        peak_below_0 = mean(dr$concave & dr$vertex < 0), peak_above_100 = mean(dr$concave & dr$vertex > 100), convex = mean(!dr$concave))
print(round(100 * sh, 1))
rec(5, "a", "Parametric bootstrap of the new stem curve (2000 draws)", "reported reduced model", nobs(bio_red_l), "Share of draws",
    sprintf("peak inside 0–100%% Fremont: %.1f%% (%d of 2000); decline across the whole gradient: %.1f%% (%d); concave with peak below 0%%: %.1f%%; not concave: %.1f%%",
            100 * sh["peak_inside"], sum(dr$peak_inside), 100 * sh["decline_whole"], sum(dr$decline_whole), 100 * sh["peak_below_0"], 100 * sh["convex"]), NA, "results/23_new_stem_curvature/bootstrap_draws.csv")
cat("Peak inside the gradient and decline across the gradient are mutually exclusive except at the boundary; together:", sum(dr$peak_inside | dr$decline_whole), "of 2000\n")
mid <- dr$vertex >= 10 & dr$vertex <= 50   # Results: "89% between 10% and 50% Fremont ancestry"
cat(sprintf("Peak between 10%% and 50%% Fremont ancestry: %.1f%% of draws (%d of 2000)\n", 100 * mean(mid), sum(mid)))

hdr("5b. Quadratic ancestry terms tested in the reported models")
inv <- list(); addq <- function(model, sample, term, string, p, source) inv[[length(inv) + 1]] <<- tibble(model, sample, term, string, P = p, source)
surv_full <- glmer(alive_01_10_2021 ~ pf50 * treatment + I(pf50^2) * treatment + pop + og_drc_mm + soil + myco_cat + (1 | geno_id), data = ds, family = binomial, control = ctl)
surv_red <- glmer(alive_01_10_2021 ~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat + (1 | geno_id), data = ds, family = binomial, control = ctl)
r <- row_z(surv_red, "I(pf50^2)", 2500); addq("Survival (binomial GLMM)", "563 cuttings", "Quadratic ancestry", r$s, r$p, "refit here (r_code/09_allocation_survival_followup.R)")
r <- row_z(surv_full, "treatmentDrought:I(pf50^2)", 2500); addq("Survival (binomial GLMM)", "563 cuttings", "Quadratic ancestry × drought", r$s, r$p, "refit here (02a_main_models.R full model)")
bf <- lmerTest::lmer(log(dryshoot_new) ~ per_fc * treatment + I(per_fc^2) * treatment + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
br <- lmerTest::lmer(log(dryshoot_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
r <- row_t(br, "I(per_fc^2)"); addq("New stem biomass (log)", "425 plants", "Quadratic ancestry", r$s, r$p, "refit here (02a_main_models.R bio_red4)")
r <- row_t(bf, "treatmentDrought:I(per_fc^2)"); addq("New stem biomass (log)", "425 plants", "Quadratic ancestry × drought", r$s, r$p, "refit here (02a_main_models.R bio_full4)")
s6 <- read_csv("results/07_allocation_appendix_S7/appendix_S7_source.csv", show_col_types = FALSE) %>% filter(term == "Quadratic ancestry")
for (i in seq_len(nrow(s6))) { addq(s6$response[i], paste(s6$n_estimated_root[i], "(estimated root mass)"), "Quadratic ancestry", s6$string_estimated_root[i], s6$p_estimated_root[i], "results/07_allocation_appendix_S7/appendix_S7_source.csv")
  addq(s6$response[i], paste(s6$n_measured_root_only[i], "(measured root mass only)"), "Quadratic ancestry", s6$string_measured_root_only[i], s6$p_measured_root_only[i], "results/07_allocation_appendix_S7/appendix_S7_source.csv") }
lr <- read_csv("results/22_allocation_log_ratios/strings.csv", show_col_types = FALSE) %>% filter(item == "a", term == "Quadratic ancestry", !grepl("leaf:root", analysis))
for (i in seq_len(nrow(lr))) addq(lr$analysis[i], paste(lr$n[i], lr$sample[i]), "Quadratic ancestry", lr$string[i], lr$p[i], "results/22_allocation_log_ratios/strings.csv (Appendix S7 log ratios)")
s10 <- read_csv("results/19_conductance_models_of_record/conductance_models.csv", show_col_types = FALSE) %>% filter(term == "Quadratic ancestry")
for (i in seq_len(nrow(s10))) addq(paste("Whole-leaf conductance:", s10$model_label[i]), paste(s10$n_plants[i], "plants"), "Quadratic ancestry", s10$string[i], s10$P[i], "results/19_conductance_models_of_record/conductance_models.csv")
ss <- read_csv("results/04_conductance_surfaces/surface_sensitivity.csv", show_col_types = FALSE) %>%
  filter(analysis == "All readings (main)", term == "I(per_fc^2)") %>% mutate(Model = str_to_sentence(surface))
for (i in seq_len(nrow(ss))) addq(paste(ss$Model[i], "conductance (plant means)"), "203 / 200 plants", "Quadratic ancestry",
  tstr(ss$beta[i], ss$t[i], ss$df[i], ss$p[i]), ss$p[i], "results/04_conductance_surfaces/surface_sensitivity.csv")
zsd <- sd(ds$emf, na.rm = TRUE)
fb <- lmerTest::lmer(log(dryshoot_new) ~ scale(amf) * treatment + scale(emf) * treatment + scale(dse) * treatment + per_fc + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
r <- row_t(fb, "I(per_fc^2)"); addq("New stem biomass, colonization model (covariate)", "77 plants", "Quadratic ancestry", r$s, r$p, "refit here (r_code/08_sample_flow_extra_tests.R)")
inv <- bind_rows(inv) %>% mutate(P_lt_0.05 = P < 0.05)
print(as.data.frame(inv %>% select(model, sample, term, string)), right = FALSE, row.names = FALSE)
readr::write_excel_csv(inv, file.path(OUT, "quadratic_inventory.csv"))
prim <- inv %>% filter(!grepl("measured root mass only", sample))
rec(5, "b", "Quadratic ancestry terms tested", "reported models", nrow(inv), "Count and terms with P < 0.05",
    sprintf("%d quadratic terms in all (%d excluding the measured-root-only repeats); P < 0.05: %s", nrow(inv), nrow(prim),
            paste(sprintf("%s [%s] (%s)", inv$model[inv$P_lt_0.05], inv$sample[inv$P_lt_0.05], inv$string[inv$P_lt_0.05]), collapse = "; ")), NA, "results/23_new_stem_curvature/quadratic_inventory.csv")
write_rec("strings.csv"); cat("\nDone.\n"); sink()
