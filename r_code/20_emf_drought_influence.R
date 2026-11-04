# ==============================================================================
# The EMF x drought association with new stem biomass (Q4).
# Reported model (r_code/08_sample_flow_extra_tests.R; r_code/09_allocation_survival_followup.R; Figure 4):
#   log(dryshoot_new) ~ amf_z * treatment + emf_z * treatment + dse_z * treatment + per_fc + I(per_fc^2)
#                       + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id),   n = 77
#   z = (x - mean) / SD over all plants with that guild scored (EMF SD 13.84 over 166 plants), as scale() in the reported fit.
# Simple slopes per 10 percentage points of EMF as 08_sample_flow_extra_tests.R (Satterthwaite df).
# Influence: Cook's distance from leave-one-plant-out refits, D_i = (b - b_(i))' V^-1 (b - b_(i)) / p with V the
# full-model covariance of the fixed effects (the influence.ME definition); leave-one-plant-out and leave-one-genotype-out
# refits of the interaction and both simple slopes. The reported model is not changed.
# Run from the project root:  Rscript r_code/20_emf_drought_influence.R
# Writes results/20_emf_drought_influence/: output.txt, strings.csv, influence_plants.csv, influence_genotypes.csv, other_responses.csv
# ==============================================================================
source("r_code/helpers/model_helpers.R")
OUT <- "results/20_emf_drought_influence"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "output.txt"), split = TRUE); run_stamp("EMF x drought influence")

zsd <- c(amf = sd(ds$amf, na.rm = TRUE), emf = sd(ds$emf, na.rm = TRUE), dse = sd(ds$dse, na.rm = TRUE))
zmu <- c(amf = mean(ds$amf, na.rm = TRUE), emf = mean(ds$emf, na.rm = TRUE), dse = mean(ds$dse, na.rm = TRUE))
zs <- function(d) d %>% mutate(amf_z = (amf - zmu["amf"]) / zsd["amf"], emf_z = (emf - zmu["emf"]) / zsd["emf"], dse_z = (dse - zmu["dse"]) / zsd["dse"])
f_rep <- log(dryshoot_new) ~ amf_z * treatment + emf_z * treatment + dse_z * treatment + per_fc + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id)
zz <- ds %>% filter(!is.na(amf), !is.na(emf), !is.na(dse), !is.na(dryshoot_new), !is.na(og_drc_mm)) %>% zs() %>%
  left_join(al %>% select(tree_idold, root_est, total_new), by = "tree_idold") %>% mutate(stem_root = dryshoot_new + root_est)
m <- lmerTest::lmer(f_rep, data = zz); stopifnot(nobs(m) == 77)
ti <- function(m) grep("^(emf_z:treatmentDrought|treatmentDrought:emf_z)$", names(fixef(m)), value = TRUE)
TI <- ti(m); stopifnot(length(TI) == 1)
K <- 10 / zsd[["emf"]]   # per SD -> per 10 percentage points
slopes <- function(m) list(ww = lin(m, setNames(K, "emf_z")), dr = lin(m, setNames(c(K, K), c("emf_z", ti(m)))))
pct_ci <- function(r) sprintf("%s%% per 10 percentage points of EMF (95%% CI %s to %s%%)", pc(100 * (exp(r$est) - 1)), pc(100 * (exp(r$est - qt(0.975, r$df) * r$se) - 1)),
                              pc(100 * (exp(r$est + qt(0.975, r$df) * r$se) - 1)))
src <- "results/20_emf_drought_influence/strings.csv"

hdr("1a. The reported 77-plant model: plants and genotypes by watering")
r0 <- row_t(m, TI); s0 <- slopes(m)
cnt <- zz %>% group_by(treatment) %>% summarise(plants = n(), genotypes = n_distinct(geno_id), .groups = "drop")
print(as.data.frame(cnt), row.names = FALSE)
cat("Genotypes in total:", n_distinct(zz$geno_id), "| in both watering treatments:", sum(table(unique(zz[, c("geno_id", "treatment")])$geno_id) == 2), "\n")
rec(1, "a", "EMF × drought, new stem (reported model)", "colonization plants with new stem mass", 77, "Plants (genotypes) by watering",
    sprintf("well-watered %d (%d genotypes); drought %d (%d genotypes); %d genotypes in all", cnt$plants[1], cnt$genotypes[1], cnt$plants[2], cnt$genotypes[2], n_distinct(zz$geno_id)), NA, src)
rec(1, "a", "EMF × drought, new stem (reported model)", "reported", 77, "EMF × drought (per SD of EMF)", r0$s, r0$p, src)
rec(1, "a", "EMF × drought, new stem (reported model)", "reported", 77, "EMF slope, well-watered", paste0(pct_ci(s0$ww), ", ", tstr(s0$ww$est, s0$ww$t, s0$ww$df, s0$ww$p), " (β per 10 points, log scale)"), s0$ww$p, src)
rec(1, "a", "EMF × drought, new stem (reported model)", "reported", 77, "EMF slope, drought", paste0(pct_ci(s0$dr), ", ", tstr(s0$dr$est, s0$dr$t, s0$dr$df, s0$dr$p), " (β per 10 points, log scale)"), s0$dr$p, src)

hdr("1b. Influence: Cook's distance and leave-one-out refits")
b0 <- fixef(m); Vinv <- solve(as.matrix(vcov(m))); p <- length(b0)
refit <- function(d) { mm <- suppressMessages(lmerTest::lmer(f_rep, data = d)); s <- slopes(mm); r <- row_t(mm, TI)
  tibble(int_b = r$b, int_t = r$t, int_df = r$df, int_p = r$p, ww_b = s$ww$est, ww_p = s$ww$p, dr_b = s$dr$est, dr_p = s$dr$p,
         cook = as.numeric(t(b0 - fixef(mm)) %*% Vinv %*% (b0 - fixef(mm))) / p, singular = isSingular(mm)) }
loo <- map_dfr(seq_len(nrow(zz)), function(i) bind_cols(zz[i, ] %>% select(tree_idold, geno_id, treatment, per_f, emf, dryshoot_new), refit(zz[-i, ]))) %>%
  arrange(desc(cook))
log <- map_dfr(unique(zz$geno_id), function(g) bind_cols(tibble(geno_id = g, n_removed = sum(zz$geno_id == g)), refit(filter(zz, geno_id != g)))) %>% arrange(int_p)
readr::write_excel_csv(loo, file.path(OUT, "influence_plants.csv")); readr::write_excel_csv(log, file.path(OUT, "influence_genotypes.csv"))
cat("Cook's distance threshold 4/n =", round(4 / nrow(zz), 3), "| plants above it:", sum(loo$cook > 4 / nrow(zz)), "\nTop 8 plants by Cook's distance:\n")
print(as.data.frame(head(loo, 8) %>% mutate(across(where(is.double), ~ signif(.x, 4)))), row.names = FALSE)
rng <- function(d, lab) {
  cat(sprintf("\n%s (%d refits):\n", lab, nrow(d)))
  for (k in list(c("int", "EMF × drought (per SD)"), c("ww", "EMF slope, well-watered (log per 10 points)"), c("dr", "EMF slope, drought (log per 10 points)"))) {
    b <- d[[paste0(k[1], "_b")]]; pv <- d[[paste0(k[1], "_p")]]
    cat(sprintf("  %-45s estimate %s to %s | P %s to %s | P < 0.05 in %d of %d\n", k[2], num(min(b)), num(max(b)), formatC(min(pv), format = "g", digits = 3),
                formatC(max(pv), format = "g", digits = 3), sum(pv < 0.05), nrow(d)))
    rec(1, "b", paste("Influence:", lab), "77-plant model", nrow(zz), k[2],
        sprintf("estimate %s to %s; P %s to %s; P < 0.05 in %d of %d refits (reported: %s)", num(min(b)), num(max(b)), sub("^P [=<] ", "", P(min(pv))), sub("^P [=<] ", "", P(max(pv))),
                sum(pv < 0.05), nrow(d), P(switch(k[1], int = r0$p, ww = s0$ww$p, dr = s0$dr$p))), max(pv), "results/20_emf_drought_influence/influence_plants.csv; influence_genotypes.csv") } }
rng(loo, "Leave one plant out"); rng(log, "Leave one genotype out")
cat("\nGenotype removals that change significance of the interaction or the drought slope:\n")
print(as.data.frame(log %>% filter(int_p >= 0.05 | dr_p >= 0.05) %>% mutate(across(where(is.double), ~ signif(.x, 4)))), row.names = FALSE)
cat("Plant removals that change significance of the interaction or the drought slope:\n")
print(as.data.frame(loo %>% filter(int_p >= 0.05 | dr_p >= 0.05) %>% mutate(across(where(is.double), ~ signif(.x, 4)))), row.names = FALSE)
cat("Refits whose well-watered slope reaches P < 0.05:", sum(loo$ww_p < 0.05), "plants;", sum(log$ww_p < 0.05), "genotypes\n")

# the Figure 4 plant at about 0.98 g new stem and 12% EMF
hi <- zz %>% filter(abs(dryshoot_new - 0.98) < 0.005, abs(emf - 12) < 1)
stopifnot(nrow(hi) == 1)
cat(sprintf("\nFigure 4 high point: %s (genotype %s, %s, %.1f%% Fremont, EMF %.1f%%, AMF %.1f%%, DSE %.1f%%, new stem %.2f g, inoculum %s, soil %s); next-largest new stem among the 77: %.2f g\n",
            hi$tree_idold, hi$geno_id, hi$treatment, hi$per_f, hi$emf, hi$amf, hi$dse, hi$dryshoot_new, hi$myco_cat, hi$soil, sort(zz$dryshoot_new, decreasing = TRUE)[2]))
lh <- loo %>% filter(tree_idold == hi$tree_idold)
cat(sprintf("Its Cook's distance: %.3f (rank %d of 77)\n", lh$cook, which(loo$tree_idold == hi$tree_idold)))
mh <- lmerTest::lmer(f_rep, data = filter(zz, tree_idold != hi$tree_idold)); rh <- row_t(mh, TI); sh <- slopes(mh)
rec(1, "b", sprintf("Without %s (Figure 4 point at %.2f g, %.1f%% EMF, drought)", hi$tree_idold, hi$dryshoot_new, hi$emf), "76 plants", nobs(mh), "EMF × drought (per SD)", rh$s, rh$p, src)
rec(1, "b", sprintf("Without %s", hi$tree_idold), "76 plants", nobs(mh), "EMF slope, well-watered", paste0(pct_ci(sh$ww), ", t = ", dec(sh$ww$t), ", df = ", dff(sh$ww$df), ", ", P(sh$ww$p)), sh$ww$p, src)
rec(1, "b", sprintf("Without %s", hi$tree_idold), "76 plants", nobs(mh), "EMF slope, drought", paste0(pct_ci(sh$dr), ", t = ", dec(sh$dr$t), ", df = ", dff(sh$dr$df), ", ", P(sh$dr$p)), sh$dr$p, src)

hdr("1c. Multiplicity")
gx <- function(g) grep(paste0("^(", g, "_z:treatmentDrought|treatmentDrought:", g, "_z)$"), names(fixef(m)), value = TRUE)
gp <- c(AMF = row_t(m, gx("amf"))$p, EMF = r0$p, DSE = row_t(m, gx("dse"))$p)
gh <- p.adjust(gp, "holm"); print(rbind(unadjusted = gp, holm = gh))
rec(1, "c", "Guild × drought, Holm across the three guilds", "reported", 77, "EMF × drought", sprintf("%s unadjusted; Holm %s (AMF Holm %s; DSE Holm %s)", P(gp["EMF"]), P(gh["EMF"]), P(gh["AMF"]), P(gh["DSE"])), gh["EMF"], src)
sp <- c(ww = s0$ww$p, dr = s0$dr$p); sh2 <- p.adjust(sp, "holm")
rec(1, "c", "EMF simple slopes, Holm within the model (2 slopes)", "reported", 77, "Drought / well-watered slope", sprintf("drought %s (Holm %s); well-watered %s (Holm %s)", P(sp["dr"]), P(sh2["dr"]), P(sp["ww"]), P(sh2["ww"])), sh2["dr"], src)
cat("Bonferroni across the 2 slopes and 3 guilds (6 tests), drought slope:", P(min(1, 6 * sp["dr"])), "\n")

hdr("1d. The same model for total new biomass and for new stem + root mass")
cat("All 77 plants were scored for colonization, so their root mass is estimated (fresh mass x predicted dry:fresh ratio);\n",
    "a measured-root-only version of these models is not possible.\n", sep = "")
oth <- list()
for (y in c("log(total_new)", "log(stem_root)", "log(dryshoot_new)")) {
  mm <- lmerTest::lmer(update(f_rep, as.formula(paste(y, "~ ."))), data = zz); r <- row_t(mm, TI); s <- slopes(mm)
  lab <- c("log(total_new)" = "Total new biomass (root mass estimated)", "log(stem_root)" = "New stem + root mass (root mass estimated)", "log(dryshoot_new)" = "New stem biomass (reported)")[[y]]
  oth[[y]] <- tibble(response = lab, n = nobs(mm), term = c("EMF × drought (per SD)", "EMF slope, well-watered", "EMF slope, drought"),
                     string = c(r$s, paste0(pct_ci(s$ww), ", t = ", dec(s$ww$t), ", df = ", dff(s$ww$df), ", ", P(s$ww$p)), paste0(pct_ci(s$dr), ", t = ", dec(s$dr$t), ", df = ", dff(s$dr$df), ", ", P(s$dr$p))),
                     p = c(r$p, s$ww$p, s$dr$p))
  if (y != "log(dryshoot_new)") for (i in 1:3) rec(1, "d", lab, "colonization plants", nobs(mm), oth[[y]]$term[i], oth[[y]]$string[i], oth[[y]]$p[i], "results/20_emf_drought_influence/other_responses.csv") }
readr::write_excel_csv(bind_rows(oth), file.path(OUT, "other_responses.csv"))

hdr("2c. EMF only (no AMF or DSE terms), every EMF-scored survivor with new stem mass")
f_emf <- log(dryshoot_new) ~ emf_z * treatment + per_fc + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id)
ze <- ds %>% filter(!is.na(emf), !is.na(dryshoot_new), !is.na(og_drc_mm)) %>% zs()
for (k in list(list("all EMF-scored survivors", ze), list("the 77 plants of the reported model", zz))) {
  mm <- lmerTest::lmer(f_emf, data = k[[2]]); r <- row_t(mm, ti(mm)); s <- slopes(mm); d <- k[[2]]
  smp <- sprintf("%s (well-watered %d, drought %d; %d genotypes)", k[[1]], sum(d$treatment == "Well-watered"), sum(d$treatment == "Drought"), n_distinct(d$geno_id))
  rec(2, "c", "EMF-only model, new stem", smp, nobs(mm), "EMF × drought (per SD)", r$s, r$p, src)
  rec(2, "c", "EMF-only model, new stem", smp, nobs(mm), "EMF slope, well-watered", paste0(pct_ci(s$ww), ", t = ", dec(s$ww$t), ", df = ", dff(s$ww$df), ", ", P(s$ww$p)), s$ww$p, src)
  rec(2, "c", "EMF-only model, new stem", smp, nobs(mm), "EMF slope, drought", paste0(pct_ci(s$dr), ", t = ", dec(s$dr$t), ", df = ", dff(s$dr$df), ", ", P(s$dr$p)), s$dr$p, src)
  cat("  singular:", isSingular(mm), "\n") }
cat("EMF-scored survivors without new stem mass or initial diameter (not in the model):", sum(!is.na(ds$emf) & ds$alive_01_10_2021 %in% 1) - nrow(ze), "\n")

hdr("9c. Narrowleaf plants only (linear ancestry only; r_code/09_allocation_survival_followup.R)")
f_nl <- log(dryshoot_new) ~ amf_z * treatment + emf_z * treatment + dse_z * treatment + per_fc + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id)
zn <- zz %>% filter(new_class == "Narrowleaf"); mn <- lmerTest::lmer(f_nl, data = zn); rn <- row_t(mn, TI); sn <- slopes(mn)
smp <- sprintf("narrowleaf plants (well-watered %d, drought %d; %d genotypes)", sum(zn$treatment == "Well-watered"), sum(zn$treatment == "Drought"), n_distinct(zn$geno_id))
rec(9, "c", "EMF × drought, new stem, narrowleaf only", smp, nobs(mn), "EMF × drought (per SD)", rn$s, rn$p, src)
rec(9, "c", "EMF × drought, new stem, narrowleaf only", smp, nobs(mn), "EMF slope, well-watered", paste0(pct_ci(sn$ww), ", t = ", dec(sn$ww$t), ", df = ", dff(sn$ww$df), ", ", P(sn$ww$p)), sn$ww$p, src)
rec(9, "c", "EMF × drought, new stem, narrowleaf only", smp, nobs(mn), "EMF slope, drought", paste0(pct_ci(sn$dr), ", t = ", dec(sn$dr$t), ", df = ", dff(sn$dr$df), ", ", P(sn$dr$p)), sn$dr$p, src)
cat("  singular:", isSingular(mn), "| residual df for 13 fixed effects:", nobs(mn) - length(fixef(mn)), "\n")

write_rec("strings.csv"); cat("\nDone.\n"); sink()
