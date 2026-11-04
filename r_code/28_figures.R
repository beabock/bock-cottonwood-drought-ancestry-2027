# ==============================================================================
# FIGURES: every figure in the manuscript and appendices except the Appendix S1 map (r_code/12_map_appendix_S1.R).
# All figures share r_code/helpers/figure_style.R (theme_bw, Arial, AJB sizes; color = watering treatment only).
#   Fig1   A survival, B leafless at harvest, along Fremont ancestry (binomial GLMMs; points = mother-tree means)
#   Fig2   A new stem biomass, B total new biomass, C stem, leaf and root fractions of total new biomass
#   Fig3   AMF, EMF and DSE colonization along Fremont ancestry (colonization LMM)
#   Fig4   new stem biomass against AMF, EMF and DSE colonization (colonization-biomass LMM)
#   FigS2  colonization by inoculum treatment (Appendix S2)
#   FigS3  colonization by watering treatment (Appendix S3)
#   AppendixS5_daypair, AppendixS5_surfaces_daypair  stomatal conductance by measurement session (Figures S5.1, S5.2)
# Models are those of record (same formulas and data as the scripts that report them; cited at each figure).
# Lines are average-marginal predictions: the design rows of the model data are averaged at each x value.
# Run from the project root (after steps 01-31):  Rscript r_code/28_figures.R
# Writes figures/<name>.pdf and .png; results/28_figures/fig2_points_at_edge.csv; results/28_figures/figure_S5_2_dropped.csv
# ==============================================================================
source("r_code/helpers/model_helpers.R")
source("r_code/helpers/figure_style.R")
library(glmmTMB); library(patchwork)
dir.create("figures", showWarnings = FALSE); dir.create("results/28_figures", showWarnings = FALSE, recursive = TRUE)
set.seed(1998)
ds <- ds %>% mutate(treatment = factor(treatment, levels = TRT_LEVELS))
grid_anc <- seq(0, 100, by = 2)
jit <- position_jitter(width = 1, height = 0, seed = 1)

# average-marginal prediction curve by watering treatment for a fitted lme4/lmerTest/glmmTMB conditional model
curve_trt <- function(beta, V, fe_terms, base, xgrid, xvar, xfun = identity, inv = identity) {
  map_dfr(xgrid, function(g) map_dfr(TRT_LEVELS, function(trt) {
    nd <- base; nd[[xvar]] <- xfun(g); nd$treatment <- factor(trt, levels = TRT_LEVELS)
    xm <- colMeans(model.matrix(fe_terms, data = nd)[, names(beta), drop = FALSE])
    tibble(x = g, treatment = trt, fit = sum(xm * beta), se = sqrt(as.numeric(t(xm) %*% V %*% xm))) })) %>%
    mutate(treatment = factor(treatment, levels = TRT_LEVELS), y = inv(fit), lower = inv(fit - 1.96 * se), upper = inv(fit + 1.96 * se)) }
curve_lmm <- function(m, d, xvar, xfun, inv = exp)
  curve_trt(fixef(m), as.matrix(vcov(m)), delete.response(terms(formula(m, fixed.only = TRUE))),
            d %>% select(og_drc_mm, myco_cat, soil, treatment), grid_anc, xvar, xfun, inv)
band <- function(cd) list(geom_ribbon(data = cd, aes(x, ymin = lower, ymax = upper, fill = treatment), alpha = BAND_ALPHA),
                          geom_line(data = cd, aes(x, y, colour = treatment, linetype = treatment), linewidth = LN_W))
pts <- function(d, x, y, size = PT_SIZE, alpha = PT_ALPHA, position = jit)
  geom_point(data = d, aes(.data[[x]], .data[[y]], fill = treatment, shape = treatment), colour = PT_EDGE, stroke = PT_STROKE,
             size = size, alpha = alpha, position = position)

# ------------------------------------------------------------------------------
# Figure 1: survival and leaflessness (models as r_code/09_allocation_survival_followup.R section 4 and r_code/24_soil_missing_stems_survival_weights.R)
# ------------------------------------------------------------------------------
surv_m <- glmer(alive_01_10_2021 ~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat + (1 | geno_id), data = ds, family = binomial, control = ctl)
sp <- ds %>% filter(if_all(c(alive_01_10_2021, pf50, treatment, pop, og_drc_mm, soil, myco_cat, geno_id), ~ !is.na(.x)))
stopifnot(nrow(sp) == nobs(surv_m), nobs(surv_m) == 563)
c1a <- curve_trt(fixef(surv_m), as.matrix(vcov(surv_m)), delete.response(terms(~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat)),
                 sp %>% select(pop, og_drc_mm, soil, myco_cat, treatment), grid_anc, "pf50", function(g) (g - 50) / 50, plogis)
lz <- ds %>% filter(alive_01_10_2021 %in% 1, !is.na(dryshoot_new)) %>% mutate(zero_leaf = as.integer(coalesce(dryleaf_weight, 0) == 0)) %>%
  filter(if_all(c(pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), ~ !is.na(.x)))
leaf_m <- glmer(as.formula(paste("zero_leaf ~ pf50 + treatment + I(pf50^2) +", cov, "+", re)), data = lz, family = binomial, control = ctl)
stopifnot(nobs(leaf_m) == 425, sum(lz$zero_leaf) == 68)
c1b <- curve_trt(fixef(leaf_m), as.matrix(vcov(leaf_m)), delete.response(terms(~ pf50 + treatment + I(pf50^2) + og_drc_mm + myco_cat + soil)),
                 lz %>% select(og_drc_mm, myco_cat, soil, treatment), grid_anc, "pf50", function(g) (g - 50) / 50, plogis)
tree_means <- function(d, y) d %>% group_by(geno_id, treatment) %>% summarise(per_f = mean(per_f), n = n(), y = mean(.data[[y]]), .groups = "drop") %>%
  mutate(per_f = per_f + ifelse(treatment == "Drought", 0.8, -0.8))     # separate the two treatments at each mother tree
g1a <- tree_means(sp, "alive_01_10_2021"); g1b <- tree_means(lz, "zero_leaf")
n_max <- max(g1a$n, g1b$n)
size_sc <- scale_size_area(max_size = 3.2, breaks = c(2, 5, 10, 15), limits = c(0, n_max), name = "Plants per mother tree")
p1 <- function(cd, g, ylab) ggplot() + band(cd) +
  geom_point(data = g, aes(per_f, y, fill = treatment, shape = treatment, size = n), colour = PT_EDGE, stroke = PT_STROKE, alpha = 0.7) +
  size_sc + trt_scales() + scale_y_continuous(limits = c(0, 1)) + labs(x = X_ANC, y = ylab) + theme_fig()
fig1 <- p1(c1a, g1a, "Survival (proportion of cuttings)") + p1(c1b, g1b, "Leafless at harvest (proportion of survivors)") +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom", legend.box = "vertical")
save_fig(fig1, "Fig1", W_2COL, 9.5 * CM)

# ------------------------------------------------------------------------------
# Figure 2: growth and allocation (models as r_code/09_allocation_survival_followup.R section 1 and r_code/22_allocation_log_ratios.R)
# Log axes have fixed limits; a point beyond a limit is drawn on the panel edge as an open symbol and listed in
# results/28_figures/fig2_points_at_edge.csv (one plant, SFF111, with 0.001 g new stem mass).
# ------------------------------------------------------------------------------
stem_m <- lmerTest::lmer(log(dryshoot_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = ds)
stem_pts <- ds %>% filter(if_all(c(dryshoot_new, per_fc, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), ~ !is.na(.x)))
stopifnot(nrow(stem_pts) == nobs(stem_m), nobs(stem_m) == 425)
a414 <- al %>% mutate(treatment = factor(treatment, levels = TRT_LEVELS))
tot_m <- lmerTest::lmer(log(total_new) ~ per_fc + treatment + I(per_fc^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = a414)
stopifnot(nobs(tot_m) == 414)
EDGE <- list()
edge_pts <- function(d, y, lims, panel) {
  d <- d %>% mutate(.y = .data[[y]], .edge = .y < lims[1] | .y > lims[2], .yp = pmin(pmax(.y, lims[1]), lims[2]))
  if (any(d$.edge)) EDGE[[length(EDGE) + 1]] <<- d %>% filter(.edge) %>%
    transmute(panel = panel, plant = tree_idold, treatment, fremont_ancestry = round(per_f, 1), value = signif(.y, 4), drawn_at = ifelse(.y < lims[1], lims[1], lims[2]))
  list(pts(filter(d, !.edge), "per_f", ".yp"),
       geom_point(data = filter(d, .edge), aes(per_f, .yp, shape = treatment), fill = NA, colour = "black", stroke = 0.6, size = PT_SIZE + 0.4, position = jit, show.legend = FALSE)) }
log_axis <- function(lims, breaks) list(scale_y_log10(breaks = breaks, labels = format(breaks, drop0trailing = TRUE, scientific = FALSE, trim = TRUE)),
                                        coord_cartesian(ylim = lims))
p2a <- ggplot() + edge_pts(stem_pts, "dryshoot_new", c(0.01, 1), "A, new stem biomass (g)") + band(curve_lmm(stem_m, stem_pts, "per_fc", function(g) g - 50)) +
  log_axis(c(0.01, 1), c(0.01, 0.03, 0.1, 0.3, 1)) + trt_scales() + labs(x = X_ANC, y = expression("New stem biomass (g, log"[10] * " scale)")) + theme_fig()
p2b <- ggplot() + edge_pts(a414, "total_new", c(0.1, 10), "B, total new biomass (g)") + band(curve_lmm(tot_m, a414, "per_fc", function(g) g - 50)) +
  log_axis(c(0.1, 10), c(0.1, 0.3, 1, 3, 10)) + trt_scales() + labs(x = X_ANC, y = expression("Total new biomass (g, log"[10] * " scale)")) + theme_fig()
# C: ordered beta GLMMs of the stem, leaf and root mass fractions; curves are the full expected value (beta mean plus the
# probabilities of exact 0 and 1 from the cutpoints psi), with 95% intervals from 1000 draws of the fixed effects
fit_ob <- function(y, d) { f <- as.formula(paste(y, "~ pf50 + treatment + I(pf50^2) +", cov, "+", re))
  for (ct in list(glmmTMBControl(), glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))) {
    m <- tryCatch(glmmTMB(f, data = d, family = ordbeta(), control = ct), error = function(e) NULL)
    if (!is.null(m) && isTRUE(m$sdr$pdHess)) return(m) }
  stop("ordbeta did not converge for ", y) }
ob_mean <- function(eta, psi) { p0 <- 1 - plogis(eta - psi[1]); p1 <- plogis(eta - psi[2]); p1 + (1 - p0 - p1) * plogis(eta) }
fr <- c(SMF = "New stem", LMF = "Leaf", RMF = "Root")
fe_ob <- delete.response(terms(as.formula(paste("y ~ pf50 + treatment + I(pf50^2) +", cov))))
c2c <- imap_dfr(set_names(names(fr)), function(nm, ...) {
  m <- fit_ob(nm, a414); b <- fixef(m)$cond; V <- vcov(m)$cond; p <- m$fit$par; psi <- unname(p[names(p) == "psi"])
  sims <- MASS::mvrnorm(1000, b, V)
  map_dfr(grid_anc, function(g) map_dfr(TRT_LEVELS, function(trt) {
    nd <- a414 %>% select(og_drc_mm, myco_cat, soil, treatment); nd$pf50 <- (g - 50) / 50; nd$treatment <- factor(trt, levels = TRT_LEVELS)
    xm <- colMeans(model.matrix(fe_ob, nd)[, names(b), drop = FALSE]); dr <- ob_mean(as.vector(sims %*% xm), psi)
    tibble(fraction = fr[[nm]], x = g, treatment = trt, y = ob_mean(sum(xm * b), psi), lower = quantile(dr, 0.025), upper = quantile(dr, 0.975)) })) }) %>%
  mutate(treatment = factor(treatment, levels = TRT_LEVELS), fraction = factor(fraction, levels = fr))
fr_pts <- a414 %>% pivot_longer(c(SMF, LMF, RMF), names_to = "fraction", values_to = "y") %>% mutate(fraction = factor(fr[fraction], levels = fr))
p2c <- ggplot() + pts(fr_pts, "per_f", "y", size = PT_SIZE * 0.8, alpha = 0.3) + band(c2c) + facet_wrap(~fraction, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) + trt_scales() + labs(x = X_ANC, y = "Fraction of total new biomass") + theme_fig()
fig2 <- ((p2a + (p2b + theme(legend.position = "none"))) / (p2c + theme(legend.position = "none"))) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A", theme = theme(legend.position = "bottom"))
save_fig(fig2, "Fig2", W_2COL, 16 * CM)
readr::write_excel_csv(bind_rows(EDGE), "results/28_figures/fig2_points_at_edge.csv")

# ------------------------------------------------------------------------------
# Figure 3: colonization by guild along Fremont ancestry (colonization LMM of r_code/21_colonization_all_scored_plants.R and 08_sample_flow_extra_tests.R;
# 78 plants with all three guilds scored). One line per guild, averaged over watering; points colored by watering.
# ------------------------------------------------------------------------------
guilds <- c(amf = "AMF", emf = "EMF", dse = "DSE")
d3 <- ds %>% select(tree_id, pop, per_f, treatment, og_drc_mm, amf, emf, dse, soil, myco_cat) %>% na.omit() %>%
  pivot_longer(names(guilds), names_to = "guild", values_to = "colonization") %>%
  mutate(guild = factor(guild, levels = names(guilds)), og_drc_mm_z = (og_drc_mm - mean(og_drc_mm)) / sd(og_drc_mm))
col_m <- lmer(colonization ~ per_f * guild * treatment + og_drc_mm_z * guild + soil + (1 | pop) + (1 | myco_cat) + (1 | tree_id), data = d3)
stopifnot(nobs(col_m) == 234)
fe3 <- delete.response(terms(colonization ~ per_f * guild * treatment + og_drc_mm_z * guild + soil))
c3 <- map_dfr(grid_anc, function(pf) map_dfr(names(guilds), function(g) {
  x <- colMeans(model.matrix(fe3, d3 %>% mutate(per_f = pf, guild = factor(g, levels = names(guilds))))[, names(fixef(col_m))])
  tibble(per_f = pf, guild = g, fit = sum(x * fixef(col_m)), se = sqrt(as.numeric(t(x) %*% as.matrix(vcov(col_m)) %*% x))) })) %>%
  mutate(guild = factor(guild, levels = names(guilds)), lower = fit - 1.96 * se, upper = fit + 1.96 * se)
fig3 <- ggplot() + pts(d3, "per_f", "colonization", position = "identity") +
  geom_ribbon(data = c3, aes(per_f, ymin = lower, ymax = upper), fill = NEUTRAL_BAND, alpha = BAND_ALPHA) +
  geom_line(data = c3, aes(per_f, fit), colour = NEUTRAL_LINE, linewidth = LN_W) +
  panel_letters(tibble(guild = factor(names(guilds), levels = names(guilds)), name = guilds), "guild", "name", 50) +
  facet_wrap(~guild, nrow = 1) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  trt_scales() + labs(x = X_ANC, y = "Colonization (%)") + theme_fig() + theme(strip.background = element_blank(), strip.text = element_blank())
save_fig(fig3, "Fig3", W_2COL, 7.5 * CM)

# ------------------------------------------------------------------------------
# Figure 4: new stem biomass against colonization (colonization-biomass LMM of r_code/02a_main_models.R and 20_emf_drought_influence.R;
# 77 plants). The drought plant with 0.98 g new stem mass, on which the drought slope depends, is circled.
# ------------------------------------------------------------------------------
m5 <- sapply(names(guilds), function(g) mean(ds[[g]], na.rm = TRUE)); s5 <- sapply(names(guilds), function(g) sd(ds[[g]], na.rm = TRUE))
b5 <- ds %>% select(tree_idold, geno_id, pop, treatment, per_fc, og_drc_mm, myco_cat, soil, amf, emf, dse, dryshoot_new) %>% na.omit() %>%
  mutate(amf_z = (amf - m5["amf"]) / s5["amf"], emf_z = (emf - m5["emf"]) / s5["emf"], dse_z = (dse - m5["dse"]) / s5["dse"])
f5 <- log(dryshoot_new) ~ (amf_z + emf_z + dse_z) * treatment + per_fc + I(per_fc^2) + og_drc_mm + myco_cat + soil
bm5 <- lmerTest::lmer(update(f5, . ~ . + (1 | pop) + (1 | pop:geno_id)), data = b5); stopifnot(nobs(bm5) == 77)
hi <- b5 %>% filter(abs(dryshoot_new - 0.98) < 0.005); stopifnot(nrow(hi) == 1)
p4 <- function(g, first) {
  cd <- curve_trt(fixef(bm5), as.matrix(vcov(bm5)), delete.response(terms(f5)), b5, seq(min(b5[[g]]), max(b5[[g]]), length.out = 50),
                  paste0(g, "_z"), function(v) (v - m5[g]) / s5[g], exp)
  ggplot() + pts(b5, g, "dryshoot_new", position = "identity") +
    geom_point(data = hi, aes(.data[[g]], dryshoot_new), shape = 1, size = 3.4, stroke = 0.5, colour = "black") + band(cd) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) + trt_scales() +
    labs(x = paste(guilds[[g]], "colonization (%)"), y = if (first) "New stem biomass (g dry mass)" else NULL) + theme_fig() }
fig4 <- p4("amf", TRUE) + p4("emf", FALSE) + p4("dse", FALSE) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
save_fig(fig4, "Fig4", W_2COL, 8 * CM)

# ------------------------------------------------------------------------------
# Figures S2 and S3: colonization by inoculum and by watering treatment (every plant scored for each guild:
# AMF and DSE 91 plants, EMF 166)
# ------------------------------------------------------------------------------
dl <- ds %>% pivot_longer(names(guilds), names_to = "guild", values_to = "colonization") %>%
  filter(!is.na(colonization), !is.na(myco_cat)) %>%
  mutate(inoculum = factor(as.character(myco_cat), levels = c("None", "EM", "AM", "Both"), labels = c("Control", "EMF", "AMF", "Dual")),
         guild = factor(guild, levels = names(guilds)))
lab_g <- tibble(guild = factor(names(guilds), levels = names(guilds)), name = guilds)
box_theme <- theme(strip.background = element_blank(), strip.text = element_blank(), panel.grid.major.x = element_blank())
s2 <- ggplot(dl, aes(inoculum, colonization)) +
  geom_boxplot(fill = NEUTRAL_FILL, colour = "grey30", linewidth = 0.3, outlier.shape = NA, width = 0.6) +
  geom_point(aes(fill = treatment, shape = treatment), position = position_jitter(width = 0.15, height = 0, seed = 1),
             colour = PT_EDGE, stroke = PT_STROKE, size = PT_SIZE, alpha = PT_ALPHA) +
  panel_letters(lab_g, "guild", "name", 2.5) + facet_wrap(~guild, nrow = 1) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  trt_scales() + labs(x = "Inoculum treatment", y = "Colonization (%)") + theme_fig() + box_theme
save_fig(s2, "FigS2", W_2COL, 7.5 * CM)
s3 <- ggplot(dl, aes(treatment, colonization)) +
  geom_boxplot(aes(fill = treatment), colour = "grey30", alpha = 0.5, linewidth = 0.3, outlier.shape = NA, width = 0.6) +
  geom_point(aes(fill = treatment, shape = treatment), position = position_jitter(width = 0.15, height = 0, seed = 1),
             colour = PT_EDGE, stroke = PT_STROKE, size = PT_SIZE, alpha = PT_ALPHA) +
  panel_letters(lab_g, "guild", "name", 1.5) + facet_wrap(~guild, nrow = 1) + scale_y_continuous(expand = expansion(mult = c(0.04, 0.14))) +
  trt_scales() + labs(x = "Watering treatment", y = "Colonization (%)") + theme_fig() + box_theme + theme(legend.position = "none")
save_fig(s3, "FigS3", W_2COL, 7 * CM)

# ------------------------------------------------------------------------------
# Figures S5.1 and S5.2: stomatal conductance, plant means by measurement session (data only; session is confounded
# with watering and ancestry). Input: results/19_conductance_models_of_record/plant_means.csv (r_code/19_conductance_models_of_record.R)
# ------------------------------------------------------------------------------
LV <- c("8/9 Sep", "17/22 Sep", "20/23 Sep", "21/26 Sep")
pm <- read_csv("results/19_conductance_models_of_record/plant_means.csv", show_col_types = FALSE) %>%
  mutate(day_pair = factor(day_pair, levels = LV), treatment = factor(treatment, levels = TRT_LEVELS))
stopifnot(nrow(pm) == 203)
n_lab <- pm %>% count(day_pair, treatment) %>% group_by(day_pair) %>%
  summarise(lab = paste(sprintf("%s %d", ifelse(treatment == "Well-watered", "WW", "D"), n), collapse = ", "), .groups = "drop")
strip_lab <- setNames(sprintf("%s (%s)", n_lab$day_pair, n_lab$lab), as.character(n_lab$day_pair))
# two lines: on one line the label was longer than the 7.5 cm panel and was clipped at the top of Figure S5.1
YLAB_GSW <- expression(atop("Stomatal conductance", "(mol·m"^"–2" * "·s"^"–1" * ", log"[10] * " scale)"))
s5_common <- function() list(geom_point(aes(fill = treatment, shape = treatment), position = position_jitter(width = 1, height = 0, seed = 1),
                                        colour = PT_EDGE, stroke = PT_STROKE, size = PT_SIZE, alpha = PT_ALPHA),
                             scale_x_continuous(limits = c(-2, 102), breaks = c(0, 50, 100)), trt_scales(), labs(x = X_ANC, y = YLAB_GSW), theme_fig())
f51 <- ggplot(pm, aes(per_f, whole_leaf)) + s5_common() +
  panel_letters(tibble(day_pair = factor(LV, levels = LV)), "day_pair") +
  facet_wrap(~day_pair, nrow = 1, labeller = as_labeller(strip_lab), drop = FALSE) + scale_y_log10(expand = expansion(mult = c(0.04, 0.14)))
save_fig(f51, "AppendixS5_daypair", W_2COL, 7.5 * CM)
sl <- pm %>% pivot_longer(c(adaxial, abaxial), names_to = "surface", values_to = "gsw") %>%
  mutate(surface = factor(surface, levels = c("adaxial", "abaxial"), labels = c("Adaxial (upper) surface", "Abaxial (lower) surface")))
dropped <- sl %>% filter(is.na(gsw) | gsw <= 0) %>% count(surface)
write.csv(dropped, "results/28_figures/figure_S5_2_dropped.csv", row.names = FALSE)
sl <- sl %>% filter(!is.na(gsw), gsw > 0)
f52 <- ggplot(sl, aes(per_f, gsw)) + s5_common() +
  panel_letters(expand_grid(surface = factor(levels(sl$surface), levels = levels(sl$surface)), day_pair = factor(LV, levels = LV)), "day_pair") +
  facet_grid(surface ~ day_pair, scales = "free_y", labeller = labeller(day_pair = as_labeller(strip_lab)), drop = FALSE) +
  scale_y_log10(expand = expansion(mult = c(0.04, 0.16)))
save_fig(f52, "AppendixS5_surfaces_daypair", W_2COL, 12.5 * CM)
cat("Done.\n")
