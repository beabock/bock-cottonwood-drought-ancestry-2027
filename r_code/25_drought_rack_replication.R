# ==============================================================================
# Drought replication when racks (6-8 plants each, separated by watering) were not recorded.
# Rack identity is not in the data (no rack, bench or position column in mega_ds.csv), so three approaches bound the effect:
#   A  df adjustment: the reported t (or z) for drought is referred to a t distribution with df = racks - 2, where
#      racks = planted cuttings per watering / plants per rack (6, 7 or 8), summed over watering treatments.
#   B  pseudo-rack aggregation: cuttings are assigned to pseudo-racks of 7 within watering treatment, and rack means are
#      analysed by weighted least squares (weights = plants per rack): response ~ drought + mean ancestry (linear and quadratic)
#      + mean initial diameter. Two assignments: (B1) random within watering, 1000 partitions, and (B2) within genotype x
#      watering (racks of up to 7 cuttings from one mother tree, so rack and genotype coincide; deterministic).
#   C  design effect: with m = 7 plants per rack and rack intraclass correlation rho, SE x sqrt(1 + 6 rho); the rho at
#      which the drought P (df = racks - 2) reaches 0.05.
# Assumptions: plants were allocated to racks without regard to genotype or ancestry (B1), and racks within a watering
# treatment are exchangeable; rack means remove within-rack replication but cannot recover a true rack variance that was
# never measured. C states how large that variance would need to be.
# Responses: survival (563 cuttings, GLMM of record), total new biomass with estimated roots (414), leafless at harvest (425).
# Run from the project root:  Rscript r_code/25_drought_rack_replication.R
# Writes results/25_drought_rack_replication/: output.txt, strings.csv, racks.csv
# ==============================================================================
source("r_code/helpers/model_helpers.R")
OUT <- "results/25_drought_rack_replication"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sink(file.path(OUT, "output.txt"), split = TRUE); run_stamp("Drought rack replication"); set.seed(1998)
src <- "results/25_drought_rack_replication/racks.csv"
cat("Columns matching rack/bench/block/position in the data:", length(grep("rack|bench|block|pos", names(ds), ignore.case = TRUE)), "\n")
planted <- table(ds$treatment); print(planted)
red <- paste("per_fc + treatment + I(per_fc^2) +", cov, "+", re)
surv <- glmer(alive_01_10_2021 ~ pf50 + treatment + I(pf50^2) + pop + og_drc_mm + soil + myco_cat + (1 | geno_id), data = ds, family = binomial, control = ctl)
tot <- lmerTest::lmer(as.formula(paste("log(total_new) ~", red)), data = al)
lz <- ds %>% filter(alive_01_10_2021 %in% 1, !is.na(dryshoot_new)) %>% mutate(zero_leaf = as.integer(coalesce(dryleaf_weight, 0) == 0)) %>%
  filter(if_all(c(pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), ~ !is.na(.x)))
leaf <- glmer(zero_leaf ~ pf50 + treatment + I(pf50^2) + og_drc_mm + myco_cat + soil + (1 | pop) + (1 | pop:geno_id), data = lz, family = binomial, control = ctl)
stopifnot(nobs(surv) == 563, nobs(tot) == 414, nobs(leaf) == 425)
mods <- list(Survival = list(m = surv, stat = "z", y = "alive_01_10_2021", ids = ds$tree_idold[ds$tree_idold %in% ds$tree_idold & !is.na(ds$alive_01_10_2021)]),
             `Total new biomass (log)` = list(m = tot, stat = "t", y = "log_total_new"), `Leafless at harvest` = list(m = leaf, stat = "z", y = "zero_leaf"))
rows <- list()
hdr("A. Reported drought test with df = number of racks - 2")
for (k in names(mods)) { m <- mods[[k]]$m; cs <- coef(summary(m))["treatmentDrought", ]; tv <- unname(cs[grep(" value$", names(cs))])
  rep_s <- if (mods[[k]]$stat == "z") zstr(cs[1], tv, cs[4]) else tstr(cs[1], tv, cs["df"], cs["Pr(>|t|)"])
  for (pp in c(6, 7, 8)) { R <- sum(round(planted / pp)); p <- 2 * pt(-abs(tv), R - 2)
    rows[[length(rows) + 1]] <- tibble(response = k, method = sprintf("A: df = racks − 2 (%d plants per rack, %d racks)", pp, R), estimate = unname(cs[1]), stat = tv, df = R - 2, P = p,
                                       string = sprintf("%s = %s, t = %s, df = %d, %s", BETA, num(cs[1]), dec(tv), R - 2, P(p)), reported = rep_s) }
  # C: design effect tipping point at 7 plants per rack
  R <- sum(round(planted / 7)); f <- function(rho) 2 * pt(-abs(tv / sqrt(1 + 6 * rho)), R - 2) - 0.05
  rho <- if (f(0) > 0) 0 else if (f(1) < 0) NA else uniroot(f, c(0, 1))$root
  rows[[length(rows) + 1]] <- tibble(response = k, method = "C: rack ICC at which drought P = 0.05 (7 per rack, df = racks − 2)", estimate = rho, stat = NA, df = R - 2, P = 0.05,
                                     string = ifelse(is.na(rho), "P < 0.05 for any ICC ≤ 1", sprintf("ICC = %.2f (design effect %.2f)", rho, 1 + 6 * rho)), reported = rep_s) }
hdr("B. Pseudo-rack means")
pr_data <- ds %>% select(tree_idold, treatment, geno_id, per_fc, og_drc_mm) %>%
  left_join(ds %>% filter(tree_idold %in% (ds %>% filter(if_all(c(alive_01_10_2021, pf50, treatment, pop, og_drc_mm, soil, myco_cat, geno_id), ~ !is.na(.x))) %>% pull(tree_idold))) %>%
              select(tree_idold, surv = alive_01_10_2021), by = "tree_idold") %>%
  left_join(al %>% transmute(tree_idold, ltot = log(total_new)), by = "tree_idold") %>% left_join(lz %>% select(tree_idold, zero_leaf), by = "tree_idold")
yv <- c(Survival = "surv", `Total new biomass (log)` = "ltot", `Leafless at harvest` = "zero_leaf")
rack_fit <- function(d) map_dfr(names(yv), function(k) { y <- yv[[k]]
  rm <- d %>% filter(!is.na(.data[[y]]), !is.na(og_drc_mm)) %>% group_by(rack, treatment) %>%
    summarise(y = mean(.data[[y]]), a = mean(per_fc), a2 = mean(per_fc^2), dia = mean(og_drc_mm), n = n(), .groups = "drop")
  fit <- lm(y ~ treatment + a + a2 + dia, data = rm, weights = n); cs <- coef(summary(fit))["treatmentDrought", ]
  tibble(response = k, racks = nrow(rm), estimate = cs[1], t = cs[3], df = fit$df.residual, P = cs[4]) })
assign_random <- function(d) d %>% group_by(treatment) %>% mutate(rack = paste(treatment, (sample(n()) - 1) %/% 7)) %>% ungroup()
assign_geno <- function(d) d %>% arrange(treatment, geno_id, tree_idold) %>% group_by(treatment, geno_id) %>% mutate(rack = paste(treatment, geno_id, (row_number() - 1) %/% 7)) %>% ungroup()
b1 <- map_dfr(1:1000, function(i) rack_fit(assign_random(pr_data)) %>% mutate(iter = i))
b2 <- rack_fit(assign_geno(pr_data))
for (k in names(yv)) { x <- filter(b1, response == k); g <- filter(b2, response == k)
  rows[[length(rows) + 1]] <- tibble(response = k, method = sprintf("B1: random pseudo-racks of 7 within watering, 1000 partitions (%d racks)", round(median(x$racks))),
    estimate = median(x$estimate), stat = median(x$t), df = median(x$df), P = median(x$P),
    string = sprintf("median %s = %s, t = %s, df = %s, %s; P < 0.05 in %.1f%% of partitions (P range %s–%s)", BETA, num(median(x$estimate)), dec(median(x$t)), dff(median(x$df)), P(median(x$P)),
                     100 * mean(x$P < 0.05), formatC(min(x$P), format = "g", digits = 2), formatC(max(x$P), format = "g", digits = 2)), reported = NA)
  rows[[length(rows) + 1]] <- tibble(response = k, method = sprintf("B2: pseudo-racks of up to 7 within genotype × watering (%d racks)", g$racks), estimate = g$estimate, stat = g$t, df = g$df, P = g$P,
    string = tstr(g$estimate, g$t, g$df, g$P), reported = NA) }
tab <- bind_rows(rows) %>% group_by(response) %>% fill(reported, .direction = "downup") %>% ungroup()
print(as.data.frame(tab %>% select(response, method, string)), right = FALSE, row.names = FALSE)
readr::write_excel_csv(tab, file.path(OUT, "racks.csv"))
cat("\nB uses the rack mean of the response (survival and leafless as proportions, total new biomass as mean log mass); estimates are on that scale, not the GLMM log-odds scale.\n")
for (i in seq_len(nrow(tab))) rec(7, "", paste("Drought effect:", tab$response[i]), tab$method[i], NA, "Drought", tab$string[i], tab$P[i], src)
for (k in names(mods)) rec(7, "", paste("Drought effect:", k), "reported model", nobs(mods[[k]]$m), "Drought (reported)", tab$reported[tab$response == k][1], NA, src)
write_rec("strings.csv"); cat("\nDone.\n"); sink()
