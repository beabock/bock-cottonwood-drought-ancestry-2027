# ==============================================================================
# Appendix S6 (was S8): biomass allocation and leaf drop; source table for the publication table.
#   Primary column  = sample (a): estimated dry root mass for colonization-scored plants (n = 414)
#   Comparison      = sample (b): colonization-scored plants excluded, measured root mass only (n = 242)
# Allocation strings are read from results/06_allocation_models/allocation_strings.csv (run
# r_code/06_allocation_models.R first). The leaf-drop model does not use root mass;
# its comparison column refits it on the unscored plants so both columns follow the same split.
# ==============================================================================
# Run from the project root (R 4.6.1):  Rscript r_code/07_allocation_appendix_S7.R
# Writes: results/07_allocation_appendix_S7/appendix_S7_source.csv (formatted into Appendix S6 by r_code/30_appendix_tables.R)
# ==============================================================================
source("r_code/helpers/utilities_data_prep.R")
library(lme4); library(lmerTest); library(tidyverse)
EN <- "–"; BETA <- "β"; TIMES <- "×"; CHI <- "χ²"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
P <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)

st <- read.csv("results/06_allocation_models/allocation_strings.csv", encoding = "UTF-8", check.names = FALSE)
resp <- c("log total new biomass" = "log total new biomass", "stem mass fraction" = "Stem mass fraction",
          "leaf mass fraction" = "Leaf mass fraction", "root mass fraction" = "Root mass fraction", "log root:shoot" = "log root:shoot")
alloc <- st %>% filter(section == "2", response %in% names(resp)) %>%
  mutate(sample = ifelse(grepl("^\\(a\\)", sample), "a", "b"), response = resp[response],
         method = ifelse(grepl("ordbeta", method), "Ordered beta GLMM", "LMM"))

# --- leaf drop: (a) = all harvested survivors (from 06_allocation_models.R), (b) = refit without colonization-scored plants ---
ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, pf50 = per_fc / 50, geno_id = paste(pop, type, geno, sep = "_"))
lz <- ds %>% filter(alive_01_10_2021 == 1, !is.na(dryshoot_new)) %>%
  mutate(zero_leaf = as.integer(coalesce(dryleaf_weight, 0) == 0), scored = !is.na(amf) | !is.na(emf) | !is.na(dse)) %>%
  filter(if_all(c(pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), ~ !is.na(.x)))
stopifnot(nrow(lz) == st$n[st$section == "4"][1])          # same sample as 06_allocation_models.R section 4
lb <- filter(lz, !scored)
cov <- "og_drc_mm + myco_cat + soil"; re <- "(1 | pop) + (1 | pop:geno_id)"
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
zf <- glmer(as.formula(paste("zero_leaf ~ pf50 * treatment + I(pf50^2) * treatment +", cov, "+", re)), data = lb, family = binomial, control = ctl)
zr <- glmer(as.formula(paste("zero_leaf ~ pf50 + treatment + I(pf50^2) +", cov, "+", re)), data = lb, family = binomial, control = ctl)
zl <- anova(zr, zf); pz <- zl$`Pr(>Chisq)`[2]; zu <- if (pz < 0.05) zf else zr
cs <- coef(summary(zu)); div <- c(pf50 = 50, "I(pf50^2)" = 2500, treatmentDrought = 1)
tl <- c(pf50 = "Linear ancestry", "I(pf50^2)" = "Quadratic ancestry", treatmentDrought = "Drought")
leaf_b <- tibble(response = "P(no leaves at harvest)", sample = "b", n = nobs(zu), method = "Binomial GLMM",
                 term = c(paste("LRT ancestry", TIMES, "drought"), tl),
                 string = c(sprintf("%s(%d) = %s, %s", CHI, zl$Df[2], dec(zl$Chisq[2]), P(pz)),
                            sapply(names(tl), function(t) sprintf("%s = %s, z = %s, %s", BETA, num(cs[t, "Estimate"] / div[t]), dec(cs[t, "z value"]), P(cs[t, "Pr(>|z|)"])))),
                 est = c(NA, cs[names(tl), "Estimate"] / div), p = c(pz, cs[names(tl), "Pr(>|z|)"]))
if (pz < 0.05) warning("Leaf-drop interaction significant in sample (b); full-model terms not tabulated")
cat("Leaf drop, sample (b): n =", nobs(zu), "| zero leaf =", sum(lb$zero_leaf), "| singular =", isSingular(zu), "\n")
leaf_a <- st %>% filter(section == "4") %>% transmute(response = "P(no leaves at harvest)", sample = "a", n, method = "Binomial GLMM", term, string, est, p)

# --- root:new-stem ratio: core LMM on log(root / new stem) in both samples, same rule as the other LMMs ---
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
rsd <- read.csv("results/06_allocation_models/allocation_data.csv") %>% select(tree_idold, scored, root_est, new_stem) %>%
  inner_join(ds %>% select(tree_idold, per_fc, treatment, og_drc_mm, myco_cat, soil, pop, geno_id), by = "tree_idold") %>% mutate(root_stem = root_est / new_stem)
stopifnot(nrow(rsd) == 414)
rs_rows <- map_dfr(c(a = "a", b = "b"), function(s) {
  d <- if (s == "a") rsd else filter(rsd, !scored)
  full <- lmerTest::lmer(as.formula(paste("log(root_stem) ~ per_fc * treatment + I(per_fc^2) * treatment +", cov, "+", re)), data = d)
  red <- lmerTest::lmer(as.formula(paste("log(root_stem) ~ per_fc + treatment + I(per_fc^2) +", cov, "+", re)), data = d)
  l <- anova(red, full, refit = TRUE); pl <- l$`Pr(>Chisq)`[2]
  if (pl < 0.05) warning("Root:new-stem interaction significant in sample ", s, "; full-model terms not tabulated")
  cr <- coef(summary(red)); k <- c(per_fc = "Linear ancestry", "I(per_fc^2)" = "Quadratic ancestry", treatmentDrought = "Drought")
  tibble(response = "log root:new-stem ratio", sample = s, n = nobs(red), method = "LMM", term = c(paste("LRT ancestry", TIMES, "drought"), unname(k)),
         string = c(sprintf("%s(%d) = %s, %s", CHI, l$Df[2], dec(l$Chisq[2]), P(pl)),
                    sapply(names(k), function(t) sprintf("%s = %s, t = %s, df = %s, %s", BETA, num(cr[t, "Estimate"]), dec(cr[t, "t value"]), dff(cr[t, "df"]), P(cr[t, "Pr(>|t|)"])))),
         est = c(NA, cr[names(k), "Estimate"]), p = c(pl, cr[names(k), "Pr(>|t|)"])) })
cat("Root:new-stem ratio: n =", paste(unique(rs_rows$n), collapse = " / "), "\n")

s8 <- bind_rows(alloc %>% select(response, sample, n, method, term, string, est, p), rs_rows, leaf_a, leaf_b) %>%
  mutate(term = sub("^LRT ancestry × drought$", paste("Ancestry", TIMES, "drought (LRT, 2 df)"), term)) %>%
  pivot_wider(names_from = sample, values_from = c(n, string, est, p)) %>%
  mutate(response = factor(response, levels = c(unname(resp), "log root:new-stem ratio", "P(no leaves at harvest)")),
         term = factor(term, levels = c(paste("Ancestry", TIMES, "drought (LRT, 2 df)"), "Linear ancestry", "Quadratic ancestry", "Drought"))) %>%
  arrange(response, term) %>% mutate(across(c(response, term), as.character),
  conclusion_changes = (p_a < 0.05) != (p_b < 0.05) | (p_a < 0.05 & p_b < 0.05 & !is.na(est_a) & sign(est_a) != sign(est_b)))
names(s8) <- sub("_a$", "_estimated_root", sub("_b$", "_measured_root_only", names(s8)))
readr::write_excel_csv(s8, "results/07_allocation_appendix_S7/appendix_S7_source.csv")
cat("Appendix S7 source written: results/07_allocation_appendix_S7/appendix_S7_source.csv |", sum(s8$conclusion_changes), "rows change conclusion\n")
