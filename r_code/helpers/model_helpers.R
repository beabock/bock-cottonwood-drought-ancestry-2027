# ==============================================================================
# Shared helpers for steps 20-29 (r_code/20_emf_drought_influence.R to r_code/29_sensitivity_checks.R).
# AJB strings (β = x, t = x, df = x, P = x; z for GLMMs; en-dash minus; capital P), the data as the
# reported models load it (r_code/helpers/utilities_data_prep.R), and a collector for each script's strings.csv rows.
# Source from the project root.
# ==============================================================================
suppressMessages({ source("r_code/helpers/utilities_data_prep.R"); library(tidyverse); library(lme4); library(lmerTest); library(emmeans) })
emm_options(lmerTest.limit = 5000, pbkrtest.limit = 5000)
EN <- "–"; BETA <- "β"; TIMES <- "×"; CHI <- "χ²"
dash <- function(s) gsub("-", EN, s, fixed = TRUE)
num <- function(x, d = 3) dash(formatC(x, digits = d, format = "fg", flag = "#"))
dec <- function(x, k = 2) dash(sprintf(paste0("%.", k, "f"), x))
dff <- function(df) if (abs(df - round(df)) < 0.05) sprintf("%d", round(df)) else sprintf("%.1f", df)
P <- function(p) if (is.na(p)) "P = NA" else if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
tstr <- function(b, t, df, p) sprintf("%s = %s, t = %s, df = %s, %s", BETA, num(b), dec(t), dff(df), P(p))
zstr <- function(b, z, p) sprintf("%s = %s, z = %s, %s", BETA, num(b), dec(z), P(p))
pc <- function(x) dash(sprintf("%.1f", x))
hdr <- function(x) cat("\n\n==============================================================================\n", x,
                       "\n==============================================================================\n", sep = "")
run_stamp <- function(name) cat(name, "run:", "|", R.version.string, "| lme4", as.character(packageVersion("lme4")),
                                "| lmerTest", as.character(packageVersion("lmerTest")), "| emmeans", as.character(packageVersion("emmeans")), "\n")
# coefficient row of a lmerTest fit as an AJB string
row_t <- function(m, t, scale = 1) { r <- coef(summary(m))[t, ]
  list(s = tstr(r["Estimate"] / scale, r["t value"], r["df"], r["Pr(>|t|)"]), b = unname(r["Estimate"]) / scale, t = unname(r["t value"]),
       df = unname(r["df"]), p = unname(r["Pr(>|t|)"]), se = unname(r["Std. Error"]) / scale) }
row_z <- function(m, t, scale = 1) { cs <- coef(summary(m)); if (is.list(cs) && "cond" %in% names(cs)) cs <- cs$cond; r <- cs[t, ]
  list(s = zstr(r[1] / scale, r[3], r[4]), b = unname(r[1]) / scale, z = unname(r[3]), p = unname(r[4])) }
# linear contrast L'beta with Satterthwaite df (lmerTest::contest)
lin <- function(m, L) { b <- fixef(m); v <- setNames(rep(0, length(b)), names(b)); v[names(L)] <- L
  ct <- lmerTest::contest(m, v, joint = FALSE); list(est = ct$Estimate, se = ct$`Std. Error`, t = ct$`t value`, df = ct$df, p = ct$`Pr(>|t|)`,
                                                       s = tstr(ct$Estimate, ct$`t value`, ct$df, ct$`Pr(>|t|)`)) }

ds <- load_and_prepare_data() %>% mutate(per_fc = per_f - 50, pf50 = per_fc / 50, geno_id = paste(pop, type, geno, sep = "_"))
cov <- "og_drc_mm + myco_cat + soil"; re <- "(1 | pop) + (1 | pop:geno_id)"
ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
# allocation data of record (r_code/06_allocation_models.R; leaf = 0 g when leafless; root mass estimated for scored plants)
al <- read_csv("results/06_allocation_models/allocation_data.csv", show_col_types = FALSE) %>%
  select(tree_idold, scored, root_est, new_stem, leaf, total_new, LMF, RMF, SMF) %>%
  inner_join(ds %>% select(tree_idold, per_f, per_fc, pf50, treatment, og_drc_mm, myco_cat, soil, pop, geno_id, new_class, wet_root, root_dryweight), by = "tree_idold")
stopifnot(nrow(al) == 414)

# strings.csv collector: section, item, analysis, sample, n, term, string, P, source (file written by the script)
REC <- list()
rec <- function(section, item, analysis, sample, n, term, string, p = NA, source = NA) {
  REC[[length(REC) + 1]] <<- tibble(section = as.character(section), item = item, analysis = analysis, sample = sample, n = as.character(n), term = term,
                                    string = string, p = p, source = source)
  cat(sprintf("  [%s%s] %s | %s | n = %s | %s: %s\n", section, item, analysis, sample, n, term, string)) }
write_rec <- function(file) { r <- bind_rows(REC); readr::write_excel_csv(r, file.path(OUT, file)); invisible(r) }
