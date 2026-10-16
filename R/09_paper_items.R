#-----------------------------------------------------------------------------------
#  09_paper_items.R  —  Create the manuscript (main text and Supplemental Material) figures and tables,
#                       under their manuscript numbers, in output/00_paper_figures_tables/
#
#   Main text: Table 1 (G1 baseline characteristics; all patients and discharge mRS 0–2 / 3–5), Table 2 (years over 15 years
#              and 15-year outcomes by discharge mRS), Figure 1 (mRS distribution at each assessment), Figure 2 (how the
#              15 years were spent by discharge mRS), Figure 3 (survival and independent survival over time),
#              Figure 4 (how the 15 years were spent by age group)
#   Supplemental Material: Tables S2–S6, Figures S1–S4 (Table S1, the list of imputation-model predictors, is written directly in the Supplemental Material)
#   Zenodo: Tables Z1–Z5 (moved from the Supplement on 2026-10-15)
#   - The object names below (s4 to s11) and the keys in paper_values.json ("Table S4", etc.) keep the numbers used before the move.
#     The correspondence with the manuscript numbers is given by PAPER_NO.
#   - Estimates come from 04_estimates.R (work/estimates_all.rds), 05_prestroke.R (estimates_prestroke02.rds), and 06_years_lived.R (estimates_years.rds).
#     Only the proportions by 4 age groups × discharge mRS and the differences between onset-period cohorts are aggregated here
#     from the imputed data (MAR) (imputation is not rerun).
#   - Tables go to Tables_main_and_supplement.xlsx (sheet = table number); values cited in the text go to its Text values sheet.
#     The same content is also saved to paper_values.json (values used in the manuscript text and tables, JSON).
#-----------------------------------------------------------------------------------

pd <- dir_create(OUT_DIR, "00_paper_figures_tables")
unlink(list.files(pd, full.names = TRUE))
co <- readRDS(file.path(WORK_DIR, "cohorts.rds")); cohorts <- co$cohorts
est_all <- readRDS(file.path(WORK_DIR, "estimates_all.rds"))
yrs <- readRDS(file.path(WORK_DIR, "estimates_years.rds"))
pre <- readRDS(file.path(WORK_DIR, "estimates_prestroke02.rds"))
H15 <- max(HORIZONS); H5 <- min(HORIZONS); VARS <- c("alive", "indep", "dep")
k15 <- names(cohorts)[sapply(cohorts, `[[`, "H") == H15]
G_OF <- c("15" = "G1", "10" = "G2", "5" = "G3")
coh <- cohorts[[k15]]; d <- coh$data; Yobs <- as.matrix(d[, coh$tp$col])
base <- d[, c("mrs_dc", "age", "male")]
e15 <- est_all %>% filter(method == "MAR", H == H15)

G3S <- c("All", "Discharge mRS 0-2", "Discharge mRS 3-5"); G6 <- paste0("Discharge mRS ", 0:5); G9 <- c(G3S, G6)
glab <- function(g) ifelse(g %in% c("All", "All patients"), "All patients", sub("Discharge mRS ", "", sub("-", "–", g)))   # "Discharge mRS" column of the tables
big <- function(x) format(x, big.mark = ",", trim = TRUE)
f1 <- function(x) formatC(x, format = "f", digits = 1)
yci <- function(e, l, u) sprintf("%s (%s–%s)", f1(e), f1(pmax(l, 0)), f1(u))                    # Years (95% CI; lower limit truncated at 0)
pick <- function(m, g, o, t, E = est_all, H = H15) E[E$method == m & E$H == H & E$group == g & E$outcome == o & E$time == t, ]
cell <- function(x) if (nrow(x)) fmt_ci(x$est, x$lcl, x$ucl) else "–"
DELTA_M <- sprintf("MNAR_delta%+.1f", DELTAS_MNAR)
rng15 <- function(g, o, E = est_all) { v <- sapply(c("MAR", DELTA_M), function(m) pick(m, g, o, H15, E)$est); sprintf("%s–%s", fmt_pct(min(v)), fmt_pct(max(v))) }
ykey <- function(g) ifelse(g == "All", "All patients", g)                                       # Group names used in 06_years_lived.R

# =================================================================================================
#  Tables
# =================================================================================================

# ---- Table 1: G1 baseline characteristics (all patients, discharge mRS 0–2, 3–5). Observed data before imputation; denominators are non-missing patients for each variable -------------
med <- function(x, dg = 0) { q <- quantile(x, c(.5, .25, .75), na.rm = TRUE); sprintf(paste0("%.", dg, "f (%.", dg, "f–%.", dg, "f)"), q[1], q[2], q[3]) }
npct <- function(x) sprintf("%s (%.1f)", big(sum(x, na.rm = TRUE)), 100 * mean(x, na.rm = TRUE))
AGE_EN <- setNames(c("<65 y", "65–74 y", "75–84 y", "≥85 y"), AGE_LABELS)
TOAST_EN <- c(CES = "Cardioembolism", LAA = "Large-artery atherosclerosis", OTS = "Other or undetermined", SVO = "Small-vessel occlusion")
stopifnot(all(levels(d$toast) %in% names(TOAST_EN)))
bg_col <- function(x) { ag <- age_group(x$age)
  c(big(nrow(x)), med(x$age), sapply(AGE_LABELS, function(a) npct(ag == a)), npct(x$male == 1), npct(x$pre_mrs >= 3), med(x$nihss), "",
    sapply(names(TOAST_EN), function(l) npct(x$toast == l)), "", npct(x$ht == 1), npct(x$dm == 1), npct(x$dl == 1), npct(x$af == 1),
    npct(x$smk == 1), npct(x$drk == 1), npct(x$ckd == 1), npct(x$prior_stroke == 1), med(x$bmi, 1), npct(x$tpa == 1), npct(x$evt == 1)) }
t1 <- data.frame(Characteristic = c("Patients, n", "Age, y", paste0("  ", AGE_EN), "Male sex", "Prestroke mRS 3–5", "NIHSS score on admission",
                                    "Stroke subtype (TOAST)", paste0("  ", TOAST_EN), "Vascular risk factors and history", "  Hypertension", "  Diabetes",
                                    "  Dyslipidemia", "  Atrial fibrillation", "  Smoking", "  Alcohol drinking", "  Chronic kidney disease",
                                    "  Prior stroke", "Body mass index, kg/m²", "Intravenous thrombolysis", "Endovascular thrombectomy"),
                 `All patients` = bg_col(d), `Discharge mRS 0–2` = bg_col(d[d$mrs_dc <= 2, ]), `Discharge mRS 3–5` = bg_col(d[d$mrs_dc >= 3, ]), check.names = FALSE)
t1_miss <- sapply(c(ckd = "ckd", prior_stroke = "prior_stroke", bmi = "bmi", tpa = "tpa", evt = "evt"), function(v) sum(is.na(d[[v]])))

# ---- Table 2: years over 15 years and 15-year outcomes by discharge mRS, applied to the recent (G3) case mix ------------------------
rec <- yrs$recent %>% filter(strata == "Age × sex (8 strata)")
rec_cell <- function(g, v) { x <- rec[rec$group == ykey(g) & rec$value == v, ]; if (!nrow(x)) "–" else sprintf("%s [%s–%s]", f1(x$est), f1(x$range_min), f1(x$range_max)) }
t2 <- do.call(rbind, lapply(G9, function(g) { y <- yrs$by_dc[yrs$by_dc$group == ykey(g), ]
  data.frame(`Discharge mRS` = glab(g), n = big(y$n),
             `Years alive` = yci(y$alive_est, y$alive_lcl, y$alive_ucl), `Years alive and independent` = yci(y$indep_est, y$indep_lcl, y$indep_ucl),
             `Years alive and dependent` = yci(y$dep_est, y$dep_lcl, y$dep_ucl),
             `Alive at 15 y, %` = cell(pick("MAR", g, "Alive", H15)), `Alive and independent at 15 y, %` = cell(pick("MAR", g, "mRS 0-2", H15)),
             `Applied to recent case mix: years independent` = rec_cell(g, "indep"), `Applied to recent case mix: independent at 15 y, %` = rec_cell(g, "p15_indep"),
             check.names = FALSE) }))

# ---- Table S2: observation status at each assessment in the 3 cohorts -----------------------------------------------------------------
chs <- cohorts[order(-sapply(cohorts, `[[`, "H"))]                    # G1 → G2 → G3
s2 <- data.frame(Time = c(time_lab_en(coh$tp$lab), "Lost to follow-up (≥1 missing mRS)"), check.names = FALSE)
for (ch in chs) { Y <- as.matrix(ch$data[, ch$tp$col]); lost <- !complete.cases(Y)
  cn <- function(what) sprintf("%s (n=%s): %s", G_OF[[as.character(ch$H)]], big(nrow(Y)), what)
  pad <- function(x, last = "") c(x, rep("–", nrow(s2) - 1 - length(x)), last)
  s2[[cn("observed alive")]] <- pad(big(colSums(!is.na(Y) & Y < 6)))
  s2[[cn("observed dead")]]  <- pad(big(colSums(Y == 6, na.rm = TRUE)))
  s2[[cn("missing, n (%)")]] <- pad(sprintf("%s (%.1f)", big(colSums(is.na(Y))), 100 * colMeans(is.na(Y))), sprintf("%s (%.1f)", big(sum(lost)), 100 * mean(lost))) }

# ---- Table S3: characteristics of complete cases and patients lost to follow-up (G1, G2, G3) -------------------------------------------------
s3 <- data.frame(Characteristic = c("Age, y", "Male sex", "Prestroke mRS", "NIHSS score on admission", paste0("Stroke subtype: ", TOAST_EN),
                                    "Hypertension", "Diabetes", "Dyslipidemia", "Atrial fibrillation", "Smoking", "Alcohol drinking",
                                    "Chronic kidney disease", "Body mass index, kg/m²", "Prior stroke", "mRS at discharge"), check.names = FALSE)
for (ch in chs) { x <- ch$data; R <- complete.cases(x[, ch$tp$col]); G <- G_OF[[as.character(ch$H)]]
  one <- function(z) c(med(z$age), npct(z$male == 1), med(z$pre_mrs), med(z$nihss), sapply(names(TOAST_EN), function(l) npct(z$toast == l)),
                       npct(z$ht == 1), npct(z$dm == 1), npct(z$dl == 1), npct(z$af == 1), npct(z$smk == 1), npct(z$drk == 1), npct(z$ckd == 1),
                       med(z$bmi, 1), npct(z$prior_stroke == 1), med(z$mrs_dc))
  s3[[sprintf("%s: complete (n=%s)", G, big(sum(R)))]] <- one(x[R, ]); s3[[sprintf("%s: lost to follow-up (n=%s)", G, big(sum(!R)))]] <- one(x[!R, ]) }

# ---- Former Table S4 (Zenodo Table Z1 in the paper): mRS 0–6 distribution at each assessment by discharge mRS in G1 (%; mean across imputations) ------------------------------------
s4 <- e15 %>% filter(group %in% c("All", G6), outcome %in% paste0("mRS ", 0:6)) %>%
  mutate(val = fmt_pct(est), Group = sprintf("%s (n = %s)", ifelse(group == "All", "All patients", group), big(n)), Time = time_lab_en(time_lab)) %>%
  select(group, time, Group, Time, outcome, val) %>% tidyr::pivot_wider(names_from = outcome, values_from = val) %>%
  arrange(match(group, c("All", G6)), time) %>% select(-group, -time) %>% rename(`mRS 6 (death)` = `mRS 6`) %>% as.data.frame(check.names = FALSE)

# ---- Former Table S5 (Zenodo Table Z2 in the paper): survival and independent survival by discharge mRS in G1 (1, 5, 10, 15 y), range under missing-data assumptions, applied to recent case mix ----------
OUTL <- c("Alive" = "Survival", "mRS 0-2" = "Independent survival (alive with mRS 0–2)")
s5 <- do.call(rbind, lapply(names(OUTL), function(o) do.call(rbind, lapply(G9, function(g)
  data.frame(Outcome = OUTL[[o]], `Discharge mRS` = glab(g), n = big(pick("MAR", g, o, 1)$n),
             `1 y` = cell(pick("MAR", g, o, 1)), `5 y` = cell(pick("MAR", g, o, 5)), `10 y` = cell(pick("MAR", g, o, 10)), `15 y` = cell(pick("MAR", g, o, 15)),
             `15 y: range under missing-data assumptions` = rng15(g, o),
             `15 y: applied to recent case mix [range]` = rec_cell(g, ifelse(o == "Alive", "p15_alive", "p15_indep")), check.names = FALSE)))))

# ---- Former Table S6 (Zenodo Table Z3 in the paper): discharge mRS × 4 age groups. A = survival and independent survival (aggregated from imputed data), B = years over 15 years -------------------
DCG <- c(list("All" = 0:5, "Discharge mRS 0-2" = 0:2, "Discharge mRS 3-5" = 3:5), setNames(as.list(0:5), G6))
add <- list()
for (dg in names(DCG)) for (a in AGE_LABELS)
  add[[sprintf("AGE4 | %s | %s", dg, a)]] <- local({ kk <- DCG[[dg]]; aa <- a; function(b) b$mrs_dc %in% kk & age_group(b$age) == aa })
for (a in AGE_LABELS) for (s in c("men", "women"))
  add[[sprintf("AGE4SEX | %s | %s", a, s)]] <- local({ aa <- a; ss <- s; function(b) age_group(b$age) == aa & b$male == (ss == "men") })
GROUPS <- c(GROUPS[setdiff(names(GROUPS), names(add))], add)
imps15 <- readRDS(file.path(WORK_DIR, sprintf("imputations_%s.rds", k15)))$imps[["delta_+0.00"]]
est_sa <- estimate_mi(imps15, base, coh$tp, groups = names(add)) %>% filter(outcome %in% c("Alive", "mRS 0-2"))
rm(imps15); gc()
gr <- group_rows(base, names(add))
miss_at <- function(g, t) { r <- gr[[g]]; m <- sum(is.na(Yobs[r, match(t, coh$tp$year)])); sprintf("%d (%.1f)", m, 100 * m / length(r)) }
s6a <- do.call(rbind, lapply(names(OUTL), function(o) do.call(rbind, lapply(names(add), function(g) {
  p <- strsplit(g, " | ", fixed = TRUE)[[1]]; x <- est_sa[est_sa$group == g & est_sa$outcome == o, ]; at <- function(t) x[x$time == t, ]
  data.frame(Outcome = OUTL[[o]], `Discharge mRS` = if (p[1] == "AGE4SEX") paste0("All patients, ", AGE_EN[[p[2]]]) else glab(p[2]),
             Subgroup = if (p[1] == "AGE4SEX") ifelse(p[3] == "men", "Men", "Women") else AGE_EN[[p[3]]], n = big(x$n[1]),
             `Missing, 5 y` = miss_at(g, 5), `Missing, 10 y` = miss_at(g, 10), `Missing, 15 y` = miss_at(g, 15),
             `1 y` = cell(at(1)), `5 y` = cell(at(5)), `10 y` = cell(at(10)), `15 y` = cell(at(15)), check.names = FALSE) }))))
ytab <- function(y, first) cbind(first, data.frame(n = big(y$n), `Years alive` = yci(y$alive_est, y$alive_lcl, y$alive_ucl),
  `Years alive and independent` = yci(y$indep_est, y$indep_lcl, y$indep_ucl), `Years alive and dependent` = yci(y$dep_est, y$dep_lcl, y$dep_ucl), check.names = FALSE))
s6b <- ytab(yrs$by_age, data.frame(`Discharge mRS` = glab(yrs$by_age$group), `Age group` = AGE_EN[yrs$by_age$age], check.names = FALSE))

# ---- Former Table S7 (Zenodo Table Z4 in the paper): by sex (unadjusted and standardized to the overall G1 age distribution) ---------------------------------------------------------
ys <- yrs$by_sex; ps <- yrs$p15_sex_std
s7 <- do.call(rbind, lapply(c("All patients", "Discharge mRS 0-2", "Discharge mRS 3-5"), function(g) do.call(rbind, lapply(c("Men", "Women", "Men − Women"), function(s) {
  u <- ys[ys$group == g & ys$sex == s & ys$type == "Unadjusted", ]; a <- ys[ys$group == g & ys$sex == s & ys$type == "Age-standardized", ]
  p <- ps[ps$group == g & ps$sex == s & ps$outcome == "mRS 0-2", ]
  dif <- s == "Men − Women"; yy <- function(x, v) if (dif) sprintf("%+.1f (%.1f to %.1f)", x[[paste0(v, "_est")]], x[[paste0(v, "_lcl")]], x[[paste0(v, "_ucl")]]) else yci(x[[paste0(v, "_est")]], x[[paste0(v, "_lcl")]], x[[paste0(v, "_ucl")]])
  data.frame(`Discharge mRS` = glab(g), Sex = ifelse(dif, "Difference (men − women)", s), n = if (dif) "" else big(u$n),
             `Years alive (unadjusted)` = yy(u, "alive"), `Years independent (unadjusted)` = yy(u, "indep"),
             `Years alive (age-standardized)` = yy(a, "alive"), `Years independent (age-standardized)` = yy(a, "indep"),
             `Years dependent (age-standardized)` = yy(a, "dep"),
             `Independent at 15 y, % (age-standardized)` = if (dif) "–" else fmt_ci(p$est, p$lcl, p$ucl), check.names = FALSE) }))))

# ---- Former Table S8 (Table S4 in the paper): patients independent before stroke (prestroke mRS 0–2). A = proportions (estimates from 05_prestroke.R), B = years ---------------------------------
ep <- pre$est; p35 <- pre$est_pre35 %>% filter(H == H15)
pj_pre <- function(g, o) { v <- sapply(c("MAR", DELTA_M), function(m)
  (std_age_sex(ep, H15, H5, outcomes = o, within = if (g == "All") NULL else g, method = m) %>% filter(time == H15))$est)
  sprintf("%s [%s–%s]", fmt_pct(v[["MAR"]]), fmt_pct(min(v)), fmt_pct(max(v))) }
s8a <- rbind(
  do.call(rbind, lapply(names(OUTL), function(o) do.call(rbind, lapply(G9, function(g)
    data.frame(Population = "Prestroke mRS 0–2", Outcome = OUTL[[o]], `Discharge mRS` = glab(g), n = big(pick("MAR", g, o, 1, ep)$n),
               `1 y` = cell(pick("MAR", g, o, 1, ep)), `5 y` = cell(pick("MAR", g, o, 5, ep)), `10 y` = cell(pick("MAR", g, o, 10, ep)), `15 y` = cell(pick("MAR", g, o, 15, ep)),
               `15 y: range under missing-data assumptions` = rng15(g, o, ep),
               `15 y: applied to recent case mix [range]` = if (g %in% G3S) pj_pre(g, o) else "–", check.names = FALSE))))),
  do.call(rbind, lapply(names(OUTL), function(o) do.call(rbind, lapply(G3S, function(g) { x <- function(t) p35[p35$group == g & p35$outcome == o & p35$time == t, ]
    data.frame(Population = "Prestroke mRS 3–5 (comparison)", Outcome = OUTL[[o]], `Discharge mRS` = glab(g), n = big(x(1)$n),
               `1 y` = cell(x(1)), `5 y` = cell(x(5)), `10 y` = cell(x(10)), `15 y` = cell(x(15)),
               `15 y: range under missing-data assumptions` = "–", `15 y: applied to recent case mix [range]` = "–", check.names = FALSE) })))))
yp <- yrs$by_prestroke %>% filter(group %in% c("All patients", "Discharge mRS 0-2", "Discharge mRS 3-5"))
s8b <- ytab(yp, data.frame(`Prestroke mRS` = sub("-", "–", sub("Prestroke mRS ", "", yp$prestroke)), `Discharge mRS` = glab(yp$group), check.names = FALSE))

# ---- Former Table S9 (Table S5 in the paper): 15-year outcomes and years independent by missing-data handling (G1) ---------------------------------------------
ysc <- yrs$by_dc_scenario; ycc <- yrs$complete_case
s9 <- do.call(rbind, lapply(G3S, function(g) {
  km <- pick("KM_obs", g, "Alive", H15); yy <- function(m, v = "indep") { x <- ysc[ysc$method == m & ysc$group == ykey(g), ]; x[[paste0(v, "_est")]] }
  ym <- ysc[ysc$method == "MAR" & ysc$group == ykey(g), ]; yc <- ycc[ycc$group == ykey(g), ]
  rows <- list(c("Alive at 15 y, %", cell(pick("MAR", g, "Alive", H15)), cell(pick("CC", g, "Alive", H15)), cell(km), sapply(DELTA_M, function(m) fmt_pct(pick(m, g, "Alive", H15)$est))),
               c("Alive and independent at 15 y, %", cell(pick("MAR", g, "mRS 0-2", H15)), cell(pick("CC", g, "mRS 0-2", H15)), "–", sapply(DELTA_M, function(m) fmt_pct(pick(m, g, "mRS 0-2", H15)$est))),
               c("Years alive", yci(ym$alive_est, ym$alive_lcl, ym$alive_ucl), yci(yc$alive_est, yc$alive_lcl, yc$alive_ucl), "–", sapply(DELTA_M, function(m) f1(yy(m, "alive")))),
               c("Years alive and independent", yci(ym$indep_est, ym$indep_lcl, ym$indep_ucl), yci(yc$indep_est, yc$indep_lcl, yc$indep_ucl), "–", sapply(DELTA_M, function(m) f1(yy(m)))))
  x <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(x) <- c("Measure", "Multiple imputation (MAR; primary)", "Complete-case analysis", "Kaplan–Meier (observed data)", sprintf("δ = %s", sub("-", "−", sprintf("%+.1f", DELTAS_MNAR))))
  cbind(`Discharge mRS` = glab(g), x) }))
n_cc <- sapply(G3S, function(g) pick("CC", g, "Alive", H15)$n)

# ---- Former Table S10 (Table S6 in the paper): characteristics of the 3 onset-period cohorts ---------------------------------------------------------------
bg3 <- function(x) c(big(nrow(x)), med(x$age), npct(x$male == 1), med(x$pre_mrs), med(x$nihss), "", sapply(names(TOAST_EN), function(l) npct(x$toast == l)), "",
                     npct(x$ht == 1), npct(x$dm == 1), npct(x$dl == 1), npct(x$af == 1), npct(x$smk == 1), npct(x$drk == 1), npct(x$ckd == 1),
                     npct(x$prior_stroke == 1), med(x$bmi, 1), "", npct(x$tpa == 1), npct(x$evt == 1), "", sapply(0:5, function(m) npct(x$mrs_dc == m)))
s10 <- data.frame(Characteristic = c("Patients, n", "Age, y", "Male sex", "Prestroke mRS", "NIHSS score on admission", "Stroke subtype (TOAST)", paste0("  ", TOAST_EN),
                                     "Vascular risk factors and history", "  Hypertension", "  Diabetes", "  Dyslipidemia", "  Atrial fibrillation", "  Smoking",
                                     "  Alcohol drinking", "  Chronic kidney disease", "  Prior stroke", "Body mass index, kg/m²", "Acute treatment",
                                     "  Intravenous thrombolysis", "  Endovascular thrombectomy", "mRS at discharge", paste0("  mRS ", 0:5)), check.names = FALSE)
for (ch in chs) s10[[sprintf("%s (%d-year follow-up)", G_OF[[as.character(ch$H)]], ch$H)]] <- bg3(ch$data)
miss3 <- sapply(c("ckd", "prior_stroke", "bmi", "tpa", "evt"), function(v) paste(sapply(chs, function(ch) sum(is.na(ch$data[[v]]))), collapse = "/"))

# ---- Former Table S11 (Zenodo Table Z5 in the paper): differences in outcomes at common time points between onset-period cohorts (per-imputation differences pooled with Rubin's rules) -----------------
props <- list()
for (k in names(cohorts)) { ch <- cohorts[[k]]; b <- ch$data[, c("mrs_dc", "age", "male")]; rows <- group_rows(b, G9)
  Yl <- readRDS(file.path(WORK_DIR, sprintf("imputations_%s.rds", k)))$imps[["delta_+0.00"]]
  for (t in intersect(c(5, 10), ch$tp$year)) for (g in G9) for (o in names(OUTL)) { j <- match(t, ch$tp$year)
    props[[length(props) + 1]] <- list(H = ch$H, t = t, g = g, o = o, n = length(rows[[g]]), p = vapply(Yl, function(Y) mean(outcome_indicator(Y[rows[[g]], j], o)), numeric(1))) }
  rm(Yl); gc() }
getp <- function(H, t, g, o) Filter(function(x) x$H == H && x$t == t && x$g == g && x$o == o, props)[[1]]
diff_rubin <- function(a, b) {               # a − b (2 independently imputed cohorts; within-imputation variance = sum of the 2 binomial variances)
  if (a$n == 0 || b$n == 0) return(c(diff = NA_real_, lcl = NA_real_, ucl = NA_real_))           # group with 0 patients
  m <- length(a$p); dd <- a$p - b$p; U <- a$p * (1 - a$p) / a$n + b$p * (1 - b$p) / b$n
  Tt <- mean(U) + (1 + 1 / m) * var(dd)
  if (Tt <= 0) {                                   # proportion 0 or 1 in every imputed dataset of both cohorts (variance 0): no interval
    warning("diff_rubin: total variance is 0; confidence interval not computed"); return(c(diff = mean(dd), lcl = NA_real_, ucl = NA_real_)) }
  lam <- (1 + 1 / m) * var(dd) / Tt; df <- max(br_df(lam, m, a$n + b$n - 1), 1)
  c(diff = mean(dd), lcl = mean(dd) - qt(0.975, df) * sqrt(Tt), ucl = mean(dd) + qt(0.975, df) * sqrt(Tt)) }
CMP <- list("G2 − G1, 5 y" = c(10, 5), "G3 − G1, 5 y" = c(5, 5), "G2 − G1, 10 y" = c(10, 10))
cdiff <- do.call(rbind, lapply(names(OUTL), function(o) do.call(rbind, lapply(G9, function(g) do.call(rbind, lapply(names(CMP), function(cn) {
  r <- diff_rubin(getp(CMP[[cn]][1], CMP[[cn]][2], g, o), getp(H15, CMP[[cn]][2], g, o)); data.frame(outcome = o, group = g, comparison = cn, t(r)) }))))))
s11 <- cdiff %>% mutate(Outcome = OUTL[outcome], `Discharge mRS` = glab(group), val = sprintf("%+.1f (%.1f to %.1f)", 100 * diff, 100 * lcl, 100 * ucl)) %>%
  select(Outcome, `Discharge mRS`, comparison, val) %>% tidyr::pivot_wider(names_from = comparison, values_from = val) %>% as.data.frame(check.names = FALSE)
minus <- function(x) gsub("-", "−", x, fixed = TRUE)
s11[, -(1:2)] <- lapply(s11[, -(1:2)], minus); s7[] <- lapply(s7, function(x) ifelse(grepl("^[+-]", x), minus(x), x))

# =================================================================================================
#  Figures
# =================================================================================================
fig <- function(name) file.path(pd, paste0(name, ".png"))

# ---- Figure 1: mRS distribution at each assessment in G1 (all patients top center; discharge mRS 0–5 in 2 columns × 3 rows) -------------------------
# Values are not printed inside bars (category values are in Table Z1)
one_panel <- function(g, show_legend = FALSE) { p <- plot_trajectory(e15, grp = g, ncol = 1, values = FALSE) + labs(x = NULL, y = NULL)
  if (!show_legend) p <- p + theme(legend.position = "none"); p }
top <- patchwork::plot_spacer() + one_panel("All", TRUE) + patchwork::plot_spacer() + patchwork::plot_layout(widths = c(0.5, 1, 0.5))
pairs <- lapply(list(0:1, 2:3, 4:5), function(k) one_panel(paste0("Discharge mRS ", k[1])) + one_panel(paste0("Discharge mRS ", k[2])))
p1 <- patchwork::wrap_elements(full = top) / pairs[[1]] / pairs[[2]] / pairs[[3]] + patchwork::plot_layout(heights = c(1.15, 1, 1, 1))
ylab <- grid::textGrob("Patients (%)", rot = 90, gp = grid::gpar(fontsize = 8, fontfamily = FIG_FONT))
p1 <- patchwork::wrap_elements(full = ylab) + patchwork::wrap_elements(full = p1) + patchwork::plot_layout(widths = c(0.03, 1))
cap <- grid::textGrob("Time after stroke onset (years unless indicated)", gp = grid::gpar(fontsize = 8, fontfamily = FIG_FONT))
save_fig(p1 / patchwork::wrap_elements(full = cap) + patchwork::plot_layout(heights = c(1, 0.025)), fig("Figure1_mRS_distribution_by_discharge_mRS"), FIG_W, 8.6)

# ---- Figures 2 and 4: how the 15 years were spent (horizontal bars of years alive and independent, alive and dependent, and after death) -------------------------------
YCOL <- c("Lived independently (mRS 0–2)" = "#1f4e79", "Lived requiring assistance (mRS 3–5)" = "#e8833a", "After death" = "#d9d9d9")
years_bars <- function(y, labels, tag = NULL, show_x = TRUE) {
  y$lab <- factor(labels, levels = rev(labels))
  long <- rbind(data.frame(lab = y$lab, state = names(YCOL)[1], yr = y$indep_est), data.frame(lab = y$lab, state = names(YCOL)[2], yr = y$dep_est),
                data.frame(lab = y$lab, state = names(YCOL)[3], yr = H15 - y$alive_est))
  long$state <- factor(long$state, levels = rev(names(YCOL)))
  ggplot(long, aes(yr, lab, fill = state)) + geom_col(width = 0.72) +
    geom_text(data = y, aes(x = indep_est / 2, y = lab, label = ifelse(indep_est >= 0.9, f1(indep_est), "")), inherit.aes = FALSE, colour = "white", size = 2.5, family = FIG_FONT) +
    geom_text(data = y, aes(x = indep_est + dep_est / 2, y = lab, label = ifelse(dep_est >= 0.9, f1(dep_est), "")), inherit.aes = FALSE, colour = "black", size = 2.5, family = FIG_FONT) +
    geom_text(data = y, aes(x = H15 + 0.2, y = lab, label = f1(alive_est)), inherit.aes = FALSE, hjust = 0, size = 2.5, family = FIG_FONT) +
    scale_fill_manual(values = YCOL, breaks = names(YCOL), name = NULL) +
    scale_x_continuous(breaks = seq(0, 15, 3), limits = c(0, H15 + 1.3), expand = c(0, 0)) +
    labs(x = if (show_x) "Years during the 15 years after stroke onset" else NULL, y = NULL, tag = tag) + theme_mrs() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
}
yd <- yrs$by_dc[match(c("All patients", G6), yrs$by_dc$group), ]
save_fig(years_bars(yd, sprintf("%s (n = %s)", c("All patients", G6), big(yd$n))), fig("Figure2_years_by_discharge_mRS"), FIG_W, 3.4)
ya <- function(g) { y <- yrs$by_age[yrs$by_age$group == g, ]; years_bars(y, sprintf("%s (n = %s)", AGE_EN[y$age], big(y$n)), tag = NULL, show_x = g == "Discharge mRS 3-5") }
ttl <- function(p, t) p + labs(tag = t)
p4 <- (ttl(ya("All patients"), "A") / ttl(ya("Discharge mRS 0-2"), "B") / ttl(ya("Discharge mRS 3-5"), "C")) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "top")
save_fig(p4, fig("Figure4_years_by_age_group"), FIG_W, 5.6)

# ---- Figure 3: survival (A) and independent survival (B) over time for each discharge mRS group 0–5 ---------------------------------------
DM_COLS <- setNames(unname(MRS_COLS[as.character(0:5)]), G6); DM_SHAPES <- setNames(c(16, 15, 17, 25, 18, 3), G6)
d3f <- e15 %>% filter(group %in% G6, outcome %in% c("Alive", "mRS 0-2")) %>% group_by(group) %>% mutate(lab = sprintf("%s (n = %s)", group, big(n[1]))) %>% ungroup()
lv <- (d3f %>% distinct(group, lab))$lab[match(G6, (d3f %>% distinct(group, lab))$group)]; d3f$lab <- factor(d3f$lab, levels = lv)
panel2 <- function(o, yl, tag, sx) ggplot(d3f %>% filter(outcome == o), aes(time, est, colour = lab, shape = lab, fill = lab)) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.3) +
  scale_colour_manual(values = setNames(DM_COLS, lv), name = NULL) + scale_fill_manual(values = setNames(DM_COLS, lv), name = NULL) + scale_shape_manual(values = setNames(DM_SHAPES, lv), name = NULL) +
  scale_x_continuous(breaks = seq(0, 15, 3), limits = c(0, 15.2), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(labels = function(x) 100 * x, limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.02)) +
  labs(x = if (sx) "Years after stroke onset" else NULL, y = yl, tag = tag) + guides(colour = guide_legend(nrow = 2, byrow = FALSE)) +
  theme_mrs() + theme(panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3))
save_fig((panel2("Alive", "Alive (%)", "A", FALSE) / panel2("mRS 0-2", "Alive with mRS 0–2 (%)", "B", TRUE)) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "top"),
         fig("Figure3_survival_and_independent_survival_by_discharge_mRS"), FIG_W, 6.6)

# ---- Figure S1: selection of analyzed patients and composition of the 3 cohorts (flowchart) -------------------------------------------
save_fig(plot_flow(co), fig("FigureS1_flowchart"), FIG_W, 4.6)

# ---- Figure S2: survival in G1 vs expected survival in the general population (left = all patients, right = discharge mRS 0–2 / 3–5; relative survival below) ----------
#   Survival = multiple imputation (MAR) estimates with 95% CI; expected survival = 06_years_lived.R (Ederer I method). Plotted values are rounded to 1 decimal place in %
GL <- c("All" = "All patients", "Discharge mRS 0-2" = "Discharge mRS 0–2", "Discharge mRS 3-5" = "Discharge mRS 3–5")
nG <- sapply(G3S, function(g) big(pick("MAR", g, "Alive", 1)$n))
r1 <- function(x) as.numeric(sprintf("%.1f", 100 * x)) / 100
lt <- e15 %>% filter(group %in% G3S, outcome == "Alive") %>% select(group, time, est, lcl, ucl) %>%
  mutate(across(c(est, lcl, ucl), ~ ifelse(time == 0, 1, .x))) %>%                    # All patients alive at discharge
  inner_join(yrs$expected_curve %>% mutate(group = ifelse(group == "All patients", "All", group)), by = c("group", "time")) %>%
  arrange(match(group, G3S), time) %>%
  mutate(rel = round(100 * est / expected, 1), obs = round(100 * est, 1) / 100, exp = round(100 * expected, 1) / 100, lcl = r1(lcl), ucl = r1(ucl),
         glab = factor(GL[group], levels = GL),
         panel = ifelse(group == "All", sprintf("All patients (n = %s)", nG[1]), sprintf("Discharge mRS 0–2 (n = %s) and 3–5 (n = %s)", nG[2], nG[3])))
GC <- setNames(c("#2b2b2a", "#2a78d6", "#eb6834"), GL)
LTY <- c("This study: multiple imputation (MAR)" = "solid", "General population (life tables, matched)" = "22")
stepd <- function(x) { r <- step_ribbon(as.data.frame(x %>% select(glab, panel, time, lcl, ucl)), "glab"); r$glab <- factor(r$glab, levels = levels(x$glab)); r }
sp <- function(x) ggplot(x, aes(time, colour = glab, fill = glab)) +
  geom_ribbon(data = stepd(x), aes(ymin = lcl, ymax = ucl), colour = NA, alpha = 0.15, outline.type = "full") +
  geom_step(aes(y = obs, linetype = names(LTY)[1]), linewidth = 0.55) + geom_step(aes(y = exp, linetype = names(LTY)[2]), linewidth = 0.5) + facet_wrap(~panel) +
  scale_colour_manual(values = GC, name = NULL) + scale_fill_manual(values = GC, guide = "none") + scale_linetype_manual(values = LTY, name = NULL) +
  scale_x_continuous(breaks = seq(0, 15, 3), expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(labels = function(x) 100 * x, limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.01)) +
  labs(x = "Years after stroke onset", y = "Survival (%)") + theme_mrs() + theme(panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3))
rs <- lt %>% filter(time %in% c(5, 10, 15)) %>% transmute(glab, time, lab = sprintf("%.1f%%", rel))
tb <- ggplot(rs, aes(factor(time), glab, label = lab)) + geom_text(size = 2.6, family = FIG_FONT) +
  scale_x_discrete(labels = function(x) paste(x, "y")) + scale_y_discrete(limits = rev(levels(rs$glab))) + labs(x = NULL, y = NULL, title = "Relative survival (observed / expected)") +
  theme_mrs() + theme(axis.line = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 8, face = "bold", family = FIG_FONT))
pS2 <- (((sp(lt %>% filter(group == "All")) + guides(colour = "none")) | (sp(lt %>% filter(group != "All")) + labs(y = NULL))) + patchwork::plot_layout(guides = "collect") &
          guides(linetype = guide_legend(nrow = 1, order = 1)) & theme(legend.position = "top", legend.box = "vertical", legend.box.just = "left")) /
  (patchwork::plot_spacer() | tb | patchwork::plot_spacer()) + patchwork::plot_layout(heights = c(1, 0.28))
pS2[[2]] <- pS2[[2]] + patchwork::plot_layout(widths = c(0.25, 1, 0.25))
save_fig(pS2, fig("FigureS2_survival_vs_general_population"), FIG_W, 4.6)

# ---- Figure S3: survival and independent survival over time by discharge mRS 0–5 × 4 age groups --------------------------------------------
OUT2 <- c("Alive" = "Alive (%)", "mRS 0-2" = "Alive with mRS 0–2 (%)")
dS3 <- est_sa %>% filter(startsWith(group, "AGE4 | ")) %>% tidyr::separate(group, c("type", "dc", "lev"), sep = " \\| ") %>% filter(dc %in% G6) %>%
  mutate(dc = factor(dc, levels = G6), lev = factor(lev, levels = AGE_LABELS, labels = c("<65 years", "65–74 years", "75–84 years", "≥85 years")), outcome = factor(OUT2[outcome], levels = OUT2))
pS3 <- ggplot(dS3, aes(time, est, colour = lev)) + geom_line(linewidth = 0.55) + facet_grid(dc ~ outcome, switch = "y") +
  scale_colour_manual(values = setNames(c("#1b8a6b", "#2a78d6", "#e39b2d", "#9b2b4b"), levels(dS3$lev)), name = NULL) +
  scale_x_continuous(breaks = seq(0, 15, 3), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(labels = function(x) 100 * x, limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0.02)) + labs(x = "Years after stroke onset", y = NULL) +
  theme_mrs() + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90, hjust = 0.5), panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3))
save_fig(pS3, fig("FigureS3_by_discharge_mRS_and_age_group"), FIG_W, 8.8)

# ---- Figure S4: survival (A) and independent survival (B) by discharge mRS 0–5 in the 3 cohorts (line type = cohort) ------------------------------
el <- est_all %>% filter(method == "MAR", group %in% G6, outcome %in% c("Alive", "mRS 0-2")) %>% mutate(G = factor(G_OF[as.character(H)], levels = c("G1", "G2", "G3")))
pan <- function(o, yl, tag, sx) ggplot(el %>% filter(outcome == o), aes(time, est, colour = group, linetype = G)) + geom_line(linewidth = 0.5) +
  scale_colour_manual(values = DM_COLS, name = NULL) + scale_linetype_manual(values = c(G1 = "solid", G2 = "22", G3 = "12"), name = NULL) + scale_x_continuous(breaks = seq(0, 15, 3)) +
  scale_y_continuous(labels = function(x) 100 * x, limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0.02)) +
  labs(x = if (sx) "Years after stroke onset" else NULL, y = yl, tag = tag) + guides(colour = guide_legend(nrow = 2)) + theme_mrs() + theme(panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3))
save_fig((pan("Alive", "Alive (%)", "A", FALSE) / pan("mRS 0-2", "Alive with mRS 0–2 (%)", "B", TRUE)) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "top"),
         fig("FigureS4_three_cohorts_by_discharge_mRS"), FIG_W, 6.8)

# =================================================================================================
#  Values cited in the text (including some not shown directly in tables or figures)
# =================================================================================================
p1f <- function(x) sprintf("%.1f", 100 * x); e <- function(g, o, t) pick("MAR", g, o, t)$est
yv <- function(g, v, tab = yrs$by_dc) f1(tab[tab$group == g, paste0(v, "_est")])
yav <- function(g, a, v) f1(yrs$by_age[yrs$by_age$group == g & yrs$by_age$age == a, paste0(v, "_est")])
ysv <- function(g, s, type, v) { x <- ys[ys$group == g & ys$sex == s & ys$type == type, ]; c(est = minus(f1(x[[paste0(v, "_est")]])), ci = minus(sprintf("%s to %s", f1(x[[paste0(v, "_lcl")]]), f1(x[[paste0(v, "_ucl")]])))) }
lost <- !complete.cases(Yobs); ex <- yrs$expected; exv <- function(g, v) f1(ex[ex$group == g, v])
r16 <- yrs$recent %>% filter(grepl("16 strata", strata)); r8 <- function(g, v, w = "est") f1(rec[rec$group == g & rec$value == v, w])
cd <- cdiff %>% filter(group == "All"); stp <- yrs$by_dc_step
V <- list(
  n_registered = big(co$flow$n[1]), n_inhosp = big(co$flow$n[2]), n_total = big(co$flow$n[4]),
  n_G1 = big(nrow(d)), n_G2 = big(nrow(chs[[2]]$data)), n_G3 = big(nrow(chs[[3]]$data)), n_G23 = big(nrow(chs[[2]]$data) + nrow(chs[[3]]$data)),
  age_med = med(d$age), male_pct = p1f(mean(d$male == 1)), nihss_med = med(d$nihss), dc02_pct = p1f(mean(d$mrs_dc <= 2)), dc35_pct = p1f(mean(d$mrs_dc >= 3)),
  dc35_age75_pct = p1f(mean(d$age[d$mrs_dc >= 3] >= 75)), dc35_pre35_pct = p1f(mean(d$pre_mrs[d$mrs_dc >= 3] >= 3)),
  women_age75_pct = p1f(mean(d$age[d$male == 0] >= 75)), men_age75_pct = p1f(mean(d$age[d$male == 1] >= 75)),
  n_lost = big(sum(lost)), lost_pct = p1f(mean(lost)), n_complete = big(sum(!lost)),
  lost_age_med = sprintf("%.0f", median(d$age[lost])), comp_age_med = sprintf("%.0f", median(d$age[!lost])),
  lost_dc_med = sprintf("%.0f", median(d$mrs_dc[lost])), comp_dc_med = sprintf("%.0f", median(d$mrs_dc[!lost])),
  ckd_miss_pct = p1f(mean(is.na(d$ckd))),
  s15_all = cell(pick("MAR", "All", "Alive", 15)), i15_all = cell(pick("MAR", "All", "mRS 0-2", 15)),
  indep_among_survivors = paste(sapply(c(5, 10, 15), function(t) p1f(e("All", "mRS 0-2", t) / e("All", "Alive", t))), collapse = " / "),
  dc02_dead15 = p1f(1 - e("Discharge mRS 0-2", "Alive", 15)), dc02_dep15 = p1f(e("Discharge mRS 0-2", "Alive", 15) - e("Discharge mRS 0-2", "mRS 0-2", 15)),
  dc02_indep = paste(sapply(c(1, 5, 10, 15), function(t) p1f(e("Discharge mRS 0-2", "mRS 0-2", t))), collapse = " / "),
  dc2_indep5 = p1f(e("Discharge mRS 2", "mRS 0-2", 5)),
  dc3_indep = paste(sapply(c(0.25, 1, 5, 15), function(t) p1f(e("Discharge mRS 3", "mRS 0-2", t))), collapse = " / "),
  dc4_indep_1_2y = paste(sapply(c(1, 2), function(t) p1f(e("Discharge mRS 4", "mRS 0-2", t))), collapse = " / "),
  dc5_indep_max = p1f(max(e15$est[e15$group == "Discharge mRS 5" & e15$outcome == "mRS 0-2" & e15$time > 0])),
  dep_share_dc3 = p1f(yrs$by_dc$dep_est[yrs$by_dc$group == "Discharge mRS 3"] / yrs$by_dc$alive_est[yrs$by_dc$group == "Discharge mRS 3"]),
  dep_share_dc4 = p1f(yrs$by_dc$dep_est[yrs$by_dc$group == "Discharge mRS 4"] / yrs$by_dc$alive_est[yrs$by_dc$group == "Discharge mRS 4"]),
  exp_years = paste(sapply(c("All patients", "Discharge mRS 0-2", "Discharge mRS 3-5"), function(g) exv(g, "expected_years")), collapse = " / "),
  exp_diff = paste(sapply(c("All patients", "Discharge mRS 0-2", "Discharge mRS 3-5"), function(g) minus(exv(g, "difference"))), collapse = " / "),
  rel_surv_dc02 = paste(sapply(c(5, 10, 15), function(t) sprintf("%.1f", lt$rel[lt$group == "Discharge mRS 0-2" & lt$time == t])), collapse = " / "),
  km15_all = cell(pick("KM_obs", "All", "Alive", 15)), cc15_all = cell(pick("CC", "All", "Alive", 15)),
  km_minus_mar = sprintf("%.1f", 100 * (pick("KM_obs", "All", "Alive", 15)$est - e("All", "Alive", 15))),
  i15_delta_range = rng15("All", "mRS 0-2"),
  yrs_indep_delta_range = { v <- ysc$indep_est[ysc$group == "All patients"]; sprintf("%s–%s", f1(min(v)), f1(max(v))) },
  yrs_step = sprintf("alive %s, independent %s, dependent %s (trapezoid: %s, %s, %s)", f1(stp$alive_est[1]), f1(stp$indep_est[1]), f1(stp$dep_est[1]), yv("All patients", "alive"), yv("All patients", "indep"), yv("All patients", "dep")),
  yrs_step_maxdiff = sprintf("%.2f", max(abs(as.matrix(stp[, paste0(VARS, "_est")]) - as.matrix(yrs$by_dc[, paste0(VARS, "_est")])))),
  cohort_diff_all_range = sprintf("%s to %s", minus(sprintf("%.1f", 100 * min(cd$diff))), minus(sprintf("%.1f", 100 * max(cd$diff)))),
  G2_dc35_surv10_diff = minus(sprintf("%.1f", 100 * cdiff$diff[cdiff$group == "Discharge mRS 3-5" & cdiff$outcome == "Alive" & cdiff$comparison == "G2 − G1, 10 y"])),
  recent16 = sprintf("years independent %s [%s–%s]; independent at 15 y %s%% [%s–%s]", f1(r16$est[r16$value == "indep"]), f1(r16$range_min[r16$value == "indep"]), f1(r16$range_max[r16$value == "indep"]),
                     f1(r16$est[r16$value == "p15_indep"]), f1(r16$range_min[r16$value == "p15_indep"]), f1(r16$range_max[r16$value == "p15_indep"])),
  recent_alive_years = r8("All patients", "alive"), recent_dep_years = r8("All patients", "dep"),
  G3_age75_pct = p1f(mean(chs[[3]]$data$age >= 75)), G1_age75_pct = p1f(mean(d$age >= 75)), G3_dc35_pct = p1f(mean(chs[[3]]$data$mrs_dc >= 3)),
  pre02_n = big(sum(d$pre_mrs <= 2)), pre35_dc02_n = big(sum(d$pre_mrs >= 3 & d$mrs_dc <= 2)),
  t1_missing = paste(names(t1_miss), t1_miss, collapse = "; "), s10_missing = paste(names(miss3), miss3, collapse = "; "),
  n_cc = paste(big(n_cc), collapse = " / "),
  w_age_G1 = paste(sprintf("%s %.1f%%", names(yrs$weights$age_G1), 100 * yrs$weights$age_G1), collapse = "; "),
  n_smallest_stratum = sprintf("%s", min(yrs$by_age$n)))
for (g in c("All patients", "Discharge mRS 0-2", "Discharge mRS 3-5", G6)) for (v in VARS) V[[paste("yrs", gsub("[^0-9A-Za-z]+", "", sub("Discharge mRS ", "dc", sub("All patients", "all", g))), v, sep = "_")]] <- yv(g, v)
for (g in c("All patients", "Discharge mRS 0-2", "Discharge mRS 3-5", G6)) for (a in AGE_LABELS) for (v in VARS)
  V[[paste("age", gsub("[^0-9A-Za-z]+", "", sub("Discharge mRS ", "dc", sub("All patients", "all", g))), gsub("[^0-9]+", "", a), v, sep = "_")]] <- yav(g, a, v)
for (g in c("All patients", "Discharge mRS 0-2", "Discharge mRS 3-5")) for (s in c("Men", "Women", "Men − Women")) for (tp in c("Unadjusted", "Age-standardized")) {
  r <- ysv(g, s, tp, "indep"); kk <- paste("sex", gsub("[^0-9A-Za-z]+", "", sub("Discharge mRS ", "dc", sub("All patients", "all", g))), c(Men = "men", Women = "women", "Men − Women" = "diff")[[s]], ifelse(tp == "Unadjusted", "crude", "std"), sep = "_")
  V[[paste0(kk, "_indep")]] <- r[["est"]]; V[[paste0(kk, "_indep_ci")]] <- r[["ci"]] }
for (g in c("All patients", "Discharge mRS 0-2", "Discharge mRS 3-5")) for (v in c("indep", "p15_indep", "p15_alive", "alive")) for (w in c("est", "range_min", "range_max"))
  V[[paste("recent", gsub("[^0-9A-Za-z]+", "", sub("Discharge mRS ", "dc", sub("All patients", "all", g))), v, w, sep = "_")]] <- r8(g, v, w)
for (g in G9) for (o in c("Alive", "mRS 0-2")) for (t in c(1, 5, 10, 15))
  V[[paste(ifelse(o == "Alive", "s", "i"), gsub("[^0-9A-Za-z]+", "", sub("Discharge mRS ", "dc", sub("All", "all", g))), t, sep = "_")]] <- p1f(e(g, o, t))
yp1 <- function(p, g, v) f1(yrs$by_prestroke[yrs$by_prestroke$prestroke == p & yrs$by_prestroke$group == g, paste0(v, "_est")])
for (v in VARS) { V[[paste0("pre02_all_", v)]] <- yp1("Prestroke mRS 0-2", "All patients", v); V[[paste0("pre02_dc35_", v)]] <- yp1("Prestroke mRS 0-2", "Discharge mRS 3-5", v)
  V[[paste0("pre02_dc02_", v)]] <- yp1("Prestroke mRS 0-2", "Discharge mRS 0-2", v) }
ci_txt <- function(l, u) sprintf("%s–%s", f1(pmax(l, 0)), f1(u))
for (g in c("All patients", "Discharge mRS 0-2", "Discharge mRS 3-5", G6)) for (v in VARS) { x <- yrs$by_dc[yrs$by_dc$group == g, ]
  V[[paste("yrs", gsub("[^0-9A-Za-z]+", "", sub("Discharge mRS ", "dc", sub("All patients", "all", g))), v, "ci", sep = "_")]] <- ci_txt(x[[paste0(v, "_lcl")]], x[[paste0(v, "_ucl")]]) }
V$exp_all <- exv("All patients", "expected_years"); V$exp_dc02 <- exv("Discharge mRS 0-2", "expected_years"); V$exp_dc35 <- exv("Discharge mRS 3-5", "expected_years")
V$expdiff_all <- f1(abs(ex$difference[ex$group == "All patients"]))
sv <- sapply(c(5, 10, 15), function(t) 100 * e("All", "mRS 0-2", t) / e("All", "Alive", t)); V$surv_indep_range <- sprintf("%.0f–%.0f", min(sv), max(sv))
V$i_dc3_3m <- p1f(e("Discharge mRS 3", "mRS 0-2", 0.25)); V$i_dc4_2 <- p1f(e("Discharge mRS 4", "mRS 0-2", 2))
V$km15_all_pt <- p1f(pick("KM_obs", "All", "Alive", 15)$est); V$cc15_all_pt <- p1f(pick("CC", "All", "Alive", 15)$est)
for (m in DELTA_M) V[[paste0("yrs_all_indep_", gsub("[^0-9]", "", m), ifelse(grepl("-", m), "m", "p"))]] <- f1(ysc$indep_est[ysc$method == m & ysc$group == "All patients"])
V$step_all <- paste(f1(stp$alive_est[1]), f1(stp$indep_est[1]), f1(stp$dep_est[1]), sep = " / ")
for (o in c("Alive", "mRS 0-2")) { x <- pick("MAR", "All", o, 15); V[[paste0(ifelse(o == "Alive", "s", "i"), "_all_15_ci")]] <- sprintf("%s–%s", fmt_pct(x$lcl), fmt_pct(x$ucl)) }
for (v in names(t1_miss)) V[[paste0("miss_G1_", v)]] <- as.character(t1_miss[[v]])
for (v in names(miss3)) V[[paste0("miss3_", v)]] <- miss3[[v]]
V$step_all_alive <- f1(stp$alive_est[1]); V$step_all_indep <- f1(stp$indep_est[1]); V$step_all_dep <- f1(stp$dep_est[1])
V$recent16_indep_est <- f1(r16$est[r16$value == "indep"]); V$recent16_p15_indep_est <- f1(r16$est[r16$value == "p15_indep"])
V$n_cc_all <- big(n_cc[[1]]); V$n_cc_dc02 <- big(n_cc[[2]]); V$n_cc_dc35 <- big(n_cc[[3]])
txt_tab <- data.frame(Key = names(V), Value = unlist(lapply(V, as.character)), row.names = NULL, check.names = FALSE)

# =================================================================================================
#  Save
# =================================================================================================
tables <- list("Table 1" = t1, "Table 2" = t2, "Table S2" = s2, "Table S3" = s3, "Table S4" = s4, "Table S5" = s5, "Table S6A" = s6a, "Table S6B" = s6b,
               "Table S7" = s7, "Table S8A" = s8a, "Table S8B" = s8b, "Table S9" = s9, "Table S10" = s10, "Table S11" = s11)
# xlsx sheet names use the manuscript numbers (Zenodo Tables Z1–Z5 are built from these objects separately and are not written here)
PAPER_NO <- c("Table 1" = "Table 1", "Table 2" = "Table 2", "Table S2" = "Table S2", "Table S3" = "Table S3",
              "Table S8A" = "Table S4A", "Table S8B" = "Table S4B", "Table S9" = "Table S5", "Table S10" = "Table S6")
paper_tables <- setNames(tables[names(PAPER_NO)], PAPER_NO)
write_xlsx_sheets(c(paper_tables, list("Text values" = txt_tab)), file.path(pd, "Tables_main_and_supplement.xlsx"),
                  notes = c("Manuscript tables (sheet name = table number) and values cited in the text (Text values). Created by 09_paper_items.R from the estimates in work/. Proportions are % (95% CI); years are mean years during the 15 years after stroke onset (95% CI).",
                            "Table 1: G1 baseline characteristics (values are n (%) or median (interquartile range); denominators are non-missing patients for each variable). Table 2: years (06_years_lived.R), 15-year proportions (04_estimates.R), applied to the recent (G3) case mix (point estimate [range over δ = −1 to +1]).",
                            "Table S4A uses estimates from 05_prestroke.R. δ columns in Table S5 are point estimates. Table S1 (imputation-model predictors) is written directly in the Supplemental Material. Zenodo Tables Z1–Z5 correspond to the former Tables S4, S5, S6, S7, and S11 (all tables, under their former numbers, are in paper_values.json; Z3A and Z5 are aggregated here from the imputed data (MAR)).",
                            "– = not calculated or outside the observation period."))
jsonlite::write_json(list(values = V, tables = lapply(tables, function(x) list(columns = names(x), rows = unname(lapply(seq_len(nrow(x)), function(i) unname(as.character(unlist(x[i, ]))))))),
                          t1_missing = as.list(t1_miss), s10_missing = as.list(miss3)),
                     file.path(pd, "paper_values.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(c("Index of manuscript figures/tables and output files (created by 09_paper_items.R)", "",
             "Tables 1, 2, S2–S6; values cited in text  : Tables_main_and_supplement.xlsx (sheet name = table number; values cited in text = Text values)",
             "Tables Z1–Z5 (Zenodo)                     : former Tables S4, S5, S6, S7, S11 (in paper_values.json under their former numbers)",
             "Table S1 (imputation-model predictors)    : written directly in the Supplemental Material (missing counts as in the Table S6 footnote)",
             "Figure 1  mRS distribution at each assessment          : Figure1_mRS_distribution_by_discharge_mRS.png / .pdf",
             "Figure 2  15 years by discharge mRS                    : Figure2_years_by_discharge_mRS.png / .pdf",
             "Figure 3  Survival and independent survival over time  : Figure3_survival_and_independent_survival_by_discharge_mRS.png / .pdf",
             "Figure 4  15 years by age group                        : Figure4_years_by_age_group.png / .pdf",
             "Figure S1 Selection of analyzed patients               : FigureS1_flowchart.png / .pdf",
             "Figure S2 Comparison with general-population expected survival : FigureS2_survival_vs_general_population.png / .pdf",
             "Figure S3 Trajectories by discharge mRS × age group    : FigureS3_by_discharge_mRS_and_age_group.png / .pdf",
             "Figure S4 Comparison of the 3 cohorts                  : FigureS4_three_cohorts_by_discharge_mRS.png / .pdf", "",
             "Estimates at each time point, all cohorts and methods  : ../estimates_all_cohorts_all_methods_long.csv (04_estimates.R)",
             "Supplemental Methods values (fit, cross-validation, convergence) : ../90_diagnostics/<cohort>/TableD2_imputed_vs_observed_mRS_MAR.xlsx, TableD3_imputation_model_fit_MAR.xlsx (04_estimates.R); FigureD1_MICE_convergence_trace.png (03_imputation.R)",
             "All estimates of years                                 : ../01_years_over_15y/TableY_years_alive_and_independent_over_15y_G1.xlsx (06_years_lived.R)",
             "Supplemental Methods values (imputation model comparison, delta = 0 regeneration) : ../90_diagnostics/TableD5_delta0_regeneration_and_model_comparison.xlsx (07_model_checks.R)",
             "Supplemental Methods values (probability of threshold order reversal)             : ../90_diagnostics/TableD6_threshold_order_reversal.xlsx (08_threshold_check.R)",
             "paper_values.json                                      : the tables above and values cited in the text (values used in the manuscript text and tables, JSON; table keys keep the former numbers S2–S11)"),
           file.path(pd, "00_index.txt"), useBytes = FALSE)
