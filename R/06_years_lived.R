#-----------------------------------------------------------------------------------
#  06_years_lived.R  —  Years alive and years lived independent over the 15 years after stroke onset (primary outcome of the paper)
#
#   Target: G1 (15-year cohort). Only summarizes the existing imputed data (work/imputations_15y.rds); imputation is not redone.
#   - Years: for each patient, the 0/1 indicators of "alive (mRS 0–5)" and "alive and independent (mRS 0–2)" at the assessment time points
#           from discharge (year 0) to 15 years (0, 0.25, 0.5, 1, …, 15 years) are integrated by the trapezoidal rule.
#           Group mean = area under the survival / independent-survival curve (restricted mean time). Years lived dependent = difference of the two.
#   - Pooling: group mean and variance (variance of patient-level values / n) per imputed dataset, pooled by Rubin's rules (Barnard–Rubin degrees of freedom).
#   - Standardization (direct method): stratum means averaged with fixed weights. Variance = Σ w² × stratum variance / stratum size (complete-data df = ∞).
#       Sex comparison = age (4 groups) distribution of all G1; application to the recent case mix = age (4 groups) × sex (8 strata) distribution of G3.
#   - Sensitivity analyses: delta-adjusted pattern-mixture model (δ = −1.0, −0.5, +0.5, +1.0), complete-case analysis, years by step function.
#   - Expected survival of the general population (Ederer I method; computed per patient from sex-, age-, and calendar-year-specific death probabilities of the abridged life tables, then averaged) and
#     expected years alive obtained by integrating it over the same assessment time points (Figure S3 of the paper).
#   Output: work/estimates_years.rds, output/01_years_over_15y/
#-----------------------------------------------------------------------------------

yd <- dir_create(OUT_DIR, "01_years_over_15y")
unlink(list.files(yd, full.names = TRUE))
cohorts <- readRDS(file.path(WORK_DIR, "cohorts.rds"))$cohorts
est_all <- readRDS(file.path(WORK_DIR, "estimates_all.rds"))
H15 <- max(HORIZONS); H5 <- min(HORIZONS)
k15 <- names(cohorts)[sapply(cohorts, `[[`, "H") == H15]
k5  <- names(cohorts)[sapply(cohorts, `[[`, "H") == H5]
coh <- cohorts[[k15]]; d <- coh$data; d3 <- cohorts[[k5]]$data
Yobs <- as.matrix(d[, coh$tp$col])
TT <- c(0, coh$tp$year)                                   # assessment time points (years; 0 = discharge)
stopifnot(max(TT) == H15)

SCEN <- c("MAR" = "delta_+0.00", setNames(mnar_key(DELTAS_MNAR), sprintf("MNAR_delta%+.1f", DELTAS_MNAR)))
imps <- readRDS(file.path(WORK_DIR, sprintf("imputations_%s.rds", k15)))$imps[SCEN]
names(imps) <- names(SCEN)
stopifnot(all(sapply(imps, length) == M_IMP), all(sapply(imps[["MAR"]], nrow) == nrow(d)))

# ---- Patient-level years (matrix of patients × imputed datasets) -----------------------------------------------------
#   method = "trapezoid": linear interpolation between adjacent assessment time points (main analysis); "step": previous state persists until the next assessment (sensitivity analysis)
integ_w <- function(tt, method) { dt <- diff(tt); if (method == "trapezoid") (c(dt, 0) + c(0, dt)) / 2 else c(dt, 0) }
years_mat <- function(Ylist, mrs_dc, method = "trapezoid") {
  w <- integ_w(TT, method)
  one <- function(Y, cut) as.vector((cbind(mrs_dc, Y) <= cut) %*% w)
  A <- vapply(Ylist, one, numeric(length(mrs_dc)), cut = 5)
  I <- vapply(Ylist, one, numeric(length(mrs_dc)), cut = 2)
  list(alive = A, indep = I, dep = A - I)
}
YM <- lapply(imps, years_mat, mrs_dc = d$mrs_dc)
YM_step <- years_mat(imps[["MAR"]], d$mrs_dc, "step")

# ---- Rubin's rules (years; pooled on the original scale) ------------------------------------------------------
rubin_pool <- function(Q, U, dfcom = Inf) {
  M <- length(Q); Qbar <- mean(Q); W <- mean(U); B <- if (M > 1) var(Q) else 0
  Tt <- W + (1 + 1 / M) * B
  if (B < 1e-12 * max(Tt, 1e-12)) nu <- dfcom else {       # complete-data df if there is no between-imputation variation
    lambda <- (1 + 1 / M) * B / Tt
    nu_old <- (M - 1) / lambda^2
    nu <- if (is.finite(dfcom)) 1 / (1 / nu_old + 1 / ((dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda))) else nu_old
  }
  q <- if (is.finite(nu)) qt(0.975, nu) else qnorm(0.975)
  c(est = Qbar, lcl = Qbar - q * sqrt(Tt), ucl = Qbar + q * sqrt(Tt))
}
col_var <- function(X) { n <- nrow(X); if (n < 2) return(rep(0, ncol(X))); (colMeans(X^2) - colMeans(X)^2) * n / (n - 1) }
VARS <- c("alive", "indep", "dep")

# Years for a group (rows: patient row numbers)
pool_years <- function(ym, rows) {
  n <- length(rows)
  r <- lapply(VARS, function(v) { X <- ym[[v]][rows, , drop = FALSE]; rubin_pool(colMeans(X), col_var(X) / n, n - 1) })
  data.frame(n = n, t(setNames(unlist(r), paste(rep(VARS, each = 3), c("est", "lcl", "ucl"), sep = "_"))))
}
# Direct standardization (rows_s: row numbers per stratum, w: stratum weights). If diff_rows_s is given, the difference of the two standardized values (rows_s − diff_rows_s)
std_QU <- function(X, rows_s, w) {
  stopifnot(length(rows_s) == length(w), abs(sum(w) - 1) < 1e-10, all(lengths(rows_s) > 0 | w == 0))
  Q <- 0; U <- 0
  for (s in which(w > 0)) { Xs <- X[rows_s[[s]], , drop = FALSE]; Q <- Q + w[s] * colMeans(Xs); U <- U + w[s]^2 * col_var(Xs) / nrow(Xs) }
  list(Q = Q, U = U)
}
std_years <- function(ym, rows_s, w, diff_rows_s = NULL) {
  r <- lapply(VARS, function(v) {
    a <- std_QU(ym[[v]], rows_s, w)
    if (!is.null(diff_rows_s)) { b <- std_QU(ym[[v]], diff_rows_s, w); a <- list(Q = a$Q - b$Q, U = a$U + b$U) }
    rubin_pool(a$Q, a$U)
  })
  data.frame(n = sum(lengths(rows_s)), t(setNames(unlist(r), paste(rep(VARS, each = 3), c("est", "lcl", "ucl"), sep = "_"))))
}

# ---- Group definitions ---------------------------------------------------------------------------------
AGE <- age_group(d$age); AGE3 <- age_group(d3$age)
SEXL <- c("Men", "Women"); sexf <- function(x) factor(ifelse(x$male == 1, "Men", "Women"), levels = SEXL)
SEX <- sexf(d); SEX3 <- sexf(d3)
DC <- c(list("All patients" = 0:5, "Discharge mRS 0-2" = 0:2, "Discharge mRS 3-5" = 3:5), setNames(as.list(0:5), paste0("Discharge mRS ", 0:5)))
DC3 <- DC[1:3]
in_dc <- function(x, g) x$mrs_dc %in% DC[[g]]
bind_groups <- function(lst) bind_rows(lapply(names(lst), function(nm) cbind(key = nm, lst[[nm]])))

# ---- 1. By discharge mRS (Table 2 and Figure 2 of the paper), step function (Supplemental Methods S5), δ scenarios (Table S7) ---------------------------------------------
by_dc <- function(ym, groups = names(DC)) bind_rows(lapply(groups, function(g) cbind(group = g, pool_years(ym, which(in_dc(d, g))))))
y_dc   <- by_dc(YM[["MAR"]])
y_step <- by_dc(YM_step)
y_scen <- bind_rows(lapply(names(SCEN), function(m) cbind(method = m, by_dc(YM[[m]]))))

# ---- 2. By age group (4 groups) (Figure 3 and Table S5 of the paper) ------------------------------------------------------------
y_age <- bind_rows(lapply(names(DC), function(g) bind_rows(lapply(AGE_LABELS, function(a)
  cbind(group = g, age = a, pool_years(YM[["MAR"]], which(in_dc(d, g) & AGE == a)))))))

# ---- 3. Sex: unadjusted, and values and male–female differences standardized to the age distribution of all G1 (Table S6 of the paper) ---------------------------------
w_age <- as.numeric(table(AGE)) / nrow(d)
rows_age <- function(sel) lapply(AGE_LABELS, function(a) which(sel & AGE == a))
y_sex <- bind_rows(lapply(names(DC3), function(g) {
  rs <- lapply(SEXL, function(s) rows_age(in_dc(d, g) & SEX == s)); names(rs) <- SEXL
  bind_rows(
    bind_rows(lapply(SEXL, function(s) cbind(group = g, sex = s, type = "Unadjusted", pool_years(YM[["MAR"]], which(in_dc(d, g) & SEX == s))))),
    bind_rows(lapply(SEXL, function(s) cbind(group = g, sex = s, type = "Age-standardized", std_years(YM[["MAR"]], rs[[s]], w_age)))),
    cbind(group = g, sex = "Men − Women", type = "Age-standardized", std_years(YM[["MAR"]], rs[["Men"]], w_age, diff_rows_s = rs[["Women"]])),
    { u <- lapply(SEXL, function(s) list(which(in_dc(d, g) & SEX == s)))
      cbind(group = g, sex = "Men − Women", type = "Unadjusted", std_years(YM[["MAR"]], u[[1]], 1, diff_rows_s = u[[2]])) })
}))
# Age-standardized survival and independent survival at 15 years (standardized per imputed dataset and pooled on the logit scale; estimate_std in 01_functions.R)
p_sex <- bind_rows(lapply(names(DC3), function(g) bind_rows(lapply(SEXL, function(s)
  cbind(group = g, sex = s, estimate_std(imps[["MAR"]], rows_age(in_dc(d, g) & SEX == s), w_age, coh$tp, times = H15)[, c("outcome", "est", "lcl", "ucl")])))))

# ---- 4. Application to the recent (G3) case mix (Table S9 of the paper): 8 strata of age (4 groups) × sex, weighted within each group by the distribution of the same group in G3 ----
#   Range = minimum to maximum of the point estimates of the main analysis and each scenario δ = −1.0 to +1.0
STRATA8 <- expand.grid(sex = c("men", "women"), age = AGE_LABELS, stringsAsFactors = FALSE)        # same order as the group names in est_all
rows8 <- function(x, sel, ag, sx) lapply(seq_len(nrow(STRATA8)), function(i) which(sel & ag == STRATA8$age[i] & sx == ifelse(STRATA8$sex[i] == "men", "Men", "Women")))
recent_years <- function(rows_s, w) {
  v <- sapply(names(SCEN), function(m) sapply(VARS, function(x) mean(std_QU(YM[[m]][[x]], rows_s, w)$Q)))
  data.frame(value = VARS, est = v[, "MAR"], range_min = apply(v, 1, min), range_max = apply(v, 1, max))
}
recent_prop <- function(strata_groups, w) {              # weighted mean of stratum-pooled proportions (est_all) (point estimates only)
  bind_rows(lapply(c("Alive", "mRS 0-2"), function(o) {
    v <- sapply(names(SCEN), function(m) {
      e <- est_all %>% filter(H == H15, method == m, outcome == o, time == H15, group %in% strata_groups)
      sum(w * e$est[match(strata_groups, e$group)]) })
    data.frame(value = ifelse(o == "Alive", "p15_alive", "p15_indep"), est = 100 * v[["MAR"]], range_min = 100 * min(v), range_max = 100 * max(v))
  }))
}
recent <- bind_rows(lapply(names(DC3), function(g) {
  r1 <- rows8(d, in_dc(d, g), AGE, SEX); n3 <- lengths(rows8(d3, d3$mrs_dc %in% DC[[g]], AGE3, SEX3)); w <- n3 / sum(n3)
  sg <- paste0(if (g == "All patients") "" else paste0(g, ", "), sprintf("Age %s, %s", STRATA8$age, STRATA8$sex))
  cbind(group = g, strata = "Age × sex (8 strata)", rbind(recent_years(r1, w), recent_prop(sg, w)))
}))
# Reference: 16 strata additionally including discharge mRS (0–2 / 3–5) (also reflects the change in the proportion with discharge mRS 3–5)
r16 <- c(rows8(d, in_dc(d, "Discharge mRS 0-2"), AGE, SEX), rows8(d, in_dc(d, "Discharge mRS 3-5"), AGE, SEX))
n16 <- lengths(c(rows8(d3, d3$mrs_dc <= 2, AGE3, SEX3), rows8(d3, d3$mrs_dc >= 3, AGE3, SEX3))); w16 <- n16 / sum(n16)
sg16 <- c(sprintf("Discharge mRS 0-2, Age %s, %s", STRATA8$age, STRATA8$sex), sprintf("Discharge mRS 3-5, Age %s, %s", STRATA8$age, STRATA8$sex))
recent <- rbind(recent, cbind(group = "All patients", strata = "Age × sex × discharge mRS 0-2/3-5 (16 strata)", rbind(recent_years(r16, w16), recent_prop(sg16, w16))))
rownames(recent) <- NULL

# ---- 5. By prestroke mRS (not included in the paper; kept for responses to reviewers) ---------------------------------------------------------------
PRE <- ifelse(d$pre_mrs <= 2, "Prestroke mRS 0-2", "Prestroke mRS 3-5")
y_pre <- bind_rows(lapply(unique(sort(PRE)), function(p) bind_rows(lapply(names(DC), function(g) {
  rows <- which(PRE == p & in_dc(d, g)); if (!length(rows)) return(NULL)
  cbind(prestroke = p, group = g, pool_years(YM[["MAR"]], rows)) }))))

# ---- 6. Years in complete-case analysis (Table S7 of the paper): patients with mRS observed at all time points. 95% CI from the t distribution -------------------------
cc <- complete.cases(Yobs)
ym_cc <- years_mat(list(Yobs[cc, , drop = FALSE]), d$mrs_dc[cc])
y_cc <- bind_rows(lapply(names(DC), function(g) {
  ix <- which(in_dc(d[cc, ], g)); n <- length(ix)
  r <- lapply(VARS, function(v) { x <- ym_cc[[v]][ix, 1]; m <- mean(x); s <- sd(x) / sqrt(n); c(m, m - qt(0.975, n - 1) * s, m + qt(0.975, n - 1) * s) })
  data.frame(group = g, n = n, t(setNames(unlist(r), paste(rep(VARS, each = 3), c("est", "lcl", "ucl"), sep = "_"))))
}))

# ---- 7. Expected years alive of the general population (Ederer I method; Figure S3 of the paper) ------------------------------------------------
#   Starting from the onset date, multiply (1 − qx) (qx = probability of dying within 1 year) for the calendar year and attained age of each follow-up year. The value at age 104 is used for ages 105 and over.
#   At 3 and 6 months, the hazard within the year is assumed constant: (1 − qx)^f (f = fraction of the year).
lt <- read.csv(LIFETABLE_PATH, fileEncoding = "UTF-8"); lt <- lt[lt$table_type == "abridged", c("year", "sex", "age", "qx")]
qx_arr <- with(lt, tapply(qx, list(age, sex, year), identity)); LT_MAXAGE <- 104; LT_YEARS <- range(lt$year)
qx_at <- function(age, male, year) qx_arr[cbind(as.character(pmin(pmax(floor(age), 0), LT_MAXAGE)), ifelse(male == 1, "male", "female"),
                                                 as.character(pmin(pmax(year, LT_YEARS[1]), LT_YEARS[2])))]
y0 <- as.numeric(format(d$onset, "%Y"))
if (min(y0) < LT_YEARS[1] || max(y0) + H15 - 1 > LT_YEARS[2]) stop("Calendar years outside the range of the life table are referenced")
S <- matrix(1, nrow(d), H15 + 1)
for (l in seq_len(H15)) S[, l + 1] <- S[, l] * (1 - qx_at(d$age + l - 1, d$male, y0 + l - 1))
Sexp <- sapply(TT, function(t) { k <- floor(t); f <- t - k; if (f > 0) S[, k + 1] * (1 - qx_at(d$age + k, d$male, y0 + k))^f else S[, k + 1] })
stopifnot(!anyNA(Sexp))
exp_one <- function(rows, lab) { Sg <- colMeans(Sexp[rows, , drop = FALSE])
  data.frame(group = lab, n = length(rows), expected_years = sum(integ_w(TT, "trapezoid") * Sg),
             S5 = 100 * Sg[TT == 5], S10 = 100 * Sg[TT == 10], S15 = 100 * Sg[TT == 15]) }
y_exp <- rbind(bind_rows(lapply(names(DC), function(g) exp_one(which(in_dc(d, g)), g))),
               bind_rows(lapply(AGE_LABELS, function(a) exp_one(which(AGE == a), paste("Age", a)))))
y_exp <- merge(y_exp, transform(y_dc, alive_years = alive_est)[, c("group", "alive_years")], by = "group", all.x = TRUE, sort = FALSE)
y_exp$difference <- y_exp$alive_years - y_exp$expected_years
# Expected survival at each assessment time point (all patients, discharge mRS 0–2 and 3–5; dashed lines in Figure S3)
exp_curve <- bind_rows(lapply(names(DC3), function(g) data.frame(group = g, time = TT, expected = colMeans(Sexp[in_dc(d, g), , drop = FALSE]))))

# ---- Check: years equal the area under the curves of survival / independent survival at each time point (averaged across imputations) ------------------------------
chk <- sapply(names(DC), function(g) { rows <- which(in_dc(d, g))
  p <- sapply(imps[["MAR"]], function(Y) colMeans(cbind(d$mrs_dc, Y)[rows, , drop = FALSE] <= 2))
  sum(integ_w(TT, "trapezoid") * rowMeans(p)) - y_dc$indep_est[y_dc$group == g] })
stopifnot(max(abs(chk)) < 1e-10)

res <- list(by_dc = y_dc, by_dc_step = y_step, by_dc_scenario = y_scen, by_age = y_age, by_sex = y_sex, p15_sex_std = p_sex,
            recent = recent, by_prestroke = y_pre, complete_case = y_cc, expected = y_exp, expected_curve = exp_curve,
            weights = list(age_G1 = setNames(w_age, AGE_LABELS)), times = TT)
saveRDS(res, file.path(WORK_DIR, "estimates_years.rds"))

# ---- Output (xlsx) ---------------------------------------------------------------------------------
write_xlsx_sheets(list("By discharge mRS" = y_dc, "By age group" = y_age, "By sex (unadj, age-std)" = y_sex, "By sex_15y proportion (age-std)" = p_sex,
                       "Recent case mix" = recent, "By prestroke mRS" = y_pre, "Complete-case analysis" = y_cc, "By δ scenario" = y_scen,
                       "Step function (sensitivity)" = y_step, "Expected years (general pop)" = y_exp),
                  file.path(yd, "TableY_years_alive_and_independent_over_15y_G1.xlsx"),
                  notes = c("Mean years during the 15 years after onset in G1 (15-year cohort) (06_years_lived.R). alive = years alive, indep = years alive and independent (mRS 0–2), dep = years alive and dependent (mRS 3–5). est = point estimate, lcl–ucl = 95% confidence interval.",
                            "Years are the 0/1 indicators of alive / alive and independent at the assessment time points from discharge (year 0) to 15 years, integrated per patient by the trapezoidal rule. Group means and variances per imputed dataset (100 datasets) are pooled by Rubin's rules (Barnard–Rubin degrees of freedom).",
                            "Age standardization by sex: weighted by the age (4 groups) distribution of all G1 (direct method). Men − Women is the male–female difference of standardized (or unadjusted) values. The 15-year proportions are standardized per imputed dataset and pooled on the logit scale.",
                            "Recent case mix: stratum-specific estimates in G1 weighted by the stratum proportions of G3 (onset 2015–2019) (point estimates). range = minimum to maximum of the point estimates of the main analysis and δ = −1.0 to +1.0. p15 = proportion at 15 years (%).",
                            "Complete-case analysis: mean in patients with mRS observed at all assessment time points (95% CI from the t distribution). δ scenarios: imputed data from the delta-adjusted pattern-mixture model. Step function: years assuming the previous state persists until the next assessment.",
                            "Expected years alive: expected survival from the abridged life tables (Ederer I method) integrated over the same assessment time points. difference = years alive − expected years alive."))
message(sprintf("  G1 all patients: years alive %.2f, years independent %.2f, years dependent %.2f", y_dc$alive_est[1], y_dc$indep_est[1], y_dc$dep_est[1]))
