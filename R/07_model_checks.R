#-----------------------------------------------------------------------------------
#  07_model_checks.R  —  Checks of the imputation model (values in Supplemental Methods: difference from δ = 0 regeneration, comparison of imputation models)
#
#  (1) Regeneration with δ = 0: the MNAR procedure (refit the model in each MAR imputed dataset and regenerate missing values) is
#      run with δ_D = δ_M = 0, and the difference from the main-analysis MI (MAR) estimates is shown.
#      Because the sensitivity analysis uses the main-analysis MI estimates as the δ = 0 reference, the difference between the two (fluctuation due to regeneration) is checked.
#      The same random-number stream as the other scenarios (seed + 2000 H + b) is used, so results are directly comparable with the δ = ±0.5 etc. estimates.
#      Scope of the comparison: survival and independent survival in all groups (all patients; discharge mRS 0–5, 0–2, 3–5; age group, sex, age group × sex),
#      and years over 15 years in the 15-year cohort (years alive, years independent, years dependent; trapezoidal rule as in 06).
#  (2) Comparison of imputation models: a single proportional odds model (mRS 0–6 with death as the highest mRS category) and
#      the two-stage model (discrete-time hazard of death + proportional odds for mRS 0–5 in survivors) are cross-validated with the same predictors and the same 5 folds,
#      comparing log loss for mRS 0–6, Brier score for death, calibration of death by preceding mRS, and
#      observed vs predicted mRS 0–5 distribution in survivors (covariates of the first imputed dataset; rows used for estimation are those with observed mRS).
#
#  Assumptions: the 3 cohorts are built by the same procedure as 02_data_prep.R (intermediate files not saved), and mice is rerun with the same
#        random seed as 03_imputation.R (the main-analysis imputed data are reproduced; agreement of estimates is checked in the output).
#  Approximate run time: about 1 hour for the 3 cohorts (slightly more aggregation now that all groups are compared). Output: output/90_diagnostics/TableD5
#-----------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
PROJECT_ROOT <- if (length(f)) normalizePath(file.path(dirname(f), ".."), winslash = "/") else normalizePath(".", winslash = "/")
setwd(PROJECT_ROOT)
source("R/00_config.R", encoding = "UTF-8")
source("R/01_functions.R", encoding = "UTF-8")
if (nzchar(Sys.getenv("MRS_CHECK_M"))) M_IMP <- as.integer(Sys.getenv("MRS_CHECK_M"))      # for test runs (normally unset)

# ---- 3 cohorts (same procedure as 02_data_prep.R; non-overlapping onset periods; not saved) ----------------------------------
d0 <- load_data()
d1 <- d0[!is.na(d0$mrs_dc) & d0$mrs_dc != 6, ]
d1 <- enforce_absorbing_death(d1)
cohorts <- lapply(HORIZONS, function(H) make_cohort(d1, H))
names(cohorts) <- sapply(cohorts, `[[`, "key")

est_saved <- readRDS(file.path(WORK_DIR, "estimates_all.rds"))
CHK_GROUPS <- names(GROUPS); CHK_OUTCOMES <- c("Alive", "mRS 0-2")   # groups (all groups) and outcomes to compare

# Group means of years over 15 years (mean across imputed datasets; same trapezoidal rule as 06_years_lived.R): group × measure
years_by_group <- function(Ylist, base, tp, groups = CHK_GROUPS) {
  TT <- c(0, tp$year); dt <- diff(TT); w <- (c(dt, 0) + c(0, dt)) / 2
  gr <- group_rows(base, groups)
  one <- function(cut) vapply(Ylist, function(Y) as.vector((cbind(base$mrs_dc, Y) <= cut) %*% w), numeric(nrow(base)))
  A <- rowMeans(one(5)); I <- rowMeans(one(2))
  bind_rows(lapply(groups, function(g) { r <- gr[[g]]
    data.frame(group = g, n = length(r), measure = c("alive", "indep", "dep"),
               value = if (length(r)) c(mean(A[r]), mean(I[r]), mean(A[r] - I[r])) else NA_real_) }))
}

# ---- Functions for (2): cross-validation of the single proportional odds model (mRS 0–6) and the two-stage model ----------------
p_all_twopart <- function(tp, Z, rows) {
  pd <- p_death_from(tp, Z, rows); pm <- p_mrs_from(tp, Z, rows)
  P <- matrix(0, length(rows), 7); P[, tp$mrs$lev + 1] <- (1 - pd) * pm; P[, 7] <- pd; P
}
fit_single_po <- function(Z, y, fit_rows) {
  cols <- usable_cols(Z, fit_rows); Zs <- Z[fit_rows, cols, drop = FALSE]; ys <- factor(y[fit_rows])
  fit <- suppressWarnings(polr(ys ~ Zs, Hess = TRUE, method = "logistic"))
  list(fit = fit, cols = sub("^Zs", "", names(coef(fit))), lev = as.numeric(levels(ys)))
}
p_all_single <- function(sp, Z, rows) {
  eta <- as.vector(Z[rows, sp$cols, drop = FALSE] %*% coef(sp$fit))
  cum <- plogis(outer(-eta, sp$fit$zeta, "+"))
  pr <- cbind(cum, 1) - cbind(0, cum)
  P <- matrix(0, length(rows), 7); P[, sp$lev + 1] <- pr; P
}

model_compare <- function(imp, coh, b = 1, K = 5) {
  cb <- complete(imp, b)
  cv <- cb[, COVARIATES]; cv$ckd <- as.numeric(as.character(cv$ckd)); cv$prior_stroke <- as.numeric(as.character(cv$prior_stroke))
  X <- build_X(cv); Yobs <- as.matrix(coh$data[, coh$tp$col]); Y <- as.matrix(cb[, coh$tp$col]); mrs_dc <- coh$data$mrs_dc
  res <- lapply(seq_len(ncol(Yobs)), function(j) {
    prev <- if (j == 1) mrs_dc else Y[, j - 1]
    fit_rows <- which(!is.na(Yobs[, j]) & prev < 6)
    Z_j <- design_Z(X, prev, Y, j)
    yj <- Yobs[fit_rows, j]; pv <- pmin(prev[fit_rows], 5)
    fold <- (seq_along(fit_rows) - 1) %% K + 1                 # same deterministic split as the earlier diagnostics (imputation_model_diag)
    P2 <- P1 <- matrix(NA_real_, length(fit_rows), 7)
    for (fo in seq_len(K)) {
      tr <- fit_rows[fold != fo]; te <- which(fold == fo)
      tp2 <- suppressWarnings(fit_twopart(Z_j, Yobs[, j], tr)); P2[te, ] <- p_all_twopart(tp2, Z_j, fit_rows[te])
      sp  <- fit_single_po(Z_j, Yobs[, j], tr);                  P1[te, ] <- p_all_single(sp, Z_j, fit_rows[te])
    }
    ll <- function(P) -mean(log(pmax(P[cbind(seq_along(yj), yj + 1)], 1e-12)))
    br <- function(P) mean((P[, 7] - (yj == 6))^2)
    alive <- yj < 6
    surv_dist <- function(P) colMeans(P[alive, 1:6, drop = FALSE] / rowSums(P[alive, 1:6, drop = FALSE]))
    obs_dist <- tabulate(yj[alive] + 1, 6) / sum(alive)
    list(
      summary = data.frame(time = coh$tp$year[j], time_lab = coh$tp$lab[j], n = length(fit_rows), deaths = sum(yj == 6),
                           logloss_single = ll(P1), logloss_twopart = ll(P2), brier_single = br(P1), brier_twopart = br(P2)),
      calib_death = bind_rows(lapply(0:5, function(k) { r <- pv == k
        data.frame(time = coh$tp$year[j], prev_mrs = k, n = sum(r), obs = mean(yj[r] == 6), pred_single = mean(P1[r, 7]), pred_twopart = mean(P2[r, 7])) })),
      calib_surv = data.frame(time = coh$tp$year[j], mrs = 0:5, n_alive = sum(alive), obs = obs_dist,
                              pred_single = surv_dist(P1), pred_twopart = surv_dist(P2))
    )
  })
  list(summary = bind_rows(lapply(res, `[[`, "summary")), calib_death = bind_rows(lapply(res, `[[`, "calib_death")),
       calib_surv = bind_rows(lapply(res, `[[`, "calib_surv")))
}

# ---- Run (in the order 15-year -> 10-year -> 5-year cohort) --------------------------------------------
out <- list()
run_keys <- rev(names(cohorts))
if (nzchar(Sys.getenv("MRS_CHECK_COHORTS"))) run_keys <- strsplit(Sys.getenv("MRS_CHECK_COHORTS"), ",")[[1]]   # for test runs
for (k in run_keys) {
  coh <- cohorts[[k]]; t0 <- Sys.time()
  message(sprintf("  [%s] rerunning mice (same seed as the main analysis)", coh$name))
  imp <- impute_mice(coh)
  Ymar <- imputed_Y(imp, coh$tp$col)
  base <- coh$data[, c("mrs_dc", "age", "male")]
  est_mar <- estimate_mi(Ymar, base, coh$tp)

  # Check that the main-analysis estimates are reproduced
  chk <- est_mar %>% inner_join(est_saved %>% filter(method == "MAR", cohort == coh$name) %>% select(group, time, outcome, est_saved = est),
                                by = c("group", "time", "outcome"))
  message(sprintf("  [%s] maximum difference from the main-analysis estimates = %.2e", coh$name, max(abs(chk$est - chk$est_saved))))

  # (1) Regeneration with δ = 0 (same random-number stream as the other scenarios)
  scen0 <- data.frame(type = "delta", delta = 0, key = "regen_0", delta_d = 0, delta_m = 0)
  Y0 <- run_mnar(imp, coh, scen0)[["regen_0"]]
  est0 <- estimate_mi(Y0, base, coh$tp, CHK_GROUPS)
  d0cmp <- est_mar %>% filter(group %in% CHK_GROUPS, outcome %in% CHK_OUTCOMES, time > 0) %>%
    select(group, time, time_lab, outcome, n, mar = est, mar_lcl = lcl, mar_ucl = ucl) %>%
    inner_join(est0 %>% select(group, time, outcome, regen0 = est, regen0_lcl = lcl, regen0_ucl = ucl), by = c("group", "time", "outcome")) %>%
    mutate(diff = regen0 - mar, cohort = coh$name, H = coh$H)
  y0cmp <- if (coh$H == max(HORIZONS)) years_by_group(Ymar, base, coh$tp) %>% rename(mar = value) %>%
    inner_join(years_by_group(Y0, base, coh$tp) %>% select(group, measure, regen0 = value), by = c("group", "measure")) %>%
    mutate(diff = regen0 - mar, cohort = coh$name, H = coh$H) else NULL

  # (2) Comparison of imputation models
  mc <- model_compare(imp, coh)
  out[[k]] <- list(delta0 = d0cmp, delta0_years = y0cmp, compare = lapply(mc, function(x) mutate(x, cohort = coh$name, H = coh$H)), reproduce_maxdiff = max(abs(chk$est - chk$est_saved)))
  saveRDS(out, file.path(WORK_DIR, "model_checks.rds"))
  message(sprintf("  [%s] done (%.1f min)", coh$name, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ---- Output (aggregate values only) ----------------------------------------------------------------
bind <- function(nm) bind_rows(lapply(out, function(x) x$compare[[nm]]))
delta0 <- bind_rows(lapply(out, `[[`, "delta0"))
delta0_years <- bind_rows(lapply(out, `[[`, "delta0_years"))
sheets <- list(
  "delta0_regen_vs_MI" = delta0,
  "delta0_regen_vs_MI_years" = delta0_years,
  "delta0_regen_max_diff" = bind_rows(
    delta0 %>% group_by(cohort, outcome) %>% summarise(max_abs_diff_pt = 100 * max(abs(diff)), .groups = "drop") %>% rename(measure = outcome),
    delta0_years %>% group_by(cohort, measure) %>% summarise(max_abs_diff_years = max(abs(diff), na.rm = TRUE), .groups = "drop")),
  "Model comparison_by time" = bind("summary"),
  "Model comparison_by cohort" = bind("summary") %>% group_by(cohort, H) %>%
    summarise(logloss_single = weighted.mean(logloss_single, n), logloss_twopart = weighted.mean(logloss_twopart, n),
              brier_single = weighted.mean(brier_single, n), brier_twopart = weighted.mean(brier_twopart, n),
              deaths = sum(deaths), n = sum(n), .groups = "drop"),                 # n is summarized last (used as weights)
  "Death calib_by prev mRS" = bind("calib_death") %>% group_by(cohort, H, prev_mrs) %>%
    summarise(obs = weighted.mean(obs, n, na.rm = TRUE), pred_single = weighted.mean(pred_single, n, na.rm = TRUE),
              pred_twopart = weighted.mean(pred_twopart, n, na.rm = TRUE), n = sum(n), .groups = "drop"),
  "Survivor mRS dist_by time" = bind("calib_surv"),
  "Survivor mRS dist_by cohort" = bind("calib_surv") %>% group_by(cohort, H, mrs) %>%
    summarise(obs = weighted.mean(obs, n_alive), pred_single = weighted.mean(pred_single, n_alive), pred_twopart = weighted.mean(pred_twopart, n_alive), .groups = "drop"),
  "Reproducibility check" = data.frame(cohort = names(out), max_abs_diff_vs_main = sapply(out, `[[`, "reproduce_maxdiff"))
)
dd <- dir_create(OUT_DIR, "90_diagnostics")
openxlsx::write.xlsx(sheets, file.path(dd, "TableD5_delta0_regeneration_and_model_comparison.xlsx"))
message("==== 07_model_checks.R done ", format(Sys.time()))
