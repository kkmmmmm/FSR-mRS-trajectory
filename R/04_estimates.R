#-----------------------------------------------------------------------------------
#  04_estimates.R  —  Estimation of mRS distribution, survival, and independent survival at each time point (3 cohorts × handling of missing data)
#
#   For each cohort (G1 = 15-year, G2 = 10-year, G3 = 5-year cohort), estimate at discharge and at each assessment time point
#     the proportion in each mRS category 0–6, survival (mRS 0–5), and independent survival (alive with mRS 0–2)
#   for all patients, by discharge mRS, by age group × sex, etc. (GROUPS in 01_functions.R).
#     MAR            : multiple imputation (main analysis). Pooled by Rubin's rules (logit scale, Barnard–Rubin degrees of freedom)
#     MNAR_delta±x.x : delta-adjusted pattern-mixture model (sensitivity analysis)
#     CC             : complete-case analysis (95% CI by the Wilson method)
#     KM_obs         : Kaplan–Meier on observed data (loss to follow-up censored at the last observed time point; survival only)
#   Also outputs imputation model diagnostics (imputed vs observed distributions, fit and cross-validation; values in Supplemental Methods).
#   Output: work/estimates_all.rds, output/estimates_all_cohorts_all_methods_long.csv, output/90_diagnostics/
#-----------------------------------------------------------------------------------

cohorts <- readRDS(file.path(WORK_DIR, "cohorts.rds"))$cohorts
SCEN <- c("MAR" = "delta_+0.00", setNames(mnar_key(DELTAS_MNAR), sprintf("MNAR_delta%+.1f", DELTAS_MNAR)))

all_est <- list()
for (k in names(cohorts)) {
  coh <- cohorts[[k]]; H <- coh$H
  Yobs <- as.matrix(coh$data[, coh$tp$col]); base <- coh$data[, c("mrs_dc", "age", "male")]   # patient characteristics used to define groups
  imp_obj <- readRDS(file.path(WORK_DIR, sprintf("imputations_%s.rds", k))); imps <- imp_obj$imps; mids <- imp_obj$mice
  fit_log <- imp_obj$fit_log; rm(imp_obj)
  message("  ", coh$name)

  est_list <- list(CC = estimate_cc(Yobs, base, coh$tp))
  for (m in names(SCEN)) est_list[[m]] <- estimate_mi(imps[[SCEN[[m]]]], base, coh$tp)
  est_list[["KM_obs"]] <- estimate_km_observed(Yobs, base, coh$tp)
  all_est[[k]] <- bind_rows(est_list, .id = "method") %>% mutate(cohort = coh$name, H = H)

  # ---- Diagnostics: mRS distribution of imputed vs observed values -------------------------------------------
  dd <- dir_create(OUT_DIR, "90_diagnostics", coh$name)
  write_xlsx_sheets(list("MAR_imputed_vs_observed" = imputed_vs_observed(imps[["delta_+0.00"]], Yobs, coh$tp)),
                    file.path(dd, "TableD2_imputed_vs_observed_mRS_MAR.xlsx"),
                    notes = c("mRS distribution (%) at each time point: observed values and values imputed into missing cells (pooled over 100 imputations).",
                              "Imputed values are based on the characteristics and preceding mRS of patients with missing data (mostly lost to follow-up), so a difference from the observed distribution is expected in itself."))

  # ---- Diagnostics: fit of the imputation model (two-stage model), calibration of the proportion dead, log of all fits ----------------
  md_k <- imputation_model_diag(mids, coh)
  write_xlsx_sheets(list("Fit diagnostics" = md_k$diag,
                         "Death calibration" = md_k$calib %>%
                           transmute(time = time_lab, prev_mrs = prev_mrs, n, deaths = deaths, `death_observed(%)` = round(100 * obs_death, 1),
                                     `death_predicted(%)` = round(100 * pred_death, 1), `death_predicted_cv(%)` = round(100 * pred_death_cv, 1),
                                     `diff_cv_predicted−observed(%points)` = round(100 * (pred_death_cv - obs_death), 1)),
                         "Log of all fits" = data.frame(stage = rownames(fit_log), fit_log, row.names = NULL)),
                    file.path(dd, "TableD3_imputation_model_fit_MAR.xlsx"),
                    notes = c("Fit diagnostics / death calibration: the imputation model for mRS at each time point (two-stage model: (1) logistic regression for death, (2) proportional odds model for mRS 0–5 in survivors) refitted with the covariates and mRS of the first imputed dataset.",
                              "obs_death = proportion dead among patients used for estimation (alive at the previous time point and observed at the current one), pred_death_fit = model-predicted proportion dead in the same patients, pred_death_cv = predicted proportion dead in 5-fold cross-validation.",
                              "brier_death_cv = cross-validated Brier score for death, logloss_mrs_cv = cross-validated log loss for mRS 0–6 (smaller is better for both).",
                              "Log of all fits: for all fits performed in MICE and in MNAR regeneration, d_ = death model, m_ = mRS model in survivors. fits = number of fits, warn = number with warnings, nonconv = number not converged, nodraw = number in which coefficients could not be drawn."))
  rm(imps, mids); gc()
}

est_all <- bind_rows(all_est)
saveRDS(est_all, file.path(WORK_DIR, "estimates_all.rds"))
write.csv(est_all, file.path(OUT_DIR, "estimates_all_cohorts_all_methods_long.csv"), row.names = FALSE, fileEncoding = "UTF-8")
