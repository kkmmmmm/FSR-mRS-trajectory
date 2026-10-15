#-----------------------------------------------------------------------------------
#  03_imputation.R  —  Multiple imputation (MAR) and delta-adjusted pattern-mixture model (MNAR)
#
#  MAR : Multivariate imputation by chained equations (MICE) with mice(), jointly imputing covariates (CKD, BMI; previous stroke if missing) and mRS at each time point.
#        - Covariates: logreg (CKD; previous stroke if missing), pmm (BMI). Predictors are all covariates and mRS at all time points
#        - mRS       : custom imputation method twopart (two-stage model). If dead (6) at the previous time point, 6 (death is absorbing).
#                      If alive at the previous time point, with Z = (X, previous mRS [categorical], mRS at each earlier time point [numeric]),
#                        (1) Death     : logit P(Y_t = 6 | Y_t-1 <= 5) = α_t + b_t'Z (discrete-time hazard per assessment interval)
#                        (2) Survivors : logit P(Y_t <= k | Y_t <= 5) = ζ_tk − β_t'Z (proportional odds model for mRS 0–5)
#                      are fitted in patients "alive at the previous time point and observed at the current one" ((2) in the survivors among them);
#                      coefficients are drawn from the approximate posterior (multivariate normal) of the estimates, and missing values are
#                      generated as death -> mRS of survivors. Imputation proceeds in time order.
#  MNAR: In each MAR imputed dataset, the above generative model is fitted, coefficients are drawn, the same δ is added to the
#        linear predictors of (1) and (2) for patients with missing data, and missing values are regenerated in time order (δ = DELTAS_MNAR).
#        After death, 6; missing cells with survival observed at a later time point are not set to death.
#        The same random-number stream is used for every scenario.
#  Output: work/imputations_<cohort>.rds (individual level; do not share), output/90_diagnostics/<cohort>/Figure D1 (convergence trace)
#-----------------------------------------------------------------------------------

co <- readRDS(file.path(WORK_DIR, "cohorts.rds"))
cohorts <- co$cohorts

for (k in names(cohorts)) {
  coh <- cohorts[[k]]
  t0 <- Sys.time()
  message(sprintf("  [%s] MAR: multiple imputation with mice (m = %d, maxit = %d)", coh$name, M_IMP, MAXIT_MICE))
  imp <- impute_mice(coh)
  scen <- mnar_scenarios(DELTAS_MNAR)
  message(sprintf("  [%s] MNAR: delta-adjusted pattern-mixture model (delta = %s)", coh$name, paste(scen$delta, collapse = ", ")))
  mn <- run_mnar(imp, coh, scen)
  imps <- c(list("delta_+0.00" = imputed_Y(imp, coh$tp$col)), mn)
  fit_log <- rbind(MICE = attr(imp, "fit_log"), MNAR = attr(mn, "fit_log"))   # log of all fits (counts of warnings, non-convergence, failed draws)
  saveRDS(list(imps = imps, mice = imp, fit_log = fit_log), file.path(WORK_DIR, sprintf("imputations_%s.rds", k)))
  message(sprintf("  [%s] done (%.1f min)", coh$name, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  # Convergence diagnostics (mice trace plots: covariates and mRS at the last time point)
  dd <- dir_create(OUT_DIR, "90_diagnostics", coh$name)
  vars <- intersect(c("ckd", "bmi", "prior_stroke", tail(coh$tp$col, 1)), names(imp$method)[imp$method != ""])
  png(file.path(dd, "FigureD1_MICE_convergence_trace.png"), width = 2400, height = 2400, res = 250)
  print(plot(imp, vars))
  dev.off()
}
