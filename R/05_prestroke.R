#-----------------------------------------------------------------------------------
#  05_prestroke.R  —  Estimates restricted to patients independent before stroke (prestroke mRS 0–2) (Table S8A of the paper)
#
#   Discharge mRS 3–5 mixes patients who newly needed assistance because of the stroke and those who needed assistance before onset.
#   The same estimation as 04_estimates.R (MAR and δ = ±0.5, ±1.0) is performed in patients with prestroke mRS 0–2.
#   - Imputations performed in all patients (03_imputation.R) are used as is and summarized over the relevant patients' rows
#     (prestroke mRS has no missing values and is a predictor in the imputation model).
#   - For contrast, estimates (MAR) in G1 patients with prestroke mRS 3–5 are also obtained.
#   - All 3 cohorts are estimated because the stratum sizes of G3 are needed to apply estimates to the recent (G3) case mix.
#   Output: work/estimates_prestroke02.rds (tables are created by 09_paper_items.R)
#-----------------------------------------------------------------------------------

PRE_MAX <- 2                                               # upper limit of prestroke mRS (0–2 = independent before stroke)
cohorts <- readRDS(file.path(WORK_DIR, "cohorts.rds"))$cohorts
SCEN <- c("MAR" = "delta_+0.00", setNames(mnar_key(DELTAS_MNAR), sprintf("MNAR_delta%+.1f", DELTAS_MNAR)))
est_s <- list(); est_pre35 <- NULL

for (k in names(cohorts)) {
  coh <- cohorts[[k]]; H <- coh$H
  sel <- coh$data$pre_mrs <= PRE_MAX
  stopifnot(!anyNA(sel))
  base <- coh$data[sel, c("mrs_dc", "age", "male")]
  imps <- readRDS(file.path(WORK_DIR, sprintf("imputations_%s.rds", k)))$imps[SCEN]; names(imps) <- names(SCEN)
  stopifnot(all(sapply(imps[["MAR"]], nrow) == nrow(coh$data)))
  message(sprintf("  %s: prestroke mRS 0–%d, %s patients (of %s)", coh$name, PRE_MAX, format(sum(sel), big.mark = ","), format(nrow(coh$data), big.mark = ",")))
  el <- lapply(imps, function(L) estimate_mi(lapply(L, function(Y) Y[sel, , drop = FALSE]), base, coh$tp))
  est_s[[k]] <- bind_rows(el, .id = "method") %>% mutate(cohort = coh$name, H = H)
  if (H == max(HORIZONS))                                  # contrast: G1 patients with prestroke mRS 3–5 (MAR; all patients, discharge mRS 0–2, 3–5)
    est_pre35 <- estimate_mi(lapply(imps[["MAR"]], function(Y) Y[!sel, , drop = FALSE]), coh$data[!sel, c("mrs_dc", "age", "male")], coh$tp,
                             groups = c("All", "Discharge mRS 0-2", "Discharge mRS 3-5")) %>% mutate(method = "MAR", cohort = coh$name, H = H)
  rm(imps); gc()
}
saveRDS(list(est = bind_rows(est_s), est_pre35 = est_pre35), file.path(WORK_DIR, "estimates_prestroke02.rds"))
