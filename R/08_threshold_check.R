#-----------------------------------------------------------------------------------
#  08_threshold_check.R  —  Check whether the order of drawn thresholds ζ can reverse in the proportional odds model for mRS 0–5 in survivors
#                           (basis for "adjacent thresholds differ by at least 1.8; probability of reversal < 10^-13" in Supplemental Methods 3)
#
#  For each cohort and assessment time point, the two-stage model (fit_twopart) is fitted to observed mRS using the covariates of the
#  first imputed dataset and the preceding mRS. The gap between adjacent thresholds ζ_{k+1} − ζ_k and the probability that the gap is negative
#  when drawn from the approximate posterior N(ζ, V), Φ(−gap / SD), are computed; the minimum gap and maximum probability per time point are output.
#  Output: output/90_diagnostics/TableD6_threshold_order_reversal.xlsx
#  Run: as part of Rscript R/run_all.R (after 03_imputation.R; takes a few minutes)
#-----------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
PROJECT_ROOT <- if (length(f)) normalizePath(file.path(dirname(f), ".."), winslash = "/") else normalizePath(".", winslash = "/")
setwd(PROJECT_ROOT)
source("R/00_config.R", encoding = "UTF-8")
source("R/01_functions.R", encoding = "UTF-8")

cohorts <- readRDS(file.path(WORK_DIR, "cohorts.rds"))$cohorts
thr <- list()
for (k in names(cohorts)) {
  coh <- cohorts[[k]]
  cb <- mice::complete(readRDS(file.path(WORK_DIR, sprintf("imputations_%s.rds", k)))$mice, 1)
  cv <- cb[, COVARIATES]; cv$ckd <- as.numeric(as.character(cv$ckd)); cv$prior_stroke <- as.numeric(as.character(cv$prior_stroke))
  X <- build_X(cv); Yobs <- as.matrix(coh$data[, coh$tp$col]); Y <- as.matrix(cb[, coh$tp$col]); mrs_dc <- coh$data$mrs_dc
  for (j in seq_len(ncol(Yobs))) {
    prev <- if (j == 1) mrs_dc else Y[, j - 1]
    fit_rows <- which(!is.na(Yobs[, j]) & prev < 6)
    assign("FIT_LOG", fit_log_new(), envir = globalenv())
    tp <- suppressWarnings(fit_twopart(design_Z(X, prev, Y, j), Yobs[, j], fit_rows))
    if (!inherits(tp$mrs$fit, "polr")) next                     # 2 or fewer survivor mRS categories (no threshold ordering issue)
    fm <- tp$mrs$fit; z <- fm$zeta; V <- vcov(fm); iz <- (length(coef(fm)) + 1):nrow(V); Vz <- V[iz, iz]
    pr <- sapply(seq_len(length(z) - 1), function(i)
      pnorm(-(z[i + 1] - z[i]) / sqrt(Vz[i, i] + Vz[i + 1, i + 1] - 2 * Vz[i, i + 1])))
    thr[[length(thr) + 1]] <- data.frame(cohort = coh$name, time = coh$tp$lab[j], n_fit = length(fit_rows),
                                          min_gap = min(diff(z)), max_p_reversal = max(pr))
  }
}
thr <- do.call(rbind, thr)
dd <- dir_create(OUT_DIR, "90_diagnostics")
openxlsx::write.xlsx(list(
  "Threshold order reversal" = thr,
  "Notes" = data.frame(Note = c(
    "Proportional odds model for mRS 0–5 in survivors for each cohort and assessment time point (fitted to observed mRS using the covariates of the first imputed dataset and the preceding mRS).",
    "min_gap = minimum gap between adjacent thresholds (logit scale). max_p_reversal = maximum probability that adjacent thresholds reverse order when drawn from the approximate posterior N(ζ, V), Φ(−gap/SD).",
    sprintf("Overall: smallest min_gap %.2f, largest max_p_reversal %.2g.", min(thr$min_gap), max(thr$max_p_reversal))))),
  file.path(dd, "TableD6_threshold_order_reversal.xlsx"))
message(sprintf("smallest min_gap = %.3f, largest max_p_reversal = %.3g", min(thr$min_gap), max(thr$max_p_reversal)))
