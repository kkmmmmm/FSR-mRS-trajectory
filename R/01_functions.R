#-----------------------------------------------------------------------------------
#  01_functions.R  —  Functions used in the analysis
#-----------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl); library(mice); library(MASS); library(survival)
  library(dplyr); library(tidyr); library(ggplot2); library(openxlsx)
  library(parallel); library(patchwork)
})

dir_create <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
fmt_pct    <- function(x, d = 1) ifelse(is.na(x), "", formatC(100 * pmax(x, 0), format = "f", digits = d))   # percentage (negative rounding errors set to 0)
fmt_ci     <- function(e, l, u, d = 1) ifelse(is.na(e), "", sprintf("%s (%s–%s)", fmt_pct(e, d), fmt_pct(l, d), fmt_pct(u, d)))

#===================================================================================
# 1. Data loading and preparation
#===================================================================================

load_data <- function(path = DATA_PATH, sheet = DATA_SHEET) {
  if (!file.exists(path)) stop("Data file not found. Set DATA_PATH in 00_config.R or the environment variable MRS_DATA_PATH.")
  raw <- as.data.frame(read_excel(path, sheet = sheet, guess_max = 1e5))
  lack <- setdiff(names(VAR_MAP), names(raw))
  if (length(lack)) stop("Required columns are missing: ", paste(lack, collapse = ", "))
  d <- raw[, names(VAR_MAP)]
  names(d) <- unname(VAR_MAP)

  d$onset <- as.Date(d$onset)
  # Prior stroke: values such as "unknown" are treated as missing (no missing values in the current data)
  d$prior_stroke <- suppressWarnings(as.numeric(ifelse(d$prior_stroke %in% c("0", "1", 0, 1), d$prior_stroke, NA)))
  num <- c("age", "male", "pre_mrs", "nihss", "ht", "dm", "dl", "af", "smk", "drk", "ckd", "bmi", "mrs_dc", TP_ALL$col)
  d[num] <- lapply(d[num], function(x) suppressWarnings(as.numeric(x)))
  d$toast <- factor(d$toast)
  for (v in names(TREATMENTS)) d[[v]] <- ifelse(d[[v]] %in% "有", 1, ifelse(d[[v]] %in% "無", 0, NA))   # yes/no -> 1/0 (blank = missing); "有" = yes, "無" = no in the source data
  d$id <- NULL   # patient ID is not needed for the analysis and is dropped immediately
  d
}

# Death (mRS 6) is an absorbing state: all time points from the first mRS 6 onward are set to 6 (corrects inconsistencies in the data)
enforce_absorbing_death <- function(d, cols = TP_ALL$col) {
  m <- as.matrix(d[, cols])
  n_fix <- 0L
  for (i in which(rowSums(m == 6, na.rm = TRUE) > 0)) {
    f <- min(which(m[i, ] == 6))
    after <- f:ncol(m)
    bad <- is.na(m[i, after]) | m[i, after] != 6
    n_fix <- n_fix + sum(bad)
    m[i, after] <- 6
  }
  d[, cols] <- m
  attr(d, "n_absorbing_fix") <- n_fix
  d
}

# Cohort with an observable period of H years: patients with onset after "extraction date − next longer horizon" and on or before "extraction date − H years"
#   The 3 cohorts do not overlap in onset period (not nested). 15-year cohort = onset <= extraction date − 15 y;
#   10-year cohort = extraction date − 15 y < onset <= extraction date − 10 y; 5-year cohort = extraction date − 10 y < onset <= extraction date − 5 y
cutoff_date <- function(H) seq(EXTRACT_DATE, by = paste0("-", H, " years"), length.out = 2)[2]
make_cohort <- function(d, H, horizons = HORIZONS) {
  cutoff <- cutoff_date(H)
  longer <- horizons[horizons > H]
  lower <- if (length(longer)) cutoff_date(min(longer)) else NA      # onset after this date (NA = no lower bound)
  tp <- TP_ALL[TP_ALL$year <= H, ]
  keep <- d$onset <= cutoff & (is.na(lower) | d$onset > lower)
  sub <- d[keep, c("onset", COVARIATES, names(TREATMENTS), tp$col)]
  rownames(sub) <- NULL
  list(H = H, cutoff = cutoff, lower = lower, tp = tp, data = sub,
       name = sprintf("%d-year cohort", H), key = sprintf("%02dy", H))
}


#===================================================================================
# 2. Multiple imputation (following program excode72.r in the Noma course materials)
#   MAR : multiple imputation by chained equations (MICE) with mice(). Covariates (CKD, BMI, prior stroke [if missing]) and mRS at each time point are imputed jointly.
#         mRS uses the custom method "twopart" (two-stage model: (1) logistic regression for death [discrete-time hazard],
#         (2) proportional-odds model for mRS 0–5 among survivors; coefficients drawn from their approximate posterior; death is absorbing)
#   MNAR: same procedure as Program 6 (pmmx) in the Noma course materials. In each MAR imputed dataset complete(imp, b), the model generating
#         missing values (the same two-stage model) is refitted, coefficients are drawn with mvrnorm, delta is added to the linear predictors of both (1) and (2), and missing values are regenerated
#===================================================================================

# Columns used for estimation: columns in which at least 5 patients (among those used for fitting) differ from the mode (drops near-constant columns)
#   Based only on value frequencies, so independent of patient order
informative_cols <- function(Zfit, min_n = 5) apply(Zfit, 2, function(z) { tb <- table(z); sum(tb) - max(tb) >= min_n })

# Log of imputation-model fits (counted only when FIT_LOG is set up on each parallel worker; does not affect estimates or random numbers)
#   For the death model (d_) and the survivor mRS model (m_): fits = number of fits, warn = number of fits with warnings,
#   nonconv = number of non-converged optimizations, nodraw = number of times coefficients could not be drawn and the point estimates were used
fit_log_new <- function() c(d_fits = 0, d_warn = 0, d_nonconv = 0, d_nodraw = 0, m_fits = 0, m_warn = 0, m_nonconv = 0, m_nodraw = 0)
fit_log_add <- function(field) if (exists("FIT_LOG", envir = globalenv())) {
  lg <- get("FIT_LOG", envir = globalenv()); lg[field] <- lg[field] + 1; assign("FIT_LOG", lg, envir = globalenv())
}

# Names of predictor columns after removing near-constant and collinear columns within the fitting rows
usable_cols <- function(Z, rows) {
  Z <- Z[, informative_cols(Z[rows, , drop = FALSE]), drop = FALSE]
  qz <- qr(cbind(1, Z[rows, , drop = FALSE]))
  colnames(Z)[sort(qz$pivot[seq_len(qz$rank)])[-1] - 1]
}

# Fit while counting warnings (warnings are muffled and do not affect estimation)
#   The "non-integer #successes" warning from weighting pseudo-observations is expected and not counted
fit_quiet <- function(expr, tag) {
  warned <- FALSE
  fit <- withCallingHandlers(tryCatch(expr, error = function(e) NULL),
                             warning = function(w) {
                               if (!grepl("non-integer #successes", conditionMessage(w))) warned <<- TRUE
                               invokeRestart("muffleWarning") })
  if (!is.null(fit)) { fit_log_add(paste0(tag, "_fits")); if (warned) fit_log_add(paste0(tag, "_warn")) }
  fit
}

# Draw from the approximate posterior of the coefficients, N(estimate, variance–covariance matrix) (point estimates if the draw fails)
draw_coef <- function(est, V, tag) {
  th <- if (!is.null(V) && all(is.finite(V))) tryCatch(mvrnorm(1, est, V), error = function(e) NULL) else NULL
  if (is.null(th)) { fit_log_add(paste0(tag, "_nodraw")); th <- est }
  th
}

# Guard against complete separation in the death model (e.g., zero deaths in a category): pseudo-observation augmentation as in mice's logreg
#   (White, Daniel & Royston. Comput Stat Data Anal 2010;54:2267-75). For each predictor, pseudo-observations with and without death are placed at
#   mean ± 0.5 SD (clipped to the observed range), with total weight w_tot patients. The mice default (number of predictors + 1) inflates
#   the death proportion at time points with few deaths, so w_tot = AUG_WEIGHT (5) was chosen based on cross-validated calibration and stability of coefficient draws
augment_binary <- function(x, y, w_tot = AUG_WEIGHT) {
  p <- ncol(x); k <- 2; nr <- 2 * p * k
  if (p == 0) return(list(x = x, y = y, w = rep(1, length(y))))
  mn <- colMeans(x); s <- apply(x, 2, sd); lo <- apply(x, 2, min); hi <- apply(x, 2, max)
  a <- matrix(mn, nr, p, byrow = TRUE)
  b <- matrix(rep(c(rep.int(c(0.5, -0.5), k), rep.int(0, nr)), length = nr * p), nr, p)
  dd <- a + b * matrix(s, nr, p, byrow = TRUE)
  dd <- pmin(pmax(dd, matrix(lo, nr, p, byrow = TRUE)), matrix(hi, nr, p, byrow = TRUE))
  colnames(dd) <- colnames(x)
  list(x = rbind(x, dd), y = c(y, rep(rep(0:1, each = 2), p)), w = c(rep(1, length(y)), rep(w_tot / nr, nr)))
}

# Fit the two-stage model (per time point t, estimated in patients fit_rows alive at the previous time point and observed at t)
#   (1) Death: logit P(Y_t = 6 | alive at t−1) = α_t + b_t'Z (discrete-time hazard per assessment interval)
#   (2) mRS 0–5 among survivors: logit P(Y_t <= k | Y_t <= 5) = ζ_tk − β_t'Z (proportional-odds model)
fit_twopart <- function(Z, y, fit_rows) {
  if (anyDuplicated(colnames(Z))) stop("Duplicated predictor column names: ", paste(unique(colnames(Z)[duplicated(colnames(Z))]), collapse = ", "))
  d <- as.integer(y[fit_rows] == 6)
  death <- NULL
  if (sum(d) > 0 && sum(d) < length(d)) {
    cd <- usable_cols(Z, fit_rows)
    au <- augment_binary(Z[fit_rows, cd, drop = FALSE], d)          # guard against complete separation (pseudo-observations)
    Zd <- au$x; da <- au$y; wa <- au$w
    fd <- fit_quiet(glm(da ~ Zd, family = binomial(), weights = wa), "d")
    if (is.null(fd)) stop("glm (death) failed")
    if (!fd$converged) fit_log_add("d_nonconv")
    death <- list(fit = fd, cols = cd)
  }
  s_rows <- fit_rows[y[fit_rows] < 6]
  if (!length(s_rows)) stop("No survivors at this time point (the survivor mRS model cannot be fitted)")
  cm <- usable_cols(Z, s_rows); Zm <- Z[s_rows, cm, drop = FALSE]; ym <- factor(y[s_rows])
  lev <- as.numeric(levels(ym))
  if (length(lev) >= 3) {                                          # usual case: proportional-odds model
    fm <- fit_quiet(polr(ym ~ Zm, Hess = TRUE, method = "logistic"), "m")
    if (is.null(fm)) stop("polr (survivors) failed")
    if (fm$convergence != 0) fit_log_add("m_nonconv")
    mrs <- list(fit = fm, beta = coef(fm), zeta = fm$zeta, V = tryCatch(vcov(fm), error = function(e) NULL),
                cols = sub("^Zm", "", names(coef(fm))), lev = lev)
  } else if (length(lev) == 2) {                                   # only 2 survivor mRS categories: polr() cannot be fitted, so the same
    fg <- fit_quiet(glm(as.integer(ym == levels(ym)[2]) ~ Zm, family = binomial()), "m")   # model (one threshold) is fitted by logistic regression
    if (is.null(fg)) stop("glm (survivors, 2 categories) failed")
    if (!fg$converged) fit_log_add("m_nonconv")
    ok <- !is.na(coef(fg))
    J <- diag(c(-1, rep(1, sum(ok) - 1)), sum(ok))[c(seq_len(sum(ok))[-1], 1), , drop = FALSE]   # (α, b) -> (β = b, ζ = −α)
    V <- tryCatch(J %*% vcov(fg)[ok, ok] %*% t(J), error = function(e) NULL)
    mrs <- list(fit = fg, beta = coef(fg)[ok][-1], zeta = setNames(-coef(fg)[[1]], paste(levels(ym), collapse = "|")), V = V,
                cols = sub("^Zm", "", names(coef(fg))[ok][-1]), lev = lev)
  } else {                                                         # only 1 survivor mRS category: fixed at that value
    mrs <- list(fit = NULL, beta = numeric(0), zeta = numeric(0), V = NULL, cols = character(0), lev = lev)
  }
  list(death = death, p_death_const = mean(d), mrs = mrs)
}

# Death probability (coefficients b: intercept + columns) and cumulative probabilities of mRS among survivors -> category probabilities
#   At a time point with zero deaths (or all deaths), no death model is fitted and the probability is fixed at the observed proportion (0 or 1);
#   adding δ then does not change the death probability (in the study data, every cohort and time point had at least 40 deaths, so this did not occur)
p_death_from <- function(tp, Z, rows, b = NULL, delta = 0) {
  if (is.null(tp$death)) return(rep(tp$p_death_const, length(rows)))
  if (is.null(b)) b <- coef(tp$death$fit)
  plogis(as.vector(cbind(1, Z[rows, tp$death$cols, drop = FALSE]) %*% b) + delta)
}
p_mrs_from <- function(tp, Z, rows, beta = NULL, zeta = NULL, delta = 0) {
  if (!length(tp$mrs$zeta)) return(matrix(1, length(rows), 1))                  # only 1 survivor mRS category
  if (is.null(beta)) { beta <- tp$mrs$beta; zeta <- tp$mrs$zeta }
  eta <- as.vector(Z[rows, tp$mrs$cols, drop = FALSE] %*% beta) + delta
  cum <- plogis(outer(-eta, zeta, "+"))                                         # P(Y <= k | alive)
  cbind(cum, 1) - cbind(0, cum)
}
draw_cat <- function(pr, lev) { u <- runif(nrow(pr)); lev[pmin(1L + rowSums(u > t(apply(pr, 1, cumsum))), length(lev))] }

# Fit the two-stage model and generate missing values by drawing coefficients from their approximate posterior (shared by MAR and MNAR)
#   Z: predictor matrix, y: outcome (numeric 0–6), fit_rows: rows used for fitting, todo: rows to generate
#   delta_d: MNAR sensitivity parameter δ_D (added to the death logit of patients with missing values); delta_m: δ_M (added to the linear predictor of survivor mRS). 0 under MAR
#   no_death: rows of todo with survival observed at a later time point (cannot be dead)
twopart_draw <- function(Z, y, fit_rows, todo, delta_d = 0, delta_m = 0, no_death = NULL) {
  tp <- fit_twopart(Z, y, fit_rows)
  # (1) Death
  bd <- if (!is.null(tp$death)) draw_coef(coef(tp$death$fit), tryCatch(vcov(tp$death$fit), error = function(e) NULL), "d")
  pd <- p_death_from(tp, Z, todo, bd, delta_d)
  if (!is.null(no_death)) pd[no_death] <- 0
  dead <- runif(length(todo)) < pd
  # (2) mRS 0–5 among survivors
  out <- rep(6L, length(todo))
  if (!length(tp$mrs$zeta)) { out[!dead] <- tp$mrs$lev; return(out) }            # only 1 survivor mRS category
  p <- length(tp$mrs$cols)
  th <- draw_coef(c(tp$mrs$beta, tp$mrs$zeta), tp$mrs$V, "m")                    # drawn even if all are dead, to keep the random-number stream
  if (any(!dead)) out[!dead] <- draw_cat(p_mrs_from(tp, Z, todo[!dead], th[seq_len(p)], sort(th[-seq_len(p)]), delta_m), tp$mrs$lev)
  out
}

# Imputation method for mice (called via mice(method = "twopart"))
#   prev: column name of the preceding mRS (passed via blots); no_death: logical vector of length n
#   x contains the covariates of Eq. (3') and mRS at all preceding time points (specified in predictorMatrix)
mice.impute.twopart <- function(y, ry, x, wy = NULL, prev, no_death = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  p1 <- x[, prev]
  out <- rep(6, sum(wy))                                   # 6 if dead at the previous time point (death is absorbing)
  alive <- p1[wy] < 6
  if (!any(alive)) return(out)
  P <- model.matrix(~ factor(pmin(p1, 5), levels = 0:5))[, -1, drop = FALSE]   # preceding mRS (categorical)
  colnames(P) <- paste0("prev", 1:5)
  Z <- cbind(x[, setdiff(colnames(x), prev), drop = FALSE], P)                  # mRS from 2 or more time points earlier kept numeric
  todo <- which(wy)[alive]
  out[alive] <- twopart_draw(Z, y, fit_rows = which(ry & p1 < 6), todo = todo,
                             no_death = if (is.null(no_death)) NULL else no_death[todo])
  out
}

# Covariate matrix of the imputation model (no intercept): used in the MNAR generating model (within mice the same terms are supplied as derived columns)
build_X <- function(cv) {
  cv$agec <- (cv$age - 75) / 10
  cv$lnihss <- log1p(cv$nihss)
  cv$bmic <- cv$bmi - 22.5
  model.matrix(~ agec + I(agec^2) + male + pre_mrs + lnihss + toast + ht + dm + dl + af + smk + drk +
                 ckd + bmic + prior_stroke + mrs_dc, data = cv)[, -1, drop = FALSE]
}

# Time of death and death indicator (for Kaplan–Meier on observed data; censored at the last observed time point if death is not observed)
death_time_status <- function(Y, years) {
  first6 <- apply(Y, 1, function(y) { w <- which(y == 6); if (length(w)) min(w) else NA })
  lastobs <- apply(Y, 1, function(y) { w <- which(!is.na(y)); if (length(w)) max(w) else NA })
  d <- as.numeric(!is.na(first6))
  t <- ifelse(d == 1, years[first6], ifelse(is.na(lastobs), 0, years[lastobs]))
  data.frame(t = t, d = d)
}

# Index of the last time point at which survival (mRS 0-5) was observed (missing cells before it cannot be death)
last_alive_index <- function(Yobs) apply(Yobs, 1, function(y) { w <- which(!is.na(y) & y < 6); if (length(w)) max(w) else 0L })

# Inputs for mice (data, methods, predictor matrix, visit sequence, blots)
mice_setup <- function(coh) {
  tp <- coh$tp$col
  dat <- coh$data[, c(COVARIATES, tp)]
  dat$ckd <- factor(dat$ckd); dat$prior_stroke <- factor(dat$prior_stroke)
  dat$agec <- (dat$age - 75) / 10; dat$agec2 <- dat$agec^2; dat$lnihss <- log1p(dat$nihss)   # no missing values
  dat$bmic <- dat$bmi - 22.5                                                                   # passively imputed from BMI
  meth <- setNames(rep("", ncol(dat)), names(dat))
  meth[c("ckd", "prior_stroke")] <- "logreg"
  meth["bmi"] <- "pmm"
  meth[tp] <- "twopart"
  meth[colSums(is.na(dat)) == 0] <- ""
  if (meth["bmi"] != "") meth["bmic"] <- "~ I(bmi - 22.5)"
  pred <- matrix(0, ncol(dat), ncol(dat), dimnames = list(names(dat), names(dat)))
  for (v in c("ckd", "bmi", "prior_stroke")) pred[v, setdiff(c(COVARIATES, tp), v)] <- 1   # covariates: all covariates and mRS at all time points
  xcov <- c("agec", "agec2", "male", "pre_mrs", "lnihss", "toast", "ht", "dm", "dl", "af", "smk", "drk",
            "ckd", "bmic", "prior_stroke", "mrs_dc")                                          # terms of Eq. (3')
  la <- last_alive_index(as.matrix(coh$data[, tp]))
  blots <- list()
  for (j in seq_along(tp)) {
    prev <- if (j == 1) "mrs_dc" else tp[j - 1]
    pred[tp[j], unique(c(xcov, tp[seq_len(j - 1)]))] <- 1                                   # mRS: predictors of Eq. (3) (all preceding time points)
    blots[[tp[j]]] <- list(prev = prev, no_death = la > j)
  }
  todo <- names(meth)[meth != ""]
  vis <- c(intersect(c("ckd", "bmi", "bmic", "prior_stroke"), todo), intersect(tp, todo))   # covariates -> mRS (in time order)
  list(data = dat, method = meth, predictorMatrix = pred, visitSequence = vis, blots = blots)
}

# MAR: multiple imputation with mice (m imputations split into chunks of 5, run in parallel, and combined with ibind)
#   Chunking and random seeds do not depend on the number of cores, so results are identical in the same computing environment (R, packages, BLAS/LAPACK)
impute_mice <- function(coh, m = M_IMP, maxit = MAXIT_MICE, seed = SEED, n_core = N_CORES, chunk_size = 5) {
  s <- mice_setup(coh)
  chunks <- split(seq_len(m), (seq_len(m) - 1) %/% chunk_size)
  cl <- makeCluster(min(n_core, length(chunks)))
  on.exit(stopCluster(cl))
  clusterEvalQ(cl, suppressPackageStartupMessages({ library(mice); library(MASS) }))
  clusterExport(cl, c(TWOPART_FUNS, "mice.impute.twopart"), envir = environment())
  res <- parLapply(cl, seq_along(chunks), function(k, s, chunks, maxit, seed) {
    assign("FIT_LOG", fit_log_new(), envir = globalenv())                      # log all fits in this chunk
    r <- mice(s$data, m = length(chunks[[k]]), maxit = maxit, method = s$method, predictorMatrix = s$predictorMatrix,
              visitSequence = s$visitSequence, blots = s$blots, seed = seed + k, printFlag = FALSE,
              eps = 0, remove.collinear = FALSE)       # collinearity is handled within twopart (mice's automatic removal disabled)
    list(imp = r, log = get("FIT_LOG", envir = globalenv()))
  }, s = s, chunks = chunks, maxit = maxit, seed = seed + 1000 * coh$H)
  fit_log <- Reduce(`+`, lapply(res, `[[`, "log"))
  res <- lapply(res, `[[`, "imp")
  out <- Reduce(ibind, res)
  # ibind drops elements for variables without missing values (imp is NULL), which makes complete() fail; rebuild preserving variable order
  imp_all <- res[[1]]$imp
  for (v in names(imp_all)) if (!is.null(imp_all[[v]])) {
    iv <- do.call(cbind, lapply(res, function(r) r$imp[[v]]))
    colnames(iv) <- seq_len(ncol(iv))
    imp_all[[v]] <- iv
  }
  out$imp <- imp_all
  attr(out, "fit_log") <- fit_log                    # log of all fits during MICE (m imputations × iterations × time points)
  out
}

# mRS matrix and covariates for each imputed dataset b
imputed_Y <- function(imp, tp_col) lapply(seq_len(imp$m), function(b) {
  Y <- as.matrix(complete(imp, b)[, tp_col]); storage.mode(Y) <- "integer"; Y
})
imputed_cov <- function(imp) lapply(seq_len(imp$m), function(b) {
  cb <- complete(imp, b)[, COVARIATES]
  cb$ckd <- as.numeric(as.character(cb$ckd)); cb$prior_stroke <- as.numeric(as.character(cb$prior_stroke))
  cb
})

# MNAR: procedure of Program 6 (pmmx) in the Noma course materials, applied to mRS in time order
#   Missing values at time point t are regenerated using the covariates of imputed dataset b and the regenerated mRS up to the preceding time point
regenerate_mnar_once <- function(Yobs, X, mrs_dc, delta_d, delta_m) {
  Y <- Yobs
  la <- last_alive_index(Yobs)
  for (j in seq_len(ncol(Yobs))) {
    prev <- if (j == 1) mrs_dc else Y[, j - 1]
    mis <- is.na(Yobs[, j])
    Y[mis & prev == 6, j] <- 6                       # 6 after death
    todo <- which(mis & prev < 6)
    if (!length(todo)) next
    Z <- design_Z(X, prev, Y, j)
    Y[todo, j] <- twopart_draw(Z, Yobs[, j], fit_rows = which(!mis & prev < 6), todo = todo,
                               delta_d = delta_d, delta_m = delta_m, no_death = la[todo] > j)
  }
  storage.mode(Y) <- "integer"
  Y
}

# Predictor matrix of the imputation model at time point j: patient characteristics X, preceding mRS (categorical), and mRS at all earlier time points (numeric)
#   At j = 1 the preceding mRS is the discharge mRS, so the numeric mrs_dc is removed from X and entered only as categorical
#   (the same predictor set as mice.impute.twopart() under MAR; entering both would change the pseudo-observation regularization relative to MAR)
design_Z <- function(X, prev, Y, j) {
  if (j == 1) X <- X[, colnames(X) != "mrs_dc", drop = FALSE]
  P <- model.matrix(~ factor(pmin(prev, 5), levels = 0:5))[, -1, drop = FALSE]
  colnames(P) <- paste0("prev", 1:5)
  L <- if (j >= 3) Y[, seq_len(j - 2), drop = FALSE] else NULL   # mRS at all time points 2 or more earlier (numeric; discharge mRS is in X)
  if (!is.null(L)) colnames(L) <- paste0("mrs_lag_", seq_len(ncol(L)))
  cbind(X, P, L)
}

# Two-stage model functions exported to parallel workers
TWOPART_FUNS <- c("twopart_draw", "fit_twopart", "augment_binary", "usable_cols", "informative_cols", "fit_quiet", "draw_coef",
                  "p_death_from", "p_mrs_from", "draw_cat", "fit_log_new", "fit_log_add", "AUG_WEIGHT")

# MNAR scenarios: for each delta (except 0), the same delta is added to the linear predictors of death and survivor mRS (δ_D = δ_M = δ)
#   key: e.g. "delta_+1.00" (MAR is "delta_+0.00")
mnar_key <- function(delta) sprintf("delta_%+.2f", delta)
mnar_scenarios <- function(deltas) {
  deltas <- setdiff(deltas, 0)
  data.frame(delta = deltas, key = mnar_key(deltas), delta_d = deltas, delta_m = deltas)
}

run_mnar <- function(imp, coh, scen, seed = SEED, n_core = N_CORES) {
  Yobs <- as.matrix(coh$data[, coh$tp$col]); mrs_dc <- coh$data$mrs_dc
  covs <- imputed_cov(imp)
  cl <- makeCluster(n_core)
  on.exit(stopCluster(cl))
  clusterEvalQ(cl, suppressPackageStartupMessages(library(MASS)))
  clusterExport(cl, c(TWOPART_FUNS, "regenerate_mnar_once", "design_Z", "build_X", "last_alive_index"), envir = environment())
  res <- parLapply(cl, seq_along(covs), function(b, covs, Yobs, mrs_dc, scen, seed_b) {
    assign("FIT_LOG", fit_log_new(), envir = globalenv())
    X <- build_X(covs[[b]])
    out <- list()
    for (i in seq_len(nrow(scen))) {
      set.seed(seed_b + b)   # common random numbers across scenarios (for smoother comparisons)
      out[[scen$key[i]]] <- regenerate_mnar_once(Yobs, X, mrs_dc, delta_d = scen$delta_d[i], delta_m = scen$delta_m[i])
    }
    list(Y = out, log = get("FIT_LOG", envir = globalenv()))
  }, covs = covs, Yobs = Yobs, mrs_dc = mrs_dc, scen = scen, seed_b = seed + 2000 * coh$H)
  fit_log <- Reduce(`+`, lapply(res, `[[`, "log"))
  res <- lapply(res, `[[`, "Y")
  scen <- names(res[[1]])
  out <- setNames(lapply(scen, function(s) lapply(res, `[[`, s)), scen)
  attr(out, "fit_log") <- fit_log                    # log of all fits in MNAR regeneration (m imputations × scenarios × time points)
  out
}

#===================================================================================
# 3. Estimation and pooling (Rubin's rules)
#===================================================================================

# Groups analyzed: functions returning group membership (TRUE/FALSE) from patient characteristics b (data frame with mrs_dc, age, male)
#   By discharge mRS, age group (<65 / 65–74 / 75–84 / >=85 years), sex, and age group × sex
age_group <- function(age) cut(age, c(-Inf, 65, 75, 85, Inf), right = FALSE, labels = AGE_LABELS)
by_dc  <- function(k) function(b) b$mrs_dc %in% k
by_age <- function(a, sex = NULL) function(b) age_group(b$age) == a & (if (is.null(sex)) TRUE else b$male == (sex == "men"))
GROUPS <- c(
  list("All" = function(b) rep(TRUE, nrow(b))),
  setNames(lapply(0:5, by_dc), paste0("Discharge mRS ", 0:5)),
  list("Discharge mRS 0-2" = by_dc(0:2), "Discharge mRS 3-5" = by_dc(3:5)),
  list("Men" = function(b) b$male == 1, "Women" = function(b) b$male == 0),
  setNames(lapply(AGE_LABELS, by_age), paste("Age", AGE_LABELS)),
  unlist(lapply(AGE_LABELS, function(a) setNames(lapply(c("men", "women"), function(s) by_age(a, s)),
                                                 sprintf("Age %s, %s", a, c("men", "women")))), recursive = FALSE),
  # Discharge mRS 0–2 / 3–5 × age group × sex (used for applying estimates to the recent case mix within each group)
  unlist(lapply(c("Discharge mRS 0-2", "Discharge mRS 3-5"), function(dg) {
    k <- if (dg == "Discharge mRS 0-2") 0:2 else 3:5
    unlist(lapply(AGE_LABELS, function(a) setNames(lapply(c("men", "women"), function(s) { f <- by_age(a, s); function(b) f(b) & b$mrs_dc %in% k }),
                                                   sprintf("%s, Age %s, %s", dg, a, c("men", "women")))), recursive = FALSE)
  }), recursive = FALSE))
group_rows <- function(b, groups = names(GROUPS)) setNames(lapply(groups, function(g) which(GROUPS[[g]](b))), groups)

# Point estimate of direct age–sex standardization: estimates in age group × sex strata (8 strata) of cohort H, weighted by the stratum composition of cohort W
#   Weighted mean of pooled stratum-specific estimates (pooled on the logit scale and back-transformed). No confidence intervals (point estimates only)
#   W = 5, H = 15: projection of 15-year-cohort stratum estimates onto the recent case mix (5-year cohort) (Table 2 of the paper; point estimates only)
#   method: estimation method ("MAR" or sensitivity analyses such as "MNAR_delta-1.0"; weights always use MAR counts)
#   within: group name when standardizing within a group (e.g. "Discharge mRS 3-5"; strata are "<group>, Age <65, men", etc.). NULL = all patients
AGE_SEX_GROUPS <- sprintf("Age %s, %s", rep(AGE_LABELS, each = 2), c("men", "women"))
std_age_sex <- function(est_all, H, W = max(HORIZONS), outcomes = c("Alive", "mRS 0-2"), within = NULL, method = "MAR") {
  strata <- if (is.null(within)) AGE_SEX_GROUPS else paste0(within, ", ", AGE_SEX_GROUPS)
  w <- est_all %>% filter(method == "MAR", H == !!W, group %in% strata, outcome == "Alive", time == 0) %>%
    transmute(group = sub(paste0("^", if (is.null(within)) "" else paste0(within, ", ")), "", group), w = n / sum(n))
  stopifnot(nrow(w) == length(AGE_SEX_GROUPS))
  est_all %>% filter(method == !!method, H == !!H, group %in% strata, outcome %in% outcomes) %>%
    mutate(group = sub(paste0("^", if (is.null(within)) "" else paste0(within, ", ")), "", group)) %>%
    inner_join(w, by = "group") %>%
    group_by(time, time_lab, outcome) %>%
    summarise(est = sum(w * est), n_strata = n(), .groups = "drop") %>%
    { stopifnot(all(.$n_strata == length(AGE_SEX_GROUPS))); . } %>%
    transmute(time, time_lab, outcome, est, H = H, weights_H = W)
}

# Direct age–sex standardization (with confidence intervals): standardize within each imputed dataset, then pool with Rubin's rules
#   rows_s: row indices per stratum (list; strata are mutually exclusive); w: stratum weights (proportions in the reference cohort; sum to 1, fixed)
#   For each imputed dataset m:  P_m = Σ_s w_s p_sm,  U_m = Σ_s w_s^2 p_sm (1 − p_sm) / n_s
#   (within an imputed dataset, strata contain different patients, so between-stratum covariance is 0; correlation induced by the
#    imputation model shared across strata is captured in the between-imputation variance B = var(P_m))
#   Point estimate and 95% CI by Rubin's rules on the logit scale (pool_logit; U_m transformed by the delta method, Barnard–Rubin df).
#   With a single stratum (w = 1), this reduces to pool_prop() (including the continuity correction)
#   For between-cohort differences, also returns var_raw = T × {p(1 − p)}², the total logit-scale variance T transformed to the
#   probability scale by the delta method (p = pooled point estimate), and its degrees of freedom df_raw (Barnard–Rubin)
br_df <- function(lam, m, n) {                         # Barnard–Rubin degrees of freedom (same formula as pool_prop())
  nu_com <- n - 1
  nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lam)
  if (lam > 0) { nu_old <- (m - 1) / lam^2; nu_old * nu_obs / (nu_old + nu_obs) } else nu_obs
}
# Rubin's rules on the logit scale (shared by pool_prop and pool_std): point estimate and 95% CI pooled on the same scale and back-transformed to a proportion (ref. 41)
#   P: proportion in each imputed dataset, U: its variance (probability scale), n: number of patients
#   Only for imputed datasets where the logit is undefined (proportion 0 or 1), the continuity-corrected proportion Pa (variance Ua) is used
#   If the proportion is 0 or 1 identically across all imputed datasets, that value is the point estimate and the 95% CI is the Wilson interval
pool_logit <- function(P, U, Pa, Ua, n) {
  m <- length(P)
  if (all(P == P[1]) && P[1] %in% c(0, 1)) {
    wl <- wilson(round(P[1] * n), n)
    return(c(est = P[1], lcl = wl[["lcl"]], ucl = wl[["ucl"]], lambda = 0))
  }
  bd <- P <= 0 | P >= 1
  pp <- ifelse(bd, Pa, P); uu <- ifelse(bd, Ua, U)
  q <- qlogis(pp); u <- uu / (pp * (1 - pp))^2                    # delta method: variance on the logit scale
  qbar <- mean(q); W <- mean(u); B <- if (m > 1) var(q) else 0
  Tt <- W + (1 + 1 / m) * B; lam <- (1 + 1 / m) * B / Tt
  ci <- plogis(qbar + c(-1, 1) * qt(0.975, max(br_df(lam, m, n), 1)) * sqrt(Tt))
  c(est = plogis(qbar), lcl = ci[1], ucl = ci[2], lambda = lam)   # lambda = (1 + 1/M)B/T (proportion of between-imputation variance; distinct from the finite-sample-corrected FMI)
}
pool_std <- function(Pm, Um, Pa, Ua, n) {
  m <- length(Pm)
  r <- pool_logit(Pm, Um, Pa, Ua, n)
  bd <- Pm <= 0 | Pm >= 1
  if (all(bd) && all(Pm == Pm[1]))                                 # 0 or 1 in all imputed datasets (variance from the Wilson interval width)
    return(c(r, var_raw = ((r[["ucl"]] - r[["lcl"]]) / (2 * qnorm(0.975)))^2, df_raw = max(n - 1, 1)))
  pp <- ifelse(bd, Pa, Pm); uu <- ifelse(bd, Ua, Um)
  q <- qlogis(pp); u <- uu / (pp * (1 - pp))^2
  Tt <- mean(u) + (1 + 1 / m) * (if (m > 1) var(q) else 0)        # total variance on the logit scale (same as pool_logit)
  p <- r[["est"]]
  c(r, var_raw = Tt * (p * (1 - p))^2, df_raw = max(br_df(r[["lambda"]], m, n), 1))   # delta method to the probability scale (for differences)
}
estimate_std <- function(Ylist, rows_s, w, tp, outcomes = c("Alive", "mRS 0-2"), times = tp$year) {
  ns <- lengths(rows_s)
  stopifnot(length(w) == length(rows_s), abs(sum(w) - 1) < 1e-10, all(ns > 0), !anyDuplicated(unlist(rows_s)))
  out <- list()
  for (t in times) {
    j <- match(t, tp$year)
    for (o in outcomes) {
      S <- matrix(vapply(Ylist, function(Y) vapply(rows_s, function(r) mean(outcome_indicator(Y[r, j], o)), numeric(1)),
                         numeric(length(rows_s))), nrow = length(rows_s))       # strata × imputed datasets
      A <- (S * ns + 0.5) / (ns + 1)
      r <- pool_std(Pm = colSums(w * S), Um = colSums(w^2 * S * (1 - S) / ns),
                    Pa = colSums(w * A), Ua = colSums(w^2 * A * (1 - A) / (ns + 1)), n = sum(ns))
      out[[length(out) + 1]] <- data.frame(time = t, time_lab = tp$lab[j], outcome = o, n = sum(ns), t(r))
    }
  }
  do.call(rbind, out)
}

# Pooling proportions: Rubin's rules on the logit scale (Barnard–Rubin df). The point estimate is also pooled on that scale and back-transformed (pool_logit)
#   mean_est = TRUE: only the point estimate is the mean proportion across imputations (each mRS 0–6 category, so the distribution sums to 100%;
#   confidence intervals are not shown in the paper)
pool_prop <- function(p, n, mean_est = FALSE) {
  if (n == 0) return(c(est = NA_real_, lcl = NA_real_, ucl = NA_real_, lambda = NA_real_))   # group with 0 patients
  pa <- (p * n + 0.5) / (n + 1)
  r <- pool_logit(p, p * (1 - p) / n, pa, pa * (1 - pa) / (n + 1), n)
  if (mean_est) r["est"] <- mean(p)
  r
}

wilson <- function(x, n) {
  if (n == 0) return(c(est = NA_real_, lcl = NA_real_, ucl = NA_real_))
  z <- qnorm(0.975); p <- x / n
  c0 <- p + z^2 / (2 * n); r <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)); dd <- 1 + z^2 / n
  c(est = p, lcl = (c0 - r) / dd, ucl = (c0 + r) / dd)
}

# Estimands: proportion with mRS k, proportion with mRS 0-2, and proportion alive (mRS < 6), per group and time point
OUTCOMES <- c(paste0("mRS ", 0:6), "mRS 0-2", "Alive")
outcome_indicator <- function(y, o) switch(o, "mRS 0-2" = y <= 2, "Alive" = y < 6, y == as.numeric(sub("mRS ", "", o)))

# Summarize multiply imputed data (base: patient characteristics [mrs_dc, age, male]; groups: groups to summarize)
estimate_mi <- function(Ylist, base, tp, groups = names(GROUPS)) {
  times <- c(0, tp$year); labs <- c("D/C", tp$lab)
  mrs_dc <- base$mrs_dc; gr <- group_rows(base, groups)
  out <- list()
  for (g in groups) {
    rows <- gr[[g]]; n <- length(rows)
    for (j in seq_along(times)) {
      for (o in OUTCOMES) {
        p <- if (j == 1) rep(mean(outcome_indicator(mrs_dc[rows], o)), length(Ylist))
             else vapply(Ylist, function(Y) mean(outcome_indicator(Y[rows, j - 1], o)), numeric(1))
        r <- if (j == 1) c(wilson(sum(outcome_indicator(mrs_dc[rows], o)), n), lambda = 0) else pool_prop(p, n, mean_est = !o %in% c("mRS 0-2", "Alive"))
        out[[length(out) + 1]] <- data.frame(group = g, time = times[j], time_lab = labs[j], outcome = o,
                                             n = n, est = r["est"], lcl = r["lcl"], ucl = r["ucl"], lambda = r["lambda"])
      }
    }
  }
  res <- do.call(rbind, out); rownames(res) <- NULL; res
}

# Fit diagnostics of the imputation model (two-stage model): refit the model at each time point using covariates and mRS of imputed dataset b
#   (fitting rows have observed mRS, so they are nearly identical across imputed datasets; b = 1 is used as representative)
#   diag : fit status per time point (convergence of the death and survivor mRS models, warnings), observed vs predicted death proportions,
#          death Brier score and mRS 0–6 log loss in 5-fold cross-validation (internal validation)
#   calib: observed, predicted (fitted), and predicted (cross-validated) death proportions until the next assessment, by time point × preceding mRS (0–5)
#   Cross-validation folds are assigned deterministically by row order (no random numbers)
imputation_model_diag <- function(imp, coh, b = 1, K = 5) {
  cb <- complete(imp, b)
  cv <- cb[, COVARIATES]; cv$ckd <- as.numeric(as.character(cv$ckd)); cv$prior_stroke <- as.numeric(as.character(cv$prior_stroke))
  X <- build_X(cv); Yobs <- as.matrix(coh$data[, coh$tp$col]); Y <- as.matrix(cb[, coh$tp$col]); mrs_dc <- coh$data$mrs_dc
  p_all <- function(tp, Z, rows) {                           # predicted probabilities of mRS 0–6 (two-stage model)
    pd <- p_death_from(tp, Z, rows); pm <- p_mrs_from(tp, Z, rows)
    P <- matrix(0, length(rows), 7); P[, tp$mrs$lev + 1] <- (1 - pd) * pm; P[, 7] <- pd; P
  }
  res <- lapply(seq_len(ncol(Yobs)), function(j) {
    prev <- if (j == 1) mrs_dc else Y[, j - 1]
    mis <- is.na(Yobs[, j]); fit_rows <- which(!mis & prev < 6); todo <- which(mis & prev < 6)
    Z_j <- design_Z(X, prev, Y, j)
    assign("FIT_LOG", fit_log_new(), envir = globalenv())                       # count warnings and non-convergence of this fit
    tp <- fit_twopart(Z_j, Yobs[, j], fit_rows)
    lg <- get("FIT_LOG", envir = globalenv()); rm("FIT_LOG", envir = globalenv())
    yj <- Yobs[fit_rows, j]; dth <- yj == 6; pv <- pmin(prev[fit_rows], 5)
    pf <- p_death_from(tp, Z_j, fit_rows)
    # 5-fold cross-validation: fit without each fold and predict for the held-out patients
    fold <- (seq_along(fit_rows) - 1) %% K + 1
    pcv <- numeric(length(fit_rows)); Pcv <- matrix(NA_real_, length(fit_rows), 7)
    for (f in seq_len(K)) {
      tr <- fit_rows[fold != f]; te <- which(fold == f)
      tpf <- suppressWarnings(fit_twopart(Z_j, Yobs[, j], tr))
      pcv[te] <- p_death_from(tpf, Z_j, fit_rows[te]); Pcv[te, ] <- p_all(tpf, Z_j, fit_rows[te])
    }
    ev <- function(h) tryCatch(min(eigen(h, symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA)
    calib <- bind_rows(lapply(0:5, function(k) { r <- pv == k
      data.frame(time = coh$tp$year[j], time_lab = coh$tp$lab[j], prev_mrs = k, n = sum(r), deaths = sum(dth[r]),
                 obs_death = if (any(r)) mean(dth[r]) else NA_real_, pred_death = if (any(r)) mean(pf[r]) else NA_real_,
                 pred_death_cv = if (any(r)) mean(pcv[r]) else NA_real_) }))
    list(diag = data.frame(time_lab = coh$tp$lab[j], n_fit = length(fit_rows), deaths = sum(dth), n_imputed = length(todo),
                           death_model = !is.null(tp$death),
                           death_converged = if (is.null(tp$death)) NA else tp$death$fit$converged,
                           mrs_converged = if (inherits(tp$mrs$fit, "polr")) tp$mrs$fit$convergence == 0 else NA, mrs_hessian_pd = if (inherits(tp$mrs$fit, "polr")) isTRUE(ev(tp$mrs$fit$Hessian) > 0) else NA,
                           warn_death = unname(lg["d_warn"]), warn_mrs = unname(lg["m_warn"]),
                           max_abs_coef_death = if (is.null(tp$death)) NA else max(abs(coef(tp$death$fit)[-1])),
                           max_se_death = if (is.null(tp$death)) NA else max(sqrt(diag(vcov(tp$death$fit)))),
                           obs_death = mean(dth), pred_death_fit = mean(pf), pred_death_cv = mean(pcv),
                           pred_death_imputed = if (length(todo)) mean(p_death_from(tp, Z_j, todo)) else NA_real_,
                           brier_death_cv = mean((pcv - dth)^2),
                           logloss_mrs_cv = -mean(log(pmax(Pcv[cbind(seq_along(yj), yj + 1)], 1e-12)))),
         calib = calib)
  })
  list(diag = bind_rows(lapply(res, `[[`, "diag")), calib = bind_rows(lapply(res, `[[`, "calib")))
}

# Complete-case analysis: restricted to patients with mRS observed at all time points up to H years (95% CI by the Wilson method)
estimate_cc <- function(Yobs, base, tp, groups = names(GROUPS)) {
  times <- c(0, tp$year); labs <- c("D/C", tp$lab)
  cc <- complete.cases(Yobs)
  mrs_dc <- base$mrs_dc; gr <- group_rows(base, groups)
  out <- list()
  for (g in groups) {
    rows <- gr[[g]]; rows <- rows[cc[rows]]
    for (j in seq_along(times)) {
      for (o in OUTCOMES) {
        y <- if (j == 1) mrs_dc[rows] else Yobs[rows, j - 1]
        r <- wilson(sum(outcome_indicator(y, o)), length(rows))
        out[[length(out) + 1]] <- data.frame(group = g, time = times[j], time_lab = labs[j], outcome = o,
                                             n = length(rows), est = r["est"], lcl = r["lcl"], ucl = r["ucl"], lambda = NA)
      }
    }
  }
  res <- do.call(rbind, out); rownames(res) <- NULL; res
}

# Standard Kaplan–Meier on observed data (loss to follow-up censored, death as the event)
estimate_km_observed <- function(Yobs, base, tp) {
  ts <- death_time_status(Yobs, tp$year)
  do.call(rbind, lapply(c("All", "Discharge mRS 0-2", "Discharge mRS 3-5"), function(g) {
    rows <- GROUPS[[g]](base)
    fit <- survfit(Surv(t, d) ~ 1, data = ts[rows, ], conf.type = "log-log")
    s <- summary(fit, times = tp$year, extend = TRUE)
    rbind(data.frame(group = g, time = 0, time_lab = "D/C", outcome = "Alive", n = sum(rows), est = 1, lcl = 1, ucl = 1, lambda = NA),
          data.frame(group = g, time = tp$year, time_lab = tp$lab, outcome = "Alive", n = s$n.risk,
                     est = s$surv, lcl = s$lower, ucl = s$upper, lambda = NA))
  }))
}

# Distribution of imputed vs observed values (diagnostics)
imputed_vs_observed <- function(Ylist, Yobs, tp) {
  do.call(rbind, lapply(seq_along(tp$col), function(j) {
    mis <- is.na(Yobs[, j])
    obs <- table(factor(Yobs[!mis, j], levels = 0:6))
    impv <- if (any(mis)) table(factor(unlist(lapply(Ylist, function(Y) Y[mis, j])), levels = 0:6)) else rep(0, 7)
    data.frame(time = tp$lab[j], type = c("Observed", "Imputed (pooled over all imputations)"),
               n = c(sum(!mis), sum(mis)),
               rbind(round(100 * prop.table(obs), 1), round(100 * prop.table(impv), 1)),
               check.names = FALSE)
  }))
}

#===================================================================================
# 4. Figures and tables
#===================================================================================

MRS_COLS <- c("0" = "#0d366b", "1" = "#3987e5", "2" = "#9ec5f4", "3" = "#f5b38f",
              "4" = "#eb6834", "5" = "#8f3512", "6" = "#2b2b2a")
MRS_TEXT <- c("0" = "white", "1" = "white", "2" = "#0b0b0b", "3" = "#0b0b0b", "4" = "#0b0b0b", "5" = "white", "6" = "white")
GRID <- "#e4e3df"

# ---- Submission format (following Stroke / AHA figure requirements) ------------------------------------------
#   - Created at final print size: two-column width FIG_W = 6.9 in (17.5 cm)
#   - Arial, 8 pt text and 7 pt tick labels (at final size). No titles inside figures (described in the figure legends)
#   - Saved as PNG (300 dpi; for Word) and PDF (vector; for submission) with the same file name
FIG_FONT <- "Arial"; FIG_W <- 6.9
theme_mrs <- function(base = 8) {
  theme_classic(base_size = base, base_family = FIG_FONT) +
    theme(axis.line = element_line(colour = "black", linewidth = 0.3),
          axis.ticks = element_line(colour = "black", linewidth = 0.3),
          axis.text = element_text(colour = "black", size = base - 1),
          axis.title = element_text(colour = "black", size = base),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", colour = "black", size = base, hjust = 0, margin = margin(0, 0, 2, 0)),
          legend.position = "top", legend.justification = "left",
          legend.title = element_text(size = base), legend.text = element_text(size = base - 1),
          legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -4, 0),
          plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
          plot.tag = element_text(face = "bold", size = base + 2, family = FIG_FONT),
          panel.spacing = unit(4, "mm"), plot.margin = margin(4, 6, 4, 4),
          plot.background = element_rect(fill = "white", colour = NA))
}

# ragg truncates long paths (over 255 bytes in UTF-8), so draw to a temporary file and then copy
#   If w, h (inches) exceed the two-column width, shrink to that width keeping the aspect ratio (font sizes stay at final pt)
save_fig <- function(p, file, w, h) {
  if (w > FIG_W) { h <- h * FIG_W / w; w <- FIG_W }
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
  ok <- file.copy(tmp, file, overwrite = TRUE); unlink(tmp)
  tmp <- tempfile(fileext = ".pdf")
  ggsave(tmp, p, width = w, height = h, device = cairo_pdf, bg = "white")
  ok <- ok && file.copy(tmp, sub("\\.png$", ".pdf", file), overwrite = TRUE); unlink(tmp)
  if (!ok) stop("Failed to save figure: ", file)
}

# Time point labels (English)
time_lab_en <- function(x) {
  x <- as.character(x)
  ifelse(x == "D/C", "Discharge", ifelse(x == "3M", "3 mo", ifelse(x == "6M", "6 mo", sub("Y$", " y", x))))
}
axis_time_lab <- function(x) { x <- as.character(x); ifelse(x == "D/C", "D/C", ifelse(x == "3M", "3 mo", ifelse(x == "6M", "6 mo", sub("Y$", "", x)))) }

# mRS distribution at each time point by discharge mRS (stacked bars)
#   grp: groups to show (default: 7 panels, all patients and discharge mRS 0–5); ncol: number of columns; values: whether to print percentages inside bars
plot_trajectory <- function(est, title = NULL, subtitle = NULL, grp = c("All", paste0("Discharge mRS ", 0:5)), ncol = 2, values = TRUE) {
  df <- est %>% filter(group %in% grp, outcome %in% paste0("mRS ", 0:6)) %>%
    mutate(mrs = factor(sub("mRS ", "", outcome), levels = 6:0),
           time_lab = factor(time_lab, levels = unique(time_lab)),
           group = factor(ifelse(group == "All", "All patients", sub("-", "–", group)), levels = ifelse(grp == "All", "All patients", sub("-", "–", grp))))
  ns <- df %>% distinct(group, n) %>% mutate(lab = sprintf("%s (n = %s)", group, prettyNum(n, big.mark = ",")))
  df$group <- factor(df$group, levels = ns$group, labels = ns$lab)
  nt <- length(unique(df$time_lab))
  ggplot(df, aes(time_lab, est, fill = mrs)) +
    geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
    geom_text(aes(label = ifelse(values & est >= ifelse(nt > 12, 0.12, 0.08), sprintf("%.0f", 100 * est), ""), colour = mrs),
              position = position_stack(vjust = 0.5), size = ifelse(nt > 12, 1.6, 1.9), family = FIG_FONT, show.legend = FALSE) +
    scale_fill_manual(values = MRS_COLS, breaks = as.character(0:6), labels = c(0:5, "6 (death)"), name = "mRS score") +
    scale_colour_manual(values = MRS_TEXT, guide = "none") +
    scale_x_discrete(labels = axis_time_lab) +
    scale_y_continuous(labels = function(x) 100 * x, expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
    facet_wrap(~group, ncol = ncol) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(x = "Time after stroke onset (years unless indicated)", y = "Patients (%)") +
    theme_mrs() + theme(axis.text.x = element_text(angle = ifelse(nt > 12, 90, 0), vjust = 0.5, hjust = ifelse(nt > 12, 1, 0.5)))
}

# Step-shaped confidence band matching geom_step
step_ribbon <- function(d, by = NULL) {
  f <- function(x) {
    x <- x[order(x$time), ]
    n <- nrow(x); if (n < 2) return(x)
    xx <- data.frame(time = c(rbind(x$time[-n], x$time[-1])), lcl = c(rbind(x$lcl[-n], x$lcl[-n])), ucl = c(rbind(x$ucl[-n], x$ucl[-n])))
    xx <- rbind(xx, data.frame(time = x$time[n], lcl = x$lcl[n], ucl = x$ucl[n]))
    xx$est <- NA
    for (v in setdiff(names(x), names(xx))) xx[[v]] <- x[[v]][1]
    xx
  }
  if (is.null(by)) f(d) else do.call(rbind, lapply(split(d, d[[by]]), f))
}

# Flowchart of the study population (STROBE style)
plot_flow <- function(co) {
  fl <- function(lab) co$flow$n[co$flow$item == lab]; ch <- co$cohorts
  f <- c(fl("Registered patients (ischemic stroke)"), fl("Excluded: in-hospital death (discharge mRS = 6)"), fl("Analyzed patients")); n_dcna <- fl("Excluded: missing discharge mRS")
  r <- range(do.call(c, co$onset_range))                   # onset period across the 3 cohorts combined
  box <- function(x0, x1, y0, y1, lab, bold = FALSE)
    list(annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1, fill = "white", colour = "black", linewidth = 0.35),
         annotate("text", x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = lab, size = 2.6, family = FIG_FONT, lineheight = 1.05,
                  fontface = if (bold) "bold" else "plain"))
  arr <- function(x0, y0, x1, y1) annotate("segment", x = x0, y = y0, xend = x1, yend = y1, linewidth = 0.35,
                                         arrow = arrow(length = unit(1.6, "mm"), type = "closed"))
  mon <- function(d) paste(month.name[as.integer(format(d, "%m"))], format(d, "%Y"))          # English month names (locale-independent)
  p <- ggplot() + xlim(0, 100) + ylim(-4, 100) + theme_void() +
    box(22, 78, 88, 99, sprintf("Patients with ischemic stroke\n(onset %s to %s)\nn = %s", mon(r[1]), mon(r[2]), format(f[1], big.mark = ","))) +
    arr(50, 88, 50, 67) +
    annotate("segment", x = 50, y = 77.5, xend = 60, yend = 77.5, linewidth = 0.35) +
    box(60, 99, 72, 83, if (n_dcna == 0) sprintf("Excluded: in-hospital death\n(discharge mRS 6)\nn = %s", format(f[2], big.mark = ","))
                        else sprintf("Excluded: in-hospital death\n(discharge mRS 6), n = %s;\nmissing discharge mRS, n = %s", format(f[2], big.mark = ","), format(n_dcna, big.mark = ","))) +
    box(22, 78, 56, 67, sprintf("Discharged alive (study population)\nn = %s", format(f[3], big.mark = ",")), bold = TRUE) +
    annotate("segment", x = 50, y = 56, xend = 50, yend = 51, linewidth = 0.35) +
    annotate("segment", x = 17, y = 51, xend = 83, yend = 51, linewidth = 0.35)
  xs <- c(17, 50, 83)
  ch <- ch[order(-sapply(ch, `[[`, "H"))]                    # oldest onset period first (15-year -> 10-year -> 5-year cohort)
  for (i in seq_along(ch)) {
    x <- ch[[i]]; n <- nrow(x$data); ncc <- sum(complete.cases(x$data[, x$tp$col]))
    r <- range(x$data$onset)
    p <- p + arr(xs[i], 51, xs[i], 45) +
      box(xs[i] - 16, xs[i] + 16, 31, 45, sprintf("%d-year cohort\n(onset %s\nto %s)\nn = %s", x$H, mon(r[1]), mon(r[2]), format(n, big.mark = ",")), bold = TRUE) +
      arr(xs[i], 31, xs[i], 25) +
      box(xs[i] - 16, xs[i] + 16, -3, 25, sprintf("Complete mRS follow-up\nto %d years\nn = %s (%.1f%%)\n\nLost to follow-up\n(≥1 missing mRS)\nn = %s (%.1f%%)",
                                                x$H, format(ncc, big.mark = ","), 100 * ncc / n, format(n - ncc, big.mark = ","), 100 * (n - ncc) / n))
  }
  p
}

#===================================================================================
# 5. Excel output
#===================================================================================

write_xlsx_sheets <- function(sheets, file, notes = NULL) {
  wb <- createWorkbook()
  hs <- createStyle(textDecoration = "bold", fgFill = "#E8EEF7", border = "bottom")
  if (!is.null(notes)) {
    addWorksheet(wb, "Notes")
    writeData(wb, "Notes", data.frame(Notes = notes), colNames = FALSE)
    setColWidths(wb, "Notes", 1, 120)
  }
  for (nm in names(sheets)) {
    sn <- substr(nm, 1, 31)
    addWorksheet(wb, sn)
    writeData(wb, sn, sheets[[nm]], headerStyle = hs)
    setColWidths(wb, sn, seq_len(ncol(sheets[[nm]])), "auto")
    freezePane(wb, sn, firstRow = TRUE)
  }
  saveWorkbook(wb, file, overwrite = TRUE)
}
