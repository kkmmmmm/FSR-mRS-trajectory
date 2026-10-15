#-----------------------------------------------------------------------------------
#  02_data_prep.R  —  Data loading, exclusions, and construction of 3 cohorts (5/10/15 years of observation; non-overlapping onset periods)
#-----------------------------------------------------------------------------------

d0 <- load_data()
n_all <- nrow(d0)

# Check target disease: the input data are assumed to be restricted to registered patients with ischemic stroke (principal diagnosis = "脳梗塞", cerebral infarction). Stop if any other is included
if (any(is.na(d0$dx) | d0$dx != "脳梗塞")) stop("Data include patients whose principal diagnosis is not cerebral infarction (or is missing): ", sum(is.na(d0$dx) | d0$dx != "脳梗塞"), " patients")

# Check input values: stop if values used for cohort assignment and integration are missing or out of range (to avoid silent exclusion/misclassification)
chk <- c("Missing onset date" = sum(is.na(d0$onset)),
         "Onset date after extraction date" = sum(d0$onset > EXTRACT_DATE, na.rm = TRUE),
         "mRS (discharge/follow-up) not an integer 0–6" = sum(!is.na(as.matrix(d0[, c("mrs_dc", TP_ALL$col)])) & !(as.matrix(d0[, c("mrs_dc", TP_ALL$col)]) %in% 0:6)),
         "Prestroke mRS not 0–5 or missing" = sum(is.na(d0$pre_mrs) | !(d0$pre_mrs %in% 0:5)),
         "Age/sex missing or out of range" = sum(is.na(d0$age) | d0$age < 0 | d0$age > 120 | !(d0$male %in% 0:1)))
if (any(chk > 0)) stop("Input data contain values that need checking: ", paste(sprintf("%s: %d", names(chk), chk)[chk > 0], collapse = "; "))

# Exclusions: in-hospital death (discharge mRS = 6) and missing discharge mRS (counted separately)
n_excl_dc6  <- sum(d0$mrs_dc == 6, na.rm = TRUE)
n_excl_dcna <- sum(is.na(d0$mrs_dc))
d1 <- d0[!is.na(d0$mrs_dc) & d0$mrs_dc != 6, ]
stopifnot(nrow(d1) == n_all - n_excl_dc6 - n_excl_dcna)

# Enforce death as an absorbing state (set all mRS after the first mRS 6 to 6)
d1 <- enforce_absorbing_death(d1)
n_fix <- attr(d1, "n_absorbing_fix")

cohorts <- lapply(HORIZONS, function(H) make_cohort(d1, H))
names(cohorts) <- sapply(cohorts, `[[`, "key")
if (sum(sapply(cohorts, function(x) nrow(x$data))) != nrow(d1))
  stop("The 3 cohorts do not sum to the analyzed patients (some patients had onset later than the extraction date minus the shortest observation period)")

# Flow data.frame (item = row label, n = number of patients); read by plot_flow() in 01_functions.R and by 09_paper_items.R
flow <- data.frame(
  item = c("Registered patients (ischemic stroke)", "Excluded: in-hospital death (discharge mRS = 6)", "Excluded: missing discharge mRS", "Analyzed patients",
          "  Cells with mRS after death set to 6 (data consistency correction)",
          sapply(cohorts, function(x) sprintf("  %s (%sonset date <= %s)", x$name,
                                              if (is.na(x$lower)) "" else paste(format(x$lower), "< "), format(x$cutoff)))),
  n = c(n_all, n_excl_dc6, n_excl_dcna, nrow(d1), n_fix, sapply(cohorts, function(x) nrow(x$data)))
)
print(flow)

onset_range <- lapply(cohorts, function(x) range(x$data$onset))

dir_create(WORK_DIR)
saveRDS(list(cohorts = cohorts, flow = flow, onset_range = onset_range), file.path(WORK_DIR, "cohorts.rds"))
