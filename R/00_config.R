#-----------------------------------------------------------------------------------
#  00_config.R  —  Analysis settings (paths, analysis parameters, variable-name mapping)
#
#  * These scripts contain no individual-level data.
#     Place the data file (data.xlsx, shared separately) in the project's data/ folder and the scripts run as is.
#     To use another location, set the environment variable MRS_DATA_PATH.
#-----------------------------------------------------------------------------------

# ---- Paths -----------------------------------------------------------------------
# Project root (one level above the folder containing run_all.R); set by run_all.R.
if (!exists("PROJECT_ROOT")) PROJECT_ROOT <- normalizePath(".", winslash = "/")

# Input data: default is <project>/data/data.xlsx (place the separately shared file under this name)
DATA_PATH  <- Sys.getenv("MRS_DATA_PATH",
                         unset = file.path(PROJECT_ROOT, "data", "data.xlsx"))
if (!file.exists(DATA_PATH)) stop("Data file not found: ", DATA_PATH,
                                 "
  Place the data file at data/data.xlsx in the project folder, or set the environment variable MRS_DATA_PATH.")
DATA_SHEET <- 1

# Japanese life tables (CSV reformatted from MHLW published values; used for expected survival of the general population)
LIFETABLE_PATH <- file.path(PROJECT_ROOT, "reference", "lifetable", "japan_lifetable_qx.csv")

OUT_DIR  <- Sys.getenv("MRS_OUT_DIR",  unset = file.path(PROJECT_ROOT, "output"))   # shareable aggregate results only
WORK_DIR <- Sys.getenv("MRS_WORK_DIR", unset = file.path(PROJECT_ROOT, "work"))     # individual-level intermediate files (do not share)

# ---- Analysis parameters ---------------------------------------------------------
EXTRACT_DATE <- as.Date("2025-07-09")   # follow-up data extraction date
HORIZONS     <- c(5, 10, 15)            # 3 cohorts defined by observable follow-up (years) (non-overlapping onset periods; make_cohort in 01_functions.R)
AGE_LABELS   <- c("<65", "65-74", "75-84", ">=85")   # age group at onset

M_IMP        <- 100                     # number of imputations (B >= 100 recommended)
MAXIT_MICE   <- 10                      # number of MICE (mice) iterations
SEED         <- 20250709                # random seed (for reproducibility)
AUG_WEIGHT   <- 5                       # total weight of pseudo-observations against complete separation in the death model (in patients; mice's logreg uses no. of predictors + 1)
N_CORES      <- max(1, min(16, parallel::detectCores() - 2))

# Sensitivity parameters for MNAR (delta-adjusted pattern-mixture model; the imputation model is the two-stage model)
#   The same delta is added to the log-odds of death and to the conditional log-odds of a "worse mRS" if alive, for patients with missing data
#   delta > 0: outcomes of missing patients worse than observed / delta < 0: better / delta = 0: MAR
DELTAS_MNAR  <- c(-1.0, -0.5, 0.5, 1.0)

# ---- Variable-name mapping (raw data column names -> analysis names) ---------------
#   Raw column names are in Japanese and must match the input data exactly; English glosses are given in comments.
VAR_MAP <- c(
  "【患者情報】臨床ID"       = "id",         # [Patient info] clinical ID
  "【診断情報】主病名"       = "dx",         # [Diagnosis] principal diagnosis
  "【診断情報】症状持続時間" = "duration",   # [Diagnosis] symptom duration
  "【入院基本情報】発症日"   = "onset",      # [Admission] onset date
  "【治療情報】TPA"                     = "tpa",   # [Treatment] intravenous tPA -- acute treatment (patient characteristics only; not in the imputation model)
  "【治療情報】超急性期血栓回収術の有無" = "evt",   # [Treatment] hyperacute endovascular thrombectomy (yes/no)
  "age"                      = "age",
  "sex(1=man)"               = "male",
  "発症前mRankin"            = "pre_mrs",    # prestroke modified Rankin Scale
  "入院時NIHSS"              = "nihss",      # NIHSS on admission
  "byokei_new"               = "toast",
  "rfht."                    = "ht",
  "rfdm."                    = "dm",
  "rfdl."                    = "dl",
  "rfaf."                    = "af",
  "rfsmk."                   = "smk",
  "rfdrk."                   = "drk",
  "rfckd."                   = "ckd",
  "bmi."                     = "bmi",
  "脳卒中既往"               = "prior_stroke",   # previous stroke (no missing values)
  "退院時mRankin"            = "mrs_dc",     # discharge modified Rankin Scale
  "mrs3mn"                   = "mrs_3m",
  "mrs6mn"                   = "mrs_6m",
  setNames(paste0("mrs_", 1:20, "y"), paste0("mrs", 1:20, "yn"))
)

# Follow-up time points (years) and column names
TP_ALL <- data.frame(col  = c("mrs_3m", "mrs_6m", paste0("mrs_", 1:20, "y")),
                     year = c(0.25, 0.5, 1:20),
                     lab  = c("3M", "6M", paste0(1:20, "Y")),
                     stringsAsFactors = FALSE)

# Covariates (display labels)
COV_LABELS <- c(
  age = "Age", male = "Sex (male)", pre_mrs = "Prestroke mRS", nihss = "NIHSS on admission",
  toast = "Stroke subtype (TOAST)", ht = "Hypertension", dm = "Diabetes mellitus", dl = "Dyslipidemia",
  af = "Atrial fibrillation", smk = "Smoking", drk = "Alcohol drinking", ckd = "Chronic kidney disease",
  bmi = "BMI", prior_stroke = "Previous stroke", mrs_dc = "Discharge mRS"
)
COVARIATES <- names(COV_LABELS)
# Acute reperfusion therapy (used only to describe patient characteristics in Table 1; not used in the imputation model or estimation)
TREATMENTS <- c(tpa = "Intravenous thrombolysis (tPA)", evt = "Endovascular thrombectomy")
