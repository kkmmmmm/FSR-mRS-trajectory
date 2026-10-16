#-----------------------------------------------------------------------------------
#  run_all.R  —  Run the full analysis
#
#  Usage (tested with R 4.6.1):
#    1) Place your data file at data/data.xlsx in the project folder
#       (or set the environment variable MRS_DATA_PATH to its path).
#       Individual-level data are not distributed with this repository.
#    2) With the project root as the working directory, run
#       Rscript R/run_all.R
#    3) Results are written to output/. Manuscript figures and tables are collected,
#       under their manuscript numbers, in output/00_paper_figures_tables/
#
#  Approximate run time: ~4 h with 16 parallel workers (03 multiple imputation: ~2.5 h, 07: ~1 h, others: ~30 min)
#-----------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grep("^--file=", args)])
PROJECT_ROOT <- if (length(f)) normalizePath(file.path(dirname(f), ".."), winslash = "/") else normalizePath(".", winslash = "/")
setwd(PROJECT_ROOT)

source("R/00_config.R", encoding = "UTF-8")
source("R/01_functions.R", encoding = "UTF-8")

steps <- c("02_data_prep.R",        # Data loading, exclusions, and construction of the 3 onset-period cohorts (G1, G2, G3)
           "03_imputation.R",       # Multiple imputation (MAR) and delta-adjusted pattern-mixture models (~2.5 h)
           "04_estimates.R",        # mRS distribution, survival, and independent survival at each time point; imputation model diagnostics
           "05_prestroke.R",        # Estimates restricted to prestroke mRS 0–2 (Table S4)
           "06_years_lived.R",      # Years alive and years independent over 15 years (primary outcome); expected years of life
           "07_model_checks.R",     # Supplemental Methods values: difference from delta = 0 regeneration, imputation model comparison (~1 h)
           "08_threshold_check.R",  # Supplemental Methods values: probability of threshold order reversal
           "09_paper_items.R")      # Manuscript figures and tables under their manuscript numbers (output/00_paper_figures_tables/)
only <- Sys.getenv("MRS_STEPS", "")          # e.g., MRS_STEPS="04,05,06,09" runs only those steps
if (nzchar(only)) steps <- steps[substr(steps, 1, 2) %in% strsplit(only, ",")[[1]]]

for (s in steps) {
  message("==== ", s, "  (", format(Sys.time()), ")")
  source(file.path("R", s), encoding = "UTF-8")
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("==== Done ", format(Sys.time()))
