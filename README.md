# FSR-mRS-trajectory

R code for the analysis of 15-year functional trajectories (modified Rankin Scale, mRS) after ischemic stroke by discharge mRS in the Fukuoka Stroke Registry (FSR; UMIN000000800).

> Manuscript under review. Citation and DOI will be added upon publication.

**This repository contains code only. No individual-level patient data are included.** The data supporting the findings are available from the corresponding author upon reasonable request from researchers trained in participant confidentiality.

## Overview

- Three cohorts are defined by period of stroke onset with 15, 10, and 5 years of potential follow-up (G1, G2, G3).
- mRS at each follow-up assessment (3 and 6 months, then yearly) is multiply imputed (m = 100) with MICE, using a two-stage model for mRS: (1) death (discrete-time logistic hazard) and (2) mRS 0–5 among survivors (proportional odds model). Death is treated as absorbing.
- Missing-not-at-random (MNAR) sensitivity analysis: delta-adjusted pattern-mixture model (delta = −1.0, −0.5, 0.5, 1.0).
- Outcomes: mRS distribution, survival, and alive-and-independent (mRS 0–2) proportions at each time point; years alive, years independent, and years dependent over 15 years; survival relative to the general Japanese population (Ederer I method).
- Comparison methods: Kaplan–Meier and complete-case analysis.

## Repository structure

```
R/
  run_all.R              Runs the whole pipeline
  00_config.R            Paths, analysis parameters, variable-name mapping
  01_functions.R         Functions (data preparation, imputation, estimation, figures)
  02_data_prep.R         Data loading, exclusions, construction of the 3 cohorts
  03_imputation.R        Multiple imputation (MAR) and delta-adjusted pattern-mixture model (MNAR)
  04_estimates.R         Estimates at each time point; imputation-model diagnostics
  05_prestroke.R         Sensitivity analysis restricted to prestroke mRS 0-2
  06_years_lived.R       Years alive / independent / dependent over 15 years; expected years
  07_model_checks.R      Reproducibility and imputation-model checks
  08_threshold_check.R   Ordering of proportional-odds thresholds
  09_paper_items.R       Tables and figures numbered as in the manuscript
  renv.lock              R package versions
reference/lifetable/
  japan_lifetable_qx.csv Age-, sex-, and year-specific death probabilities (qx) from the Japanese
                         abridged/complete life tables (Ministry of Health, Labour and Welfare)
  urls.tsv               Source URLs of the original life-table files
```

## Requirements

- R 4.6.1 (tested). Packages: readxl, mice, MASS, survival, dplyr, tidyr, ggplot2, openxlsx, patchwork (exact versions in `R/renv.lock`; restore with `renv::restore(lockfile = "R/renv.lock")`).
- Runtime: about 2–4 hours with 16 cores (multiple imputation is the main cost). The number of cores is set by `N_CORES` in `00_config.R`.

## Usage

1. Place the data file at `data/data.xlsx` (or set the environment variable `MRS_DATA_PATH`). The required columns are listed in `VAR_MAP` in `R/00_config.R`; the raw column names are in Japanese and are glossed in English in the comments.
2. From the repository root, run

   ```
   Rscript R/run_all.R
   ```

3. Aggregated results are written to `output/` (tables and figures numbered as in the manuscript in `output/00_paper_figures_tables/`). Individual-level intermediate files are written to `work/` and must not be shared. Output and work locations can be changed with `MRS_OUT_DIR` and `MRS_WORK_DIR`; a subset of steps can be run with, e.g., `MRS_STEPS="04,06,09"`.

Random seeds are fixed (`SEED` in `00_config.R`), so rerunning the pipeline on the same data reproduces the reported estimates exactly.

`data/`, `work/`, and `output/` are excluded by `.gitignore`.

## License

MIT License (see `LICENSE`). The life-table values are derived from publicly released statistics of the Ministry of Health, Labour and Welfare, Japan.
