# ADDF Plasma Biomarker Analysis

This repository implements an R targets pipeline to merge ADDF phenotype, plasma biomarkers, PET, and CSF data, then performs classification and correlation analyses with reproducible caching.

Key analyses
- Data integration: merges Excel/CSV sources; cleans, harmonizes diagnosis/race/site; derives AB42/AB40 ratio; truncates extreme biomarker values.
- ROC classification: SuperLearner AUROC for each subtype vs control and vs other dementias; CDR- and fasting-adjusted variants; sex-stratified results.
- Correlations: adjusted associations of plasma markers with PET measures (z-score, centiloid, SUVR), cognitive scores (MMSE, CDR), and CSF markers; CSF–plasma rank correlations.
- Minimal biomarker subsets: backward search to find small panels achieving near-full AUC with confidence intervals.
- Report: a Quarto HTML summary with plots and tables.

Inputs
Set environment variables to point to your data files (relative to `DATA_DIR`). See `env-sample` for names and expected files:
- `DATA_DIR`: base directory containing all inputs
- `CACHE_DIR`: directory for the targets store (cache)
- File variables such as `ADDF_FILE`, `PHENO_FILE`, `GILIPIN_FILE`, `ALL_CASES_FILE`, `TEXAL_PRELIM_FILE`, `TEXAL_UPDATE_FILE`, `ADDF_CD14_FILE`, `ADDF_CSF_FILE`, `NACC_FILE`, `NEW_WASH_COGS`

Minimum to run
- Provide all required input files and set the env vars above (copy `env-sample` to `.env` and edit, or export them in your shell).
- Use the provided Nix dev shell (recommended) or install R + packages manually.

Quick start (Nix)
1) Enter the dev shell (provides R, required packages, Quarto):
   - `nix develop`
2) Ensure environment variables are set (e.g., `cp env-sample .env` and edit paths).
3) Run the pipeline in R:
   - `R -q`
   - `targets::tar_make()`

Quick start (base R)
- R ≥ 4.3 with packages: `targets`, `tarchetypes`, `crew`, `tidyverse`, `data.table`, `readxl`, `dotenv`, `SuperLearner`, `origami`, `cvAUC`, `pROC`, `glmnet`, `ranger`, `xgboost`, `earth`, `gam`, `arm`, `rms`, `bayesplot`, `patchwork`, `gtsummary`, `kableExtra`, `visNetwork`, `qs`, `tarchetypes`, `quarto`.
- Set the same environment variables as in `env-sample`, then run `targets::tar_make()`.

Outputs
- Figures: ROC and subset-path PNGs in `plots/`.
- Report: `R/plots_and_corrs.html` (rendered by the `tar_quarto` target).
- Cache: targets store at `CACHE_DIR`.
