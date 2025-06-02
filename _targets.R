library(targets)
library(tarchetypes)
library(future)

dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
nacc_dir <- Sys.getenv("NACC_DIR")
cache_dir <- Sys.getenv("CACHE_DIR")

# set target configs
tar_config_set(store = cache_dir)

plan(multicore)

# Set target options:
tar_option_set(
  packages = c(
    "languageserver",
    "dotenv",
    "readxl",
    "performance",
    "lme4",
    "data.table",
    "SuperLearner",
    "gtsummary",
    "ggthemes",
    "vimp",
    "xgboost",
    "glmnet",
    "ranger",
    "randomForest",
    "earth",
    "gam",
    "pROC",
    "origami",
    "future",
    "future.apply",
    "logistf",
    "rms",
    "psych",
    "bayesplot",
    "knitr",
    "cardx",
    "corrr",
    "quarto",
    "arm",
    "cowplot",
    "nnet",
    "cvAUC",
    "gbm",
    "e1071",
    "class",
    "patchwork",
    "tidyverse",
    "parallel",
    "bartMachine"
  ),
  format = "qs"
)

# Run the R scripts in the R/ folder
tar_source()

# Set data table cores to 1

data.table::setDTthreads(1)

## pipeline
list(
  # declare input data
  tar_target(
    pheno_file,
    file.path(data_dir, Sys.getenv("PHENO_FILE")),
    format = "file"
  ),
  tar_target(
    gilipin_file,
    file.path(data_dir, Sys.getenv("GILIPIN_FILE")),
    format = "file"
  ),
  tar_target(
    all_cases_file,
    file.path(data_dir, Sys.getenv("ALL_CASES_FILE")),
    format = "file"
  ),
  tar_target(
    addf_file,
    file.path(data_dir, Sys.getenv("ADDF_FILE")),
    format = "file"
  ),
  tar_target(
    texas_file,
    file.path(data_dir, Sys.getenv("TEXAL_PRELIM_FILE")),
    format = "file"
  ),
  tar_target(
    texas_file_update,
    file.path(data_dir, Sys.getenv("TEXAL_UPDATE_FILE")),
    format = "file"
  ),
  tar_target(
    addf_cd14_file,
    file.path(data_dir, Sys.getenv("ADDF_CD14_FILE")),
    format = "file"
  ),
  tar_target(
    addf_csf_file,
    file.path(data_dir, Sys.getenv("ADDF_CSF_FILE")),
    format = "file"
  ),
  tar_target(
    nacc_file,
    file.path(nacc_dir, Sys.getenv("NACC_FILE")),
    format = "file"
  ),
  tar_target(
    wash_cogs_file,
    file.path(data_dir, Sys.getenv("NEW_WASH_COGS")),
    format = "file"
  ),
  # merge input data
  tar_target(
    joined,
    merge_datafiles(
      pheno_file,
      gilipin_file,
      all_cases_file,
      addf_file,
      texas_file,
      texas_file_update,
      addf_cd14_file,
      wash_cogs_file,
      addf_csf_file,
      nacc_file
    )
  ),
  # clean data
  tar_target(df, clean_data(joined)),
  # get all ROC curves
  tar_target(roc_results, all_rocs(df)),
  # get all ROC curves, including fasting status as covariate
  tar_target(roc_results_fasting, all_rocs(df, with_fasting = "yes")),
  # get combined ROC curve
  tar_target(combined_roc, get_combined_roc(roc_results)),
  # get combined ROC curve
  tar_target(combined_roc_fasting, get_combined_roc(roc_results_fasting)),
  # get sex specific ROC curves
  tar_target(sex_specific_rocs, rocs_by_sex(df)),
  # get sex specific ROC curves
  tar_target(subtypes_control, subtypes_vs_control(roc_results)),
  # MMSE correlation table
  tar_target(mmse_out, mmse_table(df)),
  # CDR correlation table
  tar_target(cdr_out, cdr_vs_markers(df)),
  # CSF AB vs plasma markers
  tar_target(csf_out, csf_vs_markers(df)),
  # PET AB vs PTAU-217
  tar_target(pet_ptau, pet_vs_ptau(df)),
  # Variable importance overall
  tar_target(vimp_full, vimp_overall(df)),
  # Variable importance for females
  tar_target(vimp_females, vimp_by_sex(df, "female")),
  # Variable importance for males
  tar_target(vimp_males, vimp_by_sex(df, "male")),
  # Marker subset for full cohort
  tar_target(subset_data, all_subset_data(df)),
  tar_target(subset_plots, all_subset_plots(subset_data)),
  # Marker subset for men
  tar_target(subset_data_men, all_subset_data(df, sex_strat = "male")),
  tar_target(
    subset_plots_men,
    all_subset_plots(subset_data_men, extra_title = " - Males")
  ),
  # Marker subset for women
  tar_target(subset_data_women, all_subset_data(df, sex_strat = "female")),
  tar_target(
    subset_plots_women,
    all_subset_plots(subset_data_women, extra_title = " - Females")
  ),

  # Marker subset for full cohort with cdr
  tar_target(subset_data_cdr, all_subset_data(df, use_cdr = TRUE)),
  tar_target(
    subset_plots_cdr,
    all_subset_plots(subset_data_cdr, use_cdr = TRUE)
  ),
  # Quarto document for results
  tar_quarto(plots_and_corrs, path = "./R/plots_and_corrs.qmd", quiet = FALSE)
)
