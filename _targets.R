library(targets)

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
    "tidyverse"
  )
)

# Run the R scripts in the R/ folder
tar_source()

dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")

## pipeline
list(
  # declare input data
  tar_target(pheno_file, file.path(data_dir, Sys.getenv("PHENO_FILE")), format = "file"),
  tar_target(gilipin_file, file.path(data_dir, Sys.getenv("GILIPIN_FILE")), format = "file"),
  tar_target(all_cases_file, file.path(data_dir, Sys.getenv("ALL_CASES_FILE")), format = "file"),
  tar_target(addf_file, file.path(data_dir, Sys.getenv("ADDF_FILE")), format = "file"),
  tar_target(texas_file, file.path(data_dir, Sys.getenv("TEXAL_PRELIM_FILE")), format = "file"),
  tar_target(texas_file_update,
    file.path(data_dir, Sys.getenv("TEXAL_UPDATE_FILE")),
    format = "file"),
  tar_target(addf_cd14_file, file.path(data_dir, Sys.getenv("ADDF_CD14_FILE")), format = "file"),
  # merge input data
  tar_target(joined, merge_datafiles(
    pheno_file, gilipin_file, all_cases_file, addf_file, texas_file, texas_file_update, addf_cd14_file
  )),
  # clean data
  tar_target(df, clean_data(joined)),
  # get all ROC curves
  tar_target(roc_results, all_rocs(df)),
  # get combined ROC curve
  tar_target(combined_roc, get_combined_roc(roc_results)),
  # get sex specific ROC curves
  tar_target(sex_specific_rocs, rocs_by_sex(df)),
  # get sex specific ROC curves
  tar_target(subtypes_control, subtypes_vs_control(roc_results)),
  # MMSE correlation table
  tar_target(mmse_out, mmse_table(df)),
  # CSF AB vs plasma markers
  tar_target(csf_out, csf_vs_markers(df))
)
