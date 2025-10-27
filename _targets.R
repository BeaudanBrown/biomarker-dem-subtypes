library(targets)
library(tarchetypes)
library(future)
library(tidyverse)
library(crew)

dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
cache_dir <- Sys.getenv("CACHE_DIR")

# set target configs
tar_config_set(store = cache_dir)

plan(multicore)

# Set up global crew controller
unlink("./logs/*", recursive = FALSE)
crew_controller_global <- crew_controller_local(
  options_local = crew_options_local(log_directory = "./logs"),
  workers = future::availableCores() - 1
)
tar_option_set(error = "null")

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
    "dbarts",
    "tmle",
    "crew"
  ),
  format = "qs",
  controller = crew_controller_global
)

# Run the R scripts in the R/ folder
tar_source()

Sys.setenv(R_DATATABLE_NUM_THREADS = 1)
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)

source("roc_targets.R")
source("subset_targets.R")

## pipeline
list(
  # declare input data
  c(source("file_targets.R")$value),
  c(source("cohort_targets.R")$value),
  c(source("pet_targets.R")$value),
  c(source("csf_targets.R")$value),
  c(source("cog_correlation_targets.R")$value),
  c(roc_targets),
  c(subset_targets),
  # clean data
  tar_target(df, clean_data(joined)),
  # add extra csf (for correlations)
  tar_target(joined_with_csf, add_extra_csf(joined, addf_csf_file)),
  # clean data
  tar_target(df_with_csf, clean_data(joined_with_csf)),
  tar_target(demo_table, demos(as.data.table(df))),
  tar_target(subtypes_control, subtypes_vs_control(roc_results)),
  # CSF marker and plasma marker partial rank order correlations
  tar_target(csf_rank_cors, csf_rank_corr(df_with_csf)),
  # Variable importance overall
  # tar_target(vimp_full, vimp_overall(df)),
  # Variable importance for females
  # tar_target(vimp_females, vimp_by_sex(df, "female")),
  # Variable importance for males
  # tar_target(vimp_males, vimp_by_sex(df, "male")),
  tar_quarto(plots_and_corrs, path = "./R/plots_and_corrs.qmd", quiet = FALSE)
)
