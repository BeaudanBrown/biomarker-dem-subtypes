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

# All TabPFN fits share one device; extra workers are an explicit opt-in.
gpu_worker_setting <- Sys.getenv("ADDF_TABPFN_GPU_WORKERS", "1")
gpu_workers <- suppressWarnings(as.integer(gpu_worker_setting))
if (
  is.na(gpu_workers) ||
    gpu_workers < 1L ||
    as.character(gpu_workers) != gpu_worker_setting
) {
  stop("ADDF_TABPFN_GPU_WORKERS must be a positive integer.", call. = FALSE)
}

# Set up crew controllers. Keep logs for progress/failure diagnosis.
crew_controller_global <- crew_controller_local(
  name = "cpu",
  options_local = crew_options_local(log_directory = "./logs"),
  workers = max(1L, future::availableCores() - 2L)
)
crew_controller_gpu <- crew_controller_local(
  name = "gpu",
  options_local = crew_options_local(log_directory = "./logs"),
  workers = gpu_workers
)
crew_controllers <- crew_controller_group(
  cpu = crew_controller_global,
  gpu = crew_controller_gpu
)

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
  controller = crew_controllers,
  resources = tar_resources(crew = tar_resources_crew(controller = "cpu")),
  error = "continue"
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
  c(source("tdp_targets.R")$value),
  c(source("tabpfn_targets.R")$value),
  if (identical(Sys.getenv("ADDF_SL_WEIGHT_PILOT"), "true")) {
    source("tabpfn_sl_weights_targets.R")$value
  } else {
    list()
  },
  c(roc_targets),
  c(subset_targets),
  # clean data
  tar_target(df, clean_data(joined)),
  # add extra csf (for correlations)
  tar_target(joined_with_csf, add_extra_csf(joined, addf_csf_file)),
  # clean data
  tar_target(df_with_csf, clean_data(joined_with_csf)),
  tar_target(descriptives, show_descriptives(df)),
  tar_target(cognitive_summary_table, cognitive_summary(df)),
  tar_target(hispanic_summary_table, hispanic_summary(df)),
  # CSF marker and plasma marker partial rank order correlations
  tar_target(csf_rank_cors, csf_rank_corr(df_with_csf)),
  # Render quarto doc
  tar_quarto(plots_and_corrs, path = "./R/plots_and_corrs.qmd", quiet = FALSE)
)
