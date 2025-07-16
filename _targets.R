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
  workers = 11
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

# Set data table cores to 1

data.table::setDTthreads(1)

subset_comparisons <- bind_rows(
  list(
    comparison = "AD_vs_Others",
    target_diagnosis = "Alzheimer's",
    comparators = list(c("Lewy bodies", "Frontotermporal"))
  ),
  list(
    comparison = "FTD_vs_Others",
    target_diagnosis = "Frontotermporal",
    comparators = list(c("Lewy bodies", "Alzheimer's"))
  ),
  list(
    comparison = "LBD_vs_Others",
    target_diagnosis = "Lewy bodies",
    comparators = list(c("Alzheimer's", "Frontotermporal"))
  )
)

## pipeline
list(
  # declare input data
  c(source("file_targets.R")$value),
  c(source("cohort_targets.R")$value),
  c(source("pet_targets.R")$value),
  c(source("cog_correlation_targets.R")$value),
  # clean data
  tar_target(df, clean_data(joined)),
  # add extra csf (for correlations)
  tar_target(joined_with_csf, add_extra_csf(joined, addf_csf_file)),
  # clean data
  tar_target(df_with_csf, clean_data(joined_with_csf)),
  # prepare ROC data
  tar_target(roc_data_prepared, prepare_roc_data(df_with_csf)),
  tar_target(
    roc_data_prepared_fasting,
    prepare_roc_data(df, with_fasting = "yes")
  ),
  # parallel ROC computations using tar_map
  tar_map(
    values = tibble(
      comparison = c(
        "AD_vs_Control",
        "FTD_vs_Control",
        "LBD_vs_Control",
        "LBD_vs_FTD",
        "AD_vs_Others",
        "LBD_vs_Others",
        "FTD_vs_Others"
      )
    ),
    tar_target(roc_result, run_single_roc(roc_data_prepared, comparison))
  ),
  # parallel ROC computations with fasting
  tar_map(
    values = tibble(
      comparison = c(
        "AD_vs_Control",
        "FTD_vs_Control",
        "LBD_vs_Control",
        "LBD_vs_FTD",
        "AD_vs_Others",
        "LBD_vs_Others",
        "FTD_vs_Others"
      )
    ),
    tar_target(
      roc_result_fasting,
      run_single_roc(roc_data_prepared_fasting, comparison)
    )
  ),
  # combine parallel results
  tar_target(
    roc_results,
    list(
      AD_vs_Control = roc_result_AD_vs_Control,
      FTD_vs_Control = roc_result_FTD_vs_Control,
      LBD_vs_Control = roc_result_LBD_vs_Control,
      LBD_vs_FTD = roc_result_LBD_vs_FTD,
      AD_vs_Others = roc_result_AD_vs_Others,
      LBD_vs_Others = roc_result_LBD_vs_Others,
      FTD_vs_Others = roc_result_FTD_vs_Others
    )
  ),
  tar_target(
    roc_results_fasting,
    list(
      AD_vs_Control = roc_result_fasting_AD_vs_Control,
      FTD_vs_Control = roc_result_fasting_FTD_vs_Control,
      LBD_vs_Control = roc_result_fasting_LBD_vs_Control,
      LBD_vs_FTD = roc_result_fasting_LBD_vs_FTD,
      AD_vs_Others = roc_result_fasting_AD_vs_Others,
      LBD_vs_Others = roc_result_fasting_LBD_vs_Others,
      FTD_vs_Others = roc_result_fasting_FTD_vs_Others
    )
  ),
  # get combined ROC curve
  tar_target(combined_roc, get_combined_roc(roc_results)),
  # get combined ROC curve
  tar_target(combined_roc_fasting, get_combined_roc(roc_results_fasting)),
  # # get sex specific ROC curves
  # tar_map(
  #   values = tibble(
  #     comparison = c(
  #       "AD_vs_Control",
  #       "FTD_vs_Control",
  #       "LBD_vs_Control",
  #       "LBD_vs_FTD",
  #       "AD_vs_Others",
  #       "LBD_vs_Others",
  #       "FTD_vs_Others"
  #     )
  #   ),
  #   tar_target(roc_result_men, run_single_roc(roc_data, comparison))
  # ),
  # get sex specific ROC curves
  # tar_target(subtypes_control, subtypes_vs_control(roc_results)),
  # CSF AB vs plasma markers
  tar_target(csf_out, csf_vs_markers(df)),
  # CSF marker and plasma marker partial rank order correlations
  tar_target(csf_rank_cors, csf_rank_corr(df_with_csf)),
  # PET AB vs PTAU-217
  # tar_target(pet_ptau, pet_vs_ptau(df)),
  # Variable importance overall
  # tar_target(vimp_full, vimp_overall(df)),
  # Variable importance for females
  # tar_target(vimp_females, vimp_by_sex(df, "female")),
  # Variable importance for males
  # tar_target(vimp_males, vimp_by_sex(df, "male")),
  # Parallel subset analysis for full cohort
  tar_map(
    values = subset_comparisons,
    names = comparison,
    tar_target(
      subset_result,
      run_single_subset(df, target_diagnosis, comparators)
    )
  ),
  # Combine parallel results for backward compatibility
  tar_target(
    subset_data,
    list(
      ad = subset_result_AD_vs_Others,
      ftd = subset_result_FTD_vs_Others,
      lbd = subset_result_LBD_vs_Others
    )
  ),
  tar_target(subset_plots, all_subset_plots(subset_data)),
  # Parallel subset analysis for men
  tar_map(
    values = subset_comparisons,
    names = comparison,
    tar_target(
      subset_result_men,
      run_single_subset(
        df,
        target_diagnosis,
        comparators,
        sex_strat = "male"
      )
    )
  ),

  # Combine parallel results for men
  tar_target(
    subset_data_men,
    list(
      ad = subset_result_men_AD_vs_Others,
      ftd = subset_result_men_FTD_vs_Others,
      lbd = subset_result_men_LBD_vs_Others
    )
  ),
  tar_target(
    subset_plots_men,
    all_subset_plots(subset_data_men, extra_title = " - Males")
  ),

  # Parallel subset analysis for women
  tar_map(
    values = subset_comparisons,
    names = comparison,
    tar_target(
      subset_result_women,
      run_single_subset(
        df,
        target_diagnosis,
        comparators,
        sex_strat = "female"
      )
    )
  ),

  # Combine parallel results for women
  tar_target(
    subset_data_women,
    list(
      ad = subset_result_women_AD_vs_Others,
      ftd = subset_result_women_FTD_vs_Others,
      lbd = subset_result_women_LBD_vs_Others
    )
  ),
  tar_target(
    subset_plots_women,
    all_subset_plots(subset_data_women, extra_title = " - Females")
  ),

  # Parallel subset analysis for CDR cohort
  tar_map(
    values = subset_comparisons,
    names = comparison,
    tar_target(
      subset_result_cdr,
      run_single_subset(df, target_diagnosis, comparators, use_cdr = TRUE)
    )
  ),

  # Combine parallel results for CDR cohort
  tar_target(
    subset_data_cdr,
    list(
      ad = subset_result_cdr_AD_vs_Others,
      ftd = subset_result_cdr_FTD_vs_Others,
      lbd = subset_result_cdr_LBD_vs_Others
    )
  ),
  tar_target(
    subset_plots_cdr,
    all_subset_plots(subset_data_cdr, use_cdr = TRUE)
  ),
  # Quarto document for results
  tar_quarto(plots_and_corrs, path = "./R/plots_and_corrs.qmd", quiet = FALSE)
)
