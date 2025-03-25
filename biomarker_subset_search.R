library(tidyverse)
library(SuperLearner)
library(gtsummary)
library(ggthemes)
library(xgboost)
library(glmnet)
library(ranger)
library(randomForest)
library(earth)
library(gam)
library(pROC)
library(origami)
library(future)
library(future.apply)
library(parallel)
library(logistf)
library(knitr)
library(cvAUC)
library(patchwork)
source("learners.R")
source("helper_functions.R")
source("subset_search_helper_functions.R")

dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
plan(multicore, workers = detectCores())
# plan(sequential)

plot_auc_steps <- function(data, outcome, extra_title = "") {
  library(ggplot2)
  library(dplyr)

  # Order by step
  data <- data %>% arrange(step)

  # Get the reference AUC value (step 0 minus 0.03)
  ref_auc <- data %>%
    filter(step == 0) %>%
    pull(auc) - 0.03

  labels <- map(data$model_vars, function(vars) {
    return(paste(vars, collapse = "\n"))
  })

  # Create the plot
  p <- ggplot(data, aes(x = step, y = auc)) +
    geom_line(size = 1, color = "steelblue") +
    geom_pointrange(aes(ymin = cil, ymax = ciu), color = "darkblue") +
    labs(
      title = paste0(outcome, " AUC Marker Subset Path", extra_title),
      x = "Markers Remaining",
      y = "AUC"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_x_continuous(
      breaks = data$step,
      labels = labels,
      minor_breaks = NULL
    ) +
    scale_y_continuous(limits = c(0.5, 1), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))

  return(p)
}

build_path <- function(full_path, ref_auc) {
  all_vars <- c(
    "mean_nfl",
    "mean_ykl",
    "mean_gfap",
    "mean_ab40",
    "mean_ab42",
    "mean_tdp",
    "mean_ptau181",
    "mean_ptau217"
  )
  # Create a named vector for variable mapping
  var_mapping <- c(
    "mean_nfl" = "NfL",
    "mean_ykl" = "YKL-40",
    "mean_gfap" = "GFAP",
    "mean_ab40" = "AB40",
    "mean_ab42" = "AB42",
    "mean_tdp" = "TDP-43",
    "mean_ptau181" = "pTau-181",
    "mean_ptau217" = "pTau-217"
  )

  # Function to get remaining variables at each step
  get_remaining_vars <- function(step, removed_vars) {
    if (step == 0) {
      return(all_vars)
    } else {
      return(setdiff(all_vars, removed_vars[1:step]))
    }
  }

  best_path <- map(full_path, function(options) {
    return(options[which.max(options$auc), ])
  })

  # Extract removed variables in order
  removed_vars_in_order <- map_chr(best_path, ~.x$removed_var)

  result1 <- best_path %>%
    map2_dfr(seq_along(.), function(tibble, index) {
      tibble %>% mutate(step = index)
    }) |>
    add_row(
      removed_var = "Full",
      auc = ref_auc$cvAUC,
      cil = ref_auc$ci[[1]],
      ciu = ref_auc$ci[[2]],
      step = 0) |>
    arrange(step) |>
    mutate(
      raw_removed_var = removed_var,
      removed_var = case_when(
        removed_var %in% names(var_mapping) ~ var_mapping[removed_var],
        TRUE ~ "Full"
      )
    ) |>
    mutate(
      # Create model_vars column with list of remaining variables at each step
      model_vars = map(step, function(s) {
        remaining_raw <- get_remaining_vars(s, removed_vars_in_order)
        remaining_display <- var_mapping[remaining_raw]
        return(remaining_display)
      })
    ) |>
    select(-raw_removed_var)  # Remove the temporary column

  return(result1)
}

read_search_data <- function() {
  df <- read_csv(file.path(data_dir, "merged_data/merged_all_vars_11_02_2025.csv"))

  # Merge diagnosis variables between Texas and Washington data
  df <-
    df |>
      mutate(Diagnosis_combined = case_when(
        Diagnosis == "normal" ~ "CO",
        Diagnosis == "AD" ~ "AD",
        Diagnosis == "FTD" ~ "FTD",
        Diagnosis == "AD/VD" ~ "AD",
        Diagnosis == "MCI" ~ "MCI",
        Diagnosis == "DLB" ~ "DLB",
        Diagnosis == "DLB, VD" ~ "DLB",
        Diagnosis == "Parkinsons" ~ "DLB",
        Final_Status == "AD" ~ "AD",
        Final_Status == "DLB" ~ "DLB",
        Final_Status == "FTD" ~ "FTD",
        Final_Status == "CO" ~ "CO",
        TRUE ~ "other"
      ))

  # Merge race
  df <- df |>
    mutate(race_combined = case_when(
      Race == "White" ~ "White",
      Race == "African American" ~ "African American",
      Race == "Asian" ~ "Asian",
      Race == "American Indian or Alaska Native" ~ "American Indian/Alaska Native",
      RACE == "White" ~ "White",
      RACE == "Black/African American" ~ "African American",
      RACE == "Asian" ~ "Asian",
      RACE == "American Indian/Alaska Native" ~ "American Indian/Alaska Native",
      is.na(Race) & is.na(RACE) ~ "Other",
      TRUE ~ "Other"
    ))

  # add site variable
  df$Site <- ifelse(!is.na(df$SEX), "Texas", "Washington")

  # exclude MCI and other
  df <- df |>
    filter(!Diagnosis_combined %in% c("MCI", "other") & !is.na(Diagnosis_combined))

  # code diagnosis as factor and set control as reference
  df <- df |>
    mutate(Diagnosis_combined = as.factor(fct_recode(Diagnosis_combined,
      "Alzheimer's" = "AD", "Control" = "CO", "Lewy bodies" = "DLB",
      "Frontotemporal" = "FTD"
    ))) |>
    mutate(Diagnosis_combined = fct_relevel(Diagnosis_combined, "Control"))

  # set sex to factor
  df$sex_combined <- as.factor(df$sex_combined)
  df$female <- ifelse(df$sex_combined == "Female", 1, 0)

  # select relevant variables
  df <-
    df |> select(
      Diagnosis_combined, age_combined, female, starts_with("mean_")
    )

  # remove CD14
  # TODO: Sensitivity
  df <- df |> select(-mean_elisa)

  # remove missing biomarkers
  df <- drop_na(df)

  # Truncate biomarkers at the 99th percentile
  df <- df |>
    mutate(across(starts_with("mean_"),
      ~ ifelse(.x > quantile(.x, 0.99), quantile(.x, 0.99), .x)))

  df <- df |>
    mutate(across(starts_with("mean_"),
      ~ ifelse(.x < quantile(.x, 0.01), quantile(.x, 0.01), .x)))

  # remove controls
  df <- df |> filter(!Diagnosis_combined == "Control")
}

df <- read_search_data()

# ad_full <- marker_subset_full(
#   df,
#   "Alzheimer's",
#   c("Lewy bodies", "Frontotemporal")
# )
# lbd_full <- marker_subset_full(
#   df,
#   "Lewy bodies",
#   c("Alzheimer's", "Frontotemporal")
# )
# ftd_full <- marker_subset_full(
#   df,
#   "Frontotemporal",
#   c("Alzheimer's", "Lewy bodies")
# )

ad_full <- read_rds("AD_full_ci.rds")
ftd_full <- read_rds("FTD_full_ci.rds")
lbd_full <- read_rds("LBD_full.rds")

ad_women_path <- build_path(ad_women$path, ad_women$reference_auc)
ad_women_plot <- plot_auc_steps(ad_women_path, "Alzheimer's")

lbd_path <- build_path(lbd_full$path, lbd_full$reference_auc)
lbd_plot <- plot_auc_steps(lbd_path, "Lewy bodies")

ftd_path <- build_path(ftd_full$path, ftd_full$reference_auc)
ftd_plot <- plot_auc_steps(ftd_path, "Frontotemporal")

plot_grid(ad_plot, lbd_plot, ftd_plot, ncol = 3)


ad_path <- build_path(ad_full$path, ad_full$reference_auc)
ad_plot <- plot_auc_steps(ad_path, "Alzheimer's")

ad_plot + ftd_plot + lbd_plot

men_and_women <- function(df) {
  men <- df |>
    filter(female == 0) |>
    select(-female)
  women <- df |>
    filter(female == 1) |>
    select(-female)

  ad_women <- marker_subset_full(
    women,
    "Alzheimer's",
    c("Lewy bodies", "Frontotemporal")
  )
  ad_women_path <- build_path(ad_women$path, ad_women$reference_auc)
  ad_women_plot <- plot_auc_steps(ad_women_path, "Alzheimer's", " - Women")

  ad_men <- marker_subset_full(
    men,
    "Alzheimer's",
    c("Lewy bodies", "Frontotemporal")
  )
  ad_men_path <- build_path(ad_men$path, ad_men$reference_auc)
  ad_men_plot <- plot_auc_steps(ad_men_path, "Alzheimer's", " - Men")

  ad_men_plot + ad_women_plot

  ftd_women <- marker_subset_full(
    women,
    "Frontotemporal",
    c("Lewy bodies", "Alzheimer's")
  )
  ftd_women_path <- build_path(ftd_women$path, ftd_women$reference_auc)
  ftd_women_plot <- plot_auc_steps(ftd_women_path, "Frontotemporal", " - Women")

  ftd_men <- marker_subset_full(
    men,
    "Frontotemporal",
    c("Lewy bodies", "Alzheimer's")
  )
  ftd_men_path <- build_path(ftd_men$path, ftd_men$reference_auc)
  ftd_men_plot <- plot_auc_steps(ftd_men_path, "Frontotemporal", " - Men")

  ftd_men_plot + ftd_women_plot

  women[women$Diagnosis_combined == "Frontotemporal", ]
  men[men$Diagnosis_combined == "Frontotemporal", ]


  lbd_women <- marker_subset_full(
    women,
    "Lewy bodies",
    c("Frontotemporal", "Alzheimer's")
  )
  lbd_women_path <- build_path(lbd_women$path, lbd_women$reference_auc)
  lbd_women_plot <- plot_auc_steps(lbd_women_path, "Lewy bodies", " - Women")

  lbd_men <- marker_subset_full(
    men,
    "Lewy bodies",
    c("Frontotemporal", "Alzheimer's")
  )
  lbd_men_path <- build_path(lbd_men$path, lbd_men$reference_auc)
  lbd_men_plot <- plot_auc_steps(lbd_men_path, "Lewy bodies", " - Men")

  lbd_men_plot + lbd_women_plot

  lbd_men <- marker_subset_full(
    men,
    "Lewy bodies",
    c("Alzheimer's", "Frontotemporal")
  )
  write_rds(lbd_men, "LBD_men.rds")
  lbd_men_path <- build_path(lbd_men$path, lbd_men$reference_auc)
  lbd_men_plot <- plot_auc_steps(lbd_men_path, "Lewy bodies")
}
