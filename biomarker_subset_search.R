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
source("learners.R")
source("helper_functions.R")
source("subset_search_helper_functions.R")

dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")
plan(multicore, workers = detectCores())

plot_auc_steps <- function(data, outcome) {
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
    # Add reference line
    geom_hline(yintercept = ref_auc, linetype = "dashed", color = "red", size = 0.8) +
    # Add annotation for reference line
    annotate("text", x = round(max(data$step) / 2), y = ref_auc - 0.002,
             label = paste("Threshold:", round(ref_auc, 3)),
             hjust = 1, color = "red", size = 3.5) +
    geom_line(size = 1, color = "steelblue") +
    geom_point(size = 3, color = "darkblue") +
    labs(
      title = paste0(outcome, " AUC Marker Subset Path"),
      x = "Markers Remaining",
      y = "AUC"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    # Use the removed_var values as the x-axis labels
    scale_x_continuous(
      breaks = data$step,
      labels = labels,
      minor_breaks = NULL
    ) +
    # Adjust y-axis to ensure the reference line is visible
    scale_y_continuous(limits = c(min(c(data$auc, ref_auc)) * 0.98, max(data$auc) * 1.02))

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
    add_row(removed_var = "Full", auc = ref_auc, step = 0) |>
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

# exclude those missing biomarkers

df <- df |> filter(!is.na(mean_ab40) & !is.na(mean_ykl))

# code diagnosis as factor and set control as reference

df <- df |>
  mutate(Diagnosis_combined = as.factor(fct_recode(Diagnosis_combined,
    "Alzheimer's" = "AD", "Control" = "CO", "Lewy bodies" = "DLB",
    "Frontotermporal" = "FTD"
  ))) |>
  mutate(Diagnosis_combined = fct_relevel(Diagnosis_combined, "Control"))

# set sex to factor

df$sex_combined <- as.factor(df$sex_combined)
df$female <- ifelse(df$sex_combined == "Female", 1, 0)

# Add ab42 to ab40 ratio

df$mean_ab42_ab40_ratio <- df$mean_ab42 / df$mean_ab40

# remove cd14

df <- df |> select(-mean_elisa)

# select relevant variables
df <-
  df |> select(
    Diagnosis_combined, age_combined, female, starts_with("mean_"),
    Site, race_combined, ends_with("_ICC")
  )

# Truncate biomarkers at the 99th percentile

df <- df |>
  mutate(across(starts_with("mean_"),
~ ifelse(.x > quantile(.x, 0.99), quantile(.x, 0.99), .x)))

df <- df |>
  mutate(across(starts_with("mean_"),
~ ifelse(.x < quantile(.x, 0.01), quantile(.x, 0.01), .x)))

# remove unneeded vars

df <- df |> select(-mean_ab42_ab40_ratio)

# remove controls

df <- df |> filter(!Diagnosis_combined == "Control")

# remove extras

df <- df |> select(-Site, -race_combined, -ends_with("_ICC"))

## Cross validated biomarker subset search

test <- get_biomarker_subset(
  df,
  "Alzheimer's",
  c("Lewy bodies", "Frontotermporal"),
  nfolds = 5
)
