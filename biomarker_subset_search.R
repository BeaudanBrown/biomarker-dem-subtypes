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

dotenv::load_dot_env()
data_dir <- Sys.getenv("DATA_DIR")

# load and prepare data

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

test <- backwards_search(df, "Alzheimer's", c("Lewy bodies", "Frontotermporal"))
