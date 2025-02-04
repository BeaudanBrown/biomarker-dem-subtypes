## Biomarker descriptions and predicting case status

library(tidyverse)
library(SuperLearner)
library(gtsummary)
library(ggthemes)
library(vimp)
library(xgboost)
library(glmnet)
library(ranger)
library(randomForest)
library(earth)
library(gam)
library(pROC)
library(origami)
source("helper_functions.R")
library(future)
library(future.apply)
library(parallel)
library(logistf)
library(rms)

# read data

df <- read_csv("../merged_data/merged_all_vars_28_01_2024.csv")

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
  df |> select(Diagnosis_combined, age_combined, female, starts_with("mean_"))

## Truncate biomarkers at the 99th percentile

df <- df |>
  mutate(across(starts_with("mean_"),
~ ifelse(.x > quantile(.x, 0.99), quantile(.x, 0.99), .x)))

df <- df |>
  mutate(across(starts_with("mean_"),
~ ifelse(.x < quantile(.x, 0.01), quantile(.x, 0.01), .x)))


#### Descriptives ####

# histograms

df |>
  select(starts_with("mean_")) |>
  mutate(across(c("mean_ykl", "mean_tdp", "mean_nfl"), log)) |>
  psych::multi.hist(global = FALSE, breaks = 30)

# density plot

df |>
  filter(!Diagnosis_combined == "Control") |>
  # remove large ab42/ab40 ratio outlier
  mutate(
    mean_ab42_ab40_ratio = ifelse(mean_ab42_ab40_ratio > 0.5, NA,
      mean_ab42_ab40_ratio),
    mean_ab40 = ifelse(mean_ab40 > 500, NA, mean_ab40),
    mean_ab42 = ifelse(mean_ab42 > 20, NA, mean_ab42)
  ) |>
  select(Diagnosis_combined, starts_with("mean_")) |>
  mutate(across(c("mean_ykl", "mean_tdp", "mean_nfl"), log)) |>
  pivot_longer(-Diagnosis_combined,
    names_to = "biomarker",
    names_prefix = "mean_"
  ) |>
  mutate(biomarker = fct_recode(biomarker,
    "AB40" = "ab40", "AB42" = "ab42","GFAP" = "gfap",
    "log(NfL)" = "nfl", "pTau181" = "ptau181", "pTau217" = "ptau217",
    "log(TDP-43)" = "tdp", "log(YKL-40)" = "ykl",
    "AB42/AB40" = "ab42_ab40_ratio"
  )) |>
  ggplot(aes(x = value, fill = Diagnosis_combined)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~biomarker, scales = "free") +
  labs(x = "Biomarker concentration", y = "Density", fill = "Dementia status") +
  scale_fill_colorblind() +
  bayesplot::theme_default() +
  theme(legend.position = "bottom")


ggsave("plots/conc_by_status.png",
  device = "png",
  width = 10, height = 10, bg = "white"
)

# table

df |>
  select(Diagnosis_combined, age_combined, female, Site, race_combined,
    starts_with("mean_")) |>
  mutate(across(c("mean_ykl", "mean_tdp", "mean_nfl"), log,
    .names = "log_{.col}"
  )) |>
  select(-mean_ykl, -mean_tdp, -mean_nfl) |>
  tbl_summary(
    by = Diagnosis_combined,
    missing = "no"
  ) |>
  add_p()


## Reliability

df |>
  select(Diagnosis_combined, ends_with("_ICC")) |>
  pivot_longer(-Diagnosis_combined,
    names_to = "biomarker"
  ) |>
  filter(!biomarker == "elisa_ICC") |>
  mutate(biomarker = fct_recode(biomarker,
    "AB40" = "ab40_ICC", "GFAP" = "gfap_ICC",
    "NfL" = "nfl_ICC", "pTau181" = "ptau181_ICC", "pTau217" = "ptau217_ICC",
    "TDP-43" = "tdp_ICC", "YKL-40" = "ykl_ICC"
  )) |>
  mutate(value = round(value, 3)) |>
  group_by(biomarker) |>
  summarise(ICC = mean(value, na.rm = T)) |>
  knitr::kable()


#### Prediction ####

# remove AB42 AB40 ratio

df <- df |> select(-mean_ab42_ab40_ratio)

# remove controls

df <- df |> filter(!Diagnosis_combined == "Control")

### ROC curves ###
# AD vs others

roc_ad <- auroc(df, "Alzheimer's", c("Frontotermporal", "Lewy bodies"))

plot_roc(roc_ad, "Alzheimer's vs other dementias")

# LB vs others

roc_lb <- auroc(df, "Lewy bodies", c("Frontotermporal", "Alzheimer's"))

plot_roc(roc_lb, "Lewy bodies vs other dementias")

# FT vs others

roc_ft <- auroc(df, "Frontotermporal", c("Lewy bodies", "Alzheimer's"))

plot_roc(roc_ft, "Frontotemporal vs other dementias")

# combined ROC plot

roc_out <- list(AD = roc_ad, LB = roc_lb, FT = roc_ft)

plot_roc_combined(roc_out, "ROC curves for dementia subtypes")

#### variable importance ####

plan(multicore, workers = 3)

vimp_out <- future_lapply(
  c("Alzheimer's", "Frontotermporal", "Lewy bodies"),
  FUN = function(.x) vimp_function(data = df, .x),
  future.seed = 1234
)

preds <- df |>
  select(age_combined, female, starts_with("mean_")) |> names()

# AD
knitr::kable(
  cbind(preds[as.numeric(vimp_out[[1]]$mat$s)],
  vimp_out[[1]]$est,
  vimp_out[[1]]$ci) |>
    as_tibble() |>
    mutate(across(c(V2, V3, V4), as.numeric)) |>
    mutate(across(c(V2, V3, V4), ~ round(.x, 3))) |>
    mutate(V2 = paste0(V2, " (", V3, ", ", V4, ")")) |>
    select(-V3, -V4)
)

# FT
knitr::kable(
  cbind(preds[as.numeric(vimp_out[[2]]$mat$s)],
  vimp_out[[2]]$est,
  vimp_out[[2]]$ci) |>
    as_tibble() |>
    mutate(across(c(V2, V3, V4), as.numeric)) |>
    mutate(across(c(V2, V3, V4), ~ round(.x, 3))) |>
    mutate(V2 = paste0(V2, " (", V3, ", ", V4, ")")) |>
    select(-V3, -V4)
)

# LB
knitr::kable(
  cbind(preds[as.numeric(vimp_out[[3]]$mat$s)],
  vimp_out[[3]]$est,
  vimp_out[[3]]$ci) |>
    as_tibble() |>
    mutate(across(c(V2, V3, V4), as.numeric)) |>
    mutate(across(c(V2, V3, V4), ~ round(.x, 3))) |>
    mutate(V2 = paste0(V2, " (", V3, ", ", V4, ")")) |>
    select(-V3, -V4)
)


### VIMP by hand ###

compare_auc <- function(data, outcome, reference, removed_var, full_auc){
  var_name <- names(data)[removed_var]
  data <- data[,-removed_var]
  auc_out <- mean(auroc(data, outcome, reference)$AUC)
  return(tibble(var = var_name, auc = auc_out, dif = full_auc - auc))
}

# AD

compare_out_ad <-
  map(2:ncol(df),
    compare_auc,
    data = df,
    outcome = "Alzheimer's",
    reference = c("Frontotermporal", "Lewy bodies"),
    full_auc = mean(roc_ad$AUC)
  )

bind_rows(compare_out_ad) |> mutate(dif = ifelse(dif < 0, 0, dif)) |>
  arrange(desc(dif))


# FTD

compare_out_ft <-
  map(2:ncol(df),
      compare_auc,
      data = df,
      outcome = "Frontotermporal",
      reference = c("Alzheimer's", "Lewy bodies"),
      full_auc = mean(roc_ft$AUC))

bind_rows(compare_out_ft) |> mutate(dif = ifelse(dif < 0, 0, dif)) |>
  arrange(desc(dif))


# LB

compare_out_lb <-
  map(2:ncol(df),
      compare_auc,
      data = df,
      outcome = "Lewy bodies",
      reference = c("Alzheimer's", "Frontotermporal"),
      full_auc = mean(roc_lb$AUC))

bind_rows(compare_out_lb) |> mutate(dif = ifelse(dif < 0, 0, dif)) |>
  arrange(desc(dif))
