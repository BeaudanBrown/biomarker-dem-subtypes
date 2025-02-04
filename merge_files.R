## Merging ADDF datafiles with biomarker, pheno, etc data

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(performance)
source("prepare_data_function.R")

#### Read data ####

### Phenotype

pheno <- read_excel("../ADDF phenotype_02222024.xlsx")

setnames(pheno, names(pheno), gsub(" ", "_", names(pheno)))

setnames(pheno, "PA_DB_UID", "PA_ID")

### match index to PA_ID

ids <- read_excel("../GILIPIN_ADDFPlasma_Manifest_Lab_20231107.xlsx")

ids <- ids |> filter(!is.na(`PA ID`))

# delete empty cols

ids <- ids[, (colSums(is.na(ids)) / nrow(ids)) < 0.75]

setnames(ids, c(
  "Index", "PA_ID", "Barcode", "Volume",
  "Destination_Plate", "Destination_Well"
))

### Case status

case <- read_excel("../ID_ALL_CASES_COMBINED.xlsx")

setnames(case, names(case), gsub(" ", "_", names(case)))

setnames(
  case,
  c("Status_at_Draw...3", "Status_at_Draw...4"),
  c("Status_at_Draw1", "Status_at_Draw2")
)

### Biomarkers

# read sheet names

marker_sheets <-
  excel_sheets("../ADDF results and master list_Jan 2025.xlsx")

ptau181 <- prepare_data(marker_sheets[2], "ptau181_")

nfl <- prepare_data(marker_sheets[3], "nfl_")

ab40 <- prepare_data(marker_sheets[4], "ab40_")

gfap <- prepare_data(marker_sheets[5], "gfap_")

ab42 <- prepare_data(marker_sheets[6], "ab42_")

ykl <- prepare_data(marker_sheets[7], "ykl_")

tdp <- prepare_data(marker_sheets[8], "tdp_")

ptau217 <- prepare_data(marker_sheets[9], "ptau217_")

elisa <- prepare_data(marker_sheets[10], "elisa_")


### Average technical replicates and calculate SD and CV
### and remove duplicates

## ptau

ptau181 <-
  filter(ptau181, !is.na(ptau181_Replicate_Conc) &
    !ptau181_Replicate_Conc == "NaN")

ptau181$ptau181_Replicate_Conc <- as.numeric(ptau181$ptau181_Replicate_Conc)

ptau181 <- ptau181 |>
  group_by(Sample_Barcode) |>
  mutate(ptau181_Replicates = 1:n()) |>
  mutate(ptau181_Replicated = ifelse(max(ptau181_Replicates) > 1, 1, 0)) |>
  ungroup()

ptau181_reps <- ptau181 |> filter(ptau181_Replicated == 1)
ptau181_noreps <- ptau181 |> filter(ptau181_Replicated == 0)

model_ptau181 <-
  lmer(ptau181_Replicate_Conc ~ 1 + (1 | Sample_Barcode),
    data = ptau181_reps
  )

ptau181_reps$mean_ptau181 <- as.numeric(predict(model_ptau181))

ptau181_noreps$mean_ptau181 <-
  if_else(ptau181_noreps$ptau181_Replicated == 1,
    NA_real_, ptau181_noreps$ptau181_Replicate_Conc
  )

ptau181 <- bind_rows(ptau181_reps, ptau181_noreps)

ptau181$ptau181_ICC <- icc(model_ptau181)$ICC_adjusted

ptau181 <- ptau181 |>
  group_by(Sample_Barcode) |>
  mutate(sd_replicates_ptau181 = sd(ptau181_Replicate_Conc, na.rm = T)) |>
  ungroup()

ptau181 <- ptau181 |>
  group_by(Sample_Barcode) |>
  mutate(
    cv_replicates_ptau181 =
      (sd_replicates_ptau181 / mean_ptau181) * 100
  ) |>
  ungroup()

ptau181 <- filter(ptau181, !duplicated(Sample_Barcode))

## nfl

nfl <-
  filter(nfl, !is.na(nfl_Replicate_Conc) &
    !nfl_Replicate_Conc == "NaN")

nfl$nfl_Replicate_Conc <- as.numeric(nfl$nfl_Replicate_Conc)

nfl <- nfl |>
  group_by(Sample_Barcode) |>
  mutate(nfl_Replicates = 1:n()) |>
  mutate(nfl_Replicated = ifelse(max(nfl_Replicates) > 1, 1, 0)) |>
  ungroup()

nfl_reps <- nfl |> filter(nfl_Replicated == 1)
nfl_noreps <- nfl |> filter(nfl_Replicated == 0)

model_nfl <-
  lmer(nfl_Replicate_Conc ~ 1 + (1 | Sample_Barcode),
    data = nfl_reps
  )

nfl_reps$mean_nfl <- as.numeric(predict(model_nfl))

nfl_noreps$mean_nfl <-
  if_else(nfl_noreps$nfl_Replicated == 1,
    NA_real_, nfl_noreps$nfl_Replicate_Conc
  )

nfl <- bind_rows(nfl_reps, nfl_noreps)

nfl$nfl_ICC <- icc(model_nfl)$ICC_adjusted

nfl <- nfl |>
  group_by(Sample_Barcode) |>
  mutate(sd_replicates_nfl = sd(nfl_Replicate_Conc, na.rm = T)) |>
  ungroup()

nfl <- nfl |>
  group_by(Sample_Barcode) |>
  mutate(
    cv_replicates_nfl =
      (sd_replicates_nfl / mean_nfl) * 100
  ) |>
  ungroup()

nfl <- filter(nfl, !duplicated(Sample_Barcode))


## ab40

ab40 <-
  filter(ab40, !is.na(ab40_Replicate_Conc) &
    !ab40_Replicate_Conc == "NaN")

ab40$ab40_Replicate_Conc <- as.numeric(ab40$ab40_Replicate_Conc)

ab40 <- ab40 |>
  group_by(Sample_Barcode) |>
  mutate(ab40_Replicates = 1:n()) |>
  mutate(ab40_Replicated = ifelse(max(ab40_Replicates) > 1, 1, 0)) |>
  ungroup()

ab40_reps <- ab40 |> filter(ab40_Replicated == 1)
ab40_noreps <- ab40 |> filter(ab40_Replicated == 0)

model_ab40 <-
  lmer(ab40_Replicate_Conc ~ 1 + (1 | Sample_Barcode),
    data = ab40_reps
  )

ab40_reps$mean_ab40 <- as.numeric(predict(model_ab40))

ab40_noreps$mean_ab40 <-
  if_else(ab40_noreps$ab40_Replicated == 1,
    NA_real_, ab40_noreps$ab40_Replicate_Conc
  )

ab40 <- bind_rows(ab40_reps, ab40_noreps)

ab40$ab40_ICC <- icc(model_ab40)$ICC_adjusted

ab40 <- ab40 |>
  group_by(Sample_Barcode) |>
  mutate(sd_replicates_ab40 = sd(ab40_Replicate_Conc, na.rm = T)) |>
  ungroup()

ab40 <- ab40 |>
  group_by(Sample_Barcode) |>
  mutate(
    cv_replicates_ab40 =
      (sd_replicates_ab40 / mean_ab40) * 100
  ) |>
  ungroup()

ab40 <- filter(ab40, !duplicated(Sample_Barcode))


## gfap

gfap <-
  filter(gfap, !is.na(gfap_Replicate_Conc) &
    !gfap_Replicate_Conc == "NaN")

gfap$gfap_Replicate_Conc <- as.numeric(gfap$gfap_Replicate_Conc)

gfap <- gfap |>
  group_by(Sample_Barcode) |>
  mutate(gfap_Replicates = 1:n()) |>
  mutate(gfap_Replicated = ifelse(max(gfap_Replicates) > 1, 1, 0)) |>
  ungroup()

gfap_reps <- gfap |> filter(gfap_Replicated == 1)
gfap_noreps <- gfap |> filter(gfap_Replicated == 0)

model_gfap <-
  lmer(gfap_Replicate_Conc ~ 1 + (1 | Sample_Barcode),
    data = gfap_reps
  )

gfap_reps$mean_gfap <- as.numeric(predict(model_gfap))

gfap_noreps$mean_gfap <-
  if_else(gfap_noreps$gfap_Replicated == 1,
    NA_real_, gfap_noreps$gfap_Replicate_Conc
  )

gfap <- bind_rows(gfap_reps, gfap_noreps)

gfap$gfap_ICC <- icc(model_gfap)$ICC_adjusted

gfap <- gfap |>
  group_by(Sample_Barcode) |>
  mutate(sd_replicates_gfap = sd(gfap_Replicate_Conc, na.rm = T)) |>
  ungroup()

gfap <- gfap |>
  group_by(Sample_Barcode) |>
  mutate(
    cv_replicates_gfap =
      (sd_replicates_gfap / mean_gfap) * 100
  ) |>
  ungroup()

gfap <- filter(gfap, !duplicated(Sample_Barcode))


## ab42

ab42 <-
  filter(ab42, !is.na(ab42_Replicate_Conc) &
    !ab42_Replicate_Conc == "NaN")

ab42$ab42_Replicate_Conc <- as.numeric(ab42$ab42_Replicate_Conc)

ab42 <- ab42 |>
  group_by(Sample_Barcode) |>
  mutate(ab42_Replicates = 1:n()) |>
  mutate(ab42_Replicated = ifelse(max(ab42_Replicates) > 1, 1, 0)) |>
  ungroup()

ab42_reps <- ab42 |> filter(ab42_Replicated == 1)
ab42_noreps <- ab42 |> filter(ab42_Replicated == 0)

model_ab42 <-
  lmer(ab42_Replicate_Conc ~ 1 + (1 | Sample_Barcode),
    data = ab42_reps
  )

ab42_reps$mean_ab42 <- as.numeric(predict(model_ab42))

ab42_noreps$mean_ab42 <-
  if_else(ab42_noreps$ab42_Replicated == 1,
    NA_real_, ab42_noreps$ab42_Replicate_Conc
  )

ab42 <- bind_rows(ab42_reps, ab42_noreps)

ab42$ab42_ICC <- icc(model_ab42)$ICC_adjusted

ab42 <- ab42 |>
  group_by(Sample_Barcode) |>
  mutate(sd_replicates_ab42 = sd(ab42_Replicate_Conc, na.rm = T)) |>
  ungroup()

ab42 <- ab42 |>
  group_by(Sample_Barcode) |>
  mutate(
    cv_replicates_ab42 =
      (sd_replicates_ab42 / mean_ab42) * 100
  ) |>
  ungroup()

ab42 <- filter(ab42, !duplicated(Sample_Barcode))


## ykl

ykl <-
  filter(ykl, !is.na(ykl_Calc_Concentration) &
    !ykl_Calc_Concentration == "NaN")

ykl$ykl_Replicate_Conc <- as.numeric(ykl$ykl_Calc_Concentration)

model_ykl <-
  lmer(ykl_Replicate_Conc ~ 1 + (1 | Sample_Barcode),
    data = ykl
  )

ykl$mean_ykl <- as.numeric(predict(model_ykl))

ykl$ykl_ICC <- icc(model_ykl)$ICC_adjusted

ykl <- ykl |>
  group_by(Sample_Barcode) |>
  mutate(sd_replicates_ykl = sd(ykl_Replicate_Conc, na.rm = T)) |>
  ungroup()

ykl <- ykl |>
  group_by(Sample_Barcode) |>
  mutate(
    cv_replicates_ykl =
      (sd_replicates_ykl / mean_ykl) * 100
  ) |>
  ungroup()

ykl <- filter(ykl, !duplicated(Sample_Barcode))


## tdp

tdp <-
  filter(tdp, !is.na(tdp_Replicate_Conc) &
    !tdp_Replicate_Conc == "NaN")

tdp$tdp_Replicate_Conc <- as.numeric(tdp$tdp_Replicate_Conc)

tdp <- filter(tdp, !is.na(tdp_Replicate_Conc))

tdp <- tdp |>
  group_by(Sample_Barcode) |>
  mutate(tdp_Replicates = 1:n()) |>
  mutate(tdp_Replicated = ifelse(max(tdp_Replicates) > 1, 1, 0)) |>
  ungroup()

tdp_reps <- tdp |> filter(tdp_Replicated == 1)
tdp_noreps <- tdp |> filter(tdp_Replicated == 0)

model_tdp <-
  lmer(tdp_Replicate_Conc ~ 1 + (1 | Sample_Barcode),
    data = tdp_reps
  )

tdp_reps$mean_tdp <- as.numeric(predict(model_tdp))

tdp_noreps$mean_tdp <-
  if_else(tdp_noreps$tdp_Replicated == 1,
    NA_real_, tdp_noreps$tdp_Replicate_Conc
  )

tdp <- bind_rows(tdp_reps, tdp_noreps)

tdp$tdp_ICC <- icc(model_tdp)$ICC_adjusted

tdp <- tdp |>
  group_by(Sample_Barcode) |>
  mutate(sd_replicates_tdp = sd(tdp_Replicate_Conc, na.rm = T)) |>
  ungroup()

tdp <- tdp |>
  group_by(Sample_Barcode) |>
  mutate(
    cv_replicates_tdp =
      (sd_replicates_tdp / mean_tdp) * 100
  ) |>
  ungroup()

tdp <- filter(tdp, !duplicated(Sample_Barcode))


## ptau217

ptau217 <-
  filter(ptau217, !is.na(ptau217_Replicate_Conc) &
    !ptau217_Replicate_Conc == "NaN")

ptau217$ptau217_Replicate_Conc <- as.numeric(ptau217$ptau217_Replicate_Conc)

ptau217 <- filter(ptau217, !is.na(ptau217_Replicate_Conc))

ptau217 <- ptau217 |>
  group_by(Sample_Barcode) |>
  mutate(ptau217_Replicates = 1:n()) |>
  mutate(ptau217_Replicated = ifelse(max(ptau217_Replicates) > 1, 1, 0)) |>
  ungroup()

ptau217_reps <- ptau217 |> filter(ptau217_Replicated == 1)
ptau217_noreps <- ptau217 |> filter(ptau217_Replicated == 0)

model_ptau217 <-
  lmer(ptau217_Replicate_Conc ~ 1 + (1 | Sample_Barcode),
    data = ptau217_reps
  )

ptau217_reps$mean_ptau217 <- as.numeric(predict(model_ptau217))

ptau217_noreps$mean_ptau217 <-
  if_else(ptau217_noreps$ptau217_Replicated == 1,
    NA_real_, ptau217_noreps$ptau217_Replicate_Conc
  )

ptau217 <- bind_rows(ptau217_reps, ptau217_noreps)

ptau217$ptau217_ICC <- icc(model_ptau217)$ICC_adjusted

ptau217 <- ptau217 |>
  group_by(Sample_Barcode) |>
  mutate(sd_replicates_ptau217 = sd(ptau217_Replicate_Conc, na.rm = T)) |>
  ungroup()

ptau217 <- ptau217 |>
  group_by(Sample_Barcode) |>
  mutate(
    cv_replicates_ptau217 =
      (sd_replicates_ptau217 / mean_ptau217) * 100
  ) |>
  ungroup()

ptau217 <- filter(ptau217, !duplicated(Sample_Barcode))


## elisa (extreme low observations removed)

elisa <-
  filter(elisa, !is.na(elisa_Mean_Concentration) &
    !elisa_Mean_Concentration == "NaN")

elisa$elisa_Replicate_Conc <- as.numeric(elisa$elisa_Mean_Concentration)

elisa <- elisa |>
  group_by(Sample_Barcode) |>
  mutate(elisa_Replicates = 1:n()) |>
  mutate(elisa_Replicated = ifelse(max(elisa_Replicates) > 1, 1, 0)) |>
  ungroup()

elisa_reps <- elisa |> filter(elisa_Replicated == 1)
elisa_noreps <- elisa |> filter(elisa_Replicated == 0)

model_elisa <-
  lmer(elisa_Replicate_Conc ~ 1 + (1 | Sample_Barcode),
    data = elisa_reps
  )

elisa_reps$mean_elisa <- as.numeric(predict(model_elisa))

elisa_noreps$mean_elisa <-
  if_else(elisa_noreps$elisa_Replicated == 1,
    NA_real_, elisa_noreps$elisa_Replicate_Conc
  )

elisa <- bind_rows(elisa_reps, elisa_noreps)

elisa$elisa_ICC <- icc(model_elisa)$ICC_adjusted

elisa <- elisa |>
  group_by(Sample_Barcode) |>
  mutate(sd_replicates_elisa = sd(elisa_Replicate_Conc, na.rm = T)) |>
  ungroup()

elisa <- elisa |>
  group_by(Sample_Barcode) |>
  mutate(
    cv_replicates_elisa =
      (sd_replicates_elisa / mean_elisa) * 100
  ) |>
  ungroup()

elisa <- filter(elisa, !duplicated(Sample_Barcode))


#### Join dataframes ####

### Join biomarker dataframes

markers <-
  ptau181 |>
  full_join(nfl, by = "Sample_Barcode") |>
  full_join(ab40, by = "Sample_Barcode") |>
  full_join(ab42, by = "Sample_Barcode") |>
  full_join(gfap, by = "Sample_Barcode") |>
  full_join(ykl, by = "Sample_Barcode") |>
  full_join(tdp, by = "Sample_Barcode") |>
  full_join(ptau217, by = "Sample_Barcode") |>
  full_join(elisa, by = "Sample_Barcode")

### Add ID variable

setnames(ids, "Barcode", "Sample_Barcode")

joined <- full_join(ids, markers, by = "Sample_Barcode")

### Add in pheno

joined <- full_join(pheno, joined, by = "PA_ID")

### Add in cases

case <- case |> select(-DrawDate, -DOB, -sex)

joined <- full_join(joined, case, by = "PA_ID")

### Add in Texas case & demo data

new_cases <-
  read_excel("../Demographic data Texas/ADDF preliminary dataset Jan 2025.xlsx")

new_cases_update <-
  read_excel("../Demographic data Texas/ADDF Jan 21 2025 update.xlsx")

new_cases_update <-
  new_cases_update |>
  select(-`Jan 21 2025 update`, -Notes)

setnames(
  new_cases,
  c(
    "ALIAS", "Is the participant Hispanic?","Fasting (hh:mm)",
    "Fasting >8 hours", "Year collected"
  ),
  c("Sample_Barcode", "Hispanic", "Fasting_time",
  "Fasting_8_hours", "Year_collected")
)

setnames(
  new_cases_update,
  c(
    "ALIAS", "Is the participant Hispanic?","Fasting (hh:mm)",
    "Fasting >8 hours", "Year collected"
  ),
  c("Sample_Barcode", "Hispanic", "Fasting_time",
  "Fasting_8_hours", "Year_collected")
)

new_cases_update$Fasting_time <- as.character(new_cases_update$Fasting_time)

new_cases <-
  new_cases |>
  filter(!Sample_Barcode %in% new_cases_update$Sample_Barcode) |>
  bind_rows(new_cases_update)

new_cases$Sample_Barcode <- as.character(new_cases$Sample_Barcode)

### Join Texas and Washington data

joined <- full_join(joined, new_cases, by = "Sample_Barcode")

# amalgamate age and sex variables

joined$age_combined <-
  ifelse(is.na(joined$Age_at_draw), joined$AGE, joined$Age_at_draw)

joined$sex_combined <-
  ifelse(is.na(joined$sex), joined$SEX, joined$sex)

### Combine data from

#### Write to csv ####

#write_csv(joined, "../merged_jan28_2025.csv")
