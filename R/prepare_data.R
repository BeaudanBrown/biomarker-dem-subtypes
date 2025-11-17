## Clean data

clean_data <- function(df, keep_cd14 = TRUE) {
  # Merge diagnosis variables between Texas and Washington data

  df <-
    df |>
    mutate(
      Diagnosis_combined = case_when(
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
      )
    )

  # Merge race

  df <- df |>
    mutate(
      race_combined = case_when(
        Race == "White" ~ "White",
        Race == "African American" ~ "African American",
        Race == "Asian" ~ "Asian",
        Race == "American Indian or Alaska Native" ~
          "American Indian/Alaska Native",
        RACE == "White" ~ "White",
        RACE == "Black/African American" ~ "African American",
        RACE == "Asian" ~ "Asian",
        RACE == "American Indian/Alaska Native" ~
          "American Indian/Alaska Native",
        is.na(Race) & is.na(RACE) ~ "Other",
        TRUE ~ "Other"
      )
    )

  # add site variable

  df$Site <- ifelse(!is.na(df$SEX), "Texas", "Washington")

  # exclude MCI and other

  df <- df |>
    filter(
      !Diagnosis_combined %in% c("MCI", "other") & !is.na(Diagnosis_combined)
    )

  # code diagnosis as factor and set control as reference

  df <- df |>
    mutate(
      Diagnosis_combined = as.factor(fct_recode(
        Diagnosis_combined,
        "Alzheimer's" = "AD",
        "Control" = "CO",
        "Lewy bodies" = "DLB",
        "Frontotermporal" = "FTD"
      ))
    ) |>
    mutate(Diagnosis_combined = fct_relevel(Diagnosis_combined, "Control"))

  # set sex to factor

  df$sex <- as.factor(df$sex)
  df$female <- ifelse(df$sex == "Female", 1, 0)

  # AB42 to AB40 ratio
  df$mean_ab42_ab40_ratio <- df$mean_ab42 / df$mean_ab40

  # remove CD14?
  if (!keep_cd14) {
    df <- select(df, -mean_elisa)
  }

  ## Truncate biomarkers at the 99th percentile

  df <- df |>
    mutate(across(
      starts_with("mean_"),
      ~ ifelse(
        .x > quantile(.x, 0.99, na.rm = TRUE),
        quantile(.x, 0.99, na.rm = TRUE),
        .x
      )
    ))

  df <- df |>
    mutate(across(
      starts_with("mean_"),
      ~ ifelse(
        .x < quantile(.x, 0.01, na.rm = TRUE),
        quantile(.x, 0.01, na.rm = TRUE),
        .x
      )
    ))

  df |> filter(complete.cases(select(df, starts_with("mean_"))))
}

## Prepare data

prepare_data <- function(sheet_name, marker_name) {
  df <-
    read_excel(
      file.path(data_dir, Sys.getenv("ADDF_FILE")),
      sheet = sheet_name,
      guess_max = 100000,
    )

  setnames(df, "Sample", "Sample_Barcode", skip_absent = TRUE)
  setnames(df, "Sample ID", "Sample_Barcode", skip_absent = TRUE)

  setnames(df, names(df), gsub(" ", "_", names(df)))
  setnames(df, names(df), gsub("\\.", "", names(df)))

  setnames(df, names(df)[-1], paste0(marker_name, names(df)[-1]))

  # remove blank rows

  df <- filter(df, !is.na(Sample_Barcode))

  # delete mostly empty cols

  df <- df[, (colSums(is.na(df)) / nrow(df)) < 0.75]

  # if sample_barcode has an asteriks, add in dilution_note

  df[[paste0(marker_name, "Dilution_note")]] <-
    ifelse(
      str_detect(df$Sample_Barcode, "\\*"),
      "indicates diluted 1:4... (see excel)",
      NA_character_
    )

  # if sample_barcode has an asteriks, remove

  df$Sample_Barcode <- as.character(gsub("\\*", "", df$Sample_Barcode))
  df
}

prepare_roc_data <- function(df) {
  df |>
    select(
      Diagnosis_combined,
      age,
      mean_elisa,
      mean_nfl,
      mean_ykl,
      mean_gfap,
      mean_ab42,
      mean_ab40,
      mean_tdp,
      mean_ptau181,
      mean_ptau217,
      female
    )
}

prepare_roc_data_ratio <- function(df) {
  df |>
    select(
      Diagnosis_combined,
      age,
      mean_elisa,
      mean_nfl,
      mean_ykl,
      mean_gfap,
      mean_ab42_ab40_ratio,
      mean_tdp,
      mean_ptau181,
      mean_ptau217,
      female
    )
}

prepare_roc_data_fasting <- function(df) {
  setDT(df)

  df[,
    fasting_combined := ifelse(
      is.na(Fasting_Status),
      Fasting_8_hours,
      Fasting_Status
    )
  ][,
    fasting_combined := ifelse(
      fasting_combined == "Unknown" |
        fasting_combined == "unknown",
      NA,
      fasting_combined
    )
  ][,
    fasting_combined := ifelse(
      fasting_combined == "non-fasted",
      "No",
      fasting_combined
    )
  ]

  df[, list(
    Diagnosis_combined,
    age,
    mean_elisa,
    mean_nfl,
    mean_ykl,
    mean_gfap,
    mean_ab42,
    mean_ab40,
    mean_tdp,
    mean_ptau181,
    mean_ptau217,
    female,
    fasting_combined
  )]
}

prepare_roc_data_cdr <- function(df) {
  df |>
    select(
      Diagnosis_combined,
      age,
      mean_elisa,
      mean_nfl,
      mean_ykl,
      mean_gfap,
      mean_ab42,
      mean_ab40,
      mean_tdp,
      mean_ptau181,
      mean_ptau217,
      female,
      cdr
    )
}
