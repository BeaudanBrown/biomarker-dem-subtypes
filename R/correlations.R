## Functions for correlations

standardise <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

prepare <- function(df) {
  washington <- df |>
    filter(Site == "Washington") |>
    mutate(
      mean_ab42_ab40_ratio = standardise(ifelse(
        mean_ab42_ab40_ratio > 0.5,
        NA,
        mean_ab42_ab40_ratio
      )),
      mean_ab40 = standardise(ifelse(mean_ab40 > 500, NA, mean_ab40)),
      mean_ab42 = standardise(ifelse(mean_ab42 > 20, NA, mean_ab42)),
      mean_ykl = standardise(log(mean_ykl)),
      mean_ptau181 = standardise(mean_ptau181),
      mean_tdp = standardise(log(mean_tdp)),
      mean_ptau217 = standardise(mean_ptau217),
      mean_elisa = standardise(mean_elisa),
      mean_gfap = standardise(mean_gfap),
      MMSE = standardise(MMSE),
      cdr = standardise(cdr),
      CSF_NfL = standardise(CSF_NfL),
      CSF_sTREM2 = standardise(CSF_sTREM2),
      CSF_VILIP1 = standardise(CSF_VILIP1),
      CSF_Ng = standardise(CSF_Ng),
      CSF_AB_Ratio = standardise(LUMIPULSE_CSF_AB42 / LUMIPULSE_CSF_AB40),
      CSF_SNAP25 = standardise(CSF_SNAP25)
    ) |>
    rename(
      "age" = "age_combined",
      "raw_suvr" = "av45/PIB_fsuvr_rsf_tot_cortmean",
      "centiloid" = "Centiloid_AV45_fSUVR_TOT_CORTMEA",
    )
  return(washington)
}


inflam_cdr <- function(df, diagnosis = NULL) {
  if (!is.null(diagnosis)) {
    df <- df[df$Diagnosis_combined %in% diagnosis, ]
  }
  washington <- prepare(df)

  get_cdr_corr <- function(marker) {
    cdr_df <- washington |>
      filter(!is.na(get(marker)) & !is.na(cdr))
    cdr_coef <- cdr_df |>
      lm(
        as.formula(paste0(
          "cdr ~ ",
          paste(marker, "age", "female", sep = " + ")
        )),
        data = _
      ) |>
      coef()

    return(sprintf("%.3f (%d)", cdr_coef[[marker]], nrow(cdr_df)))
  }

  get_mmse_corr <- function(marker) {
    mmse_df <- washington |>
      filter(!is.na(get(marker)) & !is.na(MMSE))
    mmse_coef <- mmse_df |>
      lm(
        as.formula(paste0(
          "MMSE ~ ",
          paste(marker, "age", "female", sep = " + ")
        )),
        data = _
      ) |>
      coef()

    return(sprintf("%.3f (%d)", mmse_coef[[marker]], nrow(mmse_df)))
  }
  markers <- c(
    "mean_elisa",
    "mean_nfl",
    "mean_ykl",
    "mean_gfap",
    "mean_ab42_ab40_ratio",
    "mean_ab40",
    "mean_ab42",
    "mean_tdp",
    "mean_ptau181",
    "mean_ptau217"
  )
  return(list(
    mmse = sapply(markers, get_mmse_corr),
    cdr = sapply(markers, get_cdr_corr)
  ))
}

csf_corr <- function(df, diagnosis = NULL) {
  if (!is.null(diagnosis)) {
    df <- df[df$Diagnosis_combined %in% diagnosis, ]
  }
  washington <- prepare(df) |>
    filter(!is.na(mean_ptau217) & !is.na(CSF_AB_Ratio))

  get_ptau_csf_corr <- function(csf_marker) {
    ptau_coef <- washington |>
      lm(
        as.formula(paste0(
          csf_marker,
          " ~ ",
          paste("mean_ptau217", "age", "female", sep = " + ")
        )),
        data = _
      ) |>
      coef()

    return(ptau_coef["mean_ptau217"])
  }
  csf_markers <- washington |>
    select(CSF_AB_Ratio) |>
    names()
  return(list(
    n = nrow(washington),
    coefs = sapply(csf_markers, get_ptau_csf_corr)
  ))
}

pet_corr <- function(df, diagnosis = NULL) {
  if (!is.null(diagnosis)) {
    df <- df[df$Diagnosis_combined %in% diagnosis, ]
  }
  pet_data <- prepare(df) |>
    select(
      age,
      female,
      raw_suvr,
      centiloid,
      zscore,
      mean_ptau217
    ) |>
    filter(!is.na(mean_ptau217) & !is.na(raw_suvr)) |>
    mutate(
      zscore = standardise(zscore),
      raw_suvr = standardise(raw_suvr),
      centiloid = standardise(centiloid),
      mean_ptau217 = standardise(mean_ptau217)
    )
  n <- nrow(pet_data)

  coefs <- sapply(
    c(
      "zscore",
      "raw_suvr",
      "centiloid"
    ),
    function(outcome) {
      coefs <- lm(
        as.formula(paste0(
          outcome,
          " ~ ",
          paste(
            "mean_ptau217",
            "age",
            "female",
            sep = " + "
          )
        )),
        data = pet_data
      ) |>
        coef()
      return(coefs[["mean_ptau217"]])
    }
  )
  return(list(n = n, coefs = coefs))
}

# 1. correlate each inflammatory marker (CD14, YKL, and GFAP) with CDR score and MMSE in the Wash U data only
# 2. correlate pTau217 with CSF biomarkers (which CSF biomarkers, AB ratio) in Wash U
# 3. Correlate pTau217 with AB PET measures
#   (which are the best - av45/PIB_fsuvr_rsf_tot_cortmean;
#    Centiloid_AV45_fSUVR_TOT_CORTMEA;
#    zscore;
#    Amyloid_Status Source)
# 4. Using primary analysis (SuperLearner) compare DLB and FTD and each diagnosis to controls
