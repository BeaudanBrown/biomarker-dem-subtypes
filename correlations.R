library(data.table)

standardise <- function(x) {
  return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
}

prepare <- function(df) {
  df <- df |> filter(!is.na(mean_elisa) & !is.na(mean_ab40) & !is.na(mean_ykl))

  washington <- df |>
    filter(Site == "Washington") |>
    mutate(
      mean_ab42_ab40_ratio = standardise(ifelse(mean_ab42_ab40_ratio > 0.5, NA,
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
      CSF_SNAP25 = standardise(CSF_SNAP25)
    ) |>
    setnames(
      c(
        "age_combined",
        "av45/PIB_fsuvr_rsf_tot_cortmean",
        "Centiloid_AV45_fSUVR_TOT_CORTMEA"
      ),
      c(
        "age",
        "raw_suvr",
        "centiloid"
      )
    )
  return(washington)
}


inflam_cdr <- function(df) {
  washington <- prepare(df)

  get_cdr_mmse_corr <- function(marker) {
    cdr_coef <- washington |>
      lm(
        as.formula(paste0(
          "cdr ~ ",
          paste(marker,
            "age",
            "female",
            sep = " + "
          )
        )),
        data = _
      ) |>
      coef()

    mmse_coef <- washington |>
      lm(
        as.formula(paste0(
          "MMSE ~ ",
          paste(marker,
            "age",
            "female",
            sep = " + "
          )
        )),
        data = _
      ) |>
      coef()

    return(c(cdr = cdr_coef[marker], mmse = mmse_coef[marker]))
  }
  markers <- c(
    "mean_elisa",
    "mean_ykl",
    "mean_gfap",
    "mean_ab42_ab40_ratio",
    "mean_ab40",
    "mean_ab42",
    "mean_tdp",
    "mean_ptau181",
    "mean_ptau217"
  )
  return(sapply(markers, get_cdr_mmse_corr))
}

csf_corr <- function(df) {
  washington <- prepare(df)

  get_ptau_csf_corr <- function(csf_marker) {
    ptau_coef <- washington |>
      lm(
        as.formula(paste0(
          csf_marker, " ~ ",
          paste("mean_ptau217",
            "age",
            "female",
            sep = " + "
          )
        )),
        data = _
      ) |>
      coef()

    return(ptau_coef["mean_ptau217"])
  }
  csf_markers <- washington |>
    select(starts_with("CSF_")) |>
    names()
  return(sapply(csf_markers, get_ptau_csf_corr))
}

pet_corr <- function(df) {
  pet_data <- df |>
    select(
      age_combined,
      female,
      "av45/PIB_fsuvr_rsf_tot_cortmean",
      "Centiloid_AV45_fSUVR_TOT_CORTMEA",
      zscore,
      mean_ptau217
    ) |>
    rename(
      "age" = "age_combined",
      "raw_suvr" = "av45/PIB_fsuvr_rsf_tot_cortmean",
      "centiloid" = "Centiloid_AV45_fSUVR_TOT_CORTMEA",
    ) |>
    mutate(
      zscore = standardise(zscore),
      raw_suvr = standardise(raw_suvr),
      centiloid = standardise(centiloid),
      mean_ptau217 = standardise(mean_ptau217)
    )

  sapply(c(
    "zscore",
    "raw_suvr",
    "centiloid"
  ), function(outcome) {
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
    return(coefs["mean_ptau217"])
  })
}


# 1. correlate each inflammatory marker (CD14, YKL, and GFAP) with CDR score and MMSE in the Wash U data only
# 2. correlate pTau217 with CSF biomarkers (which CSF biomarkers, AB ratio) in Wash U
# 3. Correlate pTau217 with AB PET measures
#   (which are the best - av45/PIB_fsuvr_rsf_tot_cortmean;
#    Centiloid_AV45_fSUVR_TOT_CORTMEA;
#    zscore;
#    Amyloid_Status Source)
# 4. Using primary analysis (SuperLearner) compare DLB and FTD and each diagnosis to controls
