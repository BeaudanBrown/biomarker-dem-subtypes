list(
  tar_target(
    whole_cohort,
    {
      df |> mutate(Diagnosis_combined = "Whole Cohort")
    }
  ),
  tar_target(
    all_subtypes_cohort,
    df[
      df$Diagnosis_combined %in%
        c("Alzheimer's", "Frontotermporal", "Lewy bodies"),
    ]
  ),
  tar_target(
    any_dementia_cohort,
    {
      combined_df <- df |>
        mutate(
          Diagnosis_combined = case_when(
            Diagnosis_combined == "Alzheimer's" ~ "Any Dem",
            Diagnosis_combined == "Frontotermporal" ~ "Any Dem",
            Diagnosis_combined == "Lewy bodies" ~ "Any Dem",
            TRUE ~ Diagnosis_combined
          )
        ) |>
        filter(Diagnosis_combined == "Any Dem")
      combined_df
    }
  ),
  tar_target(
    ad_cohort,
    df[df$Diagnosis_combined %in% c("Alzheimer's"), ]
  ),
  tar_target(
    ftd_cohort,
    df[df$Diagnosis_combined %in% c("Frontotermporal"), ]
  ),
  tar_target(
    lbd_cohort,
    df[df$Diagnosis_combined %in% c("Lewy bodies"), ]
  ),
  tar_target(
    all_cohorts,
    list(
      "Whole Cohort" = whole_cohort,
      "All Subtypes" = all_subtypes_cohort,
      "Any Dementia" = any_dementia_cohort,
      AD = ad_cohort,
      FTD = ftd_cohort,
      LBD = lbd_cohort
    )
  ),
  tar_target(
    standardised_cohorts,
    list(
      "Whole Cohort" = standardise_df(whole_cohort),
      "All Subtypes" = standardise_df(all_subtypes_cohort),
      "Any Dementia" = standardise_df(any_dementia_cohort),
      AD = standardise_df(ad_cohort),
      FTD = standardise_df(ftd_cohort),
      LBD = standardise_df(lbd_cohort)
    )
  ),
  tar_target(
    all_markers,
    c(
      "CD14",
      "NfL",
      "YKL-40",
      "GFAP",
      "AB42/AB40 Ratio",
      "AB40",
      "AB42",
      "TDP-43",
      "pTau-181",
      "pTau-217"
    )
  )
)
