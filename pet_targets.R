list(
  # ===================================================================
  # Z-Score
  # ===================================================================
  tar_target(
    pet_zscore_corrs,
    get_marker_corrs(standardised_cohorts, "zscore"),
    pattern = map(standardised_cohorts)
  ),
  tar_target(
    pet_zscore_corr_table,
    get_marker_table(pet_zscore_corrs)
  ),

  # ===================================================================
  # Centiloid
  # ===================================================================
  tar_target(
    pet_centiloid_corrs,
    get_marker_corrs(standardised_cohorts, "centiloid"),
    pattern = map(standardised_cohorts)
  ),
  tar_target(
    pet_centiloid_corr_table,
    get_marker_table(pet_centiloid_corrs, show_ci = TRUE)
  ),
  tar_target(
    pet_centiloid_corr_table_csv,
    write_marker_table_csv(
      pet_centiloid_corrs,
      "output/supplementary/pet_centiloid_table.csv"
    ),
    format = "file"
  ),

  # ===================================================================
  # Raw SUVR
  # ===================================================================
  tar_target(
    pet_raw_suvr_corrs,
    get_marker_corrs(standardised_cohorts, "raw_suvr"),
    pattern = map(standardised_cohorts)
  ),
  tar_target(
    pet_raw_suvr_corr_table,
    get_marker_table(pet_raw_suvr_corrs)
  )
)
