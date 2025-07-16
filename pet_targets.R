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
  tar_target(
    pet_zscore_plots,
    get_marker_plots(all_cohorts, outcome = "zscore"),
    pattern = map(all_cohorts)
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
    get_marker_table(pet_centiloid_corrs)
  ),
  tar_target(
    pet_centiloid_plots,
    get_marker_plots(all_cohorts, outcome = "centiloid"),
    pattern = map(all_cohorts)
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
  ),
  tar_target(
    pet_raw_suvr_plots,
    get_marker_plots(all_cohorts, outcome = "raw_suvr"),
    pattern = map(all_cohorts)
  )
)
