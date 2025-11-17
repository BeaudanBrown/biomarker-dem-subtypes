list(
  # ===================================================================
  # MMSE
  # ===================================================================
  tar_target(
    mmse_corrs,
    get_marker_corrs(standardised_cohorts, "MMSE"),
    pattern = map(standardised_cohorts)
  ),
  tar_target(
    mmse_corr_table,
    get_marker_table(mmse_corrs)
  ),

  # ===================================================================
  # CDR
  # ===================================================================
  tar_target(
    cdr_corrs,
    get_marker_corrs(standardised_cohorts, "cdr"),
    pattern = map(standardised_cohorts)
  ),
  tar_target(
    cdr_corr_table,
    get_marker_table(cdr_corrs)
  )
)
