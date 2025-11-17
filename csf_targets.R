list(
  tar_target(
    csf_corrs,
    get_marker_corrs(standardised_cohorts, "CSF_AB_Ratio"),
    pattern = map(standardised_cohorts)
  ),
  tar_target(
    csf_corr_table,
    get_marker_table(csf_corrs)
  )
)
