list(
  tar_target(
    csf_corrs,
    get_marker_corrs(standardised_cohorts, "CSF_AB_Ratio"),
    pattern = map(standardised_cohorts)
  ),
  tar_target(
    csf_corr_table,
    get_marker_table(csf_corrs, show_ci = TRUE)
  ),
  tar_target(
    csf_corr_table_csv,
    write_marker_table_csv(
      csf_corrs,
      "output/supplementary/csf_ab42_ab40_table.csv"
    ),
    format = "file"
  )
)
