list(
  tar_target(
    csf_corrs,
    get_marker_corrs(standardised_cohorts, "CSF_AB_Ratio"),
    pattern = map(standardised_cohorts)
  ),
  tar_target(
    csf_corr_table,
    get_marker_table(csf_corrs)
  ),
  tar_target(
    csf_plots,
    {
      # TODO: Should moved to data prep
      cohort_name <- names(all_cohorts)
      all_cohorts[[cohort_name]] <- all_cohorts[[cohort_name]] |>
        mutate(CSF_AB_Ratio = LUMIPULSE_CSF_AB42 / LUMIPULSE_CSF_AB40)
      get_marker_plots(all_cohorts, outcome = "CSF AB Ratio")
    },
    pattern = map(all_cohorts)
  )
)
