list(
  tar_target(
    pet_corrs,
    get_pet_corrs(standardised_cohorts),
    pattern = map(standardised_cohorts)
  ),
  tar_target(
    pet_corr_table,
    get_pet_table(pet_corrs)
  ),
  tar_target(
    pet_plots,
    get_pet_plots(all_cohorts),
    pattern = cross(all_cohorts)
  )
)
