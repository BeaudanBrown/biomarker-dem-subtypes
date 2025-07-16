list(
  tar_target(
    pet_corrs,
    get_pet_corrs(standardised_cohorts),
    pattern = map(standardised_cohorts)
  )
)
