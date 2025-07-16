list(
  # declare input data
  tar_target(
    pheno_file,
    file.path(data_dir, Sys.getenv("PHENO_FILE")),
    format = "file"
  ),
  tar_target(
    gilipin_file,
    file.path(data_dir, Sys.getenv("GILIPIN_FILE")),
    format = "file"
  ),
  tar_target(
    all_cases_file,
    file.path(data_dir, Sys.getenv("ALL_CASES_FILE")),
    format = "file"
  ),
  tar_target(
    addf_file,
    file.path(data_dir, Sys.getenv("ADDF_FILE")),
    format = "file"
  ),
  tar_target(
    texas_file,
    file.path(data_dir, Sys.getenv("TEXAL_PRELIM_FILE")),
    format = "file"
  ),
  tar_target(
    texas_file_update,
    file.path(data_dir, Sys.getenv("TEXAL_UPDATE_FILE")),
    format = "file"
  ),
  tar_target(
    addf_cd14_file,
    file.path(data_dir, Sys.getenv("ADDF_CD14_FILE")),
    format = "file"
  ),
  tar_target(
    addf_csf_file,
    file.path(data_dir, Sys.getenv("ADDF_CSF_FILE")),
    format = "file"
  ),
  tar_target(
    nacc_file,
    file.path(data_dir, Sys.getenv("NACC_FILE")),
    format = "file"
  ),
  tar_target(
    wash_cogs_file,
    file.path(data_dir, Sys.getenv("NEW_WASH_COGS")),
    format = "file"
  ),
  # merge input data
  tar_target(
    joined,
    merge_datafiles(
      pheno_file,
      gilipin_file,
      all_cases_file,
      addf_file,
      texas_file,
      texas_file_update,
      addf_cd14_file,
      wash_cogs_file,
      addf_csf_file,
      nacc_file
    )
  )
)
