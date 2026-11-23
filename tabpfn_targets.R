# Exploratory branch in the existing store. When crew is enabled, only these
# fitting targets use the GPU pool. Each fold retains its isolated subprocess.
list(
  tar_target(
    tabpfn_config,
    tabpfn_default_config(),
    cue = tar_cue(mode = "always"),
    packages = character()
  ),
  tar_target(
    tabpfn_adapter_source,
    c("R/tabpfn_prediction.R", "scripts/tabpfn_fold.R"),
    format = "file",
    packages = character()
  ),
  tar_target(
    tabpfn_cached_sl,
    list(
      results = list(
        AD = roc_result_AD_vs_Control,
        FTD = roc_result_FTD_vs_Control,
        LBD = roc_result_LBD_vs_Control
      )
    ),
    packages = character()
  ),
  tar_map(
    values = dplyr::mutate(
      dplyr::filter(
        roc_defs,
        comparison %in% c("AD_vs_Control", "FTD_vs_Control", "LBD_vs_Control")
      ),
      seed_offset = c(100L, 200L, 300L)
    ),
    names = comparison,
    tar_target(
      tabpfn_roc_result,
      tabpfn_auroc(
        roc_data_prepared,
        target_diagnosis,
        comparators,
        utils::modifyList(
          tabpfn_config,
          list(seed = tabpfn_config$seed + seed_offset)
        ),
        tabpfn_adapter_source
      ),
      resources = tar_resources(crew = tar_resources_crew(controller = "gpu")),
      packages = character()
    )
  ),
  # Existing control target names and seeds above remain unchanged.
  tar_map(
    values = dplyr::mutate(
      dplyr::filter(
        roc_defs,
        comparison %in%
          c(
            "AD_vs_Others",
            "LBD_vs_Others",
            "FTD_vs_Others",
            "LBD_vs_FTD"
          )
      ),
      seed_offset = c(400L, 500L, 600L, 700L)
    ),
    names = comparison,
    tar_target(
      tabpfn_roc_result,
      tabpfn_auroc(
        roc_data_prepared,
        target_diagnosis,
        comparators,
        utils::modifyList(
          tabpfn_config,
          list(seed = tabpfn_config$seed + seed_offset)
        ),
        tabpfn_adapter_source
      ),
      resources = tar_resources(crew = tar_resources_crew(controller = "gpu")),
      packages = character()
    )
  ),
  tar_target(
    tabpfn_roc_results,
    list(
      AD = tabpfn_roc_result_AD_vs_Control,
      FTD = tabpfn_roc_result_FTD_vs_Control,
      LBD = tabpfn_roc_result_LBD_vs_Control
    ),
    packages = character()
  ),
  tar_target(
    tabpfn_comparison,
    tabpfn_comparison_plot(tabpfn_roc_results, tabpfn_cached_sl),
    packages = character()
  ),
  tar_target(
    tabpfn_comparison_files,
    tabpfn_export_comparison(tabpfn_comparison),
    format = "file",
    packages = character()
  ),
  tar_target(
    tabpfn_cached_sl_others,
    list(
      results = list(
        AD = roc_result_AD_vs_Others,
        FTD = roc_result_FTD_vs_Others,
        LBD = roc_result_LBD_vs_Others
      )
    ),
    packages = character()
  ),
  tar_target(
    tabpfn_roc_results_others,
    list(
      AD = tabpfn_roc_result_AD_vs_Others,
      FTD = tabpfn_roc_result_FTD_vs_Others,
      LBD = tabpfn_roc_result_LBD_vs_Others
    ),
    packages = character()
  ),
  tar_target(
    tabpfn_comparison_others,
    tabpfn_comparison_plot(
      tabpfn_roc_results_others,
      tabpfn_cached_sl_others,
      "other dementias",
      "Dementia subtypes versus other dementias"
    ),
    packages = character()
  ),
  tar_target(
    tabpfn_comparison_lbd_ftd,
    tabpfn_comparison_plot(
      list(LBD = tabpfn_roc_result_LBD_vs_FTD),
      list(results = list(LBD = roc_result_LBD_vs_FTD)),
      "FTD",
      "LBD/PD versus FTD"
    ),
    packages = character()
  ),
  tar_target(
    tabpfn_comparison_others_files,
    tabpfn_export_comparison(
      tabpfn_comparison_others,
      stem = "subtypes_vs_others"
    ),
    format = "file",
    packages = character()
  ),
  tar_target(
    tabpfn_comparison_lbd_ftd_files,
    tabpfn_export_comparison(tabpfn_comparison_lbd_ftd, stem = "lbd_vs_ftd"),
    format = "file",
    packages = character()
  ),
  tar_target(
    tabpfn_comparison_pngs,
    {
      files <- c(
        tabpfn_comparison_files,
        tabpfn_comparison_others_files,
        tabpfn_comparison_lbd_ftd_files
      )
      files[grepl("\\.png$", files)]
    },
    format = "file",
    packages = character()
  )
)
