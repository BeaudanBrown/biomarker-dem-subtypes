roc_defs <- bind_rows(
  list(
    comparison = "AD_vs_Control",
    target_diagnosis = "Alzheimer's",
    comparators = list(c("Control")),
    title = "Alzheimer's vs control"
  ),
  list(
    comparison = "FTD_vs_Control",
    target_diagnosis = "Frontotermporal",
    comparators = list(c("Control")),
    title = "Frontotermporal vs control"
  ),
  list(
    comparison = "LBD_vs_Control",
    target_diagnosis = "Lewy bodies",
    comparators = list(c("Control")),
    title = "Lewy bodies vs control"
  ),
  list(
    comparison = "LBD_vs_FTD",
    target_diagnosis = "Lewy bodies",
    comparators = list(c("Frontotermporal")),
    title = "Lewy bodies vs frontotemporal"
  ),
  list(
    comparison = "AD_vs_Others",
    target_diagnosis = "Alzheimer's",
    comparators = list(c("Frontotermporal", "Lewy bodies")),
    title = "Alzheimer's vs other dementias"
  ),
  list(
    comparison = "LBD_vs_Others",
    target_diagnosis = "Lewy bodies",
    comparators = list(c("Frontotermporal", "Alzheimer's")),
    title = "Lewy bodies vs other dementias"
  ),
  list(
    comparison = "FTD_vs_Others",
    target_diagnosis = "Frontotermporal",
    comparators = list(c("Lewy bodies", "Alzheimer's")),
    title = "Frontotemporal vs other dementias"
  )
)

roc_targets <- list(
  tar_target(roc_data_prepared, prepare_roc_data(df_with_csf)),
  tar_target(
    roc_data_prepared_fasting,
    prepare_roc_data(df, with_fasting = "yes")
  ),

  # ===================================================================
  # Primary analysis
  # ===================================================================
  tar_map(
    values = roc_defs,
    names = comparison,
    tar_target(
      roc_result,
      auroc(roc_data_prepared, target_diagnosis, comparators)
    ),
    tar_target(
      roc_plot,
      plot_roc(roc_result, title)
    ),
    tar_target(
      roc_merged,
      list(roc_result = roc_result, roc_plot = roc_plot)
    )
  ),
  tar_target(
    roc_results,
    list(
      AD_vs_Control = roc_merged_AD_vs_Control,
      FTD_vs_Control = roc_merged_FTD_vs_Control,
      LBD_vs_Control = roc_merged_LBD_vs_Control,
      LBD_vs_FTD = roc_merged_LBD_vs_FTD,
      AD_vs_Others = roc_merged_AD_vs_Others,
      LBD_vs_Others = roc_merged_LBD_vs_Others,
      FTD_vs_Others = roc_merged_FTD_vs_Others
    )
  ),
  tar_target(combined_roc, get_combined_roc(roc_results)),

  # ===================================================================
  # Fasting Adjusted
  # ===================================================================
  tar_map(
    values = roc_defs,
    names = comparison,
    tar_target(
      roc_fasting_result,
      auroc(roc_data_prepared_fasting, target_diagnosis, comparators)
    ),
    tar_target(
      roc_fasting_plot,
      plot_roc(roc_fasting_result, title)
    ),
    tar_target(
      roc_fasting_merged,
      list(roc_result = roc_fasting_result, roc_plot = roc_fasting_plot)
    )
  ),
  tar_target(
    roc_results_fasting,
    list(
      AD_vs_Control = roc_fasting_merged_AD_vs_Control,
      FTD_vs_Control = roc_fasting_merged_FTD_vs_Control,
      LBD_vs_Control = roc_fasting_merged_LBD_vs_Control,
      LBD_vs_FTD = roc_fasting_merged_LBD_vs_FTD,
      AD_vs_Others = roc_fasting_merged_AD_vs_Others,
      LBD_vs_Others = roc_fasting_merged_LBD_vs_Others,
      FTD_vs_Others = roc_fasting_merged_FTD_vs_Others
    )
  ),
  tar_target(combined_roc_fasting, get_combined_roc(roc_results_fasting)),

  # ===================================================================
  # By Sex
  # Men
  # ===================================================================
  tar_map(
    values = roc_defs,
    names = comparison,
    tar_target(
      roc_men_result,
      auroc(
        roc_data_prepared |> filter(female == 0) |> select(-female),
        target_diagnosis,
        comparators
      )
    ),
    tar_target(
      roc_men_plot,
      plot_roc(roc_men_result, title)
    ),
    tar_target(
      roc_men_merged,
      list(roc_result = roc_men_result, roc_plot = roc_men_plot)
    )
  ),
  tar_target(
    roc_results_men,
    list(
      AD_vs_Control = roc_men_merged_AD_vs_Control,
      FTD_vs_Control = roc_men_merged_FTD_vs_Control,
      LBD_vs_Control = roc_men_merged_LBD_vs_Control,
      LBD_vs_FTD = roc_men_merged_LBD_vs_FTD,
      AD_vs_Others = roc_men_merged_AD_vs_Others,
      LBD_vs_Others = roc_men_merged_LBD_vs_Others,
      FTD_vs_Others = roc_men_merged_FTD_vs_Others
    )
  ),

  # ===================================================================
  # Women
  # ===================================================================
  tar_map(
    values = roc_defs,
    names = comparison,
    tar_target(
      roc_women_result,
      auroc(
        roc_data_prepared |> filter(female == 1) |> select(-female),
        target_diagnosis,
        comparators
      )
    ),
    tar_target(
      roc_women_plot,
      plot_roc(roc_women_result, title)
    ),
    tar_target(
      roc_women_merged,
      list(roc_result = roc_women_result, roc_plot = roc_women_plot)
    )
  ),
  tar_target(
    roc_results_women,
    list(
      AD_vs_Control = roc_women_merged_AD_vs_Control,
      FTD_vs_Control = roc_women_merged_FTD_vs_Control,
      LBD_vs_Control = roc_women_merged_LBD_vs_Control,
      LBD_vs_FTD = roc_women_merged_LBD_vs_FTD,
      AD_vs_Others = roc_women_merged_AD_vs_Others,
      LBD_vs_Others = roc_women_merged_LBD_vs_Others,
      FTD_vs_Others = roc_women_merged_FTD_vs_Others
    )
  ),
  tar_target(
    sex_rocs,
    list(
      men = roc_results_men,
      women = roc_results_women
    )
  ),
  tar_target(
    sex_specific_rocs,
    rocs_by_sex(df, sex_rocs)
  )
)
