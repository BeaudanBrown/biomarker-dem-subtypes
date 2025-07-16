subset_comparisons <- roc_defs |> filter(str_ends(comparison, "Others"))
subset_targets <- list(
  tar_map(
    values = subset_comparisons,
    names = comparison,
    tar_target(
      subset_result_full,
      run_single_subset(
        df,
        target_diagnosis,
        comparators
      )
    ),
    tar_target(
      subset_path_full,
      build_path(subset_result_full)
    ),
    tar_target(
      subset_plot_full,
      plot_auc_steps(subset_path_full, comparison)
    )
  ),
  tar_target(
    subset_plots_full,
    list(
      ad = subset_plot_full_AD_vs_Others,
      ftd = subset_plot_full_FTD_vs_Others,
      lbd = subset_plot_full_LBD_vs_Others
    )
  ),
  tar_map(
    values = subset_comparisons,
    names = comparison,
    tar_target(
      subset_result_men,
      run_single_subset(
        df,
        target_diagnosis,
        comparators,
        sex_strat = "male"
      )
    ),
    tar_target(
      subset_path_men,
      build_path(subset_result_men)
    ),
    tar_target(
      subset_plot_men,
      plot_auc_steps(
        subset_path_men,
        target_diagnosis,
        extra_title = " - Men"
      )
    )
  ),
  tar_target(
    subset_plots_men,
    list(
      ad = subset_plot_men_AD_vs_Others,
      ftd = subset_plot_men_FTD_vs_Others,
      lbd = subset_plot_men_LBD_vs_Others
    )
  ),
  tar_map(
    values = subset_comparisons,
    names = comparison,
    tar_target(
      subset_result_women,
      run_single_subset(
        df,
        target_diagnosis,
        comparators,
        sex_strat = "female"
      )
    ),
    tar_target(
      subset_path_women,
      build_path(subset_result_women)
    ),
    tar_target(
      subset_plot_women,
      plot_auc_steps(
        subset_path_women,
        target_diagnosis,
        extra_title = " - Women"
      )
    )
  ),
  tar_target(
    subset_plots_women,
    list(
      ad = subset_plot_women_AD_vs_Others,
      ftd = subset_plot_women_FTD_vs_Others,
      lbd = subset_plot_women_LBD_vs_Others
    )
  ),
  tar_map(
    values = subset_comparisons,
    names = comparison,
    tar_target(
      subset_result_cdr,
      run_single_subset(
        df,
        target_diagnosis,
        comparators,
        use_cdr = TRUE
      )
    ),
    tar_target(
      subset_path_cdr,
      build_path(subset_result_cdr)
    ),
    tar_target(
      subset_plot_cdr,
      plot_auc_steps(
        subset_path_cdr,
        target_diagnosis,
        extra_title = " - CDR Adjusted"
      )
    )
  ),
  tar_target(
    subset_plots_cdr,
    list(
      ad = subset_plot_cdr_AD_vs_Others,
      ftd = subset_plot_cdr_FTD_vs_Others,
      lbd = subset_plot_cdr_LBD_vs_Others
    )
  )
)
