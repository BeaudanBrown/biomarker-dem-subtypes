subset_comparisons <- roc_defs |> filter(str_ends(comparison, "Others"))
subset_comparisons_ftd <- roc_defs |>
  filter(str_ends(comparison, "Others") | str_ends(comparison, "FTD"))
subset_targets <- list(
  tar_map(
    values = subset_comparisons_ftd,
    names = comparison,
    tar_target(
      subset_result_full,
      find_minimal_biomarker_subset(
        df,
        target_diagnosis,
        reference_group = comparators
      )
    ),
    tar_target(
      subset_path_full,
      build_path(subset_result_full)
    ),
    tar_target(
      subset_plot_full,
      plot_auc_steps(subset_path_full, title)
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
  tar_target(subset_plot_combined, combine_subset_plots(subset_plots_full)),
  tar_target(subset_plot_men_combined, combine_subset_plots(subset_plots_men)),
  tar_target(
    subset_plot_women_combined,
    combine_subset_plots(subset_plots_women)
  ),
  tar_map(
    values = subset_comparisons,
    names = comparison,
    tar_target(
      subset_result_men,
      find_minimal_biomarker_subset(
        df,
        target_diagnosis,
        reference_group = comparators,
        stratify_by_sex = "male"
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
        title,
        use_sex = FALSE,
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
      find_minimal_biomarker_subset(
        df,
        target_diagnosis,
        reference_group = comparators,
        stratify_by_sex = "female"
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
        title,
        use_sex = FALSE,
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
      find_minimal_biomarker_subset(
        df,
        target_diagnosis,
        reference_group = comparators,
        adjust_for_cdr = TRUE
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
        title,
        use_cdr = TRUE,
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
