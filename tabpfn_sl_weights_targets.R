# Opt-in AD vs Control pilot, independent of the historical ROC fit targets.
list(
  tar_target(ad_sl_prepared, tabpfn_sl_prepare(roc_data_prepared)),
  tar_target(ad_sl_fold, seq_len(10L), iteration = "vector"),
  tar_target(
    ad_sl_baseline_fold,
    tabpfn_sl_baseline_fold(ad_sl_prepared, ad_sl_fold, SL.library),
    pattern = map(ad_sl_fold),
    iteration = "list"
  ),
  tar_target(
    ad_sl_baseline_weights,
    tabpfn_sl_weight_table(ad_sl_baseline_fold)
  ),
  tar_target(
    ad_sl_baseline_weights_file,
    tabpfn_sl_export_weights(ad_sl_baseline_weights),
    format = "file"
  ),
  tar_target(
    ad_sl_augmented_fold,
    tabpfn_sl_augmented_fold(
      ad_sl_prepared,
      ad_sl_baseline_fold,
      tabpfn_config,
      tabpfn_adapter_source
    ),
    pattern = map(ad_sl_baseline_fold),
    iteration = "list",
    resources = tar_resources(crew = tar_resources_crew(controller = "gpu"))
  ),
  tar_target(
    ad_sl_augmented_weights,
    tabpfn_sl_weight_table(ad_sl_augmented_fold)
  ),
  tar_target(
    ad_sl_paired_weights_files,
    tabpfn_sl_export_weights(ad_sl_baseline_weights, ad_sl_augmented_weights),
    format = "file"
  )
)
