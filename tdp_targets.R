# Independent exploratory branch: tar_make(names = tdp43_report).
# No changes to existing ROC, PET/CSF regression or descriptive targets.
list(
  tar_target(tdp43_data, tdp43_prepare(df_with_csf)),

  # Diagnostic distributions, followed by the six exploratory analyses.
  tar_target(
    tdp43_subtype,
    tdp43_comparison(tdp43_data, "Diagnosis", "Diagnostic group", "subtype", by_diagnosis = FALSE)
  ),
  tar_target(tdp43_age, tdp43_scatter(tdp43_data, "age", "Age (years)", "age")),
  tar_target(tdp43_sex, tdp43_comparison(tdp43_data, "Sex", "Sex", "sex")),
  tar_target(tdp43_biomarker_correlations, tdp43_correlations(tdp43_data)),
  tar_target(tdp43_cdr, tdp43_comparison(tdp43_data, "cdr", "Global CDR", "cdr")),
  tar_target(tdp43_pet, tdp43_scatter(tdp43_data, "centiloid", "Amyloid PET (centiloids)", "pet_centiloid")),
  tar_target(tdp43_csf, tdp43_scatter(tdp43_data, "csf_ab42_ab40", "CSF Aβ42/40", "csf_ab42_ab40")),

  # Track all PNG/CSV files so deleted exports are regenerated.
  tar_target(
    tdp43_exports,
    tdp43_export(list(tdp43_subtype, tdp43_age, tdp43_sex, tdp43_biomarker_correlations, tdp43_cdr, tdp43_pet, tdp43_csf)),
    format = "file"
  ),
  # Track presentation sources explicitly, including the browser-only controls.
  tar_target(
    tdp43_widget_sources,
    c("R/tdp_widgets.R", "R/tdp_exploration.R", "R/tdp_controls.js"),
    format = "file"
  ),
  tar_quarto(
    tdp43_report_render,
    path = "R/tdp_exploration.qmd",
    quiet = FALSE
  ),
  tar_target(tdp43_report, tdp43_publish_report(tdp43_report_render), format = "file")
)
