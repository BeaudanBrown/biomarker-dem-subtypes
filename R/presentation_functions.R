## Presentation functions

dem_names <- c(
  "Alzheimer's",
  "Frontotermporal",
  "Lewy bodies"
)

## combined roc curve

get_combined_roc <- function(roc_results) {
  roc_combined <- list(
    AD = roc_results$AD_vs_Others[["roc_result"]],
    LBD = roc_results$LBD_vs_Others[["roc_result"]],
    FTD = roc_results$FTD_vs_Others[["roc_result"]]
  )

  outcome_label_map <- c(
    "AD" = "AD vs others (AUC: %s)",
    "FTD" = "FTD vs others (AUC: %s)",
    "LBD" = "LBD/PD vs others (AUC: %s)"
  )

  plot_roc_combined(
    roc_list = roc_combined,
    title_text = "ROC curves for dementia subtypes",
    label_map = outcome_label_map
  )
}

combine_subtype_control <- function(subtype_rocs, control_rocs) {
  subtype_rocs <- subtype_rocs +
    ggtitle("A") +
    theme(
      legend.spacing.y = unit(1, "cm"),
      legend.key.height = unit(1, "cm"),
      plot.title = element_text(face = "bold")
    )
  control_rocs <- control_rocs +
    ggtitle("B") +
    theme(
      legend.spacing.y = unit(1, "cm"),
      legend.key.height = unit(1, "cm"),
      plot.title = element_text(face = "bold")
    )
  subtype_rocs / control_rocs
}

## ROCs by sex

rocs_by_sex <- function(df, sex_roc_results) {
  roc_df <-
    df |>
    select(
      Diagnosis_combined,
      age,
      mean_elisa,
      mean_nfl,
      mean_ykl,
      mean_gfap,
      mean_ab40,
      mean_ab42,
      mean_tdp,
      mean_ptau181,
      mean_ptau217,
      female
    )

  sex_label_map <- c(
    "MEN" = "Men (AUC: %s)",
    "WOMEN" = "Women (AUC: %s)"
  )

  #### AD ####
  ad_plot <- plot_roc_combined(
    roc_list = list(
      MEN = sex_roc_results$men$AD_vs_Others$roc_result,
      WOMEN = sex_roc_results$women$AD_vs_Others$roc_result
    ),
    title_text = "AD vs other dementias",
    label_map = sex_label_map
  )

  #### LBD ####
  lbd_plot <- plot_roc_combined(
    roc_list = list(
      MEN = sex_roc_results$men$LBD_vs_Others$roc_result,
      WOMEN = sex_roc_results$women$LBD_vs_Others$roc_result
    ),
    title_text = "LBD vs other dementias",
    label_map = sex_label_map
  )

  #### FTD ####
  ftd_plot <- plot_roc_combined(
    roc_list = list(
      MEN = sex_roc_results$men$FTD_vs_Others$roc_result,
      WOMEN = sex_roc_results$women$FTD_vs_Others$roc_result
    ),
    title_text = "FTD vs other dementias",
    label_map = sex_label_map
  )

  gender_diag_table <- roc_df |>
    mutate(Gender = ifelse(female == 1, "Women", "Men")) |>
    group_by(Gender, Diagnosis_combined) |>
    summarise(Count = n(), .groups = "drop") |>
    pivot_wider(names_from = Gender, values_from = Count, values_fill = 0) |>
    mutate(
      Total = Men + Women,
      Men = sprintf("%d (%.1f%%)", Men, Men / Total * 100),
      Women = sprintf("%d (%.1f%%)", Women, Women / Total * 100)
    ) |>
    select(Diagnosis_combined, Men, Women, Total)

  # If you want a prettier table for publication/reports
  gender_diag_table <- kable(
    gender_diag_table,
    caption = "Distribution of men and women across diagnostic groups",
    col.names = c("Diagnosis", "Men N (%)", "Women N (%)", "Total")
  )

  list(
    gender_diag_table = gender_diag_table,
    ad_plot = ad_plot,
    lbd_plot = lbd_plot,
    ftd_plot = ftd_plot
  )
}

combine_sex_roc_plots <- function(sex_specific_plots) {
  ad <- sex_specific_plots[["ad_plot"]] + ggtitle("AD Vs Other Dementias")
  lbd <- sex_specific_plots[["lbd_plot"]] +
    ggtitle("LBD/PD Vs Other Dementias")
  ftd <- sex_specific_plots[["ftd_plot"]] + ggtitle("FTD Vs Other Dementias")
  ad / lbd / ftd
}

## Dementia subtypes vs control

subtypes_vs_control <- function(roc_results) {
  roc_combined <- list(
    AD = roc_results$AD_vs_Control[["roc_result"]],
    LBD = roc_results$LBD_vs_Control[["roc_result"]],
    FTD = roc_results$FTD_vs_Control[["roc_result"]]
  )
  outcome_label_map <- c(
    "AD" = "AD vs control (AUC: %s)",
    "FTD" = "FTD vs control (AUC: %s)",
    "LBD" = "LBD/PD vs control (AUC: %s)"
  )
  plot_roc_combined(
    roc_list = roc_combined,
    title_text = "ROC curves for dementia subtypes",
    label_map = outcome_label_map
  )
}

## Descriptives
show_descriptives <- function(df) {
  # density plot
  density_plot <- df |>
    filter(!Diagnosis_combined == "Control") |>
    # remove large ab42/ab40 ratio outlier
    mutate(
      mean_ab42_ab40_ratio = ifelse(
        mean_ab42_ab40_ratio > 0.5,
        NA,
        mean_ab42_ab40_ratio
      ),
      mean_ab40 = ifelse(mean_ab40 > 500, NA, mean_ab40),
      mean_ab42 = ifelse(mean_ab42 > 20, NA, mean_ab42)
    ) |>
    select(Diagnosis_combined, starts_with("mean_")) |>
    mutate(across(c("mean_ykl", "mean_tdp", "mean_nfl"), log)) |>
    pivot_longer(
      -Diagnosis_combined,
      names_to = "biomarker",
      names_prefix = "mean_"
    ) |>
    mutate(
      biomarker = fct_recode(
        biomarker,
        "AB40" = "ab40",
        "AB42" = "ab42",
        "GFAP" = "gfap",
        "log(NfL)" = "nfl",
        "pTau181" = "ptau181",
        "pTau217" = "ptau217",
        "log(TDP-43)" = "tdp",
        "log(YKL-40)" = "ykl",
        "AB42/AB40" = "ab42_ab40_ratio"
      )
    ) |>
    ggplot(aes(x = value, fill = Diagnosis_combined)) +
    geom_density(alpha = 0.7, size = 0.3) +
    facet_wrap(~biomarker, scales = "free", ncol = 3) +
    labs(
      x = "Biomarker concentration",
      y = "Density",
      fill = "Diagnosis"
    ) +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey95", size = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      panel.spacing = unit(0.5, "cm")
    )

  ## Reliability

  reliability <- df |>
    select(Diagnosis_combined, ends_with("_ICC")) |>
    pivot_longer(-Diagnosis_combined, names_to = "biomarker") |>
    filter(!biomarker == "elisa_ICC") |>
    mutate(
      biomarker = fct_recode(
        biomarker,
        "AB40" = "ab40_ICC",
        "GFAP" = "gfap_ICC",
        "NfL" = "nfl_ICC",
        "pTau181" = "ptau181_ICC",
        "pTau217" = "ptau217_ICC",
        "TDP-43" = "tdp_ICC",
        "YKL-40" = "ykl_ICC"
      )
    ) |>
    mutate(value = round(value, 3)) |>
    group_by(biomarker) |>
    summarise(ICC = mean(value, na.rm = TRUE)) |>
    knitr::kable()

  ## Demo table

  dt <- setDT(df)

  dt[,
    race_combined := ifelse(
      race_combined == "American Indian/Alaska Native",
      "Other",
      race_combined
    )
  ]
  demos <- dt[, .(
    Diagnosis_combined,
    Site,
    age,
    female,
    race_combined,
    APOE,
    mean_ptau181,
    mean_nfl,
    mean_ab40,
    mean_ab42,
    mean_ab42_ab40_ratio,
    mean_gfap,
    mean_ykl,
    mean_tdp,
    mean_ptau217,
    mean_elisa
  )] |>
    gtsummary::tbl_summary(
      by = Diagnosis_combined,
      missing = "ifany"
    )

  return(list(
    density_plot = density_plot,
    reliability = reliability,
    demos = demos
  ))
}
