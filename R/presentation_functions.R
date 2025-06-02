## Presentation functions

markers <- c(
  "CD14",
  "NfL",
  "YKL-40",
  "GFAP",
  "AB42/AB40 Ratio",
  "AB40",
  "AB42",
  "TDP-43",
  "pTau-181",
  "pTau-217"
)

dem_names <- c(
  "Alzheimer's",
  "Frontotermporal",
  "Lewy bodies"
)


## Get all ROC curves

all_rocs <- function(df, with_fasting = "no") {
  if (with_fasting == "yes") {
    setDT(df)

    df[,
      fasting_combined := ifelse(
        is.na(Fasting_Status),
        Fasting_8_hours,
        Fasting_Status
      )
    ][,
      fasting_combined := ifelse(
        fasting_combined == "Unknown" |
          fasting_combined == "unknown",
        NA,
        fasting_combined
      )
    ][,
      fasting_combined := ifelse(
        fasting_combined == "non-fasted",
        "No",
        fasting_combined
      )
    ] |>
      as_tibble()

    roc_df <-
      df |>
      select(
        Diagnosis_combined,
        age_combined,
        mean_elisa,
        mean_nfl,
        mean_ykl,
        mean_gfap,
        mean_ab40,
        mean_ab42,
        mean_tdp,
        mean_ptau181,
        mean_ptau217,
        female,
        fasting_combined
      )
  } else {
    roc_df <-
      df |>
      select(
        Diagnosis_combined,
        age_combined,
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
  }

  get_all_rocs(roc_df)
}

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
    "LBD" = "LBD vs others (AUC: %s)"
  )

  plot_roc_combined(
    roc_list = roc_combined,
    title_text = "ROC curves for dementia subtypes",
    label_map = outcome_label_map
  )
}

## ROCs by sex

rocs_by_sex <- function(df) {
  roc_df <-
    df |>
    select(
      Diagnosis_combined,
      age_combined,
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

  sex_rec_results <- get_all_sex_rocs(roc_df)
  sex_label_map <- c(
    "MEN" = "Men (AUC: %s)",
    "WOMEN" = "Women (AUC: %s)"
  )

  #### AD ####
  ad_plot <- plot_roc_combined(
    roc_list = list(
      MEN = sex_rec_results$men$AD_vs_Others$roc_result,
      WOMEN = sex_rec_results$women$AD_vs_Others$roc_result
    ),
    title_text = "AD vs other dementias",
    label_map = sex_label_map
  )

  #### LBD ####
  lbd_plot <- plot_roc_combined(
    roc_list = list(
      MEN = sex_rec_results$men$LBD_vs_Others$roc_result,
      WOMEN = sex_rec_results$women$LBD_vs_Others$roc_result
    ),
    title_text = "LBD vs other dementias",
    label_map = sex_label_map
  )

  #### FTD ####
  ftd_plot <- plot_roc_combined(
    roc_list = list(
      MEN = sex_rec_results$men$FTD_vs_Others$roc_result,
      WOMEN = sex_rec_results$women$FTD_vs_Others$roc_result
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

  return(list(gender_diag_table, ad_plot, lbd_plot, ftd_plot))
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
    "LBD" = "LBD vs control (AUC: %s)"
  )
  plot_roc_combined(
    roc_list = roc_combined,
    title_text = "ROC curves for dementia subtypes",
    label_map = outcome_label_map
  )
}

## NOT USED CURRENTLY

print_pet_table <- function(pet_output, caption) {
  df <- as.data.frame(t(pet_output$coefs))
  rownames(df) <- NULL
  df <- df |>
    rename(
      "Z-score" = "zscore.mean_ptau217",
      "Raw SUVR" = "raw_suvr.mean_ptau217",
      "Centiloid" = "centiloid.mean_ptau217"
    )

  knitr::kable(
    df,
    digits = 3,
    caption = paste0(caption, " (n = ", pet_output$n, ")")
  )
}

## MMSE table

mmse_table <- function(df) {
  whole_inflam <- inflam_cdr(df)
  Any_inflam <- inflam_cdr(
    df,
    diagnosis = c("Alzheimer's", "Frontotermporal", "Lewy bodies")
  )
  AD_inflam <- inflam_cdr(df, diagnosis = c("Alzheimer's"))
  FTD_inflam <- inflam_cdr(df, diagnosis = c("Frontotermporal"))
  LBD_inflam <- inflam_cdr(df, diagnosis = c("Lewy bodies"))

  mmse_list <- list(
    "All (n)" = as.data.frame(whole_inflam$mmse),
    "Any Dem (n)" = as.data.frame(Any_inflam$mmse),
    "AD (n)" = as.data.frame(AD_inflam$mmse),
    "LBD (n)" = as.data.frame(LBD_inflam$mmse),
    "FTD (n)" = as.data.frame(FTD_inflam$mmse)
  )

  mmse <- bind_cols(mmse_list) |>
    setNames(names(mmse_list)) |>
    rownames_to_column("Marker") |>
    relocate(Marker) |>
    mutate(
      Marker = case_when(
        Marker == "mean_elisa" ~ "CD14",
        Marker == "mean_nfl" ~ "NfL",
        Marker == "mean_ykl" ~ "YKL-40",
        Marker == "mean_gfap" ~ "GFAP",
        Marker == "mean_ab42_ab40_ratio" ~ "AB42/AB40 Ratio",
        Marker == "mean_ab40" ~ "AB40",
        Marker == "mean_ab42" ~ "AB42",
        Marker == "mean_tdp" ~ "TDP-43",
        Marker == "mean_ptau181" ~ "pTau-181",
        Marker == "mean_ptau217" ~ "pTau-217",
        TRUE ~ Marker
      )
    )
  table <-
    knitr::kable(
      mmse,
      digits = 3,
      caption = "Inflamatory Markers and MMSE (age/sex adjusted)"
    )

  plots <- list()

  for (marker in markers) {
    plots[[marker]] <-
      df |>
      rename(
        "CD14" = "mean_elisa",
        "NfL" = "mean_nfl",
        "YKL-40" = "mean_ykl",
        "GFAP" = "mean_gfap",
        "AB42/AB40 Ratio" = "mean_ab42_ab40_ratio",
        "AB40" = "mean_ab40",
        "AB42" = "mean_ab42",
        "TDP-43" = "mean_tdp",
        "pTau-181" = "mean_ptau181",
        "pTau-217" = "mean_ptau217"
      ) |>
      filter(!is.na(get(marker)) & !is.na(MMSE)) |>
      filter(Site == "Washington") |>
      ggplot(aes(x = get(marker), y = MMSE, color = Diagnosis_combined)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(color = "Diagnosis") +
      xlab(marker) +
      ggtitle(paste0(marker, " vs MMSE (not adjusted)"))
  }

  return(list(table, plots))
}

## CDR vs markers

cdr_vs_markers <- function(df) {
  whole_inflam <- inflam_cdr(df)
  Any_inflam <- inflam_cdr(
    df,
    diagnosis = c("Alzheimer's", "Frontotermporal", "Lewy bodies")
  )
  AD_inflam <- inflam_cdr(df, diagnosis = c("Alzheimer's"))
  FTD_inflam <- inflam_cdr(df, diagnosis = c("Frontotermporal"))
  LBD_inflam <- inflam_cdr(df, diagnosis = c("Lewy bodies"))

  cdr_list <- list(
    "All" = as.data.frame(whole_inflam$cdr),
    "Any Dem" = as.data.frame(Any_inflam$cdr),
    "AD" = as.data.frame(AD_inflam$cdr),
    "LBD" = as.data.frame(LBD_inflam$cdr),
    "FTD" = as.data.frame(FTD_inflam$cdr)
  )

  cdr <- bind_cols(cdr_list) |>
    setNames(names(cdr_list)) |>
    rownames_to_column("Marker") |>
    relocate(Marker) |>
    mutate(
      Marker = case_when(
        Marker == "mean_elisa" ~ "CD14",
        Marker == "mean_nfl" ~ "NfL",
        Marker == "mean_ykl" ~ "YKL-40",
        Marker == "mean_gfap" ~ "GFAP",
        Marker == "mean_ab42_ab40_ratio" ~ "AB42/AB40 Ratio",
        Marker == "mean_ab40" ~ "AB40",
        Marker == "mean_ab42" ~ "AB42",
        Marker == "mean_tdp" ~ "TDP-43",
        Marker == "mean_ptau181" ~ "pTau-181",
        Marker == "mean_ptau217" ~ "pTau-217",
        TRUE ~ Marker
      )
    )
  cdr_table <-
    knitr::kable(
      cdr,
      digits = 3,
      caption = "Inflamatory Markers and CDR (age/sex adjusted)"
    )

  plots <- list()

  for (marker in markers) {
    plots[[marker]] <- df |>
      rename(
        "CD14" = "mean_elisa",
        "NfL" = "mean_nfl",
        "YKL-40" = "mean_ykl",
        "GFAP" = "mean_gfap",
        "AB42/AB40 Ratio" = "mean_ab42_ab40_ratio",
        "AB40" = "mean_ab40",
        "AB42" = "mean_ab42",
        "TDP-43" = "mean_tdp",
        "pTau-181" = "mean_ptau181",
        "pTau-217" = "mean_ptau217"
      ) |>
      filter(!is.na(get(marker)) & !is.na(cdr)) |>
      filter(Site == "Washington") |>
      ggplot(aes(x = get(marker), y = cdr, color = Diagnosis_combined)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(color = "Diagnosis") +
      xlab(marker) +
      ggtitle(paste0(marker, " vs CDR (not adjusted)"))
  }

  return(list(cdr_table = cdr_table, plots = plots))
}


## CSF AB vs plasma markers plot

csf_vs_markers <- function(df) {
  cohort_info <- list(
    "All" = NULL,
    "Any Dem" = c("Alzheimer's", "Frontotermporal", "Lewy bodies"),
    "Alzheimer's" = c("Alzheimer's"),
    "Frontotermporal" = c("Frontotermporal"),
    "Lewy bodies" = c("Lewy bodies")
  )

  # Loop over each cohort, run csf_corr(), extract the AB ratio + n,
  # and record the cohort name.
  results <- lapply(names(cohort_info), function(cohort) {
    diagnosis <- cohort_info[[cohort]]
    csf_out <- csf_corr(df, diagnosis = diagnosis)
    tmp_df <- as.data.frame(t(csf_out$coefs))
    ab_ratio <- round(tmp_df$CSF_AB_Ratio.mean_ptau217, 3)
    list(Cohort = cohort, "AB42/AB40 Ratio" = ab_ratio, n = csf_out$n)
  })
  final_df <- do.call(rbind, results)
  knitr::kable(
    final_df,
    caption = "pTau-217 and CSF AB-Ratio (age/sex adjusted)"
  )

  df |>
    mutate(CSF_AB_Ratio = LUMIPULSE_CSF_AB42 / LUMIPULSE_CSF_AB40) |>
    filter(!is.na(mean_ptau217) & !is.na(CSF_AB_Ratio)) |>
    filter(Site == "Washington") |>
    ggplot(aes(
      x = mean_ptau217,
      y = CSF_AB_Ratio,
      color = Diagnosis_combined
    )) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      color = "Diagnosis",
      title = "CSF AB-Ratio vs pTau-217 (not age/sex adjusted)",
      x = "Mean pTau-217",
      y = "CSF AB-Ratio"
    )
}


## PET vs PTAU217

pet_vs_ptau <- function(df) {
  whole_pet <- pet_corr(df)
  Any_pet <- pet_corr(
    df,
    diagnosis = c("Alzheimer's", "Frontotermporal", "Lewy bodies")
  )
  AD_pet <- pet_corr(df, diagnosis = c("Alzheimer's"))
  FTD_pet <- pet_corr(df, diagnosis = c("Frontotermporal"))
  LBD_pet <- pet_corr(df, diagnosis = c("Lewy bodies"))

  pet <- bind_cols(
    as.data.frame(whole_pet$coefs),
    as.data.frame(Any_pet$coefs),
    as.data.frame(AD_pet$coefs),
    as.data.frame(FTD_pet$coefs),
    as.data.frame(LBD_pet$coefs)
  ) |>
    setNames(c(
      paste0("All (n = ", whole_pet$n, ")"),
      paste0("Any Dem (n = ", Any_pet$n, ")"),
      paste0("AD (n = ", AD_pet$n, ")"),
      paste0("LBD (n = ", FTD_pet$n, ")"),
      paste0("FTD (n = ", LBD_pet$n, ")")
    )) |>
    rownames_to_column("Measure") |>
    mutate(
      Measure = case_when(
        Measure == "zscore" ~ "Z-score",
        Measure == "raw_suvr" ~ "Raw SUVR",
        Measure == "centiloid" ~ "Centiloid",
        TRUE ~ Measure
      )
    ) |>
    relocate(Measure)

  pet_table <-
    knitr::kable(
      pet,
      digits = 3,
      caption = "PET Measures and pTau-217 (age/sex adjusted)"
    )

  plot1 <- df |>
    rename(
      "raw_suvr" = "av45/PIB_fsuvr_rsf_tot_cortmean",
      "centiloid" = "Centiloid_AV45_fSUVR_TOT_CORTMEA",
    ) |>
    filter(!is.na(mean_ptau217) & !is.na(raw_suvr)) |>
    filter(Site == "Washington") |>
    ggplot(aes(x = mean_ptau217, y = raw_suvr, color = Diagnosis_combined)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      color = "Diagnosis",
      title = "Raw SUVR AB-Ratio vs pTau-217 (not age/sex adjusted)",
      x = "Mean pTau-217",
      y = "SUVR"
    )

  plot2 <- df |>
    rename(
      "raw_suvr" = "av45/PIB_fsuvr_rsf_tot_cortmean",
      "centiloid" = "Centiloid_AV45_fSUVR_TOT_CORTMEA",
    ) |>
    filter(!is.na(mean_ptau217) & !is.na(centiloid)) |>
    filter(Site == "Washington") |>
    ggplot(aes(x = mean_ptau217, y = centiloid, color = Diagnosis_combined)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      color = "Diagnosis",
      title = "Centiloid AB-Ratio vs pTau-217 (not age/sex adjusted)",
      x = "Mean pTau-217",
      y = "Centiloid"
    )

  return(list(pet_table, plot1, plot2))
}

## VIMP overall

vimp_overall <- function(df) {
  vimp_data <- df |>
    filter(Diagnosis_combined != "Control") |>
    select(
      Diagnosis_combined,
      age_combined,
      mean_elisa,
      mean_nfl,
      mean_ykl,
      mean_gfap,
      mean_ab40,
      mean_ab42,
      mean_tdp,
      mean_ptau181,
      mean_ptau217
    ) |>
    drop_na()

  vimp_full <- lapply(
    dem_names,
    function(outcome) {
      vimp_function_par(
        data = vimp_data,
        outcome_subtype = outcome
      )$mat |>
        select(c(s, est, test))
    }
  ) |>
    setNames(dem_names)

  preds <- vimp_data |>
    select(-Diagnosis_combined) |>
    names()

  tables <- list()

  for (dem_name in dem_names) {
    tables[[dem_name]] <- kable(
      cbind(
        preds[as.numeric(vimp_full[[dem_name]]$s)],
        vimp_full[[dem_name]]$est,
        vimp_full[[dem_name]]$test
      ) |>
        as_tibble() |>
        setNames(c(
          "var_name",
          "imp",
          "test"
        )) |>
        mutate(imp = round(as.numeric(imp), 3)) |>
        mutate(
          imp = paste0(imp, ifelse(test, "*", ""))
        ) |>
        select(-test) |>
        mutate(
          var_name = case_when(
            var_name == "mean_elisa" ~ "CD14",
            var_name == "mean_nfl" ~ "NfL",
            var_name == "mean_ykl" ~ "YKL-40",
            var_name == "mean_gfap" ~ "GFAP",
            var_name == "mean_ab42_ab40_ratio" ~ "AB42/AB40 Ratio",
            var_name == "mean_ab40" ~ "AB40",
            var_name == "mean_ab42" ~ "AB42",
            var_name == "mean_tdp" ~ "TDP-43",
            var_name == "mean_ptau181" ~ "pTau-181",
            var_name == "mean_ptau217" ~ "pTau-217",
            TRUE ~ var_name
          )
        ),
      col.names = c("Predictor", "Importance"),
      caption = paste0(
        "VImp for ",
        dem_name,
        " vs Other Dem",
        " (n = ",
        nrow(filter(vimp_data, Diagnosis_combined == dem_name)),
        ")"
      )
    )
  }

  return(tables)
}

## VIMP by sex

vimp_by_sex <- function(df, sex) {
  plan(multicore, workers = detectCores())

  vimp_data <- df |>
    filter(
      female == ifelse(sex == "female", 1, 0) & Diagnosis_combined != "Control"
    ) |>
    select(
      Diagnosis_combined,
      age_combined,
      mean_elisa,
      mean_nfl,
      mean_ykl,
      mean_gfap,
      mean_ab40,
      mean_ab42,
      mean_tdp,
      mean_ptau181,
      mean_ptau217
    ) |>
    drop_na()

  vimp <- lapply(
    dem_names,
    function(outcome) {
      vimp_function_par(
        data = vimp_data,
        outcome_subtype = outcome,
        stratification = ifelse(sex == "female", "female", "male")
      )$mat |>
        select(c(s, est, test))
    }
  ) |>
    setNames(dem_names)

  preds <- vimp_data |>
    select(-Diagnosis_combined) |>
    names()

  tables <- list()

  for (dem_name in dem_names) {
    tables[[dem_name]] <-
      kable(
        cbind(
          preds[as.numeric(vimp[[dem_name]]$s)],
          vimp[[dem_name]]$est,
          vimp[[dem_name]]$test
        ) |>
          as_tibble() |>
          setNames(c(
            "var_name",
            "imp",
            "test"
          )) |>
          mutate(imp = round(as.numeric(imp), 3)) |>
          mutate(
            imp = paste0(imp, ifelse(test, "*", ""))
          ) |>
          select(-test) |>
          mutate(
            var_name = case_when(
              var_name == "mean_elisa" ~ "CD14",
              var_name == "mean_nfl" ~ "NfL",
              var_name == "mean_ykl" ~ "YKL-40",
              var_name == "mean_gfap" ~ "GFAP",
              var_name == "mean_ab42_ab40_ratio" ~ "AB42/AB40 Ratio",
              var_name == "mean_ab40" ~ "AB40",
              var_name == "mean_ab42" ~ "AB42",
              var_name == "mean_tdp" ~ "TDP-43",
              var_name == "mean_ptau181" ~ "pTau-181",
              var_name == "mean_ptau217" ~ "pTau-217",
              TRUE ~ var_name
            )
          ),
        col.names = c("Predictor", "Importance"),
        caption = paste0(
          "VImp for ",
          dem_name,
          ifelse(
            sex == "female",
            " vs Other Dem in Females",
            " vs Other Dem in Males"
          ),
          " (n = ",
          nrow(filter(vimp_data, Diagnosis_combined == dem_name)),
          ")"
        )
      )
  }

  return(tables)
}
