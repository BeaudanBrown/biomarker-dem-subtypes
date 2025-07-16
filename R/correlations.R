## Functions for correlations

standardise <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

standardise_df <- function(df) {
  df |>
    mutate(
      mean_ab42_ab40_ratio = standardise(ifelse(
        mean_ab42_ab40_ratio > 0.5,
        NA,
        mean_ab42_ab40_ratio
      )),
      mean_ab40 = standardise(ifelse(mean_ab40 > 500, NA, mean_ab40)),
      mean_ab42 = standardise(ifelse(mean_ab42 > 20, NA, mean_ab42)),
      mean_ykl = standardise(log(mean_ykl)),
      mean_ptau181 = standardise(mean_ptau181),
      mean_tdp = standardise(log(mean_tdp)),
      mean_ptau217 = standardise(mean_ptau217),
      mean_elisa = standardise(mean_elisa),
      mean_gfap = standardise(mean_gfap),
      MMSE = standardise(MMSE),
      cdr = standardise(cdr),
      CSF_NfL = standardise(CSF_NfL),
      CSF_sTREM2 = standardise(CSF_sTREM2),
      CSF_VILIP1 = standardise(CSF_VILIP1),
      CSF_Ng = standardise(CSF_Ng),
      CSF_AB_Ratio = standardise(LUMIPULSE_CSF_AB42 / LUMIPULSE_CSF_AB40),
      CSF_SNAP25 = standardise(CSF_SNAP25),
      zscore = standardise(zscore),
      raw_suvr = standardise(raw_suvr),
      centiloid = standardise(centiloid)
    )
}

get_adjusted_corr <- function(cohort, outcome, predictor) {
  format_estimate <- function(x) sprintf("%.3f", round(x, 3))
  format_pval <- function(p) {
    if (is.na(p)) {
      return(NA)
    }
    if (p < 0.05) {
      base <- "*"
    } else {
      base <- ""
    }
    if (p < 0.001) {
      paste0(base, formatC(p, format = "e", digits = 2))
    } else {
      paste0(base, sprintf("%.3f", round(p, 3)))
    }
  }

  cohort_name <- names(cohort)
  df <- cohort[[cohort_name]]
  filtered_df <- df |>
    filter(!is.na(get(outcome)) & !is.na(get(predictor)))
  model <- lm(
    as.formula(paste0(
      outcome,
      " ~ ",
      paste(predictor, "age", "female", sep = " + ")
    )),
    data = filtered_df
  )
  summary_model <- summary(model)
  coef_predictor <- coef(summary_model)[predictor, "Estimate"]
  pval_predictor <- coef(summary_model)[predictor, "Pr(>|t|)"]

  tibble(
    cohort = cohort_name,
    outcome = outcome,
    predictor = predictor,
    estimate = format_estimate(coef_predictor),
    p = format_pval(pval_predictor),
    n = nrow(filtered_df)
  )
}

make_beta_table <- function(
  data,
  id_col,
  id_label
) {
  idq <- enquo(id_col)
  id_nm <- quo_name(idq)

  data <- data |>
    filter(cohort != "All Subtypes") |>
    group_by(cohort) |>
    mutate(cohort = glue::glue("{cohort} (n={max(n)})")) |>
    ungroup()

  cohorts <- data |>
    distinct(cohort) |>
    pull(cohort)
  nc <- length(cohorts)

  df_wide <- data |>
    pivot_wider(
      id_cols = !!idq,
      names_from = cohort,
      values_from = c(estimate, p),
      names_glue = "{cohort}_{.value}"
    )

  col_order <- c(
    id_nm,
    unlist(lapply(cohorts, function(co) paste0(co, "_", c("estimate", "p"))))
  )
  df_wide <- df_wide[, col_order]

  bottom_header <- c(id_label, rep(c("Beta", "P-value"), times = nc))
  top_header <- c(" " = 1, setNames(rep(2, nc), cohorts))

  df_wide |>
    kableExtra::kable(
      format = "html",
      booktabs = TRUE,
      col.names = bottom_header
    ) |>
    kableExtra::add_header_above(top_header) |>
    kableExtra::kable_styling(latex_options = "striped", full_width = FALSE)
}

get_marker_table <- function(corr_data) {
  corr_data <- corr_data |>
    filter(cohort != "All Subtypes") |>
    group_by(cohort) |>
    mutate(
      cohort = glue::glue("{cohort} (n={max(n)})")
    ) |>
    ungroup()

  cohorts <- corr_data |>
    dplyr::distinct(cohort) |>
    dplyr::pull(cohort)
  nc <- length(cohorts)

  df_wide <- corr_data |>
    tidyr::pivot_wider(
      id_cols = predictor,
      names_from = cohort,
      values_from = c(estimate, p),
      names_glue = "{cohort}_{.value}"
    )

  col_order <- c(
    "predictor",
    unlist(
      lapply(cohorts, function(x) paste0(x, "_", c("estimate", "p")))
    )
  )
  df_wide <- df_wide[, col_order]

  bottom_header <- c("Predictor", rep(c("Beta", "P‐value"), times = nc))
  top_header <- c(" " = 1, setNames(rep(2, nc), cohorts))

  df_wide |>
    kableExtra::kable(
      format = "html",
      booktabs = TRUE,
      col.names = bottom_header
    ) |>
    kableExtra::add_header_above(top_header) |>
    kableExtra::kable_styling(latex_options = "striped", full_width = FALSE)
}

get_marker_corrs <- function(cohort, outcome) {
  markers <- c(
    "mean_elisa",
    "mean_nfl",
    "mean_ykl",
    "mean_gfap",
    "mean_ab42_ab40_ratio",
    "mean_ab40",
    "mean_ab42",
    "mean_tdp",
    "mean_ptau181",
    "mean_ptau217"
  )
  bind_rows(lapply(
    markers,
    function(marker) {
      get_adjusted_corr(
        cohort = cohort,
        outcome = outcome,
        predictor = marker
      )
    }
  ))
}

csf_corr <- function(df, diagnosis = NULL) {
  if (!is.null(diagnosis)) {
    df <- df[df$Diagnosis_combined %in% diagnosis, ]
  }
  washington <- standardise_df(df) |>
    filter(!is.na(mean_ptau217) & !is.na(CSF_AB_Ratio))

  get_ptau_csf_corr <- function(csf_marker) {
    ptau_coef <- washington |>
      lm(
        as.formula(paste0(
          csf_marker,
          " ~ ",
          paste("mean_ptau217", "age", "female", sep = " + ")
        )),
        data = _
      ) |>
      coef()

    ptau_coef["mean_ptau217"]
  }
  csf_markers <- washington |>
    select(CSF_AB_Ratio) |>
    names()
  return(list(
    n = nrow(washington),
    coefs = sapply(csf_markers, get_ptau_csf_corr)
  ))
}

csf_rank_corr <- function(df, diagnosis = NULL) {
  if (!is.null(diagnosis)) {
    df <- df[df$Diagnosis_combined %in% diagnosis, ]
  }

  get_ptau_csf_corr <- function(df, csf_marker, plasma_marker) {
    # complete observations
    df <- filter(
      df,
      !is.na(.data[[csf_marker]]) & !is.na(.data[[plasma_marker]])
    )

    corr <- cor(
      df[[csf_marker]],
      df[[plasma_marker]],
      use = "complete.obs",
      method = "spearman"
    )

    corr <- data.frame(
      csf_marker = csf_marker,
      plasma_marker = plasma_marker,
      corr = corr,
      n = nrow(df)
    )

    return(corr)
  }

  grid <- expand_grid(
    csf_marker = c(
      "LUMIPULSE_CSF_AB42",
      "LUMIPULSE_CSF_tTau",
      "LUMIPULSE_CSF_pTau",
      "LUMIPULSE_CSF_pTau_AB42"
    ),
    plasma_marker = c(
      "mean_ptau181",
      "mean_ab42",
      "mean_nfl",
      "mean_ab40",
      "mean_gfap",
      "mean_ykl",
      "mean_tdp",
      "mean_ptau217",
      "mean_elisa"
    )
  )

  out <- pmap(
    list(
      grid$csf_marker,
      grid$plasma_marker
    ),
    get_ptau_csf_corr,
    df = df
  )

  out <- bind_rows(out) |>
    as_tibble() |>
    arrange(csf_marker) |>
    mutate(
      csf_marker = gsub("LUMIPULSE_", "", csf_marker),
      plasma_marker = gsub("mean_", "", plasma_marker)
    )

  return(out)
}

get_pet_corrs <- function(cohort) {
  measures <- c(
    "zscore",
    "raw_suvr",
    "centiloid"
  )
  bind_rows(lapply(
    measures,
    function(outcome) {
      get_adjusted_corr(cohort, outcome, predictor = "mean_ptau217")
    }
  ))
}

get_csf_corrs <- function(cohort) {
  get_adjusted_corr(
    cohort,
    outcome = "CSF_AB_Ratio",
    predictor = "mean_ptau217"
  )
}
