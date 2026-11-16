# Descriptive summaries only: use the same df population as Table 1.
# Percentages describe the diagnostic distribution of observed scores per panel.
cognitive_summary <- function(df) {
  diagnoses <- c("Control", "Alzheimer's", "Frontotermporal", "Lewy bodies")
  diagnosis_labels <- c("Controls", "AD", "FTD", "LBD/PD")
  panels <- list(
    All = df,
    Males = dplyr::filter(df, female == 0),
    Females = dplyr::filter(df, female == 1)
  )

  summarise_score <- function(x, denominator) {
    observed <- x[!is.na(x)]
    n <- length(observed)
    summary <- if (n == 0L) {
      "—"
    } else {
      quartiles <- stats::quantile(observed, probs = c(0.25, 0.5, 0.75))
      sprintf("%.1f (%.1f, %.1f)", quartiles[2], quartiles[1], quartiles[3])
    }
    availability <- if (denominator == 0L) {
      "0 (—)"
    } else {
      sprintf("%d (%.1f%%)", n, 100 * n / denominator)
    }
    c(summary, availability)
  }

  purrr::imap_dfr(panels, function(panel, panel_name) {
    panel <- dplyr::filter(panel, Diagnosis_combined %in% diagnoses)
    cdr_denominator <- sum(!is.na(panel$cdr))
    mmse_denominator <- sum(!is.na(panel$MMSE))
    result <- tibble::tibble(
      Panel = panel_name,
      Summary = c(
        "CDR median (Q1, Q3)", "CDR available n (%)",
        "MMSE median (Q1, Q3)", "MMSE available n (%)"
      )
    )
    for (i in seq_along(diagnoses)) {
      group <- dplyr::filter(panel, Diagnosis_combined == diagnoses[i])
      result[[diagnosis_labels[i]]] <- c(
        summarise_score(group$cdr, cdr_denominator),
        summarise_score(group$MMSE, mmse_denominator)
      )
    }
    result
  })
}
