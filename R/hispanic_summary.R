# Describe Hispanic ethnicity without changing the analysis cohort.
# Percentages within each diagnosis use Yes + No + Unknown as the denominator.
# Unknown/Unsure are combined; missing/blank values are excluded.
hispanic_summary <- function(df) {
  if (!"Hispanic" %in% names(df)) {
    stop("Hispanic column is required for the ethnicity summary.")
  }
  diagnoses <- c("Control", "Alzheimer's", "Frontotermporal", "Lewy bodies")
  diagnosis_labels <- c("Controls", "AD", "FTD", "LBD/PD")
  responses <- tolower(trimws(as.character(df$Hispanic)))
  if (any(!is.na(responses) & !responses %in% c("", "yes", "no", "unknown", "unsure"))) {
    stop("Unexpected Hispanic response encoding; confirm the coding before summarising.")
  }
  df <- dplyr::mutate(
    df,
    hispanic_response = dplyr::case_when(
      responses == "yes" ~ "Yes",
      responses == "no" ~ "No",
      responses %in% c("unknown", "unsure") ~ "Unknown",
      TRUE ~ NA_character_
    )
  )
  df <- dplyr::filter(df, Diagnosis_combined %in% diagnoses)
  categories <- c("Yes", "No", "Unknown")

  result <- tibble::tibble(
    Summary = paste("Hispanic", categories, "n (%)")
  )
  for (i in seq_along(diagnoses)) {
    group <- dplyr::filter(df, Diagnosis_combined == diagnoses[i])
    denominator <- sum(group$hispanic_response %in% categories)
    result[[diagnosis_labels[i]]] <- vapply(seq_along(categories), function(j) {
      n <- sum(group$hispanic_response == categories[j], na.rm = TRUE)
      if (denominator == 0L) {
        "0 (—)"
      } else {
        sprintf("%d (%.1f%%)", n, 100 * n / denominator)
      }
    }, character(1))
  }
  result
}
