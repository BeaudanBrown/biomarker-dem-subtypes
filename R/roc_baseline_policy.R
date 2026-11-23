# Freeze historical primary ROC results (control and differential diagnoses).
# Separate opt-ins preserve the original control policy.
# A cue alone is insufficient: targets can execute when metadata is absent or
# the last attempt errored. The guard must run before any fitting call.
roc_refit_controls_allowed <- function(variable = "ADDF_REFIT_SL_CONTROLS") {
  value <- tolower(trimws(Sys.getenv(variable, "false")))
  if (!value %in% c("true", "false")) {
    stop(variable, " must be true or false.", call. = FALSE)
  }
  identical(value, "true")
}

roc_guard_historical_controls <- function(
  refit_allowed,
  variable = "ADDF_REFIT_SL_CONTROLS"
) {
  if (!isTRUE(refit_allowed)) {
    stop(
      "Historical SuperLearner results are frozen. ",
      "A cached result/metadata is unavailable or previously errored. ",
      "Restore the existing cache; no fitting was attempted. ",
      "Only if deliberately recomputing, set ",
      variable,
      "=true.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
