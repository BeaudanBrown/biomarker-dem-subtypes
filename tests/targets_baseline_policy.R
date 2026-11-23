# USER-RUN ONLY: isolated temporary targets stores, synthetic values, no models/data.
# From the repository root: source("tests/targets_baseline_policy.R")
source("R/roc_baseline_policy.R")

run_baseline_checks <- function() {
  expect_error <- function(expression) {
    error <- tryCatch(
      {
        force(expression)
        NULL
      },
      error = identity
    )
    stopifnot(inherits(error, "error"))
  }
  withr::with_envvar(c(ADDF_REFIT_SL_CONTROLS = NA), {
    stopifnot(!roc_refit_controls_allowed())
    expect_error(roc_guard_historical_controls(FALSE))
  })
  withr::with_envvar(c(ADDF_REFIT_SL_CONTROLS = "true"), {
    stopifnot(roc_refit_controls_allowed(), roc_guard_historical_controls(TRUE))
  })
  withr::with_envvar(c(ADDF_REFIT_SL_CONTROLS = "maybe"), {
    expect_error(roc_refit_controls_allowed())
  })

  directory <- tempfile("targets-baseline-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  script <- file.path(directory, "_targets.R")
  store <- file.path(directory, "store")
  marker <- file.path(directory, "fit-was-called")
  quoted <- function(x) encodeString(x, quote = '"')
  header <- c(
    "library(targets)",
    "tar_option_set(format = 'rds', packages = character(), error = 'stop')",
    paste0("source(", quoted(normalizePath("R/roc_baseline_policy.R")), ")")
  )
  make <- function()
    targets::tar_make(
      script = script,
      store = store,
      callr_function = NULL,
      use_crew = FALSE,
      reporter = "silent"
    )
  writeLines(
    c(
      header,
      "list(tar_target(baseline, 42L), tar_target(comparison, list(value = baseline)))"
    ),
    script
  )
  make()
  stopifnot(identical(
    targets::tar_read_raw("comparison", store = store),
    list(value = 42L)
  ))

  frozen <- c(
    header,
    "list(",
    "tar_target(baseline, { roc_guard_historical_controls(FALSE);",
    paste0(
      "writeLines('FIT', ",
      quoted(marker),
      "); 43L }, cue = tar_cue(mode = 'never')),"
    ),
    "tar_target(comparison, list(value = baseline)))"
  )
  writeLines(frozen, script)
  make() # Changed command/dependencies must not rebuild the successful baseline.
  stopifnot(
    !file.exists(marker),
    targets::tar_read_raw("baseline", store = store) == 42L
  )

  # Existing metadata with missing output must not cause a replacement fit.
  targets::tar_delete(names = baseline, store = store)
  # Force a downstream load; an unchanged, already-cached comparison could skip.
  writeLines(
    sub(
      "list(value = baseline)",
      "list(value = baseline + 0L)",
      frozen,
      fixed = TRUE
    ),
    script
  )
  expect_error(make())
  stopifnot(!file.exists(marker))

  # A new store has no metadata. The first failed guard records an error; the
  # second attempt therefore also exercises the previous-error exception.
  store <- file.path(directory, "empty-store")
  expect_error(make())
  expect_error(make())
  stopifnot(!file.exists(marker))

  writeLines(
    c(
      header,
      "list(",
      "tar_target(baseline, { roc_guard_historical_controls(TRUE);",
      paste0(
        "writeLines('FIT', ",
        quoted(marker),
        "); 43L }, cue = tar_cue(mode = 'thorough')),"
      ),
      "tar_target(comparison, list(value = baseline)))"
    ),
    script
  )
  make() # Explicit opt-in is the only route to this synthetic fitting stand-in.
  stopifnot(
    file.exists(marker),
    targets::tar_read_raw("comparison", store = store)$value == 43L
  )

  # Inspect the real target declarations without sourcing _targets.R or data.
  manifest_script <- file.path(directory, "manifest.R")
  writeLines(
    c(
      "library(targets); library(tarchetypes); library(dplyr)",
      "tar_option_set(format = 'rds', packages = character())",
      paste0("source(", quoted(normalizePath("R/roc_baseline_policy.R")), ")"),
      paste0("source(", quoted(normalizePath("roc_targets.R")), ")"),
      paste0(
        "c(roc_targets, source(",
        quoted(normalizePath("tabpfn_targets.R")),
        ")$value)"
      )
    ),
    manifest_script
  )
  withr::with_envvar(
    c(ADDF_REFIT_SL_CONTROLS = "false", ADDF_REFIT_SL_SUBTYPES = "false"),
    {
      manifest <- targets::tar_manifest(
        script = manifest_script,
        fields = c(name, command, cue_mode),
        callr_function = NULL
      )
      stopifnot(!anyDuplicated(manifest$name))
      baseline_names <- paste0(
        "roc_result_",
        c("AD", "FTD", "LBD"),
        "_vs_Control"
      )
      subtype_names <- paste0(
        "roc_result_",
        c("AD_vs_Others", "FTD_vs_Others", "LBD_vs_Others", "LBD_vs_FTD")
      )
      stopifnot(
        all(c(baseline_names, subtype_names) %in% manifest$name),
        setequal(
          manifest$name[manifest$cue_mode == "never"],
          c(baseline_names, subtype_names)
        ),
        all(
          paste0("tabpfn_", c(baseline_names, subtype_names)) %in% manifest$name
        ),
        "tabpfn_comparison_pngs" %in% manifest$name,
        !any(grepl("^tabpfn_report", manifest$name))
      )
      aggregate <- manifest$command[manifest$name == "tabpfn_cached_sl"]
      stopifnot(
        length(aggregate) == 1L,
        !grepl("tar_read", aggregate),
        all(vapply(
          baseline_names,
          grepl,
          logical(1),
          x = aggregate,
          fixed = TRUE
        ))
      )
    }
  )
  withr::with_envvar(
    c(ADDF_REFIT_SL_CONTROLS = "true", ADDF_REFIT_SL_SUBTYPES = "true"),
    {
      manifest <- targets::tar_manifest(
        script = manifest_script,
        fields = c(name, cue_mode),
        callr_function = NULL
      )
      stopifnot(all(
        manifest$cue_mode[grepl("^roc_result_", manifest$name)] == "thorough"
      ))
    }
  )
}
run_baseline_checks()
message("Synthetic frozen-baseline and normal-dependency checks passed.")
