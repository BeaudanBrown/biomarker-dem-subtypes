# Exploratory prediction branch. Existing SuperLearner functions are unchanged.
tabpfn_default_config <- function() {
  list(
    folds = 10L,
    members = 4L,
    seed = 20261116L,
    threads = 3L,
    timeout = 600,
    device = "cuda:0",
    checkpoint = path.expand(
      "~/.local/share/tabpfn/models/tabpfn-v3.5-20260909.safetensors"
    ),
    sha256 = "ece4d67eadfea42eb0e610df5189bea60cb7f31073d81e9c7a019b76eacf0be3",
    rscript = Sys.getenv("ADDF_TABPFN_RSCRIPT")
  )
}

tabpfn_prepare_comparison <- function(data, outcome, reference) {
  # Match auroc(): restrict diagnoses, complete cases, then standardize once
  # across the comparison BEFORE outer CV. This limitation is intentional.
  row_ids <- which(data$Diagnosis_combined %in% c(outcome, reference))
  data <- data[row_ids, , drop = FALSE]
  complete <- stats::complete.cases(data)
  data <- data[complete, , drop = FALSE]
  row_ids <- row_ids[complete]
  y <- as.integer(data$Diagnosis_combined == outcome)
  x <- as.data.frame(data[,
    setdiff(names(data), "Diagnosis_combined"),
    drop = FALSE
  ])
  for (name in names(x)) {
    if (!is.numeric(x[[name]]) || any(!is.finite(x[[name]]))) {
      stop("Primary TabPFN predictors must be finite numeric columns: ", name)
    }
    if (length(unique(x[[name]])) > 5L) {
      x[[name]] <- (x[[name]] - mean(x[[name]])) / stats::sd(x[[name]])
    }
  }
  list(x = x, y = y, row_ids = row_ids)
}

tabpfn_stratified_folds <- function(y, folds, seed) {
  if (
    length(folds) != 1L ||
      folds < 2L ||
      folds != as.integer(folds) ||
      !identical(sort(unique(y)), 0:1) ||
      any(table(y) < folds)
  ) {
    stop("Each class must have at least one participant per validation fold.")
  }
  # Explicit, reproducible new folds; no claim to match historical SL folds.
  withr::with_seed(seed, {
    assignment <- integer(length(y))
    for (label in 0:1) {
      indices <- which(y == label)
      assignment[indices] <- sample(rep(
        seq_len(folds),
        length.out = length(indices)
      ))
    }
    assignment
  })
}

tabpfn_roc_values <- function(labels, predictions) {
  if (
    length(labels) != length(predictions) ||
      !length(labels) ||
      anyNA(labels) ||
      !all(labels %in% 0:1) ||
      length(unique(labels)) != 2L ||
      any(!is.finite(predictions)) ||
      any(predictions < 0 | predictions > 1)
  ) {
    stop("Invalid held-out binary predictions.")
  }
  rates <- vapply(
    seq(0, 1, length.out = 100),
    function(threshold) {
      predicted <- predictions >= threshold
      c(TPR = mean(predicted[labels == 1]), FPR = mean(predicted[labels == 0]))
    },
    numeric(2)
  )
  # Fix orientation: larger probabilities always mean the target diagnosis.
  auc <- as.numeric(pROC::auc(pROC::roc(
    labels,
    predictions,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )))
  list(
    labels = labels,
    preds = predictions,
    TPR = rates["TPR", ],
    FPR = rates["FPR", ],
    AUC = auc
  )
}

tabpfn_run_fold <- function(x, y, newx, config, seed, adapter_source) {
  if (!nzchar(config$rscript) || !file.exists(config$rscript)) {
    stop("GPU R runtime missing: enter nix develop .#tabpfn-gpu first.")
  }
  # Separate R installations: preserve the existing analysis R/packages and
  # let the GPU wrapper select its own pinned R/Python libraries.
  directory <- tempfile("tabpfn-fold-")
  dir.create(directory, mode = "0700")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  input <- file.path(directory, "input.rds")
  output <- file.path(directory, "output.rds")
  saveRDS(list(x = x, y = y, newx = newx, config = config, seed = seed), input)
  Sys.chmod(input, "0600")
  sources <- normalizePath(adapter_source, mustWork = TRUE)
  processx::run(
    config$rscript,
    c(
      "--vanilla",
      sources[basename(sources) == "tabpfn_fold.R"],
      input,
      output,
      sources[basename(sources) == "tabpfn_prediction.R"]
    ),
    timeout = config$timeout,
    error_on_status = TRUE
  )
  # Python/GPU objects never cross the process/cache boundary. A fresh process
  # per fold releases model memory even after an error or timeout.
  readRDS(output)
}

tabpfn_auroc <- function(
  data,
  outcome,
  reference,
  config,
  adapter_source,
  predict_fold = tabpfn_run_fold
) {
  prepared <- tabpfn_prepare_comparison(data, outcome, reference)
  fold_ids <- tabpfn_stratified_folds(prepared$y, config$folds, config$seed)
  results <- vector("list", config$folds)
  for (fold in seq_len(config$folds)) {
    training <- which(fold_ids != fold)
    heldout <- which(fold_ids == fold)
    seed <- config$seed + fold
    start <- proc.time()[["elapsed"]]
    fit <- predict_fold(
      prepared$x[training, , drop = FALSE],
      prepared$y[training],
      prepared$x[heldout, , drop = FALSE],
      config,
      seed,
      adapter_source
    )
    values <- tabpfn_roc_values(prepared$y[heldout], fit$predictions)
    results[[fold]] <- c(
      values,
      list(
        fold = fold,
        seed = seed,
        row_ids = prepared$row_ids[heldout],
        training_row_ids = prepared$row_ids[training],
        seconds = proc.time()[["elapsed"]] - start,
        diagnostics = fit$diagnostics
      )
    )
  }
  list(
    aucs = list(results = results),
    # IDs are positions in roc_data_prepared, NOT participant identifiers.
    row_ids = prepared$row_ids,
    fold_ids = fold_ids,
    outcome = outcome,
    reference = reference,
    config = config,
    preprocessing = "Existing upstream truncation; complete cases and whole-comparison scaling before CV"
  )
}

tabpfn_curve_data <- function(results, learner) {
  labels <- c(AD = "AD", FTD = "FTD", LBD = "LBD/PD")
  purrr::imap_dfr(results, function(result, comparison) {
    folds <- result$aucs$results
    if (!length(folds)) stop("No cached ROC folds for ", comparison)
    data.frame(
      learner = learner,
      subtype = unname(labels[comparison]),
      TPR = Reduce(`+`, lapply(folds, `[[`, "TPR")) / length(folds),
      FPR = Reduce(`+`, lapply(folds, `[[`, "FPR")) / length(folds),
      AUC = mean(vapply(folds, function(fold) mean(fold$AUC), numeric(1)))
    )
  })
}

tabpfn_comparison_plot <- function(
  tabpfn_results,
  cached_sl,
  reference_label = "control",
  title = "Dementia subtypes versus control"
) {
  curves <- rbind(
    tabpfn_curve_data(cached_sl$results, "SuperLearner (cached)"),
    tabpfn_curve_data(tabpfn_results, "TabPFN")
  )
  palette <- c(AD = "#0072B2", FTD = "#D55E00", `LBD/PD` = "#009E73")
  panel <- function(learner) {
    data <- curves[curves$learner == learner, ]
    aucs <- unique(data[, c("subtype", "AUC")])
    legend <- stats::setNames(
      sprintf("%s vs %s (AUC: %.2f)", aucs$subtype, reference_label, aucs$AUC),
      aucs$subtype
    )
    ggplot2::ggplot(data, ggplot2::aes(FPR, TPR, colour = subtype)) +
      ggplot2::geom_path(linewidth = 1) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        alpha = 0.5
      ) +
      ggplot2::scale_colour_manual(
        values = palette,
        breaks = names(palette)[names(palette) %in% aucs$subtype],
        labels = legend[names(palette)[names(palette) %in% aucs$subtype]]
      ) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(
        title = learner,
        x = "False Positive Rate (1 - Specificity)",
        y = "True Positive Rate (Sensitivity)",
        colour = NULL
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))
  }
  plot <- patchwork::wrap_plots(
    panel("SuperLearner (cached)"),
    panel("TabPFN"),
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = title,
      caption = "Exploratory comparison: existing preprocessing before CV; historical SuperLearner and new TabPFN folds need not match. AUC = mean fold AUC."
    )
  list(plot = plot, curves = curves)
}

tabpfn_export_comparison <- function(
  comparison,
  directory = file.path("output", "tabpfn"),
  stem = "subtypes_vs_control"
) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  png <- file.path(directory, paste0(stem, "_comparison.png"))
  csv <- file.path(directory, paste0(stem, "_curves.csv"))
  ggplot2::ggsave(
    png,
    comparison$plot,
    width = 14,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  utils::write.csv(comparison$curves, csv, row.names = FALSE)
  c(png, csv)
}
