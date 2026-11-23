# This pilot matches the exploratory branch's pre-CV preprocessing.
# It does not read or update any historical SuperLearner results.
tabpfn_sl_prepare <- function(data, seed = 20261216L) {
  prepared <- tabpfn_prepare_comparison(data, "Alzheimer's", "Control")
  prepared$outer <- tabpfn_stratified_folds(prepared$y, 10L, seed)
  prepared$seed <- seed
  prepared
}

tabpfn_sl_baseline_fold <- function(
  prepared,
  fold,
  library,
  inner_folds = 20L
) {
  train <- which(prepared$outer != fold)
  holdout <- which(prepared$outer == fold)
  fit <- withr::with_seed(
    prepared$seed + fold,
    SuperLearner::SuperLearner(
      Y = prepared$y[train],
      X = prepared$x[train, , drop = FALSE],
      newX = prepared$x[holdout, , drop = FALSE],
      family = stats::binomial(),
      SL.library = library,
      method = "method.AUC",
      cvControl = list(V = inner_folds, stratifyCV = TRUE),
      control = list(saveFitLibrary = FALSE),
      env = environment(tabpfn_sl_baseline_fold)
    )
  )
  failed <- as.logical(fit$errorsInCVLibrary | fit$errorsInLibrary) |
    !apply(is.finite(fit$Z), 2, all) |
    !apply(is.finite(fit$library.predict), 2, all)
  if (any(!is.finite(fit$coef)) || any(fit$coef[failed] != 0)) {
    stop("Invalid baseline weights in outer fold ", fold)
  }
  list(
    fold = fold,
    train = train,
    holdout = holdout,
    inner = fit$validRows,
    oof = fit$Z,
    heldout = fit$library.predict,
    weights = fit$coef,
    failed = failed,
    auc = tabpfn_roc_values(prepared$y[holdout], as.numeric(fit$SL.predict))$AUC
  )
}

tabpfn_sl_augmented_fold <- function(
  prepared,
  baseline,
  config,
  adapter_source,
  predict_fold = tabpfn_run_fold
) {
  train <- baseline$train
  holdout <- baseline$holdout
  y <- prepared$y[train]
  inner <- baseline$inner
  rows <- unlist(inner, use.names = FALSE)
  if (length(rows) != length(train) || !setequal(rows, seq_along(train))) {
    stop("Inner folds must predict each outer training row exactly once.")
  }
  predict_rows <- function(training, validation, step) {
    probs <- predict_fold(
      prepared$x[training, , drop = FALSE],
      prepared$y[training],
      prepared$x[validation, , drop = FALSE],
      config,
      prepared$seed + 100L * baseline$fold + step,
      adapter_source
    )$predictions
    if (
      length(probs) != length(validation) ||
        any(!is.finite(probs)) ||
        any(probs < 0 | probs > 1)
    ) {
      stop(
        "Invalid TabPFN probabilities in outer fold ",
        baseline$fold,
        ", step ",
        step
      )
    }
    probs
  }
  oof <- rep(NA_real_, length(train))
  for (step in seq_along(inner)) {
    validation <- inner[[step]]
    oof[validation] <- predict_rows(train[-validation], train[validation], step)
  }
  # Exclude flagged learners rather than letting the optimizer assign weight
  # to a placeholder prediction. Preserve their slots in the exported table.
  keep <- !baseline$failed
  learners <- c(names(baseline$weights), "TabPFN")
  active <- c(keep, TRUE)
  Z <- cbind(baseline$oof[, keep, drop = FALSE], TabPFN = oof)
  method <- SuperLearner::method.AUC()
  control <- SuperLearner::SuperLearner.control()
  meta <- method$computeCoef(
    Z = Z,
    Y = y,
    libraryNames = learners[active],
    obsWeights = rep(1, length(y)),
    control = control,
    verbose = FALSE,
    errorsInLibrary = rep(FALSE, ncol(Z))
  )
  if (
    any(!is.finite(meta$coef)) ||
      any(meta$coef < 0) ||
      abs(sum(meta$coef) - 1) > 1e-8
  ) {
    stop("Invalid augmented weights in outer fold ", baseline$fold)
  }
  # Score the untouched outer holdout only after learning the weights.
  heldout <- cbind(
    baseline$heldout[, keep, drop = FALSE],
    TabPFN = predict_rows(train, holdout, length(inner) + 1L)
  )
  predictions <- as.numeric(method$computePred(heldout, meta$coef, control))
  weights <- stats::setNames(numeric(length(learners)), learners)
  weights[active] <- meta$coef
  list(
    fold = baseline$fold,
    weights = weights,
    failed = !active,
    auc = tabpfn_roc_values(prepared$y[holdout], predictions)$AUC
  )
}

tabpfn_sl_weight_table <- function(folds) {
  do.call(
    rbind,
    lapply(folds, function(fit) {
      data.frame(
        fold = fit$fold,
        learner = names(fit$weights),
        weight = unname(fit$weights),
        failed = fit$failed,
        auc = fit$auc
      )
    })
  )
}

tabpfn_sl_plot_weights <- function(weights) {
  ensembles <- c("Original library", "With TabPFN")
  folds <- sort(unique(weights$fold))
  learners <- sort(unique(as.character(weights$learner)))
  learners <- c(setdiff(learners, "TabPFN"), intersect("TabPFN", learners))
  palette <- stats::setNames(
    grDevices::hcl.colors(length(learners), "Dark 3"),
    learners
  )
  weights$learner <- factor(weights$learner, levels = learners)
  # Offset the two stacks around the fold tick. Stacking (not dodging) then
  # applies independently at each offset; missing/zero weights take no space.
  weights$bar_x <- match(weights$fold, folds) +
    c(-0.20, 0.20)[match(weights$ensemble, ensembles)]

  ggplot2::ggplot(weights, ggplot2::aes(bar_x, weight, fill = learner)) +
    ggplot2::geom_col(
      width = 0.36,
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::scale_x_continuous(breaks = seq_along(folds), labels = folds) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_percent()
    ) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    # Read the legend top to bottom, like the vertical stacks.
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, reverse = TRUE)) +
    ggplot2::labs(
      title = "AD vs Control: SuperLearner weights by outer fold",
      subtitle = "Within each fold: original library on the left, with TabPFN on the right.",
      caption = "Zero-weight models have no visible segment. See the paired weights CSV for fit failures.",
      x = "Outer fold",
      y = "Ensemble weight",
      fill = "Model"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}

tabpfn_sl_export_weights <- function(baseline, augmented = NULL) {
  directory <- file.path("output", "tabpfn")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  baseline$ensemble <- "Original library"
  if (is.null(augmented)) {
    file <- file.path(directory, "ad_vs_control_sl_baseline_weights.csv")
    utils::write.csv(baseline, file, row.names = FALSE)
    return(file)
  }
  augmented$ensemble <- "With TabPFN"
  weights <- rbind(baseline, augmented)
  weights$ensemble <- factor(
    weights$ensemble,
    levels = c("Original library", "With TabPFN")
  )
  auc <- merge(
    unique(baseline[c("fold", "auc")]),
    unique(augmented[c("fold", "auc")]),
    by = "fold",
    suffixes = c("_baseline", "_augmented")
  )
  auc$auc_delta <- auc$auc_augmented - auc$auc_baseline
  tabpfn <- augmented$weight[augmented$learner == "TabPFN"]
  summary <- data.frame(
    tabpfn_weight_median = stats::median(tabpfn),
    tabpfn_weight_min = min(tabpfn),
    tabpfn_weight_max = max(tabpfn),
    tabpfn_nonzero_folds = sum(tabpfn > 0),
    folds = length(tabpfn),
    mean_paired_auc_delta = mean(auc$auc_delta)
  )
  files <- file.path(
    directory,
    paste0(
      "ad_vs_control_sl_",
      c(
        "paired_weights.csv",
        "paired_auc.csv",
        "tabpfn_summary.csv",
        "paired_weights.png"
      )
    )
  )
  utils::write.csv(weights, files[1], row.names = FALSE)
  utils::write.csv(auc, files[2], row.names = FALSE)
  utils::write.csv(summary, files[3], row.names = FALSE)
  plot <- tabpfn_sl_plot_weights(weights)
  ggplot2::ggsave(files[4], plot, width = 11, height = 6, dpi = 150)
  files
}
