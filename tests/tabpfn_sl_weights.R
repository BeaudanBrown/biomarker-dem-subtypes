# Run from the repository root in the analysis R environment. No GPU or raw data.
source("R/tabpfn_analysis.R")
source("R/tabpfn_sl_weights.R")

local({
  y <- rep(0:1, each = 50L)
  prepared <- tabpfn_sl_prepare(data.frame(
    Diagnosis_combined = ifelse(y == 1, "Alzheimer's", "Control"),
    marker = sin(seq_along(y)) + y,
    age = cos(seq_along(y))
  ))
  baseline <- tabpfn_sl_baseline_fold(prepared, 1L, "SL.glm", inner_folds = 3L)
  reference <- withr::with_seed(
    prepared$seed + 1L,
    SuperLearner::SuperLearner(
      Y = prepared$y[baseline$train],
      X = prepared$x[baseline$train, , drop = FALSE],
      family = stats::binomial(),
      SL.library = "SL.glm",
      method = "method.AUC",
      cvControl = list(V = 3L, stratifyCV = TRUE)
    )
  )
  expected <- predict(
    reference,
    prepared$x[baseline$holdout, , drop = FALSE]
  )$pred
  stopifnot(
    isTRUE(all.equal(reference$coef, baseline$weights)),
    isTRUE(all.equal(as.numeric(expected), as.numeric(baseline$heldout))),
    identical(reference$validRows, baseline$inner)
  )

  calls <- 0L
  mock <- function(x, y, newx, config, seed, adapter_source) {
    calls <<- calls + 1L
    if (calls <= length(baseline$inner)) {
      validation <- baseline$inner[[calls]]
      train <- baseline$train[-validation]
      heldout <- baseline$train[validation]
    } else {
      train <- baseline$train
      heldout <- baseline$holdout
    }
    stopifnot(
      identical(x, prepared$x[train, , drop = FALSE]),
      identical(y, prepared$y[train]),
      identical(newx, prepared$x[heldout, , drop = FALSE]),
      !any(rownames(x) %in% rownames(newx))
    )
    list(predictions = stats::plogis(newx$marker))
  }
  augmented <- tabpfn_sl_augmented_fold(prepared, baseline, NULL, NULL, mock)
  weights <- tabpfn_sl_weight_table(list(augmented))
  stopifnot(
    calls == 4L,
    identical(weights$learner, c(names(baseline$weights), "TabPFN")),
    all(!weights$failed),
    abs(sum(weights$weight) - 1) < 1e-8,
    augmented$auc >= 0,
    augmented$auc <= 1
  )
})
