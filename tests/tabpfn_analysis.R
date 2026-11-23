# USER-RUN ONLY: synthetic, offline, no TabPFN/Python/GPU imports or project data.
# Run from the repository root: source("tests/tabpfn_analysis.R")
source("R/tabpfn_analysis.R")

expect_error <- function(expression) {
  failed <- tryCatch(
    {
      force(expression)
      FALSE
    },
    error = function(e) TRUE
  )
  stopifnot(failed)
}

fixture <- data.frame(
  Diagnosis_combined = rep(
    c("Alzheimer's", "Control", "Lewy bodies"),
    each = 20
  ),
  age = seq_len(60) + 30,
  mean_nfl = seq_len(60) / 10,
  female = rep(0:1, 30)
)
# Complete-case exclusions match the existing ROC path. Other diagnoses and
# row identifiers must never enter the learner's features.
fixture$age[2] <- NA_real_
prepared <- tabpfn_prepare_comparison(fixture, "Alzheimer's", "Control")
stopifnot(
  nrow(prepared$x) == 39L,
  identical(prepared$row_ids, c(1L, 3:40)),
  identical(names(prepared$x), c("age", "mean_nfl", "female")),
  abs(mean(prepared$x$age)) < 1e-12,
  abs(stats::sd(prepared$x$mean_nfl) - 1) < 1e-12,
  all(prepared$x$female %in% 0:1)
)
set.seed(9)
before <- .Random.seed
folds <- tabpfn_stratified_folds(prepared$y, 10L, 100L)
stopifnot(
  identical(before, .Random.seed),
  identical(folds, tabpfn_stratified_folds(prepared$y, 10L, 100L)),
  all(table(prepared$y, folds) >= 1L)
)
expect_error(tabpfn_stratified_folds(c(0L, 1L), 10L, 100L))
expect_error(tabpfn_stratified_folds(rep(0L, 20), 10L, 100L))
expect_error(tabpfn_stratified_folds(prepared$y, 2.5, 100L))

perfect <- tabpfn_roc_values(c(0L, 1L), c(0.1, 0.9))
reversed <- tabpfn_roc_values(c(0L, 1L), c(0.9, 0.1))
stopifnot(perfect$AUC == 1, reversed$AUC == 0, length(perfect$TPR) == 100L)
expect_error(tabpfn_roc_values(c(0L, 1L), c(0.2, Inf)))
expect_error(tabpfn_roc_values(c(0L, 1L), c(0.2, 1.1)))
expect_error(tabpfn_roc_values(c(0L, 1L), 0.2))

calls <- 0L
fake_predict <- function(x, y, newx, config, seed, adapter_source) {
  calls <<- calls + 1L
  stopifnot(
    nrow(x) + nrow(newx) == 39L,
    identical(names(x), names(newx)),
    length(y) == nrow(x),
    length(unique(y)) == 2L
  )
  list(predictions = rep(0.5, nrow(newx)), diagnostics = list(fake = TRUE))
}
config <- tabpfn_default_config()
result <- tabpfn_auroc(
  fixture,
  "Alzheimer's",
  "Control",
  config,
  "unused",
  predict_fold = fake_predict
)
stopifnot(calls == 10L, length(result$aucs$results) == 10L)
heldout <- unlist(lapply(result$aucs$results, `[[`, "row_ids"))
stopifnot(identical(sort(heldout), prepared$row_ids), !anyDuplicated(heldout))
for (fold in result$aucs$results) {
  stopifnot(
    !length(intersect(fold$row_ids, fold$training_row_ids)),
    length(fold$preds) == length(fold$labels),
    length(fold$preds) == length(fold$row_ids),
    fold$AUC == 0.5,
    isTRUE(fold$diagnostics$fake)
  )
}
results <- list(AD = result, FTD = result, LBD = result)
curves <- tabpfn_curve_data(results, "test")
stopifnot(
  nrow(curves) == 300L,
  all(curves$AUC == 0.5),
  setequal(unique(curves$subtype), c("AD", "FTD", "LBD/PD"))
)
plot <- tabpfn_comparison_plot(results, list(results = results))
stopifnot(inherits(plot$plot, "patchwork"), nrow(plot$curves) == 600L)

others_plot <- tabpfn_comparison_plot(
  results,
  list(results = results),
  "other dementias",
  "Differential diagnosis"
)
pair_plot <- tabpfn_comparison_plot(
  list(LBD = result),
  list(results = list(LBD = result)),
  "FTD",
  "LBD/PD versus FTD"
)
stopifnot(
  nrow(others_plot$curves) == 600L,
  nrow(pair_plot$curves) == 200L,
  all(pair_plot$curves$subtype == "LBD/PD")
)

# Frozen historical baseline/dependency behavior is tested separately in
# tests/targets_baseline_policy.R using temporary synthetic targets stores.
message("Synthetic TabPFN analysis checks passed (no model inference).")
