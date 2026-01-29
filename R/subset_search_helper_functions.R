## Function for cross-validated biomarker subset identification
run_single_subset <- function(
  data,
  outcome,
  reference,
  sex_strat = "",
  use_cdr = FALSE
) {
  if (sex_strat != "") {
    data <- data |>
      filter(female == ifelse(sex_strat == "female", 1, 0)) |>
      select(-female)
    vars <- c("Diagnosis_combined", "age")
  } else {
    vars <- c("Diagnosis_combined", "age", "female")
  }
  if (use_cdr) {
    vars <- append(vars, c("cdr"))
  }

  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ] |>
    select(all_of(vars), starts_with("mean_")) |>
    select(-mean_ab42_ab40_ratio) |>
    drop_na()

  backwards_search(
    data,
    outcome = outcome,
    reference = reference,
    threshold = 1
  )
}

## Backwards search for best minimal subset of biomarkers

backwards_search <- function(data, outcome, reference, threshold = 0.03) {
  # reference (full model) AUC

  auc_results <- auroc(data, outcome, reference)
  preds <- lapply(auc_results$aucs$results, function(a) a$preds)
  labels <- lapply(auc_results$aucs$results, function(a) a$labels)
  reference_auc <- ci.cvAUC(preds, labels)

  # start backwards search

  total_preds <- data |>
    select(starts_with("mean")) |>
    names() |>
    length()

  preds_removed <- 1
  smallest_auc_drop <- 0
  aucs <- list()
  result <- list()
  full_data <- data

  while (preds_removed <= total_preds && smallest_auc_drop < threshold) {
    indices <-
      which(
        !names(data) %in%
          c("age", "female", "Diagnosis_combined", "cdr")
      )

    get_submodel_auc <- function(i) {
      trim <- select(data, -names(data)[i])
      # auc <- mean(auroc(trim, outcome, reference)$AUC)
      auc_results <- auroc(trim, outcome, reference)
      preds <- lapply(auc_results$aucs$results, function(a) a$preds)
      labels <- lapply(auc_results$aucs$results, function(a) a$labels)
      auc <- ci.cvAUC(preds, labels)
      cil <- auc$ci[[1]]
      ciu <- auc$ci[[2]]
      return(tibble(
        removed_var = names(data)[i],
        auc = auc$cvAUC,
        cil = cil,
        ciu = ciu
      ))
    }

    auc_out <- map(indices, get_submodel_auc)

    aucs[[preds_removed]] <- bind_rows(auc_out)

    # variable to remove before next step
    best_auc_var <- filter(aucs[[preds_removed]], auc == max(auc))

    to_remove <- best_auc_var$removed_var

    data <- select(data, -all_of(to_remove))

    result[[preds_removed]] <- list(vars = colnames(data), auc = best_auc_var)

    smallest_auc_drop <- reference_auc$cvAUC - best_auc_var$auc

    preds_removed <- preds_removed + 1
  }

  # superlearner fit for full model

  X <- select(full_data, -Diagnosis_combined)
  Y <- ifelse(full_data$Diagnosis_combined == outcome, 1, 0)

  reference_model <-
    SuperLearner(
      Y = Y,
      X = X,
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = 10, stratifyCV = TRUE)
    )

  # superlearner fit for best subset model

  best_biomarkers <- aucs[[length(aucs)]]$removed_var

  covars <- data |>
    select(-starts_with("mean"), -Diagnosis_combined) |>
    names()
  X <- select(X, all_of(covars), all_of(best_biomarkers))

  subset_model <-
    SuperLearner(
      Y = Y,
      X = X,
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = 10, stratifyCV = TRUE)
    )

  # return full model superlearner, subset model superlearner, and best subset
  return(list(
    path = aucs,
    reference_auc = reference_auc,
    reference_model = reference_model,
    subset_model = subset_model,
    n = nrow(data[data$Diagnosis_combined == outcome, ])
  ))
}
