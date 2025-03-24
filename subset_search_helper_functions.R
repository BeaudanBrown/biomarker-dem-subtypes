### Function for cross-validated biomarker subset identification

get_fold_stats <- function(data, outcome, reference, nfolds = 5, stuff) {
  # filter to outcome and reference

  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]

  data <- drop_na(data)

  set.seed(1234)

  Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)

  folds <- make_folds(nrow(data),
    fold_fun = folds_vfold,
    V = nfolds,
    strata_ids = Y
  )

  for (i in 1:length(folds)) {
    path <- stuff$path[[i]]
    subset_model <- stuff$subset_model[[i]]
    reference_model <- stuff$reference_model[[i]]

    reference_auc <- mean(auroc(data, outcome, reference)$AUC)

    train_auc <- get_ref_sub_AUC(data = data[folds[[i]]$training_set, ],
      reference_model = reference_model,
      subset_model = subset_model,
      best_biomarkers = path[[length(path)]]$removed_var,
      outcome = outcome)
    valid_auc <- get_ref_sub_AUC(data = data[folds[[i]]$validation_set, ],
      reference_model = reference_model,
      subset_model = subset_model,
      best_biomarkers = path[[length(path)]]$removed_var,
      outcome = outcome)
    print(mean(train_auc$reference_auc$AUC))
    print(mean(valid_auc$reference_auc$AUC))
    print("~~~~~~~~~~")
  }
}

get_biomarker_subset <- function(data, outcome, reference, nfolds = 5) {
  # filter to outcome and reference

  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]

  data <- drop_na(data)

  set.seed(1234)

  Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)

  folds <- make_folds(nrow(data),
    fold_fun = folds_vfold,
    V = nfolds,
    strata_ids = Y
  )

  out <-
    cross_validate(
      folds = folds,
      cv_fun = cv_biomarker_subset,
      data = data,
      outcome = outcome,
      reference = reference
    )

  return(out)
}

cv_biomarker_subset <-
  function(fold,
           data,
           outcome,
           reference) {
    train_data <- as_tibble(training(data))
    valid_data <- as_tibble(validation(data))

    out <-
      backwards_search(
        train_data,
        outcome = outcome,
        reference = reference
      )

    out$ref_sub_auc <- get_ref_sub_AUC(data = valid_data,
                                       reference_model = out$reference_model,
                                       subset_model = out$subset_model,
                                       best_biomarkers = out$path[[length(out$path)]]$removed_var,
                                       outcome = outcome)

    return(out)
}

marker_subset_full <- function(data, outcome, reference) {
  set.seed(1234)
  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]
  out <-
    backwards_search(
      data,
      outcome = outcome,
      reference = reference,
      threshold = 1
    )
  return(out)
}

## Backwards search for best minimal subset of biomarkers

backwards_search <- function(data, outcome, reference, threshold = 0.03) {

  # reference (full model) AUC

  reference_auc <- mean(auroc(data, outcome, reference)$AUC)

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

  while (preds_removed < total_preds && smallest_auc_drop < threshold) {
    indices <-
      which(!names(data) %in% c("age_combined", "female", "Diagnosis_combined"))

    get_submodel_auc <- function(i) {
      trim <- select(data, -names(data)[i])
      auc <- mean(auroc(trim, outcome, reference)$AUC)
      return(tibble(removed_var = names(data)[i], auc = auc))
    }

    auc_out <- map(indices, get_submodel_auc)

    aucs[[preds_removed]] <- bind_rows(auc_out)

    # variable to remove before next step
    best_auc_var <- filter(aucs[[preds_removed]], auc == max(auc))

    to_remove <- best_auc_var$removed_var

    data <- select(data, -all_of(to_remove))

    result[[preds_removed]] <- list(vars = colnames(data), auc = best_auc_var)

    smallest_auc_drop <- reference_auc - best_auc_var$auc

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
    subset_model = subset_model
  ))
}

get_ref_sub_AUC <- function(data, reference_model, subset_model, best_biomarkers, outcome) {
  all_markers <- names(select(data, -Diagnosis_combined))
  best_biomarkers <- c("age_combined", "female", best_biomarkers)

  get_auc <- function(data, predictors, model) {
    X <- select(data, all_of(predictors))
    Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)
    data <- select(data, -Diagnosis_combined)
    preds <- as.numeric(predict(model, newdata = X)$pred)
    roc_data <- tibble(true_labels = Y, predicted_probs = preds)
    # Get unique thresholds (predicted probabilities)
    thresholds <- seq(0, 1, length.out = 100)

    # Calculate TPR and FPR for each threshold
    roc_data <- data.frame(
      t(sapply(thresholds, calculate_tpr_fpr, data = roc_data))
    )
    colnames(roc_data) <- c("TPR", "FPR")

    roc_data$AUC <- as.numeric(roc(Y, preds)$auc)

    return(list(
      TPR = roc_data$TPR,
      FPR = roc_data$FPR,
      AUC = roc_data$AUC
    ))
  }

  # reference model AUC

  reference_auc <-
    get_auc(
      data,
      all_markers,
      reference_model)

  # subset model AUC

  subset_auc <-
    get_auc(
      data,
      best_biomarkers,
      subset_model)

  return(list(
    reference_auc = reference_auc,
    subset_auc = subset_auc
  ))
}
