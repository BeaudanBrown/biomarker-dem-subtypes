### Function for cross-validated biomarker subset identification

get_biomarker_subset <- function(data, outcome, reference, nfolds = 5) {
  # filter to outcome and reference

  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]

  data <- drop_na(data)

  data$Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)

  data <- data |> select(-Diagnosis_combined)

  set.seed(1234)

  folds <- make_folds(nrow(data),
    fold_fun = folds_vfold,
    V = nfolds,
    strata_ids = data$Y
  )

  # plan(multicore, workers = as.numeric(Sys.getenv("N_CORES")))

  out <-
    cross_validate(
      folds = folds,
      cv_fun = cv_biomarker_subset,
      data = data,
      nfolds = nfolds
    )

  return(out)
}

cv_biomarker_subset <- function(
    fold,
    data,
    outcome = outcome,
    reference = reference) {
  train_data <- as_tibble(training(data))
  valid_data <- as_tibble(validation(data))

  out <- backwards_search(train_data, outcome = outcome, reference = reference)

  aucs_out <- get_AUC(valid_data, out)

  return(aucs_out)
}

## Backwards search for best minimal subset of biomarkers

backwards_search <- function(data, outcome, reference) {
  X <- select(data, -Y)
  Y <- data$Y

  reference_model <-
    SuperLearner(
      Y = Y,
      X = X,
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = 10, stratifyCV = TRUE)
    )

  # AUC

  reference_auc <- mean(auroc(data, outcome, reference)$AUC)

  # start backwards search

  total_preds <- length(
    which(!names(data) %in% c("age_combined", "female", "Diagnosis_combined"))
  )

  preds_removed <- 1
  smallest_auc_drop <- 0
  aucs <- list()

  while (preds_removed < total_preds & smallest_auc_drop < 0.03) {
    indices <-
      which(!names(data) %in% c("age_combined", "female", "Diagnosis_combined"))

    get_submodel_auc <- function(i) {
      trim <- select(data, -names(data)[i])
      auc <- mean(auroc(trim, outcome, reference)$AUC)
      return(tibble(removed_var = names(data)[i], auc = auc))
    }

    auc_out <- map(indices, get_submodel_auc)

    aucs[[preds_removed]] <- bind_rows(auc_out)

    print(aucs[[preds_removed]])

    # variable to remove before next step

    to_remove <- filter(aucs[[preds_removed]], auc == max(auc))$removed_var

    data <- select(data, -all_of(to_remove))

    smallest_auc_drop <-
      reference_auc - filter(aucs[[preds_removed]], auc == max(auc))$auc

    preds_removed <- preds_removed + 1
  }

  # Fit superlearner for best subset model

  best_biomarkers <- aucs[[length(aucs)]]$removed_var

  X <- select(data, age_combined, female, all_of(best_biomarkers))

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
    best_biomarkers = best_biomarkers,
    reference_model = reference_model,
    subset_model = subset_model
  ))
}

get_AUC <- function(data, out) {

  get_auc <- function(data, out, model) {
    X <- select(data, -Y)
    Y <- data$Y
    if (model == subset_model) {
      X <- select(X, age_combined, female, all_of(out$best_biomarkers))
    }
    preds <- as.numeric(predict(out$model, newdata = X)$pred)
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
      out,
      reference_model)

  # subset model AUC

  subset_auc <-
    get_auc(
      data,
      out,
      subset_model)

  return(list(
    reference_auc = reference_auc,
    subset_auc = subset_auc,
    best_biomarkers = out$best_biomarkers
  ))
}
