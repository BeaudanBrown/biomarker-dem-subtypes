## Function for cross-validated biomarker subset identification

get_fold_stats <- function(data, outcome, reference, nfolds = 5, stuff) {
  # filter to outcome and reference

  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]
  data <- drop_na(data)

  set.seed(Sys.getenv("SEED"))

  Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)

  folds <- make_folds(
    nrow(data),
    fold_fun = folds_vfold,
    V = nfolds,
    strata_ids = Y
  )

  for (i in seq_along(folds)) {
    path <- stuff$path[[i]]
    subset_model <- stuff$subset_model[[i]]
    reference_model <- stuff$reference_model[[i]]

    reference_auc <- mean(auroc(data, outcome, reference)$AUC)

    train_auc <- get_ref_sub_AUC(
      data = data[folds[[i]]$training_set, ],
      reference_model = reference_model,
      subset_model = subset_model,
      best_biomarkers = path[[length(path)]]$removed_var,
      outcome = outcome
    )
    valid_auc <- get_ref_sub_AUC(
      data = data[folds[[i]]$validation_set, ],
      reference_model = reference_model,
      subset_model = subset_model,
      best_biomarkers = path[[length(path)]]$removed_var,
      outcome = outcome
    )
    print(mean(train_auc$reference_auc$AUC))
    print(mean(valid_auc$reference_auc$AUC))
    print("~~~~~~~~~~")
  }
}

get_biomarker_subset <- function(data, outcome, reference, nfolds = 5) {
  # filter to outcome and reference

  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]

  data <- drop_na(data)

  set.seed(Sys.getenv("SEED"))

  Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)

  folds <- make_folds(
    nrow(data),
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
  function(fold, data, outcome, reference) {
    train_data <- as_tibble(training(data))
    valid_data <- as_tibble(validation(data))

    out <-
      backwards_search(
        train_data,
        outcome = outcome,
        reference = reference
      )

    out$ref_sub_auc <- get_ref_sub_AUC(
      data = valid_data,
      reference_model = out$reference_model,
      subset_model = out$subset_model,
      best_biomarkers = out$path[[length(out$path)]]$removed_var,
      outcome = outcome
    )

    return(out)
  }

all_subset_data <- function(data, sex_strat = "", use_cdr = FALSE) {
  ad <- marker_subset(
    data,
    "Alzheimer's",
    c("Lewy bodies", "Frontotermporal"),
    sex_strat = sex_strat,
    use_cdr = use_cdr
  )
  ftd <- marker_subset(
    data,
    "Frontotermporal",
    c("Lewy bodies", "Alzheimer's"),
    sex_strat = sex_strat,
    use_cdr = use_cdr
  )
  lbd <- marker_subset(
    data,
    "Lewy bodies",
    c("Alzheimer's", "Frontotermporal"),
    sex_strat = sex_strat,
    use_cdr = use_cdr
  )
  list(ad = ad, ftd = ftd, lbd = lbd)
}

all_subset_plots <- function(data, extra_title = "", use_cdr = FALSE) {
  ad_title <- paste0(extra_title, " - n = ", data$ad$n)
  ad_path <- build_path(data$ad$path, data$ad$reference_auc)
  ad_plot <- plot_auc_steps(ad_path, "Alzheimer's", ad_title, use_cdr = use_cdr)

  ftd_title <- paste0(extra_title, " - n = ", data$ftd$n)
  ftd_path <- build_path(data$ftd$path, data$ftd$reference_auc)
  ftd_plot <- plot_auc_steps(
    ftd_path,
    "Frontotermporal",
    ftd_title,
    use_cdr = use_cdr
  )

  lbd_title <- paste0(extra_title, " - n = ", data$lbd$n)
  lbd_path <- build_path(data$lbd$path, data$lbd$reference_auc)
  lbd_plot <- plot_auc_steps(
    lbd_path,
    "Lewy bodies",
    lbd_title,
    use_cdr = use_cdr
  )
  list(ad = ad_plot, ftd = ftd_plot, lbd = lbd_plot)
}


run_single_subset <- function(
  data,
  outcome,
  reference,
  sex_strat = "",
  use_cdr = FALSE
) {
  set.seed(Sys.getenv("SEED"))
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
    select(-mean_ab42, -mean_ab40) |>
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

get_ref_sub_AUC <- function(
  data,
  reference_model,
  subset_model,
  best_biomarkers,
  outcome
) {
  all_markers <- names(select(data, -Diagnosis_combined))
  best_biomarkers <- c("age", "female", best_biomarkers)

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
      reference_model
    )

  # subset model AUC

  subset_auc <-
    get_auc(
      data,
      best_biomarkers,
      subset_model
    )

  return(list(
    reference_auc = reference_auc,
    subset_auc = subset_auc
  ))
}

## Helper functions for parallel subset analysis

prepare_subset_data <- function(data, sex_strat = "", use_cdr = FALSE) {
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

  data |>
    select(all_of(vars), starts_with("mean_")) |>
    select(-mean_ab42_ab40_ratio) |>
    drop_na()
}
