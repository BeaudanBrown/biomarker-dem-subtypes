### Functions for biomarker descriptions

auroc <- function(data, outcome, reference, nfolds = 10) {
  # filter to outcome and reference
  set.seed(Sys.getenv("SEED"))

  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]

  data <- drop_na(data)

  data$Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)

  data <- data |> select(-Diagnosis_combined)

  # standardise continuous predictors and outcomes
  contvars <- names(data)[sapply(data, function(col) length(unique(col)) > 5)]
  std <- function(x) {
    (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
  }
  data <- mutate(data, across(all_of(contvars), std))

  # CV
  folds <- make_folds(
    nrow(data),
    fold_fun = folds_vfold,
    V = nfolds,
    strata_ids = data$Y
  )

  out <-
    cross_validate(
      folds = folds,
      cv_fun = cv_auroc,
      data = data,
      nfolds = nfolds
    )

  return(list(aucs = out))
}


cv_auroc <- function(fold, data, nfolds) {
  train_data <- as_tibble(training(data))
  valid_data <- as_tibble(validation(data))

  out <- predict_status(train_data, nfolds)

  roc_out <- get_roc(valid_data, out)

  return(list(results = roc_out))
}

predict_status <- function(data, nfolds) {
  train_Y <- data$Y

  train_X <- data |> select(-Y)

  # fit superlearner

  if (ncol(train_X) == 1) {
    fit <- SuperLearner(
      Y = train_Y,
      X = train_X,
      family = binomial(),
      method = "method.AUC",
      SL.library = univariate_library,
      cvControl = list(V = 20, stratifyCV = TRUE)
    )
  } else {
    fit <- SuperLearner(
      Y = train_Y,
      X = train_X,
      family = binomial(),
      method = "method.AUC",
      SL.library = SL.library,
      cvControl = list(V = 20, stratifyCV = TRUE)
    )
  }

  print(fit$coef)

  pred_vars <- names(train_X)

  return(list(fit = fit, pred_vars = pred_vars))
}

get_roc <- function(data, out) {
  valid_X <- data |> select(-Y)
  valid_Y <- data$Y
  preds <- as.numeric(predict(out$fit, newdata = valid_X)$pred)
  roc_data <- tibble(true_labels = valid_Y, predicted_probs = preds)
  # Get unique thresholds (predicted probabilities)
  thresholds <- seq(0, 1, length.out = 100)

  # Calculate TPR and FPR for each threshold
  roc_data <- data.frame(
    t(sapply(thresholds, calculate_tpr_fpr, data = roc_data))
  )
  colnames(roc_data) <- c("TPR", "FPR")

  roc_data$AUC <- as.numeric(roc(valid_Y, preds)$auc)

  return(list(
    labels = valid_Y,
    preds = preds,
    TPR = roc_data$TPR,
    FPR = roc_data$FPR,
    AUC = roc_data$AUC
  ))
}

calculate_tpr_fpr <- function(threshold, data) {
  predictions <- ifelse(data$predicted_probs >= threshold, 1, 0)
  TP <- sum(predictions == 1 & data$true_labels == 1)
  TN <- sum(predictions == 0 & data$true_labels == 0)
  FP <- sum(predictions == 1 & data$true_labels == 0)
  FN <- sum(predictions == 0 & data$true_labels == 1)

  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)

  return(c(TPR, FPR))
}


plot_roc <- function(roc_data, title_text, nfolds = 5) {
  AUC <- mean(map_vec(roc_data$aucs$results, function(a) mean(a$AUC)))
  TPR <- Reduce(`+`, lapply(roc_data$aucs$results, function(a) a$TPR)) / nfolds
  FPR <- Reduce(`+`, lapply(roc_data$aucs$results, function(a) a$FPR)) / nfolds

  plot <- tibble(
    TPR = TPR,
    FPR = FPR
  ) |>
    mutate(comparison = paste0("LBD vs FTD (AUC: ", round(AUC, 2), ")")) |>
    ggplot(aes(x = FPR, y = TPR, colour = comparison)) +
    geom_path(size = 1) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "black",
      alpha = 0.5
    ) +
    labs(
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      title = title_text,
      colour = NULL
    ) +
    bayesplot::theme_default() +
    theme(legend.position = "bottom")

  ggsave(
    paste0("plots/", title_text, ".png"),
    plot = plot,
    width = 10,
    height = 10,
    bg = "white"
  )

  return(plot)
}


plot_roc_combined <- function(roc_list, title_text, label_map, nfolds = 5) {
  # A helper function to extract the summary results from a single ROC
  extract_results <- function(roc_data) {
    AUC <- mean(map_vec(roc_data$aucs$results, function(a) mean(a$AUC)))
    TPR <- Reduce(`+`, lapply(roc_data$aucs$results, function(a) a$TPR)) /
      nfolds
    FPR <- Reduce(`+`, lapply(roc_data$aucs$results, function(a) a$FPR)) /
      nfolds

    data <- tibble(
      TPR = TPR,
      FPR = FPR,
      AUC = AUC
    )

    return(data)
  }

  # Process each ROC element in the list
  combined_data <- map(roc_list, extract_results)
  combined_data <- bind_rows(combined_data, .id = "comparison")

  # For each group calculate the average AUC and build a new label
  # The new label is created by taking the sprintf template from label_map
  groups <- unique(combined_data[["comparison"]])
  new_levels <- sapply(
    groups,
    function(g) {
      auc_val <- round(
        mean(combined_data[combined_data[["comparison"]] == g, "AUC"][[1]]),
        2
      )
      sprintf(label_map[[g]], auc_val)
    },
    simplify = TRUE,
    USE.NAMES = FALSE
  )

  # Convert the grouping variable into a factor with the new labels
  combined_data[["comparison"]] <- factor(
    combined_data[["comparison"]],
    levels = groups
  )
  levels(combined_data[["comparison"]]) <- new_levels

  # Create the ROC plot
  plot <- ggplot(
    combined_data,
    aes_string(x = "FPR", y = "TPR", colour = "comparison")
  ) +
    geom_path(size = 1) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "black",
      alpha = 0.5
    ) +
    labs(
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      title = title_text,
      colour = NULL
    ) +
    bayesplot::theme_default() +
    scale_color_colorblind() +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2))

  # Save the plot and return it
  ggsave(
    filename = paste0("plots/", title_text, ".png"),
    plot = plot,
    width = 13,
    height = 13,
    bg = "white"
  )

  return(plot)
}

#### VIMP ####

vimp_function_par <- function(data, outcome_subtype, stratification = "") {
  # Remove rows with NA
  data <- drop_na(data)
  # Create response variable
  Y <- ifelse(data$Diagnosis_combined == outcome_subtype, 1, 0)
  # Remove the outcome column from predictors
  X <- select(data, -Diagnosis_combined)
  # Parallelize over columns in X; note that each parallel call will then run cv_vim
  vimp_list <- future_lapply(
    seq_len(ncol(X)),
    function(i) {
      # Optionally: Instead of print, you could write to a log file to avoid console jumbles.
      print("Processing predictor ", i)
      message("Processing predictor ", i)
      cv_vim(
        Y = Y,
        X = X,
        indx = i,
        type = "auc",
        V = 5,
        run_regression = TRUE,
        SL.library = SL.library,
        sample_splitting = FALSE,
        stratified = TRUE,
        cvControl = list(V = 20, stratifyCV = TRUE),
        family = binomial()
      )
    },
    future.seed = 1234
  )
  # Merge results from each predictor (assuming merge_vim is associative)
  out <- Reduce(function(x, y) merge_vim(x, y), vimp_list)
  return(out)
}

vimp_function <- function(data, outcome_subtype) {
  # create dataset
  data <- drop_na(data)

  Y <- ifelse(data$Diagnosis_combined == outcome_subtype, 1, 0)

  X <- select(data, -Diagnosis_combined)

  for (i in seq_len(ncol(X))) {
    print(i)
    vimp1 <- cv_vim(
      Y = Y,
      X = X,
      indx = i,
      type = "auc",
      V = 10,
      run_regression = TRUE,
      SL.library = SL.library,
      sample_splitting = FALSE,
      stratified = TRUE,
      cvControl = list(V = 20, stratifyCV = TRUE),
      family = binomial()
    )
    if (i == 1) {
      out <- vimp1
    } else {
      out <- merge_vim(out, vimp1)
    }
  }

  return(out)
}
