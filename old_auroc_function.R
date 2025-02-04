auroc <- function(data, outcome, reference) {
  # filter to outcome and reference

  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]
  # create dataset

  data <-
    data |>
    select(Diagnosis_combined, age_combined, female, starts_with("mean_"))

  data <- drop_na(data)

  data$Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)

  data <- data |> select(-Diagnosis_combined)

  folds <- make_folds(data, strata_ids = data$Y)

  out <-
    cross_validate(
      folds = folds,
      cv_fun = cv_auroc,
      data = data
    )

  return(out)
}

cv_auroc <- function(fold, data) {
  train_data <- as_tibble(training(data))
  valid_data <- as_tibble(validation(data))

  out_basic <-
    predict_status(train_data,
      include_CD14 = FALSE,
      include_YKL = FALSE,
      include_GFAP = FALSE
    )

  out_complete <-
    predict_status(train_data,
      include_CD14 = TRUE,
      include_YKL = TRUE,
      include_GFAP = TRUE
    )

  roc_out <-
    get_roc(
      valid_data,
      out_basic,
      out_complete
    )

  return(roc_out)
}

predict_status <-
  function(data, include_CD14, include_YKL, include_GFAP) {
    train_Y <- data$Y

    train_X <- data |> select(-Y)

    if (isFALSE(include_CD14)) {
      train_X <- select(train_X, -mean_elisa)
    }

    if (isFALSE(include_YKL)) {
      train_X <- select(train_X, -mean_ykl)
    }

    if (isFALSE(include_GFAP)) {
      train_X <- select(train_X, -mean_gfap)
    }

    # fit superlearner

    fit <- SuperLearner(
      Y = train_Y,
      X = train_X,
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = 10, stratifyCV = TRUE)
    )

    pred_vars <- names(train_X)

    return(list(fit = fit, pred_vars = pred_vars))
  }


get_roc <- function(data, out_basic, out_complete) {
  auroc_data <- function(data, out) {
    valid_X <- data |> select(all_of(out$pred_vars))
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
      TPR = roc_data$TPR,
      FPR = roc_data$FPR,
      AUC = roc_data$AUC
    ))
  }

  auroc_basic <- auroc_data(data, out_basic)
  auroc_complete <- auroc_data(data, out_complete)

  return(
    list(
      basic_TPR = auroc_basic$TPR,
      basic_FPR = auroc_basic$FPR,
      basic_AUC = auroc_basic$AUC,
      complete_TPR = auroc_complete$TPR,
      complete_FPR = auroc_complete$FPR,
      complete_AUC = auroc_complete$AUC
    )
  )
}


plot_roc <- function(roc_data, title_text) {
  basic_AUC <- round(mean(roc_data$basic_AUC), 3)
  basic_data <-
    tibble(
      TPR = roc_data$basic_TPR,
      FPR = roc_data$basic_FPR
    ) |>
    mutate(threshold = rep(1:(length(roc_data$basic_TPR)/10), 10)) |>
    group_by(threshold) |>
    summarise(
      TPR = mean(TPR),
      FPR = mean(FPR)
    ) |>
    mutate(Model =
      paste0("No CD14, YKL-40, or GFAP", " (AUC: ", basic_AUC, ")"), )

  complete_AUC <- round(mean(roc_data$complete_AUC), 3)
  complete_data <-
    tibble(
      TPR = roc_data$complete_TPR,
      FPR = roc_data$complete_FPR
    ) |>
    mutate(threshold = rep(1:(length(roc_data$basic_TPR)/10), 10)) |>
    group_by(threshold) |>
    summarise(
      TPR = mean(TPR),
      FPR = mean(FPR)
    ) |>
    mutate(
      Model =
        paste0("With CD14, YKL-40, and GFAP", " (AUC: ", complete_AUC, ")"),
    )

  plot_data <- bind_rows(basic_data, complete_data)

  plot <- ggplot(plot_data, aes(x = FPR, y = TPR, colour = Model)) +
    geom_path(size = 1) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "black", alpha = 0.5
    ) +
    labs(
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      title = title_text
    ) +
    bayesplot::theme_default() +
    theme(legend.position = "bottom")

  ggsave(paste0("plots/", title_text, ".png"),
    plot = plot,
    width = 10, height = 10
  )

  return(plot)
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
