### Functions for biomarker descriptions

source("learners.R")

### AUROC function

auroc <- function(data, outcome, reference) {
  # filter to outcome and reference

  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]

  data <- drop_na(data)

  data$Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)

  data <- data |> select(-Diagnosis_combined)

  set.seed(1234)

  folds <- make_folds(nrow(data),
    fold_fun = folds_vfold,
    V = 10,
    strata_ids = data$Y
  )

  plan(multicore, workers = as.numeric(Sys.getenv("N_CORES")))

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

  out <- predict_status(train_data)

  roc_out <- get_roc(valid_data, out)

  return(roc_out)
}

predict_status <-
  function(data) {
    train_Y <- data$Y

    train_X <- data |> select(-Y)

    # fit superlearner

    fit <- SuperLearner(
      Y = train_Y,
      X = train_X,
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = 10, stratifyCV = TRUE)
    )

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


plot_roc <- function(roc_data, title_text, nfolds = 10) {

  AUC <- round(mean(roc_data$AUC), 2)

  data <-
    tibble(
      TPR = roc_data$TPR,
      FPR = roc_data$FPR
    ) |>
    mutate(threshold = rep(1:(length(roc_data$TPR) / nfolds), nfolds)) |>
    group_by(threshold) |>
    summarise(
      TPR = mean(TPR),
      FPR = mean(FPR)
    )

  plot <- ggplot(data, aes(x = FPR, y = TPR)) +
  geom_path(size = 1) +
  geom_abline(
    slope = 1, intercept = 0,
    linetype = "dashed", color = "black", alpha = 0.5
  ) +
  annotate(
    geom = "text", x = 0.1, y = 0.95,
    label = paste("AUC =", AUC),
    size = 6, fontface = "bold"
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
    width = 10, height = 10, bg = "white"
  )

  return(plot)
}

plot_roc_combined <- function(roc_list, title_text, nfolds = 10) {

  extract_results <- function(roc_data) {
    AUC <- round(mean(roc_data$AUC), 2)

    data <-
      tibble(
        TPR = roc_data$TPR,
        FPR = roc_data$FPR
      ) |>
      mutate(threshold = rep(1:(length(roc_data$TPR) / nfolds), nfolds)) |>
      group_by(threshold) |>
      summarise(
        TPR = mean(TPR),
        FPR = mean(FPR)
      )
    data$AUC <- AUC
    return(data)
  }

  combined_data <- map(roc_list, extract_results)
  combined_data <- bind_rows(combined_data, .id = "Outcome")
  AD_AUC <- round(mean(combined_data[combined_data$Outcome == "AD", ]$AUC), 2)
  LB_AUC <- round(mean(combined_data[combined_data$Outcome == "LB", ]$AUC), 2)
  FT_AUC <- round(mean(combined_data[combined_data$Outcome == "FT", ]$AUC), 2)
  combined_data$Outcome <- as.factor(combined_data$Outcome)
  levels(combined_data$Outcome) <-
    c(
      paste0("AD vs other dementias (AUC: ", AD_AUC, ")"),
      paste0("FTD vs other dementias (AUC: ", FT_AUC, ")"),
      paste0("LBD vs other dementias (AUC: ", LB_AUC, ")")
    )


  plot <- ggplot(combined_data, aes(x = FPR, y = TPR, colour = Outcome)) +
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
  scale_color_colorblind() +
  theme(legend.position = "bottom")


  ggsave(paste0("plots/", title_text, ".png"),
    plot = plot,
    width = 13, height = 13, bg = "white"
  )

  return(plot)
}

#### VIMP ####

vimp_function <- function(data, outcome_subtype) {
  # create dataset

  data <-
    data |>
    select(Diagnosis_combined, age_combined, female, starts_with("mean_"))

  data <- drop_na(data)

  Y <- ifelse(data$Diagnosis_combined == outcome_subtype, 1, 0)

  X <- select(data, -Diagnosis_combined)

  if (file.exists(paste0("vimp_", sub("'", "", outcome_subtype), ".rds"))) {
    out <- read_rds(paste0("vimp_", sub("'", "", outcome_subtype), ".rds"))
  } else {

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
        cvControl = list(V = 10, stratifyCV = TRUE),
        family = binomial()
      )
      if (i == 1) {
        out <- vimp1
      } else {
        out <- merge_vim(out, vimp1)
      }
    }
  }

  write_rds(out, paste0("vimp_", sub("'", "", outcome_subtype), ".rds"))

  return(out)
}
