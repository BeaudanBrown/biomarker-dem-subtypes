### Functions for biomarker descriptions

source("learners.R")

### AUROC function

auroc <- function(data, outcome, reference, nfolds = 10) {
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

  plan(multicore, workers = as.numeric(Sys.getenv("N_CORES")))

  out <-
    cross_validate(
      folds = folds,
      cv_fun = cv_auroc,
      data = data,
      nfolds = nfolds
    )

  return(out)
}


cv_auroc <- function(fold, data, nfolds) {
  train_data <- as_tibble(training(data))
  valid_data <- as_tibble(validation(data))

  out <- predict_status(train_data, nfolds)

  roc_out <- get_roc(valid_data, out)

  return(roc_out)
}

predict_status <- function(data, nfolds) {
    train_Y <- data$Y

    train_X <- data |> select(-Y)

    # fit superlearner

    fit <- SuperLearner(
      Y = train_Y,
      X = train_X,
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = nfolds, stratifyCV = TRUE)
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

plot_roc_combined <- function(roc_list, title_text, label_map, nfolds = 10) {

  # A helper function to extract the summary results from a single ROC
  extract_results <- function(roc_data) {
    AUC <- round(mean(roc_data$AUC), 2)
    data <-
      tibble(
        TPR = roc_data$TPR,
        FPR = roc_data$FPR
      ) |>
      mutate(
        threshold = rep(1:(length(roc_data$TPR) / nfolds), nfolds)
      ) |>
      group_by(threshold) |>
      summarise(
        TPR = mean(TPR),
        FPR = mean(FPR),
        .groups = "drop"
      )
    data$AUC <- AUC
    return(data)
  }

  # Process each ROC element in the list
  combined_data <- map(roc_list, extract_results)
  combined_data <- bind_rows(combined_data, .id = "comparison")

  # For each group calculate the average AUC and build a new label
  # The new label is created by taking the sprintf template from label_map
  groups <- unique(combined_data[["comparison"]])
  new_levels <- sapply(groups, function(g) {
    auc_val <- round(mean(combined_data[combined_data[["comparison"]] == g, "AUC"][[1]]), 2)
    sprintf(label_map[[g]], auc_val)
  }, simplify = TRUE, USE.NAMES = FALSE)

  # Convert the grouping variable into a factor with the new labels
  combined_data[["comparison"]] <- factor(combined_data[["comparison"]], levels = groups)
  levels(combined_data[["comparison"]]) <- new_levels

  # Create the ROC plot
  plot <- ggplot(combined_data, aes_string(x = "FPR", y = "TPR", colour = "comparison")) +
    geom_path(size = 1) +
    geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "black", alpha = 0.5) +
    labs(
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      title = title_text
    ) +
    bayesplot::theme_default() +
    scale_color_colorblind() +
    theme(legend.position = "bottom")

  # Save the plot and return it
  ggsave(filename = paste0("plots/", title_text, ".png"),
                  plot = plot,
                  width = 13, height = 13, bg = "white")

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
  # Define output file name (be sure to clean the outcome_subtype if needed)
  outfile <- paste0(
    "vimp_",
    sub("'", "", outcome_subtype),
    ifelse(nzchar(stratification), paste0("_", stratification), ""),
    ".rds"
  )
  # Check if already computed
  if (file.exists(outfile)) {
    return(read_rds(outfile))
  }
  # Parallelize over columns in X; note that each parallel call will then run cv_vim
  library(future.apply)
  vimp_list <- future_lapply(seq_len(ncol(X)), function(i) {
    # Optionally: Instead of print, you could write to a log file to avoid console jumbles.
    message("Processing predictor ", i)
    cv_vim(
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
  }, future.seed = 1234)
  # Merge results from each predictor (assuming merge_vim is associative)
  out <- Reduce(function(x, y) merge_vim(x, y), vimp_list)
  # Save the results for the current outcome subtype
  write_rds(out, outfile)
  return(out)
}

vimp_function <- function(data, outcome_subtype) {
  # create dataset
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
