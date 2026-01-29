## Helper: Compute cross-validated AUC with confidence intervals
compute_cv_auc_with_ci <- function(data, outcome, reference_group) {
  auc_results <- auroc(data, outcome, reference_group)
  preds <- lapply(auc_results$aucs$results, function(a) a$preds)
  labels <- lapply(auc_results$aucs$results, function(a) a$labels)
  ci.cvAUC(preds, labels)
}

## Function for cross-validated biomarker subset identification
find_minimal_biomarker_subset <- function(
  data,
  outcome,
  reference_group,
  stratify_by_sex = NULL,
  adjust_for_cdr = FALSE
) {
  if (!is.null(stratify_by_sex) && stratify_by_sex != "") {
    sex_filter_value <- ifelse(stratify_by_sex == "female", 1, 0)
    data <- data |>
      filter(female == sex_filter_value) |>
      select(-female)
    covariates <- c("Diagnosis_combined", "age")
  } else {
    covariates <- c("Diagnosis_combined", "age", "female")
  }

  if (adjust_for_cdr) {
    covariates <- append(covariates, c("cdr"))
  }

  data <- data[data$Diagnosis_combined %in% c(outcome, reference_group), ] |>
    select(all_of(covariates), starts_with("mean_")) |>
    select(-mean_ab42_ab40_ratio) |>
    drop_na()

  backward_elimination_search(
    data,
    outcome = outcome,
    reference_group = reference_group,
    max_auc_drop = 1
  )
}

## Backwards search for best minimal subset of biomarkers
backward_elimination_search <- function(
  data,
  outcome,
  reference_group,
  max_auc_drop = 0.03
) {
  # Full model AUC
  full_model_auc <- compute_cv_auc_with_ci(data, outcome, reference_group)

  # Start backward search
  n_biomarkers <- data |>
    select(starts_with("mean")) |>
    names() |>
    length()

  n_removed <- 1
  current_auc_drop <- 0
  step_aucs <- list()
  original_data <- data

  COVARIATE_COLUMNS <- c("age", "female", "Diagnosis_combined", "cdr")

  while (n_removed <= n_biomarkers && current_auc_drop < max_auc_drop) {
    biomarker_indices <- which(!names(data) %in% COVARIATE_COLUMNS)

    compute_leave_one_out_auc <- function(i) {
      reduced_data <- select(data, -names(data)[i])

      auc <- compute_cv_auc_with_ci(reduced_data, outcome, reference_group)

      cil <- auc$ci[[1]]
      ciu <- auc$ci[[2]]

      return(tibble(
        removed_var = names(data)[i],
        auc = auc$cvAUC,
        cil = cil,
        ciu = ciu
      ))
    }

    iteration_results <- map(biomarker_indices, compute_leave_one_out_auc)

    step_aucs[[n_removed]] <- bind_rows(iteration_results)

    # Variable to remove before next step
    best_removal_candidate <- filter(step_aucs[[n_removed]], auc == max(auc))

    biomarker_to_remove <- best_removal_candidate$removed_var

    data <- select(data, -all_of(biomarker_to_remove))

    # Calculate drop relative to full model
    current_auc_drop <- full_model_auc$cvAUC - best_removal_candidate$auc

    n_removed <- n_removed + 1
  }

  # SuperLearner fit for full model
  X <- select(original_data, -Diagnosis_combined)
  Y <- ifelse(original_data$Diagnosis_combined == outcome, 1, 0)

  full_model_fit <-
    SuperLearner(
      Y = Y,
      X = X,
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = 10, stratifyCV = TRUE)
    )

  # SuperLearner fit for best subset model
  # The last element of step_aucs contains the biomarkers that were candidates for removal
  # in the final iteration - these are the ones we keep (selected subset)
  selected_biomarkers <- step_aucs[[length(step_aucs)]]$removed_var

  covariates <- original_data |>
    select(-starts_with("mean"), -Diagnosis_combined) |>
    names()

  X_subset <- select(X, all_of(covariates), all_of(selected_biomarkers))

  subset_model_fit <-
    SuperLearner(
      Y = Y,
      X = X_subset,
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = 10, stratifyCV = TRUE)
    )

  # Return results
  return(list(
    elimination_path = step_aucs,
    full_model_auc = full_model_auc,
    full_model_fit = full_model_fit,
    subset_model_fit = subset_model_fit,
    selected_biomarkers = selected_biomarkers,
    n = nrow(original_data[original_data$Diagnosis_combined == outcome, ])
  ))
}
