source("biomarker_subset_search.R")

build_path_forward <- function(full_path, ref_auc) {
  all_vars <- c(
    "mean_elisa",
    "mean_nfl",
    "mean_ykl",
    "mean_gfap",
    "mean_ab40",
    "mean_ab42",
    "mean_tdp",
    "mean_ptau181",
    "mean_ptau217"
  )

  # Create a named vector for variable mapping
  var_mapping <- c(
    "mean_elisa" = "CD14",
    "mean_nfl" = "NfL",
    "mean_ykl" = "YKL-40",
    "mean_gfap" = "GFAP",
    "mean_ab40" = "AB40",
    "mean_ab42" = "AB42",
    "mean_tdp" = "TDP-43",
    "mean_ptau181" = "pTau-181",
    "mean_ptau217" = "pTau-217"
  )

  # Function to get included variables at each step
  get_included_vars <- function(step, added_vars) {
    if (step == 0) {
      return(character(0))  # Start with no variables
    } else {
      return(added_vars[1:step])
    }
  }

  best_path <- map(full_path, function(options) {
    return(options[which.max(options$auc), ])
  })

  # Extract added variables in order
  added_vars_in_order <- map_chr(best_path, ~.x$added_var)

  base_auc <- best_path[[1]]$auc

  result1 <- best_path %>%
    map2_dfr(seq_along(.), function(tibble, index) {
      tibble %>% mutate(step = index)
    }) |>
    arrange(step) |>
    mutate(
      raw_added_var = added_var,
      added_var = case_when(
        added_var %in% names(var_mapping) ~ var_mapping[added_var],
        TRUE ~ "Base"
      )
    ) |>
    mutate(
      # Create model_vars column with list of included variables at each step
      model_vars = map(step, function(s) {
        included_raw <- get_included_vars(s, added_vars_in_order)
        included_display <- var_mapping[included_raw]
        return(included_display)
      })
    ) |>
    # Rename columns to reflect forward search
    select(-raw_added_var)  # Remove the temporary column

  return(result1)
}

marker_subset_forward <- function(data, outcome, reference) {
  set.seed(1234)
  plan(multicore, workers = detectCores())
  data <- data[data$Diagnosis_combined %in% c(outcome, reference), ]
  out <-
    forward_search(
      data,
      outcome = outcome,
      reference = reference
    )
  return(out)
}

## Forwards search
forward_search <- function(data, outcome, reference) {
  # Start with only demographic variables (age_combined, female)
  base_vars <- c("age_combined", "female", "Diagnosis_combined")
  start_data <- data %>% select(all_of(base_vars))

  # All potential predictor variables (those starting with "mean")
  all_pred_vars <- data %>%
    select(starts_with("mean")) %>%
    names()

  total_preds <- length(all_pred_vars)

  # Initialize
  preds_added <- 1
  aucs <- list()
  result <- list()

  # Calculate AUC for the base model
  auc_results <- auroc(start_data, outcome, reference)
  preds <- lapply(auc_results$aucs$results, function(a) a$preds)
  labels <- lapply(auc_results$aucs$results, function(a) a$labels)
  base_auc <- ci.cvAUC(preds, labels)

  current_data <- start_data
  cil <- base_auc$ci[[1]]
  ciu <- base_auc$ci[[2]]
  aucs[[preds_added]] <- tibble(
    added_var = "",
    auc = base_auc$cvAUC,
    cil = cil,
    ciu = ciu
  )

  while (preds_added <= total_preds) {

    # Variables not yet in the model
    remaining_vars <- setdiff(all_pred_vars, names(current_data))

    get_addmodel_auc <- function(var) {
      # Add this variable to the current model
      add_data <- bind_cols(current_data, data %>% select(all_of(var)))

      # Calculate AUC
      auc_results <- auroc(add_data, outcome, reference)
      preds <- lapply(auc_results$aucs$results, function(a) a$preds)
      labels <- lapply(auc_results$aucs$results, function(a) a$labels)
      auc <- ci.cvAUC(preds, labels)
      cil <- auc$ci[[1]]
      ciu <- auc$ci[[2]]

      return(tibble(
        added_var = var,
        auc = auc$cvAUC,
        cil = cil,
        ciu = ciu
      ))
    }

    # Calculate AUC for each potential addition
    auc_out <- map(remaining_vars, get_addmodel_auc)

    # Record this step's results
    preds_added <- preds_added + 1
    aucs[[preds_added]] <- bind_rows(auc_out)

    # Find the best variable to add
    best_auc_var <- filter(aucs[[preds_added]], auc == max(auc))
    to_add <- best_auc_var$added_var

    # Add the best variable to the current model
    current_data <- bind_cols(
      current_data, 
      data %>% select(all_of(to_add))
    )

    result[[preds_added]] <- list(vars = colnames(current_data), auc = best_auc_var)
  }

  # Calculate reference model (full model) AUC
  auc_results <- auroc(data, outcome, reference)
  preds <- lapply(auc_results$aucs$results, function(a) a$preds)
  labels <- lapply(auc_results$aucs$results, function(a) a$labels)
  reference_auc <- ci.cvAUC(preds, labels)

  # SuperLearner fit for full model
  X <- select(data, -Diagnosis_combined)
  Y <- ifelse(data$Diagnosis_combined == outcome, 1, 0)

  reference_model <-
    SuperLearner(
      Y = Y,
      X = X,
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = 10, stratifyCV = TRUE)
    )

  # SuperLearner fit for model at the last step (all variables added)
  final_model <-
    SuperLearner(
      Y = Y,
      X = X,  # Same as full model
      family = binomial(),
      SL.library = SL.library,
      cvControl = list(V = 10, stratifyCV = TRUE)
    )

  # Return results in the same format as backward search
  return(list(
    path = aucs,
    reference_auc = reference_auc,
    reference_model = reference_model,
    subset_model = final_model
  ))
}
