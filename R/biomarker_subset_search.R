
plot_auc_steps <- function(data, outcome, extra_title = "", use_cdr = TRUE) {
  # Order by step
  data <- data %>% arrange(step)

  # Get the removal order from the data
  removed_vars_order <- data %>%
    filter(removed_var != "Full") %>%
    arrange(step) %>%
    pull(removed_var)

  # Create a fixed list of all variables, ordered by removal (reverse order)
  # Last removed at the top, first removed at the bottom, with covariates at the very bottom
  all_vars <- c(
    rev(removed_vars_order),  # Biomarkers in reverse removal order
    "Age",                    # Covariates at the bottom
    "Sex"
  )
  if (use_cdr) {
    all_vars <- c(all_vars, "CDR")
  }

  create_step_labels <- function(step) {
    if (step == 0) {
      # At step 0, all biomarkers are active
      active_vars <- all_vars
    } else {
      # Variables removed up to this step
      removed_up_to_step <- data %>%
        filter(step > 0, step <= !!step) %>%
        pull(removed_var)

      # Active variables
      active_vars <- setdiff(all_vars, removed_up_to_step)
    }

    # Join active variables with newlines
    paste(active_vars, collapse = "\n")
  }

  # Create labels for each step
  labels <- map_chr(data$step, create_step_labels)

  # Create the plot
  p <- ggplot(data, aes(x = step, y = auc)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_pointrange(aes(ymin = cil, ymax = ciu), color = "darkblue") +
    labs(
      title = paste0(outcome, " Marker Path", extra_title),
      x = "Markers Remaining",
      y = "AUC"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_x_continuous(
      breaks = data$step,
      labels = labels,
      minor_breaks = NULL
    ) +
    scale_y_continuous(limits = c(0.5, 1), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))

  return(p)
}

build_path <- function(full_path, ref_auc) {
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

  # Function to get remaining variables at each step
  get_remaining_vars <- function(step, removed_vars) {
    if (step == 0) {
      return(all_vars)
    } else if (step == length(all_vars)) {
      return(NA)
    } else {
      return(setdiff(all_vars, removed_vars[1:step]))
    }
  }

  best_path <- map(full_path, function(options) {
    return(options[which.max(options$auc), ])
  })

  # Extract removed variables in order
  removed_vars_in_order <- map_chr(best_path, ~.x$removed_var)

  result1 <- best_path %>%
    map2_dfr(seq_along(.), function(tibble, index) {
      tibble %>% mutate(step = index)
    }) |>
    add_row(
      removed_var = "Full",
      auc = ref_auc$cvAUC,
      cil = ref_auc$ci[[1]],
      ciu = ref_auc$ci[[2]],
      step = 0) |>
    arrange(step) |>
    mutate(
      raw_removed_var = removed_var,
      removed_var = case_when(
        removed_var %in% names(var_mapping) ~ var_mapping[removed_var],
        TRUE ~ "Full"
      )
    ) |>
    mutate(
      # Create model_vars column with list of remaining variables at each step
      model_vars = map(step, function(s) {
        remaining_raw <- get_remaining_vars(s, removed_vars_in_order)
        remaining_display <- var_mapping[remaining_raw]
        return(remaining_display)
      })
    ) |>
    select(-raw_removed_var)  # Remove the temporary column

  return(result1)
}
