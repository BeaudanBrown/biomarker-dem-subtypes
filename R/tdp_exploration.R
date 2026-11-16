# TDP-43 exploratory report helpers. Existing analysis functions are unchanged.
# These functions run only when the user executes the corresponding targets.

tdp43_diagnoses <- function() {
  c("Control", "Alzheimer's", "Frontotermporal", "Lewy bodies")
}

tdp43_labels <- function() c("Controls", "AD", "FTD", "LBD/PD")

tdp43_palette <- function() {
  stats::setNames(c("#666666", "#0072B2", "#D55E00", "#009E73"), tdp43_labels())
}

tdp43_cohort_palette <- function() {
  c(All = "#222222", Overall = "#222222", tdp43_palette())
}

tdp43_markers <- function() {
  c(
    CD14 = "mean_elisa", NfL = "mean_nfl", `YKL-40` = "mean_ykl",
    GFAP = "mean_gfap", `Aβ42` = "mean_ab42", `Aβ40` = "mean_ab40",
    `TDP-43` = "mean_tdp", `p-tau181` = "mean_ptau181",
    `p-tau217` = "mean_ptau217", `Aβ42/40 (derived)` = "mean_ab42_ab40_ratio"
  )
}

tdp43_check_values <- function(x, label, positive = FALSE) {
  # Ordinary NA is allowed. NaN, infinity and invalid observed values are not.
  if (!is.numeric(x) || any(is.nan(x)) ||
      any(!is.na(x) & !is.finite(x)) ||
      (positive && any(x <= 0, na.rm = TRUE))) {
    stop(paste("Invalid observed values in", label), call. = FALSE)
  }
  invisible(NULL)
}

tdp43_prepare <- function(df) {
  required <- c(
    "Diagnosis_combined", "age", "female", "cdr", "centiloid",
    unname(tdp43_markers()), "LUMIPULSE_CSF_AB42", "LUMIPULSE_CSF_AB40"
  )
  if (!all(required %in% names(df))) {
    stop("Required columns are missing for the TDP-43 report.", call. = FALSE)
  }
  df <- tibble::as_tibble(df) |>
    dplyr::select(dplyr::all_of(required))
  if (any(is.na(df$Diagnosis_combined)) ||
      any(!df$Diagnosis_combined %in% tdp43_diagnoses())) {
    stop("Unexpected or missing diagnosis in the TDP-43 cohort.", call. = FALSE)
  }
  for (column in setdiff(required, "Diagnosis_combined")) {
    tdp43_check_values(
      df[[column]], column,
      positive = column %in% c("mean_tdp", "LUMIPULSE_CSF_AB42", "LUMIPULSE_CSF_AB40")
    )
  }
  if (any(!is.na(df$female) & !df$female %in% c(0, 1))) {
    stop("Unexpected sex encoding in the TDP-43 cohort.", call. = FALSE)
  }
  df <- df |>
    dplyr::mutate(
      Diagnosis = factor(Diagnosis_combined, levels = tdp43_diagnoses(), labels = tdp43_labels()),
      Sex = factor(female, levels = c(0, 1), labels = c("Males", "Females")),
      csf_ab42_ab40 = LUMIPULSE_CSF_AB42 / LUMIPULSE_CSF_AB40
    )
  tdp43_check_values(df$csf_ab42_ab40, "CSF Aβ42/40", positive = TRUE)
  df
}

tdp43_scatter <- function(df, predictor, label, key) {
  points <- df |>
    dplyr::filter(!is.na(mean_tdp), !is.na(.data[[predictor]])) |>
    dplyr::mutate(x = .data[[predictor]]) |>
    dplyr::select(Diagnosis, x, mean_tdp)
  counts <- points |>
    dplyr::count(Diagnosis, .drop = FALSE, name = "n") |>
    dplyr::mutate(Diagnosis = as.character(Diagnosis))
  counts <- dplyr::bind_rows(tibble::tibble(Diagnosis = "Total", n = nrow(points)), counts)
  caption <- paste(paste0(counts$Diagnosis, ": n=", counts$n), collapse = "; ")

  # Independent, unadjusted linear regressions on log(TDP-43).
  # No custom sparse-group checks: fitting/prediction errors surface naturally.
  all_trends <- purrr::map_dfr(c("All", tdp43_labels()), function(diagnosis) {
    group <- if (diagnosis == "All") points else
      dplyr::filter(points, as.character(Diagnosis) == diagnosis)
    fit <- stats::lm(log(mean_tdp) ~ x, data = group)
    grid <- data.frame(x = seq(min(group$x), max(group$x), length.out = 100))
    predicted <- stats::predict(fit, newdata = grid, interval = "confidence", level = 0.95)
    tibble::tibble(
      Diagnosis = factor(diagnosis, levels = c("All", tdp43_labels())),
      x = grid$x, fit = exp(predicted[, "fit"]),
      lower = exp(predicted[, "lwr"]), upper = exp(predicted[, "upr"])
    )
  })

  trends <- dplyr::filter(all_trends, as.character(Diagnosis) != "All")
  pooled_trends <- dplyr::filter(all_trends, as.character(Diagnosis) == "All")

  base <- ggplot2::ggplot(points, ggplot2::aes(x = x, y = mean_tdp, colour = Diagnosis)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.7) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = tdp43_palette(), drop = FALSE) +
    ggplot2::labs(x = label, y = "Plasma TDP-43 (concentration; log scale)", caption = caption) +
    cowplot::theme_cowplot(font_size = 12) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0))
  with_trends <- base +
    ggplot2::geom_ribbon(
      data = trends,
      ggplot2::aes(x = x, ymin = lower, ymax = upper, fill = Diagnosis),
      inherit.aes = FALSE, alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = trends, ggplot2::aes(x = x, y = fit, colour = Diagnosis),
      inherit.aes = FALSE, linewidth = 0.8
    ) +
    ggplot2::scale_fill_manual(values = tdp43_palette(), drop = FALSE, guide = "none")
  # Build the pooled view without a diagnosis colour scale: removing only the
  # colour mapping from `base` leaves an unused manual scale and emits warnings.
  pooled_base <- ggplot2::ggplot(points, ggplot2::aes(x = x, y = mean_tdp)) +
    ggplot2::geom_point(colour = unname(tdp43_cohort_palette()["All"]), alpha = 0.65, size = 1.7) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = label, y = "Plasma TDP-43 (concentration; log scale)",
      caption = paste0("All (pooled): n=", nrow(points))
    ) +
    cowplot::theme_cowplot(font_size = 12) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
  pooled_with_trends <- pooled_base +
    ggplot2::geom_ribbon(
      data = pooled_trends,
      ggplot2::aes(x = x, ymin = lower, ymax = upper),
      inherit.aes = FALSE, fill = unname(tdp43_cohort_palette()["All"]), alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = pooled_trends, ggplot2::aes(x = x, y = fit),
      inherit.aes = FALSE, colour = unname(tdp43_cohort_palette()["All"]), linewidth = 0.8
    )
  list(
    plots = stats::setNames(
      list(base, with_trends, pooled_base, pooled_with_trends),
      c(key, paste0(key, "_trends"), paste0(key, "_all"), paste0(key, "_all_trends"))
    ),
    tables = stats::setNames(list(counts), paste0(key, "_counts")),
    interactive = list(trends = all_trends, label = label)
  )
}

tdp43_box_fill_palette <- function(variable, groups) {
  if (variable == "Sex") return(c(Males = "#56B4E9", Females = "#CC79A7"))
  if (variable == "cdr") {
    return(stats::setNames(
      grDevices::colorRampPalette(c("#DEEBF7", "#08519C"))(length(groups)), groups
    ))
  }
  tdp43_cohort_palette()
}

tdp43_comparison <- function(df, variable, label, key, by_diagnosis = TRUE) {
  observed <- df |>
    dplyr::filter(!is.na(mean_tdp), !is.na(.data[[variable]]))
  if (variable == "cdr") {
    observed <- dplyr::mutate(observed, Group = factor(cdr, levels = sort(unique(cdr))))
  } else {
    observed <- dplyr::mutate(observed, Group = .data[[variable]])
  }
  # Define fills once over all observed levels so panel selection never
  # changes the colour assigned to a sex or CDR category.
  fill_palette <- tdp43_box_fill_palette(variable, levels(observed$Group))
  fill_alpha <- if (variable %in% c("Sex", "cdr")) 0.65 else 0.25
  panelled <- dplyr::mutate(observed, Panel = "Overall")
  if (by_diagnosis) {
    panelled <- dplyr::bind_rows(
      panelled,
      dplyr::mutate(observed, Panel = as.character(Diagnosis))
    )
  }
  panelled <- panelled |>
    dplyr::mutate(Panel = factor(Panel, levels = c("Overall", tdp43_labels()))) |>
    dplyr::select(Panel, Group, mean_tdp) |>
    dplyr::mutate(ColourGroup = if (by_diagnosis) as.character(Panel) else as.character(Group))
  summary <- panelled |>
    dplyr::group_by(Panel, Group) |>
    dplyr::summarise(
      n = dplyr::n(),
      median = stats::median(mean_tdp),
      Q1 = stats::quantile(mean_tdp, 0.25),
      Q3 = stats::quantile(mean_tdp, 0.75),
      .groups = "drop"
    )
  plot <- ggplot2::ggplot(panelled, ggplot2::aes(x = Group, y = mean_tdp, colour = ColourGroup)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = Group), outlier.shape = NA, alpha = fill_alpha) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 42),
      alpha = 0.45, size = 1.3
    ) +
    ggplot2::scale_colour_manual(values = tdp43_cohort_palette(), guide = "none") +
    ggplot2::scale_fill_manual(values = fill_palette, guide = "none", drop = FALSE) +
    # Calculate box statistics in concentration units, then transform the view.
    ggplot2::coord_trans(y = "log10") +
    ggplot2::labs(x = label, y = "Plasma TDP-43 (concentration; log scale)") +
    cowplot::theme_cowplot(font_size = 12)
  if (by_diagnosis) {
    plot <- plot + ggplot2::facet_wrap(~Panel, ncol = 3, drop = FALSE)
  }
  list(
    plots = stats::setNames(list(plot), key),
    tables = stats::setNames(list(summary), key),
    interactive = list(fill_palette = fill_palette, fill_alpha = fill_alpha)
  )
}

tdp43_correlations <- function(df) {
  markers <- tdp43_markers()
  panels <- c(list(Overall = df), stats::setNames(
    lapply(tdp43_labels(), function(diagnosis) {
      dplyr::filter(df, as.character(Diagnosis) == diagnosis)
    }), tdp43_labels()
  ))
  keys <- c("overall", "controls", "ad", "ftd", "lbd_pd")
  plots <- list()
  tables <- list()
  for (i in seq_along(panels)) {
    panel <- panels[[i]]
    correlations <- matrix(NA_real_, length(markers), length(markers), dimnames = list(names(markers), names(markers)))
    counts <- matrix(0L, length(markers), length(markers), dimnames = dimnames(correlations))
    for (j in seq_along(markers)) {
      for (k in seq_len(j)) {
        x <- panel[[markers[j]]]
        y <- panel[[markers[k]]]
        complete <- !is.na(x) & !is.na(y)
        if (sum(complete) < 2L || length(unique(x[complete])) < 2L || length(unique(y[complete])) < 2L) {
          stop(paste("Undefined Spearman correlation:", names(panels)[i], names(markers)[j], names(markers)[k]), call. = FALSE)
        }
        correlations[j, k] <- correlations[k, j] <- stats::cor(x[complete], y[complete], method = "spearman")
        counts[j, k] <- counts[k, j] <- sum(complete)
      }
    }
    long <- as.data.frame(as.table(correlations), responseName = "rho")
    plot <- ggplot2::ggplot(long, ggplot2::aes(x = Var1, y = Var2, fill = rho)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 3) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "Spearman ρ") +
      ggplot2::coord_equal() +
      ggplot2::labs(x = NULL, y = NULL, title = names(panels)[i]) +
      cowplot::theme_cowplot(font_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    key <- paste0("correlations_", keys[i])
    plots[[key]] <- plot
    tables[[key]] <- tibble::rownames_to_column(as.data.frame(correlations), "Biomarker")
    tables[[paste0(key, "_n")]] <- tibble::rownames_to_column(as.data.frame(counts), "Biomarker")
  }
  list(plots = plots, tables = tables)
}

tdp43_export <- function(outputs) {
  directory <- file.path("output", "tdp43")
  dir.create(file.path(directory, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(directory, "tables"), recursive = TRUE, showWarnings = FALSE)
  paths <- character()
  for (output in outputs) {
    for (key in names(output$plots)) {
      path <- file.path(directory, "figures", paste0(key, ".png"))
      ggplot2::ggsave(path, output$plots[[key]], width = 11, height = 8, dpi = 300, bg = "white")
      paths <- c(paths, path)
    }
    for (key in names(output$tables)) {
      path <- file.path(directory, "tables", paste0(key, ".csv"))
      data.table::fwrite(output$tables[[key]], path)
      paths <- c(paths, path)
    }
  }
  paths
}

tdp43_publish_report <- function(rendered_files) {
  # Explicit dependency on the tar_quarto file target; publish only the HTML.
  source <- rendered_files[basename(rendered_files) == "tdp43_report.html"]
  if (length(source) != 1L) stop("Expected one rendered TDP-43 HTML report.")
  destination <- file.path("output", "tdp43", "tdp43_report.html")
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, destination, overwrite = TRUE)) stop("Could not publish TDP-43 report.")
  destination
}
