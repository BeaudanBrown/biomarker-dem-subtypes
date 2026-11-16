# Browser-only presentation helpers. They reuse precomputed models/summaries.
# Cowplot supplies visual styling; boxes use explicit per-group summaries.

tdp43_tag_traces <- function(traces, panel, kind) {
  lapply(traces, function(trace) {
    trace$meta <- list(panel = panel, kind = kind)
    trace$showlegend <- FALSE
    trace$visible <- kind != "trend" && panel != "All"
    trace
  })
}

tdp43_widget_finish <- function(widget, traces, kind, script, counts = NULL) {
  widget$x$data <- traces
  widget$x$layout$showlegend <- FALSE
  widget$x$layout$title <- NULL
  # Preserve the full-cohort ranges when toggling groups/layers.
  widget$x$layout$xaxis$autorange <- FALSE
  widget$x$layout$yaxis$autorange <- FALSE
  widget$x$layout$uirevision <- "tdp43"
  widget <- plotly::config(widget, displaylogo = FALSE, responsive = TRUE)
  htmlwidgets::onRender(
    widget,
    paste(readLines(script, warn = FALSE), collapse = "\n"),
    data = list(kind = kind, counts = counts)
  )
}

tdp43_scatter_widget <- function(output, script) {
  base <- output$plots[[1]] + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
  # Include pooled and subgroup confidence bands in the fixed axis ranges.
  range_plot <- output$plots[[2]] + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) +
    ggplot2::geom_blank(
      data = output$interactive$trends, ggplot2::aes(x = x, y = lower), inherit.aes = FALSE
    ) +
    ggplot2::geom_blank(
      data = output$interactive$trends, ggplot2::aes(x = x, y = upper), inherit.aes = FALSE
    )
  widget <- plotly::ggplotly(range_plot)
  traces <- list()
  for (diagnosis in c("All", tdp43_labels())) {
    points <- base
    if (diagnosis == "All") {
      points <- output$plots[[3]] + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
    } else {
      points$data <- dplyr::filter(base$data, as.character(Diagnosis) == diagnosis)
    }
    points$data$hover <- paste0(
      diagnosis, "<br>", output$interactive$label, ": ", signif(points$data$x, 5),
      "<br>TDP-43: ", signif(points$data$mean_tdp, 5)
    )
    points <- points + ggplot2::aes(text = hover)
    points_widget <- plotly::ggplotly(points, tooltip = "text")
    traces <- c(traces, tdp43_tag_traces(points_widget$x$data, diagnosis, "points"))
    trends <- dplyr::filter(output$interactive$trends, as.character(Diagnosis) == diagnosis)
    colour <- unname(tdp43_cohort_palette()[diagnosis])
    trends$hover <- paste0(
      diagnosis, "<br>", output$interactive$label, ": ", signif(trends$x, 5),
      "<br>Fitted TDP-43: ", signif(trends$fit, 5),
      "<br>95% CI: ", signif(trends$lower, 5), "–", signif(trends$upper, 5)
    )
    trend_plot <- ggplot2::ggplot(trends, ggplot2::aes(x = x, y = fit, text = hover, group = 1)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = colour, alpha = 0.15) +
      ggplot2::geom_line(colour = colour, linewidth = 0.8) +
      ggplot2::scale_y_log10() +
      cowplot::theme_cowplot(font_size = 12)
    trend_widget <- plotly::ggplotly(trend_plot, tooltip = "text")
    traces <- c(traces, tdp43_tag_traces(trend_widget$x$data, diagnosis, "trend"))
  }
  counts <- output$tables[[1]]
  counts <- counts[counts$Diagnosis != "Total", ]
  tdp43_widget_finish(
    widget, traces, "scatter", script,
    as.list(stats::setNames(counts$n, counts$Diagnosis))
  )
}

tdp43_box_statistics <- function(data) {
  # Compute each box independently from its observed concentrations. Whiskers
  # are the extreme observations within 1.5 IQR, not global plot limits.
  data |>
    dplyr::group_by(Group) |>
    dplyr::summarise(
      n = dplyr::n(),
      Q1 = stats::quantile(mean_tdp, 0.25),
      median = stats::median(mean_tdp),
      Q3 = stats::quantile(mean_tdp, 0.75),
      lowerfence = min(mean_tdp[mean_tdp >= Q1 - 1.5 * (Q3 - Q1)]),
      upperfence = max(mean_tdp[mean_tdp <= Q3 + 1.5 * (Q3 - Q1)]),
      .groups = "drop"
    )
}

tdp43_box_trace <- function(
  statistics, group_levels, panel,
  colour = unname(tdp43_cohort_palette()[panel]),
  fill_colour = colour, fill_alpha = 0.25
) {
  # Use Plotly's explicit-summary API. No synthetic y samples or automatic
  # ggplotly box conversion are allowed to redefine extrema or quartiles.
  array <- function(x) unname(as.list(x))
  list(
    type = "box",
    x = array(match(as.character(statistics$Group), group_levels)),
    q1 = array(statistics$Q1),
    median = array(statistics$median),
    q3 = array(statistics$Q3),
    lowerfence = array(statistics$lowerfence),
    upperfence = array(statistics$upperfence),
    boxpoints = FALSE,
    boxmean = FALSE,
    width = 0.75,
    fillcolor = {
      rgb <- grDevices::col2rgb(fill_colour)
      sprintf("rgba(%d,%d,%d,%g)", rgb[1], rgb[2], rgb[3], fill_alpha)
    },
    line = list(color = colour),
    name = panel,
    meta = list(panel = panel, kind = "box"),
    visible = identical(panel, "Overall"),
    showlegend = FALSE
  )
}

tdp43_box_widget <- function(output, script) {
  # ggplotly supplies the cowplot layout and jittered points only. Box geometry
  # and hover statistics come from explicit, original-unit summaries below.
  base <- output$plots[[1]] + ggplot2::facet_null() + ggplot2::coord_cartesian() +
    ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
  base <- base + ggplot2::scale_x_discrete(limits = levels(base$data$Group), drop = FALSE)
  log_range <- range(log10(base$data$mean_tdp))
  padding <- max(diff(log_range) * 0.05, 0.05)
  panels <- unique(as.character(base$data$Panel))
  traces <- list()
  widget <- NULL
  for (panel in panels) {
    panel_plot <- base
    panel_plot$data <- dplyr::filter(base$data, as.character(Panel) == panel)
    converted <- plotly::ggplotly(panel_plot, tooltip = c("x", "y"))
    if (is.null(widget)) widget <- converted
    statistics <- tdp43_box_statistics(panel_plot$data)
    for (i in seq_len(nrow(statistics))) {
      group <- as.character(statistics$Group[i])
      colour_group <- if (group %in% tdp43_labels()) group else panel
      traces[[length(traces) + 1L]] <- tdp43_box_trace(
        statistics[i, ], levels(base$data$Group), panel,
        colour = unname(tdp43_cohort_palette()[colour_group]),
        fill_colour = unname(output$interactive$fill_palette[group]),
        fill_alpha = output$interactive$fill_alpha
      )
    }
    for (trace in converted$x$data) {
      is_points <- identical(trace$type, "scatter") &&
        !is.null(trace$mode) && grepl("markers", trace$mode)
      if (!is_points) next
      trace$meta <- list(panel = panel, kind = "points")
      trace$visible <- identical(panel, "Overall")
      trace$showlegend <- FALSE
      traces[[length(traces) + 1L]] <- trace
    }
  }
  # Each category now has its own coloured box trace; retain category centres.
  widget$x$layout$boxmode <- "overlay"
  widget$x$layout$yaxis$type <- "log"
  widget$x$layout$yaxis$range <- c(log_range[1] - padding, log_range[2] + padding)
  widget$x$layout$yaxis$tickmode <- "auto"
  widget$x$layout$yaxis$tickvals <- NULL
  widget$x$layout$yaxis$ticktext <- NULL
  widget$x$layout$yaxis$tick0 <- NULL
  widget$x$layout$yaxis$dtick <- NULL
  widget$x$layout$yaxis$tickformat <- NULL
  widget$x$layout$yaxis$hoverformat <- ".5g"
  tdp43_widget_finish(widget, traces, "box", script)
}

tdp43_matrix_widget <- function(output, script) {
  panels <- c("Overall", tdp43_labels())
  traces <- list()
  widget <- NULL
  for (i in seq_along(output$plots)) {
    key <- names(output$plots)[i]
    panel_plot <- output$plots[[i]]
    counts <- as.matrix(output$tables[[paste0(key, "_n")]][, -1])
    panel_plot$data$hover <- paste0(
      panel_plot$data$Var1, " / ", panel_plot$data$Var2,
      "<br>Spearman ρ: ", sprintf("%.2f", panel_plot$data$rho),
      "<br>n: ", as.vector(counts)
    )
    panel_plot <- panel_plot + ggplot2::aes(text = hover) +
      ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL)
    converted <- plotly::ggplotly(panel_plot, tooltip = "text")
    if (is.null(widget)) widget <- converted
    for (trace in converted$x$data) {
      trace$meta <- list(panel = panels[i], kind = "matrix")
      trace$visible <- i == 1L
      trace$showlegend <- FALSE
      traces[[length(traces) + 1L]] <- trace
    }
  }
  tdp43_widget_finish(widget, traces, "matrix", script)
}

tdp43_panel_select <- function(panels) {
  htmltools::tags$label(
    "Cohort: ",
    htmltools::tags$select(
      `data-tdp-panel-select` = "",
      lapply(panels, function(panel) htmltools::tags$option(
        value = panel, style = paste0("color:", tdp43_cohort_palette()[panel]), panel
      ))
    )
  )
}

tdp43_scatter_view <- function(output, script) {
  controls <- htmltools::tags$fieldset(
    htmltools::tags$legend("Show cohorts"),
    htmltools::tags$label(
      style = paste0("color:", tdp43_cohort_palette()["All"]),
      htmltools::tags$input(type = "checkbox", `data-tdp-cohort` = "All"), " All (pooled)"
    ),
    lapply(tdp43_labels(), function(diagnosis) {
      htmltools::tags$label(
        style = paste0("color:", tdp43_palette()[diagnosis]),
        htmltools::tags$input(type = "checkbox", `data-tdp-cohort` = diagnosis, checked = "checked"),
        paste0(" ", diagnosis)
      )
    })
  )
  htmltools::tags$div(
    class = "tdp43-widget", controls,
    tdp43_scatter_widget(output, script),
    htmltools::tags$p(`data-tdp-counts` = "", `aria-live` = "polite")
  )
}

tdp43_summary_html <- function(summary) {
  display <- summary |>
    dplyr::mutate(`Median (Q1, Q3)` = sprintf("%.1f (%.1f, %.1f)", median, Q1, Q3)) |>
    dplyr::select(Group, n, `Median (Q1, Q3)`)
  htmltools::HTML(as.character(knitr::kable(display, format = "html", row.names = FALSE)))
}

tdp43_box_view <- function(output, script) {
  summary <- output$tables[[1]]
  panels <- unique(as.character(summary$Panel))
  htmltools::tags$div(
    class = "tdp43-widget",
    if (length(panels) > 1L) tdp43_panel_select(panels),
    htmltools::tags$label(
      htmltools::tags$input(type = "checkbox", `data-tdp-points` = "", checked = "checked"),
      " Show individual points"
    ),
    tdp43_box_widget(output, script),
    lapply(panels, function(panel) {
      htmltools::tags$div(
        `data-tdp-table-panel` = panel,
        style = if (panel == "Overall") "" else "display:none",
        tdp43_summary_html(dplyr::filter(summary, as.character(Panel) == panel))
      )
    })
  )
}

tdp43_matrix_view <- function(output, script) {
  panels <- c("Overall", tdp43_labels())
  htmltools::tags$div(
    class = "tdp43-widget", tdp43_panel_select(panels),
    tdp43_matrix_widget(output, script),
    lapply(seq_along(panels), function(i) {
      key <- names(output$plots)[i]
      htmltools::tags$details(
        `data-tdp-table-panel` = panels[i],
        style = if (i == 1L) "" else "display:none",
        htmltools::tags$summary("Pairwise sample sizes"),
        htmltools::HTML(as.character(knitr::kable(
          output$tables[[paste0(key, "_n")]], format = "html", row.names = FALSE
        )))
      )
    })
  )
}
