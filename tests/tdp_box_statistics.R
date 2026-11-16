# USER-RUN ONLY: synthetic regression checks; no project data or R pipeline.
# From the repository root: source("tests/tdp_box_statistics.R")
source("R/tdp_exploration.R")
source("R/tdp_widgets.R")

fixture <- data.frame(
  Group = factor(rep(c("Controls", "AD"), each = 5), levels = c("Controls", "AD")),
  mean_tdp = c(1, 2, 3, 4, 50, 10, 20, 30, 40, 500)
)
statistics <- tdp43_box_statistics(fixture)
stopifnot(
  identical(as.character(statistics$Group), c("Controls", "AD")),
  all(statistics$Q1 == c(2, 20)),
  all(statistics$median == c(3, 30)),
  all(statistics$Q3 == c(4, 40)),
  all(statistics$lowerfence == c(1, 10)),
  all(statistics$upperfence == c(4, 40))
)
trace <- tdp43_box_trace(statistics, levels(fixture$Group), "Overall")
stopifnot(
  identical(trace$type, "box"),
  is.null(trace$y), # No reconstructed samples from automatic conversion.
  all(unlist(trace$x) == c(1, 2)),
  all(unlist(trace$lowerfence) == c(1, 10)),
  all(unlist(trace$upperfence) == c(4, 40)),
  identical(trace$boxpoints, FALSE)
)
# One group retains its original category position and explicit summary arrays.
ad <- tdp43_box_statistics(fixture[fixture$Group == "AD", ])
ad_trace <- tdp43_box_trace(ad, levels(fixture$Group), "AD")
stopifnot(
  identical(ad_trace$x, list(2L)),
  identical(ad_trace$median, list(30)),
  identical(ad_trace$upperfence, list(40)),
  identical(ad_trace$visible, FALSE)
)
# Direct diagnostic boxes and sex/CDR panels share the scatterplot palette.
for (cohort in tdp43_labels()) {
  colour <- unname(tdp43_cohort_palette()[cohort])
  coloured <- tdp43_box_trace(ad, levels(fixture$Group), cohort, colour)
  stopifnot(identical(coloured$line$color, unname(tdp43_palette()[cohort])))
}
stopifnot(tdp43_cohort_palette()["All"] == tdp43_cohort_palette()["Overall"])
sex_fills <- tdp43_box_fill_palette("Sex", c("Males", "Females"))
stopifnot(sex_fills["Males"] != sex_fills["Females"])
cdr_fills <- tdp43_box_fill_palette("cdr", c("0", "0.5", "1", "2", "3"))
stopifnot(
  identical(names(cdr_fills), c("0", "0.5", "1", "2", "3")),
  all(diff(colSums(grDevices::col2rgb(cdr_fills))) < 0)
)
filled <- tdp43_box_trace(
  ad, levels(fixture$Group), "AD",
  fill_colour = unname(sex_fills["Males"]), fill_alpha = 0.65
)
stopifnot(
  identical(filled$fillcolor, "rgba(86,180,233,0.65)"),
  identical(filled$line$color, unname(tdp43_palette()["AD"]))
)
message("Synthetic box statistics, trace and colour checks passed.")
