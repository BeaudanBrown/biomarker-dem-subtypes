# TDP-43 exploratory report

## Organisation

- `tdp_targets.R`: independent `tdp43_*` target branch, registered once in `_targets.R`.
- `R/tdp_exploration.R`: preparation, plotting, summaries and exports.
- `R/tdp_exploration.qmd`: HTML presentation only; analyses live in targets.
- `R/tdp_widgets.R`: cowplot-to-Plotly conversion and HTML controls.
- `R/tdp_controls.js`: display-only browser state; no fitting or data retrieval.
- `tests/tdp_controls.test.cjs`: synthetic JavaScript control tests.

Existing analyses, supplementary regression tables and the existing report are
unchanged. This branch uses `df_with_csf` with its existing upstream cleaning.
The supplementary-table code assessment remains in `REVIEWER_RESPONSES.md`,
not in the exploratory report.

## Run (user only)

Use the same configured R environment and targets store as the existing pipeline.
The Nix environment now includes `plotly`, `htmlwidgets` and `htmltools`; re-enter
the development shell after updating it. No Shiny server is required:


```r
targets::tar_make(names = tdp43_report)
```

This builds only the report and its dependencies, not unrelated ROC analyses.
Missing or outdated upstream dependencies may still need to run.

Outputs:

```text
output/tdp43/
  tdp43_report.html
  figures/                 # PNGs, including both scatterplot trend variants
  tables/                  # Numeric CSVs without display rounding
```

`tar_quarto` also tracks the intermediate `R/tdp43_report.html`; the final file
target publishes a self-contained copy to the output directory. Both locations
contain the same report. Keep the intermediate file if avoiding a rerender matters.

To build only the PNG/CSV exports:

```r
targets::tar_make(names = tdp43_exports)
```

Individual analyses are `tdp43_subtype`, `tdp43_age`, `tdp43_sex`,
`tdp43_biomarker_correlations`, `tdp43_cdr`, `tdp43_pet` and `tdp43_csf`.
Their values contain named `plots` and `tables` lists.

## Agreed behaviour

- Seven numbered sections in a continuous, self-contained HTML document with a
  table of contents; no introductory methods text or supplementary-table audit.
  Interactive plots omit internal titles/subtitles/captions; section headings,
  controls and scatterplot sample counts remain. Warning text is not printed
  in the HTML; errors still stop rendering.
- All ggplots use `cowplot::theme_cowplot()`, including those converted to Plotly.
  Cohort colours are shared by scatterplots, box outlines, individual points
  and selector options: Controls grey, AD blue, FTD orange, LBD/PD green,
  and pooled All/Overall dark grey. The direct diagnostic-group plot also fills
  each box by diagnosis. Sex plots use blue male fills and pink-purple female
  fills; CDR fills progress from light to dark blue with increasing observed
  CDR level. Fills are defined once across the full report cohort, so selecting
  a diagnostic panel cannot change a level's shade. Sex/CDR box outlines and
  points retain the cohort colour. The HTML and PNGs use the same palettes.
  Correlation heatmap cells retain their diverging correlation-value scale.
- A single toggle, off by default, controls precomputed pooled and
  diagnosis-specific log-linear trends and 95% confidence bands on all three
  scatterplots. All (pooled) is an independent checkbox, initially off: it shows
  every participant in one colour with a separately fitted pooled trend. It
  never selects or deselects the diagnostic-group checkboxes. Diagnostic views
  start enabled and can be hidden independently. Pooled and subgroup layers
  may overlap; displayed totals count each participant once. Hidden layers
  stay hidden when trends are enabled. Axis ranges include all bands and stay
  fixed during selection; interactive zoom is still available.
- Scatterplots show the visible total and diagnosis-specific n below each
  widget. Complete counts remain in PNG captions and CSVs. Axes use
  unstandardised age, centiloids and CSF Aβ42/40.
- TDP-43 axes are logarithmic but retain concentration values, not z-scores.
  No unverified physical unit conversion is applied.
- A diagnostic-group boxplot compares Controls, AD, FTD and LBD/PD directly,
  with individual points and n, median (Q1, Q3). It exports `subtype.png` and
  `subtype.csv`, with no p-values.
- Sex/CDR boxplots have an Overall/diagnosis selector and matching summary
  tables with n and median (Q1, Q3). All boxplots have a point-visibility toggle;
  the direct diagnostic-group boxplot always retains all diagnostic groups.
  Box statistics are calculated in original concentration units before the
  logarithmic display transform. Interactive boxplots use a native Plotly log
  axis with untransformed trace values, so point and box-summary hovers report
  concentrations rather than log10 values. Each box is explicitly supplied
  with its own quartiles and 1.5-IQR whisker endpoints; ggplotly's automatically
  converted box traces are discarded. Whisker endpoints are not necessarily
  the observed sample minimum and maximum when outliers are present.
  Static PNGs also transform only
  the display, after calculating box statistics.
  Plotly's built-in box-trace points are disabled explicitly, so turning off
  individual points also removes outlier dots without hiding boxes or whiskers.
  No p-values. CDR is categorical.
- Spearman matrices include the nine biomarkers plus the labelled derived plasma
  Aβ42/40 ratio, overall and by diagnosis. Matrices use pairwise-complete values,
  fixed order and colour limits. A cohort selector switches between Overall
  (a pooled calculation) and diagnosis-specific matrices. Hover shows the
  biomarker pair, coefficient and pairwise n. Matching pairwise n tables are
  expandable in HTML and exported alongside full-precision coefficients.
- Ordinary missing measurements are omitted only for the relevant analysis.
  Invalid observed values cause errors rather than silent replacement/filtering.
  Positive TDP-43 and CSF ratio components are required. Missing sex excludes a
  participant from sex comparisons, not unrelated analyses.
- Undefined correlations fail clearly. No custom sparse-group checks are added
  to trend fitting; model/prediction errors surface naturally.
- HTML displays correlations to two decimals and concentration summaries to one;
  CSVs retain numeric precision. Counts are integers.

## Sharing and exports

The self-contained HTML embeds plotted participant-level coordinates and
cohort labels, without participant identifiers. It is intended for the agreed
colleague sharing; recipients can extract the plotted values. All widget
libraries and controls are embedded for offline use. No models are refitted
when controls change. PNG/CSV exports always contain the complete results,
including diagnosis-coloured and pooled scatterplots with and without trends,
regardless of browser selections. Pooled PNGs use `_all.png` and
`_all_trends.png` suffixes.

## Verification status

No raw data have been inspected and no R or Quarto execution has been performed
by the assistant. Synthetic browser-state tests can be run without R or data:

```sh
node --test tests/tdp_controls.test.cjs
```

A separate synthetic R test is provided for the user to run (not executed by
this assistant): `source("tests/tdp_box_statistics.R")`. It checks distinct
per-group quartiles/whiskers and explicit Plotly box traces without project data.

These tests do not establish R-to-Plotly conversion or rendered HTML correctness.
The user should render the report, open the copied HTML offline and check cohort
checkboxes, the global trend toggle, point toggles, selectors, hover values,
matching summary tables, fixed axes and complete exports. Do not use an R-backed
language server or test runner for assistant-side verification.
