# Reviewer response drafts

## Figure 1: model construction

### Response

All comparisons used the same nine biomarkers plus age and sex, and the same
SuperLearner candidate algorithms. However, models were trained separately
for each comparison, with model parameters and ensemble weights re-estimated
within each cross-validation training fold. No single diagnostic cut-off was
selected: ROC curves show performance across thresholds. The following text
is proposed for the Methods and figure legend.

### Proposed Methods clarification

Separate binary SuperLearner classifiers were trained for each diagnostic
comparison. All primary models used CD14, NfL, YKL-40, GFAP, Aβ42, Aβ40,
TDP-43, p-tau181 and p-tau217, together with age and sex. The candidate library
comprised logistic regression, penalised logistic regression, Bayesian
logistic regression, Bayesian additive regression trees, a generalised
additive model, multivariate adaptive regression splines, random forest,
and two gradient-boosted-tree configurations. Ensemble weights were estimated
using AUC optimisation (`method.AUC`) with 20-fold stratified internal
cross-validation. Performance was evaluated using 10-fold stratified outer
cross-validation. Candidate models and ensemble weights were fitted anew
within each outer training fold, independently for each comparison.

ROC curves were constructed by averaging held-out-fold sensitivity and
false-positive rate at each of 100 thresholds from 0 to 1. Displayed AUCs
are the arithmetic means of the outer-fold AUCs, calculated from held-out
predictions. No single diagnostic cut-off was selected.

### Proposed Figure 1 legend clarification

Each curve represents a separately trained binary classifier using the same
nine biomarkers, age and sex, and the same SuperLearner candidate library.
Model parameters and ensemble weights were re-estimated within each training
fold for each comparison. Curves show mean held-out-fold sensitivity and
false-positive rate across thresholds; labels report mean outer-fold AUC.
No single diagnostic cut-off was selected.

### Scope note

This revision clarifies the existing analysis; it does not alter model fitting
or preprocessing. In the current code, standardisation occurs before outer
fold creation and biomarker truncation occurs upstream. The clarification
should not be read as a claim that preprocessing is fitted within training
folds. The quoted manuscript AUCs have not been independently verified here.

## CDR/MMSE descriptive table

A separate table is specified using the same `df` population as the existing
Table 1 demographic summaries. Panels are All, Males and Females; columns are
Controls, AD, FTD and LBD/PD, without an overall column. Each panel reports
CDR and MMSE median (Q1, Q3), plus available n (%). Percentage denominators
are the participants with an observed score across the four diagnostic groups
in the same sex panel, separately for CDR and MMSE. Thus each percentage row
sums to approximately 100%, allowing for rounding, and describes diagnostic
composition rather than within-group score availability. Scores, quartiles and percentages use one
decimal place; counts are integers. No p-values are included.

Current CDR source precedence is retained. Race and ethnicity revisions are
deferred. No additional range/missing-code validation step is included, as
requested. These changes have not been executed in R; generated tables remain
for the user to produce and inspect.

## Supplementary Tables 5 and 6: code assessment (outside TDP-43 report)

Status: source inspection only; numerical results have not been verified.
These notes do not authorise changes to the existing analyses. Verification
and corrections remain separate from the planned exploratory report.

### Pipeline mapping

- Plasma biomarkers versus amyloid PET centiloids: `pet_centiloid_corrs`
  and `pet_centiloid_corr_table` in `pet_targets.R`.
- Plasma biomarkers versus CSF Aβ42/40: `csf_corrs` and `csf_corr_table`
  in `csf_targets.R`.
- Both use `get_adjusted_corr()` and `get_marker_table()` in
  `R/correlations.R`.

### Interpretation and labelling

The models regress the standardised outcome on one standardised plasma
biomarker, age and sex. TDP-43 and YKL-40 are log-transformed before
standardisation. The extracted estimate and p-value correspond to the
biomarker term in `summary(lm())`.

The manuscript label “R” is not appropriate: these are adjusted standardised
regression coefficients (β), not Pearson, Spearman or partial correlations.
Recommended wording is “Standardised β” and “age- and sex-adjusted linear
regression”, with transformations documented. Actual partial correlations
would require a separate calculation, not relabelling the estimates.
Standardised multiple-regression coefficients need not lie between −1 and 1.

### Reporting defects to address separately

- P-values are formatted to two decimal places before the `<0.001` check.
  For example, 0.004 can become `0.00` and be labelled `<0.001`. Apply the
  threshold to the original numeric value before formatting.
- Reported n excludes missing outcome/predictor values but does not account
  for further exclusions of missing age or sex by `lm()`. Use actual model n.
- Cohort headers display the maximum n across biomarkers, concealing
  model-specific sample sizes. Preserve and report model-specific n.
- Coefficients are rounded to three decimals and subsequently two, allowing
  double-rounding discrepancies. Retain numeric full precision until display.

### Methodological caveats, not confirmed numerical errors

- Standardisation occurs within each cohort before model-specific exclusions;
  SD units therefore need not correspond to the exact fitted sample.
- Pooled models do not adjust for diagnosis or site; pooled associations can
  reflect differences between groups rather than within-group relationships.
- Ordinary linear-model p-values rely on appropriate residual behaviour and
  independent observations. These assumptions cannot be verified from source.
- No multiple-testing correction is applied; exploratory interpretation should
  acknowledge this.
- Upstream CSF handling can supplement Aβ42 from another source and convert
  censored measurements to boundary values. Assay compatibility and censoring
  implications require assessment rather than assumptions.

Recommendation: correct reporting separately, then have the user regenerate
outputs and compare them with the published tables. No R execution or raw-data
inspection by the assistant is permitted.

