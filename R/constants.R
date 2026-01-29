## Central constants and mappings for biomarker analysis

# Biomarker variable names to display names
BIOMARKER_NAMES <- c(
  "mean_elisa" = "CD14",
  "mean_nfl" = "NfL",
  "mean_ykl" = "YKL-40",
  "mean_gfap" = "GFAP",
  "mean_ab40" = "AB40",
  "mean_ab42" = "AB42",
  "mean_tdp" = "TDP-43",
  "mean_ptau181" = "pTau-181",
  "mean_ptau217" = "pTau-217",
  "mean_ab42_ab40_ratio" = "AB42/AB40 Ratio",
  "cdr" = "CDR",
  "zscore" = "PET Z-Score",
  "centiloid" = "PET Centiloid",
  "raw_suvr" = "PET Raw SUVR"
)

# Helper function to rename biomarker variables to display names
rename_biomarkers <- function(x) {
  ifelse(x %in% names(BIOMARKER_NAMES), BIOMARKER_NAMES[x], x)
}

# Base columns for ROC analysis
BASE_ROC_COLUMNS <- c(
  "Diagnosis_combined",
  "age",
  "mean_elisa",
  "mean_nfl",
  "mean_ykl",
  "mean_gfap",
  "mean_ab40",
  "mean_ab42",
  "mean_tdp",
  "mean_ptau181",
  "mean_ptau217",
  "female"
)
