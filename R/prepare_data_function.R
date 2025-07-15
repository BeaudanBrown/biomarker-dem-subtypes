## helper function

prepare_data <- function(sheet_name, marker_name) {
  df <-
    read_excel(
      file.path(data_dir, Sys.getenv("ADDF_FILE")),
      sheet = sheet_name,
      guess_max = 100000,
    )

  setnames(df, "Sample", "Sample_Barcode", skip_absent = TRUE)
  setnames(df, "Sample ID", "Sample_Barcode", skip_absent = TRUE)

  setnames(df, names(df), gsub(" ", "_", names(df)))
  setnames(df, names(df), gsub("\\.", "", names(df)))

  setnames(df, names(df)[-1], paste0(marker_name, names(df)[-1]))

  # remove blank rows

  df <- filter(df, !is.na(Sample_Barcode))

  # delete mostly empty cols

  df <- df[, (colSums(is.na(df)) / nrow(df)) < 0.75]

  # if sample_barcode has an asteriks, add in dilution_note

  df[[paste0(marker_name, "Dilution_note")]] <-
    ifelse(
      str_detect(df$Sample_Barcode, "\\*"),
      "indicates diluted 1:4... (see excel)",
      NA_character_
    )

  # if sample_barcode has an asteriks, remove

  df$Sample_Barcode <- as.character(gsub("\\*", "", df$Sample_Barcode))
  df
}
