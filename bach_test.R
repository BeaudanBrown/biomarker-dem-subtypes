library(readxl)
library(data.table)

load_data <- function() {
  bach <- as.data.table(read_excel("./data/BACH_spreadsheet.xlsx"))
  setnames(
    bach,
    c(
      "STUDY ID",
      "BARCODE",
      "DOB",
      "GENDER",
      "TYPE",
      "DATE COLLECTED",
      "PROCESSING TIME",
      "BOX #",
      "BOX POSITION",
      "RACK/DRAWER/ POSITION",
      "BOX LOCATION",
      "FREEZER LOCATION",
      "FREEZER ID",
      "COMMENTS"
    ),
    c(
      "id",
      "barcode",
      "dob",
      "gender",
      "type",
      "date",
      "processing_time",
      "box_no",
      "box_pos",
      "rack_pos",
      "box_loc",
      "freezer_loc",
      "freezer_id",
      "comments"
    )
  )
  bach <- bach[
    ,
    c(
      "id",
      "barcode",
      "dob",
      "gender",
      "type",
      "date",
      "processing_time",
      "box_no",
      "box_pos",
      "rack_pos",
      "box_loc",
      "freezer_loc",
      "freezer_id",
      "comments"
    )
  ]
  # ONE RECORD IS LABELLED AS AND####
  bach[, bach_id := gsub(".*(?:BACH_?|ANB)([0-9]+).*", "\\1", id, perl = TRUE)]
  bach <- bach[complete.cases(bach_id), ]
  bach <- bach[complete.cases(type), ]
  bach <- bach[nchar(bach_id) <= 4]
  bach[, bach_id := sprintf("%04d", as.numeric(bach_id))]
  bach
}


bach <- load_data()
# BACH Test
table(bach$type)
bach[, type := toupper(trimws(type))]
bach[, type := gsub("^[WHOLE\\s+BLOOD|W/BLD|WBLD].*$", "WHOLE_BLOOD", type)]
bach[, type := gsub("(SERUM|CSF|PLASMA)[-_].*", "\\1", type)]
table(bach[, c(type, bach_id)])

counts_dt <- bach[, .(count = .N), by = .(bach_id, type)]
bach[is.na(type), ]
summary_dt <- counts_dt[, .(mean_count = mean(count),
                            min_count = min(count),
                            max_count = max(count),
                            range = max(count) - min(count),
                            n_participants = .N),
                        by = type]

summary_dt <- bach[, .(mean = mean(value), 
                     min = min(value), 
                     max = max(value),
                     range = max(value) - min(value),
                     n = .N),
                 by = .(bach_id, type)]
bach$type
names(bach)
