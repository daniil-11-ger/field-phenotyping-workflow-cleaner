# 01_import_fieldbook.R
# Import and consolidate Field Book CSV files

source("R/00_config.R")
load_required_packages()

read_fieldbook_csv <- function(file_path) {
  readr::read_csv(
    file_path,
    na = missing_value_codes,
    show_col_types = FALSE,
    guess_max = 10000
  ) |>
    janitor::clean_names() |>
    mutate(source_file = basename(file_path))
}

import_fieldbook_data <- function(raw_dir = project_paths$raw_dir) {
  files <- list.files(raw_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No CSV files found in data/raw. Add Field Book exports first.")
  }
  purrr::map_dfr(files, read_fieldbook_csv)
}

standardize_fieldbook_columns <- function(df) {
  # Expected minimum columns: uuid, entry, block, row, column, trait, value.
  # Rename common Field Book variants if present.
  df |>
    rename_with(~ str_replace_all(.x, "range", "row")) |>
    mutate(
      uuid = as.character(coalesce(!!!rlang::syms(intersect(c("uuid", "unique_id", "plant_id"), names(.))))),
      entry = as.character(coalesce(!!!rlang::syms(intersect(c("entry", "accession", "line", "genotype"), names(.)))))
    )
}

raw_data <- import_fieldbook_data()
raw_data <- standardize_fieldbook_columns(raw_data)

readr::write_csv(raw_data, file.path(project_paths$processed_dir, "01_consolidated_raw.csv"))
message("Saved: data/processed/01_consolidated_raw.csv")
