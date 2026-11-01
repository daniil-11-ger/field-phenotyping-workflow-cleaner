# 02_clean_transform.R
# Clean, type-convert, and transform long Field Book data to wide plant-level data

source("R/00_config.R")
load_required_packages()

input_file <- file.path(project_paths$processed_dir, "01_consolidated_raw.csv")
if (!file.exists(input_file)) stop("Run R/01_import_fieldbook.R first.")

df <- readr::read_csv(input_file, show_col_types = FALSE)

required_cols <- c("uuid", "entry", "block", "row", "column", "trait", "value")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

clean_long <- df |>
  mutate(
    trait = trait |>
      stringr::str_to_lower() |>
      stringr::str_replace_all("[^a-z0-9]+", "_") |>
      stringr::str_replace_all("_$", ""),
    value_chr = as.character(value),
    value_num = suppressWarnings(as.numeric(value_chr)),
    block = as.factor(block),
    row = as.integer(row),
    column = as.integer(column),
    entry = as.factor(entry),
    uuid = as.character(uuid)
  ) |>
  distinct(uuid, entry, block, row, column, trait, .keep_all = TRUE)

wide_data <- clean_long |>
  select(uuid, entry, block, row, column, trait, value_num) |>
  tidyr::pivot_wider(
    names_from = trait,
    values_from = value_num,
    values_fn = list(value_num = ~ mean(.x, na.rm = TRUE))
  ) |>
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA_real_, .x)))

readr::write_csv(clean_long, file.path(project_paths$processed_dir, "02_clean_long.csv"))
readr::write_csv(wide_data, file.path(project_paths$processed_dir, "02_clean_wide.csv"))
message("Saved cleaned long and wide datasets.")
