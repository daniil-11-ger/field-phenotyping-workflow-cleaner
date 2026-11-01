# 03_quality_control.R
# Basic quality-control checks and summary tables

source("R/00_config.R")
load_required_packages()

input_file <- file.path(project_paths$processed_dir, "02_clean_wide.csv")
if (!file.exists(input_file)) stop("Run R/02_clean_transform.R first.")

df <- readr::read_csv(input_file, show_col_types = FALSE)

qc_missing <- df |>
  summarise(across(everything(), ~ sum(is.na(.x)))) |>
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "missing_n") |>
  mutate(missing_pct = round(100 * missing_n / nrow(df), 2)) |>
  arrange(desc(missing_pct))

qc_duplicates <- df |>
  count(uuid, name = "n") |>
  filter(n > 1)

numeric_summary <- df |>
  select(where(is.numeric)) |>
  summarise(across(everything(), list(
    n = ~ sum(!is.na(.x)),
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    min = ~ min(.x, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE),
    cv_pct = ~ 100 * sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)
  ), .names = "{.col}_{.fn}")) |>
  tidyr::pivot_longer(everything(), names_to = "metric", values_to = "value")

readr::write_csv(qc_missing, file.path(project_paths$tables_dir, "qc_missing_values.csv"))
readr::write_csv(qc_duplicates, file.path(project_paths$tables_dir, "qc_duplicate_uuid.csv"))
readr::write_csv(numeric_summary, file.path(project_paths$tables_dir, "qc_numeric_summary.csv"))

p <- qc_missing |>
  filter(!variable %in% c("uuid", "entry")) |>
  slice_max(missing_pct, n = 20) |>
  ggplot(aes(x = reorder(variable, missing_pct), y = missing_pct)) +
  geom_col() +
  coord_flip() +
  labs(x = "Variable", y = "Missing values, %", title = "Missing data by variable") +
  theme_minimal()

ggsave(file.path(project_paths$figures_dir, "missing_values_by_variable.png"), p, width = 8, height = 5, dpi = 300)
message("Saved QC tables and missing-value figure.")
