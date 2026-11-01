# 05_mixed_models.R
# Prepare and fit basic linear mixed models for breeding-value-oriented analysis

source("R/00_config.R")
load_required_packages()

input_file <- file.path(project_paths$processed_dir, "03_cleaned_outliers_removed.csv")
if (!file.exists(input_file)) stop("Run R/04_outlier_detection.R first.")

df <- readr::read_csv(input_file, show_col_types = FALSE) |>
  mutate(
    entry = as.factor(entry),
    block = as.factor(block),
    row = as.factor(row),
    column = as.factor(column)
  )

fit_lmm_for_trait <- function(data, trait) {
  if (!trait %in% names(data)) return(NULL)
  model_data <- data |>
    select(entry, block, row, column, all_of(trait)) |>
    filter(!is.na(.data[[trait]]))
  if (nrow(model_data) < 20 || n_distinct(model_data$entry) < 2) return(NULL)

  formula <- as.formula(paste0(trait, " ~ block + (1|entry)"))
  model <- lme4::lmer(formula, data = model_data, REML = TRUE)

  variance_components <- as.data.frame(VarCorr(model)) |>
    mutate(trait = trait)

  entry_blups <- ranef(model)$entry |>
    tibble::rownames_to_column("entry") |>
    rename(blup = `(Intercept)`) |>
    mutate(trait = trait)

  list(model = model, variance_components = variance_components, entry_blups = entry_blups)
}

traits_to_model <- intersect(c("plant_height", "internode_length", "branch_number", "capsule_count", "days_to_flowering"), names(df))
model_results <- purrr::map(traits_to_model, ~ fit_lmm_for_trait(df, .x))
model_results <- model_results[!purrr::map_lgl(model_results, is.null)]

variance_table <- purrr::map_dfr(model_results, "variance_components")
blup_table <- purrr::map_dfr(model_results, "entry_blups")

readr::write_csv(variance_table, file.path(project_paths$tables_dir, "lmm_variance_components.csv"))
readr::write_csv(blup_table, file.path(project_paths$tables_dir, "entry_blups.csv"))
message("Saved LMM variance components and entry BLUPs.")
