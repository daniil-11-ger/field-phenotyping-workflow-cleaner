# 00_config.R
# Project configuration for sesame phenotyping workflow

required_packages <- c(
  "tidyverse", "janitor", "lubridate", "readr", "stringr",
  "lme4", "broom.mixed", "ggplot2", "here", "yaml"
)

install_missing_packages <- function(pkgs = required_packages) {
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing) > 0) install.packages(missing, dependencies = TRUE)
}

load_required_packages <- function(pkgs = required_packages) {
  invisible(lapply(pkgs, library, character.only = TRUE))
}

project_paths <- list(
  raw_dir = "data/raw",
  processed_dir = "data/processed",
  metadata_dir = "metadata",
  figures_dir = "outputs/figures",
  tables_dir = "outputs/tables"
)

biological_limits <- tibble::tribble(
  ~trait, ~min_value, ~max_value, ~unit, ~action,
  "plant_height", 0, 250, "cm", "flag_or_na",
  "internode_length", 0, 50, "cm", "flag_or_na",
  "branch_number", 0, 80, "count", "flag_or_na",
  "capsule_count", 0, 1000, "count", "flag_or_na",
  "days_to_flowering", 15, 120, "days", "flag_or_na",
  "days_to_maturity", 40, 200, "days", "flag_or_na"
)

missing_value_codes <- c("", "NA", "N/A", "na", "null", "NULL", "-990", "-999", "9999")
