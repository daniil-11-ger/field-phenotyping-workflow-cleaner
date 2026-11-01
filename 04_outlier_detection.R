# 04_outlier_detection.R
# Detect outliers using biological limits and IQR rules

source("R/00_config.R")
load_required_packages()

input_file <- file.path(project_paths$processed_dir, "02_clean_wide.csv")
if (!file.exists(input_file)) stop("Run R/02_clean_transform.R first.")

df <- readr::read_csv(input_file, show_col_types = FALSE)

flag_iqr_outliers <- function(x, k = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  x < (q1 - k * iqr) | x > (q3 + k * iqr)
}

cleaned <- df
outlier_log <- list()

for (tr in biological_limits$trait) {
  if (!tr %in% names(cleaned)) next
  lim <- biological_limits |> filter(trait == tr)
  biological_flag <- cleaned[[tr]] < lim$min_value | cleaned[[tr]] > lim$max_value
  iqr_flag <- flag_iqr_outliers(cleaned[[tr]])
  combined_flag <- biological_flag | iqr_flag
  combined_flag[is.na(combined_flag)] <- FALSE

  outlier_log[[tr]] <- tibble(
    uuid = cleaned$uuid,
    entry = cleaned$entry,
    block = cleaned$block,
    row = cleaned$row,
    column = cleaned$column,
    trait = tr,
    value = cleaned[[tr]],
    biological_outlier = biological_flag,
    iqr_outlier = iqr_flag,
    outlier = combined_flag
  ) |>
    filter(outlier)

  cleaned[[paste0(tr, "_raw")]] <- cleaned[[tr]]
  cleaned[[tr]][combined_flag] <- NA_real_
}

outlier_table <- dplyr::bind_rows(outlier_log)
outlier_summary <- outlier_table |>
  count(trait, biological_outlier, iqr_outlier, name = "n") |>
  arrange(trait)

readr::write_csv(cleaned, file.path(project_paths$processed_dir, "03_cleaned_outliers_removed.csv"))
readr::write_csv(outlier_table, file.path(project_paths$tables_dir, "outlier_log.csv"))
readr::write_csv(outlier_summary, file.path(project_paths$tables_dir, "outlier_summary.csv"))

for (tr in intersect(biological_limits$trait, names(df))) {
  p <- cleaned |>
    ggplot(aes(x = block, y = .data[[tr]])) +
    geom_boxplot(na.rm = TRUE) +
    labs(x = "Block", y = tr, title = paste("Boxplot by block:", tr)) +
    theme_minimal()
  ggsave(file.path(project_paths$figures_dir, paste0("boxplot_", tr, "_by_block.png")), p, width = 7, height = 5, dpi = 300)
}

message("Saved cleaned dataset, outlier log, and boxplots.")
