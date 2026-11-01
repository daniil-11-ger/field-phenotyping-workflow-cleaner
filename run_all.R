# run_all.R
# Execute the complete reproducible workflow from raw Field Book CSV files to outputs

steps <- c(
  "R/01_import_fieldbook.R",
  "R/02_clean_transform.R",
  "R/03_quality_control.R",
  "R/04_outlier_detection.R",
  "R/05_mixed_models.R",
  "R/06_dna_extraction_summary.R"
)

for (script in steps) {
  message("Running: ", script)
  source(script)
}

message("Workflow complete. Check data/processed, outputs/tables, and outputs/figures.")
