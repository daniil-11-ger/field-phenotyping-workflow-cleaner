# Sesame Field Phenotyping and Bioinformatics Workflow

This repository contains a reproducible technical workflow for cleaning, standardizing, visualizing, and modeling field phenotyping data collected for sesame (*Sesamum indicum* L.). The workflow is designed for Field Book CSV exports and prepares data for downstream breeding models such as BLUP/GBLUP and future phenotype-genotype integration.

## Project structure

```text
sesame_workflow_project/
├── R/
│   ├── 00_config.R
│   ├── 01_import_fieldbook.R
│   ├── 02_clean_transform.R
│   ├── 03_quality_control.R
│   ├── 04_outlier_detection.R
│   ├── 05_mixed_models.R
│   ├── 06_dna_extraction_summary.R
│   └── run_all.R
├── data/
│   ├── raw/
│   │   └── example_fieldbook_export.csv
│   └── processed/
├── metadata/
│   ├── data_dictionary.csv
│   └── biological_limits.csv
├── outputs/
│   ├── figures/
│   └── tables/
└── docs/
```

## Required software

- R version 4.2 or newer
- RStudio is recommended
- Required R packages: `tidyverse`, `janitor`, `lubridate`, `readr`, `stringr`, `lme4`, `broom.mixed`, `ggplot2`, `here`, `yaml`

The script `R/00_config.R` contains a helper function for installing missing packages.

## Input data format

The expected input is one or more CSV files exported from Field Book and placed in `data/raw/`.

Minimum required columns:

- `uuid`: unique plant/sample identifier
- `entry`: accession, line, or genotype name
- `block`: replication/block
- `row`: field row
- `column`: within-row coordinate or plant position
- `trait`: measured trait name
- `value`: measured value

Optional columns include `person`, `timestamp`, `device_name`, and `notes`.

## How to run

Open the project folder in RStudio and run:

```r
source("R/run_all.R")
```

Or run individual steps:

```r
source("R/01_import_fieldbook.R")
source("R/02_clean_transform.R")
source("R/03_quality_control.R")
source("R/04_outlier_detection.R")
source("R/05_mixed_models.R")
source("R/06_dna_extraction_summary.R")
```

## Workflow stages

1. **Import and consolidation**: reads all Field Book CSV files and combines them into one master file.
2. **Cleaning and transformation**: standardizes names, converts trait values to numeric format, removes duplicated records, and transforms data from long to wide format.
3. **Quality control**: calculates missing values, duplicate identifiers, and descriptive statistics.
4. **Outlier detection**: applies biological limits and IQR-based statistical detection.
5. **Mixed models**: fits simple linear mixed models with block as a fixed effect and genotype/entry as a random effect.
6. **DNA extraction summary**: summarizes HMW DNA extraction protocol comparison.

## Main outputs

Processed data:

- `data/processed/01_consolidated_raw.csv`
- `data/processed/02_clean_long.csv`
- `data/processed/02_clean_wide.csv`
- `data/processed/03_cleaned_outliers_removed.csv`

Tables:

- `outputs/tables/qc_missing_values.csv`
- `outputs/tables/qc_duplicate_uuid.csv`
- `outputs/tables/qc_numeric_summary.csv`
- `outputs/tables/outlier_log.csv`
- `outputs/tables/outlier_summary.csv`
- `outputs/tables/lmm_variance_components.csv`
- `outputs/tables/entry_blups.csv`
- `outputs/tables/dna_extraction_summary.csv`

Figures:

- `outputs/figures/missing_values_by_variable.png`
- `outputs/figures/boxplot_<trait>_by_block.png`
- `outputs/figures/dna_extraction_yield.png`

## Notes for article appendix

This package supports the technical reproducibility of the manuscript by documenting the complete data-processing path from raw field records to cleaned phenotype tables and preliminary breeding-model inputs. Real raw data can be substituted for the example CSV files without changing the analytical logic, provided that the required columns are present.

## Interactive dashboard

The project now includes an interactive web dashboard:

```r
source("shiny_app/packages.R")
shiny::runApp("shiny_app")
```

Main files:

- `shiny_app/app.R` — full R Shiny application.
- `shiny_app/www/app.css` — visual styling.
- `docs/interactive_app_guide_ru.md` — Russian guide for local launch and deployment.
- `deploy/shinyapps_publish.R` — shinyapps.io deployment template.
- `Dockerfile` and `deploy/docker_run.sh` — Docker deployment.
- `app_static_demo/index.html` — lightweight browser-only demo.

The Shiny app allows CSV upload, quality control, outlier detection, spatial field visualization, simple LMM/BLUP demonstration, DNA extraction comparison, and export of cleaned datasets.
