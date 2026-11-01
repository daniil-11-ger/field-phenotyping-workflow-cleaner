# FIELD CLEANER | Phenotype QC Studio

This is the enhanced interactive version of the sesame phenotyping workflow.

Main improvements:

- modern dashboard title and visual structure;
- automatic Field Book column detection;
- configurable column mapping;
- one-click full workflow execution;
- processing progress modal and persistent percent indicator;
- configurable biological ranges for each trait;
- QC summary with overall error map;
- full error log and top-problematic-trait table;
- trait explorer, spatial QC, basic BLUP demo, DNA module;
- presentation story panel based on the provided slide content.

## Run in RStudio

```r
setwd("D:/Windows/Рабочий стол/стат/стат 3 конф коз готовая стат и заявк (готовое прилодение)/код ии приложение/4 field_cleaner_studio/field_cleaner_studio")
shiny::runApp("shiny_app", launch.browser = TRUE)
```

If packages are missing:

```r
source("shiny_app/packages.R")
```

Then run again:

```r
shiny::runApp("shiny_app", launch.browser = TRUE)
```
