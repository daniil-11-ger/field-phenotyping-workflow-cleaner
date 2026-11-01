# FIELD CLEANER Studio Plus

This version improves the interactive UI for the sesame phenotype QC workflow:

- responsive buttons and labels that fit without browser zooming out;
- configurable project cards for crop, location, scale, design and goal;
- visible upload/status cards and workflow progress bar;
- quick example result cards for QC, trait plots, spatial map and error log;
- more modern blue/teal palette and Bootstrap icons;
- robust project path detection for demo dataset loading.

Run in RStudio from the project folder:

```r
shiny::runApp("shiny_app", launch.browser = TRUE)
```

If packages are missing:

```r
source("shiny_app/packages.R")
```
