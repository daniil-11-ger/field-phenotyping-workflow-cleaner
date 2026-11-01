# FIELD CLEANER Studio: interactive Shiny dashboard for field phenotyping QC
# Designed for Field Book CSV exports and sesame breeding workflow technical appendix

options(shiny.maxRequestSize = 200 * 1024^2)
options(shiny.reactlog = FALSE)

required_packages <- c(
  "shiny", "bslib", "bsicons", "DT", "plotly", "dplyr", "tidyr", "readr", "janitor",
  "ggplot2", "stringr", "lubridate", "purrr", "lme4", "broom.mixed", "tibble"
)

missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them with: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

start_dir <- normalizePath(getwd(), mustWork = TRUE)
if (basename(start_dir) == "shiny_app") {
  app_dir <- start_dir
  project_dir <- dirname(start_dir)
} else if (dir.exists(file.path(start_dir, "shiny_app"))) {
  project_dir <- start_dir
  app_dir <- file.path(project_dir, "shiny_app")
} else {
  app_dir <- start_dir
  project_dir <- dirname(start_dir)
}
example_path <- file.path(project_dir, "data", "raw", "example_fieldbook_export.csv")
limits_path <- file.path(project_dir, "metadata", "biological_limits.csv")

missing_value_codes <- c("", "NA", "N/A", "na", "null", "NULL", "-990", "-999", "9999")

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}


slide_facts <- list(
  title = "Интеграция полевых фенотипических данных и биоинформатического анализа",
  crop = "Sesamum indicum L. / кунжут, F4",
  location = "KAUST, Thuwal, Saudi Arabia",
  scale = ">5000 observations, 15 variables, 20 stable lines",
  design = "RCBD, 3 blocks",
  problem = "Garbage in - garbage out: raw field noise can distort BLUP/GWAS outputs",
  target = "Clean, model-ready phenotype matrix for LMM, BLUP, GBLUP, GWAS"
)

default_limits <- tibble::tribble(
  ~trait, ~min_value, ~max_value, ~unit, ~action,
  "plant_height", 0, 250, "cm", "flag_or_na",
  "internode_length", 0, 50, "cm", "flag_or_na",
  "branch_number", 0, 80, "count", "flag_or_na",
  "capsule_count", 0, 1000, "count", "flag_or_na",
  "days_to_flowering", 15, 120, "days", "flag_or_na",
  "days_to_maturity", 40, 200, "days", "flag_or_na",
  "maturity_date", NA_real_, NA_real_, "date", "keep",
  "lodging_status", NA_real_, NA_real_, "category", "keep",
  "disease", NA_real_, NA_real_, "category", "keep"
)

synonyms <- list(
  uuid = c("uuid", "unique_id", "sample_id", "plant_uuid", "id"),
  entry = c("entry", "line", "genotype", "accession", "pi", "cultivar"),
  block = c("block", "rep", "replicate", "replication"),
  row = c("row", "range", "field_row"),
  column = c("column", "col", "plant_position", "position"),
  plant_id = c("plant_id", "plantid", "plant", "plant_number"),
  trait = c("trait", "variable", "trait_name", "observation_variable"),
  value = c("value", "trait_value", "measurement", "score"),
  timestamp = c("timestamp", "time", "date_time", "created_at"),
  person = c("person", "operator", "collector", "observer"),
  notes = c("notes", "note", "comment", "comments", "observation"),
  disease = c("disease", "phyllody", "phyllody_disease", "pest", "health_status")
)

clean_col <- function(x) janitor::make_clean_names(x)

read_fieldbook <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, na = missing_value_codes, guess_max = 10000) |>
    janitor::clean_names() |>
    mutate(across(where(is.character), ~ stringr::str_squish(.x)))
}

detect_col <- function(names_vec, field) {
  s <- synonyms[[field]]
  if (is.null(s)) return("")
  hit <- intersect(s, names_vec)
  if (length(hit) > 0) hit[[1]] else ""
}

is_blank <- function(x) is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")

choose_vec <- function(cols) c("Not available" = "", stats::setNames(cols, cols))

numericish_cols <- function(df) {
  cols <- names(df)
  cols[purrr::map_lgl(df, ~ {
    if (is.numeric(.x) || is.integer(.x)) return(TRUE)
    xx <- suppressWarnings(as.numeric(as.character(.x)))
    mean(!is.na(xx)) >= 0.60
  })]
}

prepare_long <- function(df, mapping, selected_traits = NULL) {
  df <- janitor::clean_names(df)
  map_value <- function(field) {
    val <- mapping[[field]]
    if (is_blank(val) || !val %in% names(df)) NA_character_ else val
  }

  uuid_col <- map_value("uuid")
  entry_col <- map_value("entry")
  block_col <- map_value("block")
  row_col <- map_value("row")
  column_col <- map_value("column")
  plant_id_col <- map_value("plant_id")
  trait_col <- map_value("trait")
  value_col <- map_value("value")
  timestamp_col <- map_value("timestamp")
  person_col <- map_value("person")
  notes_col <- map_value("notes")
  disease_col <- map_value("disease")

  add_meta <- function(tbl) {
    tbl |>
      mutate(
        uuid = if (!is.na(uuid_col)) as.character(df[[uuid_col]]) else if (!is.na(plant_id_col)) as.character(df[[plant_id_col]]) else paste0("plant_", row_number()),
        entry = if (!is.na(entry_col)) as.character(df[[entry_col]]) else "unknown_entry",
        block = if (!is.na(block_col)) as.character(df[[block_col]]) else NA_character_,
        row = if (!is.na(row_col)) suppressWarnings(as.integer(df[[row_col]])) else NA_integer_,
        column = if (!is.na(column_col)) suppressWarnings(as.integer(df[[column_col]])) else NA_integer_,
        plant_id = if (!is.na(plant_id_col)) as.character(df[[plant_id_col]]) else uuid,
        timestamp = if (!is.na(timestamp_col)) suppressWarnings(lubridate::ymd_hms(df[[timestamp_col]])) else as.POSIXct(NA),
        person = if (!is.na(person_col)) as.character(df[[person_col]]) else NA_character_,
        notes = if (!is.na(notes_col)) as.character(df[[notes_col]]) else NA_character_,
        disease_status = if (!is.na(disease_col)) as.character(df[[disease_col]]) else NA_character_
      )
  }

  if (!is.na(trait_col) && !is.na(value_col)) {
    long <- tibble(row_index = seq_len(nrow(df)), trait = as.character(df[[trait_col]]), value_raw = as.character(df[[value_col]])) |>
      add_meta()
  } else {
    meta_cols <- c(uuid_col, entry_col, block_col, row_col, column_col, plant_id_col, timestamp_col, person_col, notes_col, disease_col)
    meta_cols <- meta_cols[!is.na(meta_cols)]
    trait_candidates <- selected_traits
    if (is.null(trait_candidates) || length(trait_candidates) == 0) {
      trait_candidates <- setdiff(numericish_cols(df), meta_cols)
    }
    validate_traits <- intersect(trait_candidates, names(df))
    if (length(validate_traits) == 0) stop("No trait columns selected or detected. Select numeric trait columns in the Control Center.")

    long <- df |>
      mutate(row_index = row_number()) |>
      tidyr::pivot_longer(cols = all_of(validate_traits), names_to = "trait", values_to = "value_raw", values_transform = list(value_raw = as.character)) |>
      select(row_index, trait, value_raw, everything())
    # Rebuild metadata from pivoted rows safely
    long <- long |>
      mutate(
        uuid = if (!is.na(uuid_col)) as.character(.data[[uuid_col]]) else if (!is.na(plant_id_col)) as.character(.data[[plant_id_col]]) else paste0("plant_", row_index),
        entry = if (!is.na(entry_col)) as.character(.data[[entry_col]]) else "unknown_entry",
        block = if (!is.na(block_col)) as.character(.data[[block_col]]) else NA_character_,
        row = if (!is.na(row_col)) suppressWarnings(as.integer(.data[[row_col]])) else NA_integer_,
        column = if (!is.na(column_col)) suppressWarnings(as.integer(.data[[column_col]])) else NA_integer_,
        plant_id = if (!is.na(plant_id_col)) as.character(.data[[plant_id_col]]) else uuid,
        timestamp = if (!is.na(timestamp_col)) suppressWarnings(lubridate::ymd_hms(.data[[timestamp_col]])) else as.POSIXct(NA),
        person = if (!is.na(person_col)) as.character(.data[[person_col]]) else NA_character_,
        notes = if (!is.na(notes_col)) as.character(.data[[notes_col]]) else NA_character_,
        disease_status = if (!is.na(disease_col)) as.character(.data[[disease_col]]) else NA_character_
      ) |>
      select(row_index, uuid, plant_id, entry, block, row, column, trait, value_raw, timestamp, person, notes, disease_status)
  }

  long |>
    mutate(
      trait = janitor::make_clean_names(trait),
      value_num = suppressWarnings(as.numeric(value_raw)),
      block = as.factor(block),
      value_is_missing = is.na(value_raw) | value_raw %in% missing_value_codes,
      value_is_non_numeric = !value_is_missing & is.na(value_num),
      disease_flag = str_detect(str_to_lower(coalesce(disease_status, notes, "")), "yes|phyllody|disease|pest|pilostethus|pandurus|бол|вред|филлод")
    )
}

make_wide <- function(cleaned_long) {
  id_cols <- intersect(c("uuid", "plant_id", "entry", "block", "row", "column", "disease_status"), names(cleaned_long))
  cleaned_long |>
    group_by(across(all_of(c(id_cols, "trait")))) |>
    summarise(value_clean = dplyr::first(na.omit(value_clean)), .groups = "drop") |>
    tidyr::pivot_wider(names_from = trait, values_from = value_clean)
}

make_limits_for_traits <- function(long_df, uploaded_limits = NULL) {
  observed <- long_df |>
    filter(!is.na(value_num)) |>
    group_by(trait) |>
    summarise(
      obs_min = suppressWarnings(min(value_num, na.rm = TRUE)),
      obs_max = suppressWarnings(max(value_num, na.rm = TRUE)),
      q01 = suppressWarnings(quantile(value_num, 0.01, na.rm = TRUE)),
      q99 = suppressWarnings(quantile(value_num, 0.99, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(
      min_value = floor(if_else(is.finite(q01), q01, obs_min)),
      max_value = ceiling(if_else(is.finite(q99), q99, obs_max)),
      unit = "auto",
      action = "flag_or_na"
    ) |>
    select(trait, min_value, max_value, unit, action)

  limits <- default_limits
  if (!is.null(uploaded_limits)) {
    limits <- uploaded_limits |> janitor::clean_names() |> select(any_of(c("trait", "min_value", "max_value", "unit", "action")))
  } else if (file.exists(limits_path)) {
    limits <- readr::read_csv(limits_path, show_col_types = FALSE) |> janitor::clean_names()
  }

  observed |>
    rows_update(limits |> mutate(trait = janitor::make_clean_names(trait)), by = "trait", unmatched = "ignore") |>
    rows_insert(limits |> mutate(trait = janitor::make_clean_names(trait)), by = "trait", conflict = "ignore") |>
    filter(trait %in% unique(long_df$trait)) |>
    mutate(
      min_value = suppressWarnings(as.numeric(min_value)),
      max_value = suppressWarnings(as.numeric(max_value)),
      unit = if_else(is.na(unit) | unit == "", "auto", unit),
      action = if_else(is.na(action) | action == "", "flag_or_na", action)
    )
}

flag_range_outliers <- function(long_df, limits) {
  long_df |>
    left_join(limits, by = "trait") |>
    mutate(
      range_outlier = !is.na(value_num) & !is.na(min_value) & !is.na(max_value) &
        (value_num < min_value | value_num > max_value),
      range_reason = case_when(
        value_is_missing ~ "missing_value",
        value_is_non_numeric ~ "non_numeric_value",
        range_outlier ~ "outside_biological_limits",
        TRUE ~ "ok"
      )
    )
}

flag_iqr_outliers <- function(long_df) {
  stats <- long_df |>
    filter(!is.na(value_num)) |>
    group_by(trait) |>
    summarise(
      q1 = quantile(value_num, 0.25, na.rm = TRUE),
      q3 = quantile(value_num, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      lower_iqr = q1 - 1.5 * iqr,
      upper_iqr = q3 + 1.5 * iqr,
      .groups = "drop"
    ) |>
    mutate(lower_iqr = if_else(iqr == 0, -Inf, lower_iqr), upper_iqr = if_else(iqr == 0, Inf, upper_iqr))

  long_df |>
    left_join(stats, by = "trait") |>
    mutate(
      iqr_outlier = !is.na(value_num) & !is.na(lower_iqr) & !is.na(upper_iqr) &
        (value_num < lower_iqr | value_num > upper_iqr),
      iqr_reason = if_else(iqr_outlier, "outside_iqr_limits", "ok")
    )
}

clean_long <- function(long_df, limits, exclude_disease_from_yield = TRUE) {
  flagged <- long_df |>
    flag_range_outliers(limits) |>
    flag_iqr_outliers()

  yield_like <- str_detect(trait, "capsule|yield|height|branch|internode|biomass")

  flagged |>
    mutate(
      disease_yield_flag = isTRUE(exclude_disease_from_yield) & disease_flag & yield_like,
      value_clean = if_else(range_outlier | iqr_outlier | disease_yield_flag | value_is_non_numeric, NA_real_, value_num),
      final_status = case_when(
        value_is_missing ~ "missing",
        value_is_non_numeric ~ "non_numeric",
        range_outlier ~ "range_error",
        iqr_outlier ~ "statistical_outlier",
        disease_yield_flag ~ "biological_confounder",
        TRUE ~ "kept"
      )
    )
}

summarise_qc <- function(cleaned_long) {
  cleaned_long |>
    group_by(trait) |>
    summarise(
      records = n(),
      plants = n_distinct(uuid),
      missing = sum(final_status == "missing", na.rm = TRUE),
      non_numeric = sum(final_status == "non_numeric", na.rm = TRUE),
      range_errors = sum(final_status == "range_error", na.rm = TRUE),
      statistical_outliers = sum(final_status == "statistical_outlier", na.rm = TRUE),
      biological_confounders = sum(final_status == "biological_confounder", na.rm = TRUE),
      kept = sum(final_status == "kept", na.rm = TRUE),
      flagged_total = records - kept,
      flagged_pct = round(100 * flagged_total / records, 2),
      mean_clean = round(mean(value_clean, na.rm = TRUE), 2),
      sd_clean = round(sd(value_clean, na.rm = TRUE), 2),
      cv_clean_pct = round(100 * sd(value_clean, na.rm = TRUE) / mean(value_clean, na.rm = TRUE), 2),
      .groups = "drop"
    ) |>
    arrange(desc(flagged_pct), desc(records))
}

make_error_long <- function(qc) {
  qc |>
    select(trait, missing, non_numeric, range_errors, statistical_outliers, biological_confounders) |>
    pivot_longer(-trait, names_to = "error_type", values_to = "count") |>
    mutate(error_type = factor(error_type, levels = c("missing", "non_numeric", "range_errors", "statistical_outliers", "biological_confounders")))
}

fit_simple_lmm <- function(wide_df, trait_name) {
  if (!trait_name %in% names(wide_df)) return(tibble(message = "Trait not found"))
  if (!all(c("entry", "block") %in% names(wide_df))) return(tibble(message = "Columns entry and block are required for LMM"))

  model_df <- wide_df |>
    select(entry, block, all_of(trait_name)) |>
    rename(y = all_of(trait_name)) |>
    filter(!is.na(y), !is.na(entry), !is.na(block))

  if (nrow(model_df) < 5 || length(unique(model_df$entry)) < 2) return(tibble(message = "Not enough observations for model"))

  fit <- tryCatch(lme4::lmer(y ~ block + (1 | entry), data = model_df, REML = TRUE), error = function(e) e)
  if (inherits(fit, "error")) return(tibble(message = fit$message))

  lme4::ranef(fit)$entry |>
    tibble::rownames_to_column("entry") |>
    rename(blup = `(Intercept)`) |>
    arrange(desc(blup))
}

dna_table <- tibble::tribble(
  ~protocol, ~sample, ~dsdna_ng_ul, ~ratio_260_280, ~ratio_260_230, ~integrity, ~interpretation,
  "CTAB-PK", "Leaf 1", 250.3, 1.83, 2.01, "HMW > 50 kb", "Optimal yield and purity",
  "CTAB-PK", "Leaf 2", 215.1, 1.79, 2.15, "HMW > 50 kb", "Optimal yield and purity",
  "Qiagen Genomics-Tip", "Leaf 1", 38.7, 1.45, 1.05, "Medium / partial fragmentation", "Lower yield; impurities present",
  "Qiagen Genomics-Tip", "Leaf 2", 42.4, 1.40, 0.98, "Medium / partial fragmentation", "Lower yield; impurities present"
)

explain_card <- function(title, body, bullets = NULL) {
  card(
    class = "help-card",
    card_header(tags$span(class = "section-icon", "i"), title),
    p(body),
    if (!is.null(bullets)) tags$ul(lapply(bullets, tags$li))
  )
}

metric_box <- function(title, value, icon, class = "metric-teal") {
  value_box(title = title, value = value, showcase = bsicons::bs_icon(icon), class = class)
}

nav_label <- function(icon, text) {
  span(class = "nav-label", bsicons::bs_icon(icon), span(text))
}

quick_result_card <- function(id, icon, title, text) {
  div(class = "quick-card",
      div(class = "quick-icon", bsicons::bs_icon(icon)),
      h4(title),
      p(text),
      actionButton(id, "Показать пример", class = "btn btn-outline-primary quick-btn"))
}

ui <- page_navbar(
  title = div(class = "brand-title", span("FIELD CLEANER"), tags$small("Phenotype QC Studio")),
  theme = bs_theme(version = 5, bootswatch = "minty", primary = "#2563eb", secondary = "#14b8a6", success = "#10b981", base_font = font_google("Inter")),
  header = tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "app.css")),
  window_title = "FIELD CLEANER | Phenotype QC Studio",

  nav_panel(
    nav_label("speedometer2", "Overview"),
    div(class = "hero",
        div(class = "hero-text",
            h1("FIELD CLEANER"),
            h3("Интерактивная очистка полевых фенотипических данных"),
            p("От Field Book CSV до model-ready phenotype matrix для BLUP, GWAS и геномной селекции."),
            div(class = "hero-tags", span("Field Book"), span("R Workflow"), span("QC"), span("Outliers"), span("Spatial effects"), span("BLUP-ready"))
        ),
        div(class = "hero-flow",
            div(class = "flow-node", "Raw Field Data"), div(class = "flow-arrow", "→"),
            div(class = "flow-node", "QC Engine"), div(class = "flow-arrow", "→"),
            div(class = "flow-node", "Model-Ready Data")
        )
    ),
    uiOutput("project_metric_boxes"),
    layout_columns(
      explain_card("Зачем это нужно", "Перспективы геномной селекции снижаются, если исходные полевые данные содержат ошибки, выбросы и биологический шум.", c("человеческий фактор", "экологический шум", "стресс 38-40 C+", "болезни и вредители")),
      explain_card("Что делает приложение", "Приложение автоматически определяет структуру CSV, предлагает колонки, проверяет диапазоны, считает ошибки и показывает, где именно данные требуют внимания.", c("автоопределение столбцов", "настройка биологических диапазонов", "общая карта ошибок", "экспорт очищенной матрицы")),
      col_widths = c(6, 6)
    ),
    card(card_header("Workflow map"),
         div(class = "workflow-strip",
             div("1. Upload"), div("2. Auto map"), div("3. Limits"), div("4. Clean"), div("5. QC"), div("6. Export")
         ))
  ),

  nav_panel(
    nav_label("sliders", "Project Setup"),
    layout_sidebar(
      sidebar = sidebar(
        width = 390,
        h4("Настраиваемая карточка проекта"),
        textInput("project_crop", "Crop / культура", value = slide_facts$crop),
        textInput("project_location", "Location / локация", value = slide_facts$location),
        textInput("project_scale", "Scale / масштаб", value = slide_facts$scale),
        textInput("project_design", "Design / дизайн опыта", value = slide_facts$design),
        textAreaInput("project_goal", "Goal / цель", value = slide_facts$target, rows = 3),
        actionButton("apply_project", "Apply project card", class = "btn-primary modern-btn")
      ),
      layout_columns(
        explain_card("Для чего эта вкладка", "Эти поля формируют карточки на главном экране и подписи в отчете. Их можно менять под другую культуру, другую площадку, другой год или другой эксперимент."),
        card(card_header("Live project preview"), uiOutput("project_preview")),
        col_widths = c(5, 7)
      ),
      card(card_header("Что можно менять"),
           div(class = "settings-grid",
               div(class = "setting-pill", bsicons::bs_icon("flower1"), "Crop / generation"),
               div(class = "setting-pill", bsicons::bs_icon("geo-alt"), "Location / site"),
               div(class = "setting-pill", bsicons::bs_icon("database"), "Scale / variables"),
               div(class = "setting-pill", bsicons::bs_icon("diagram-3"), "Design / blocks"),
               div(class = "setting-pill", bsicons::bs_icon("bullseye"), "Goal / model target")
           ))
    )
  ),

  nav_panel(
    nav_label("cloud-upload", "Control Center"),
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        h4("1. Загрузка данных"),
        fileInput("fieldbook_file", "Upload Field Book CSV", accept = c(".csv"), buttonLabel = "Browse...", placeholder = "No file selected"),
        checkboxInput("use_demo", "Use demo dataset if no file uploaded", TRUE),
        actionButton("quick_demo", "Load demo + run", class = "btn-success modern-btn"),
        helpText("После загрузки файла появится статус: rows, columns, detected format."),
        hr(),
        h4("2. Авто-настройка колонок"),
        uiOutput("mapping_ui"),
        hr(),
        h4("3. Признаки для анализа"),
        checkboxInput("auto_select_traits", "Analyze all detected traits", TRUE),
        uiOutput("trait_selector_ui"),
        checkboxInput("exclude_disease", "Exclude diseased plants from yield-like traits", TRUE),
        hr(),
        h4("4. Запуск"),
        actionButton("run", "Run workflow", class = "btn-primary btn-lg modern-btn"),
        helpText("Прогресс показывается в браузере и в Console. После 100% можно открывать все вкладки."),
        downloadButton("download_clean", "Clean matrix"),
        downloadButton("download_qc", "QC summary")
      ),
      layout_columns(
        card(card_header("Upload and processing status"),
             uiOutput("data_status_cards"),
             uiOutput("progress_bar_ui"),
             textOutput("processing_status"),
             DTOutput("run_log")),
        card(card_header("Raw data preview"), DTOutput("raw_preview")),
        col_widths = c(4, 8)
      ),
      explain_card("Что происходит на этой вкладке", "Здесь приложение загружает CSV, автоматически предлагает соответствие колонок и запускает весь workflow одним нажатием. После запуска остальные вкладки используют уже подготовленные результаты, а не пересчитывают всё заново.")
    )
  ),

  nav_panel(
    nav_label("magic", "Column Autopilot"),
    layout_columns(
      explain_card("Автоматическое определение", "Приложение ищет стандартные поля Field Book: uuid, entry, block, row, column, trait, value, timestamp, notes. Если названия отличаются, их можно вручную выбрать в Control Center."),
      card(card_header("Detected columns"), DTOutput("detected_columns")),
      col_widths = c(4, 8)
    ),
    card(card_header("Trait preview"), DTOutput("trait_preview"))
  ),

  nav_panel(
    nav_label("rulers", "Limits Builder"),
    layout_sidebar(
      sidebar = sidebar(
        selectInput("limit_trait", "Trait", choices = NULL),
        numericInput("limit_min", "Minimum allowed value", value = NA),
        numericInput("limit_max", "Maximum allowed value", value = NA),
        textInput("limit_unit", "Unit", value = "auto"),
        actionButton("save_limit", "Save range for this trait", class = "btn-primary"),
        fileInput("limits_file", "Optional: upload biological_limits.csv", accept = c(".csv"))
      ),
      explain_card("Для чего нужны диапазоны", "Биологические диапазоны помогают отделить невозможные значения от реальной вариабельности. Например, -990 часто является кодом ошибки, а 85 для internode length может быть ошибкой десятичной точки."),
      card(card_header("Current limits table"), DTOutput("limits_table"))
    )
  ),

  nav_panel(
    nav_label("stars", "Example Results"),
    explain_card("Быстрые примеры результатов", "Если полный анализ занимает время, можно нажимать на карточки ниже и быстро увидеть, какие типы результатов приложение строит после обработки данных."),
    div(class = "quick-grid",
        quick_result_card("example_qc", "clipboard2-check", "QC summary", "Сводка ошибок по признакам: пропуски, диапазоны, выбросы."),
        quick_result_card("example_trait", "bar-chart-line", "Trait plot", "Красивый график распределения выбранного признака."),
        quick_result_card("example_spatial", "grid-3x3-gap", "Field map", "Карта поля по row и column с проблемными зонами."),
        quick_result_card("example_errors", "bug", "Error log", "Таблица записей, которые требуют проверки.")
    ),
    card(card_header(textOutput("example_title")), uiOutput("example_output"))
  ),

  nav_panel(
    nav_label("clipboard2-pulse", "QC Summary"),
    layout_columns(
      metric_box("Records", textOutput("n_records"), "database"),
      metric_box("Plants", textOutput("n_plants"), "flower1"),
      metric_box("Traits", textOutput("n_traits"), "bar-chart"),
      metric_box("Flagged", textOutput("n_flagged"), "exclamation-triangle", "metric-amber"),
      col_widths = c(3, 3, 3, 3)
    ),
    layout_columns(
      card(card_header("Overall error map"), plotOutput("error_overview_plot", height = "420px")),
      card(card_header("Error types"), plotOutput("error_type_plot", height = "420px")),
      col_widths = c(7, 5)
    ),
    explain_card("Как читать QC", "Красные зоны на общей карте показывают признаки, где больше всего проблем: пропуски, нечисловые значения, выход за диапазон, статистические выбросы или биологические конфаундеры."),
    card(card_header("QC table by trait"), DTOutput("qc_table"))
  ),

  nav_panel(
    nav_label("bar-chart-line", "Trait Explorer"),
    layout_sidebar(
      sidebar = sidebar(
        selectInput("trait_outlier", "Trait", choices = NULL),
        checkboxInput("show_removed", "Show removed / flagged observations", TRUE),
        checkboxInput("use_plotly_trait", "Interactive plot", TRUE)
      ),
      explain_card("Что показывает окно", "Здесь можно выбрать конкретный признак и увидеть распределение значений, статусы записей и конкретные растения с ошибками."),
      card(card_header("Trait distribution and flags"), uiOutput("trait_plot_ui")),
      card(card_header("Flagged observations for selected trait"), DTOutput("outlier_table"))
    )
  ),

  nav_panel(
    nav_label("grid-3x3-gap", "Spatial QC"),
    layout_sidebar(
      sidebar = sidebar(
        selectInput("trait_spatial", "Trait", choices = NULL),
        checkboxInput("show_spatial_flags", "Show flags instead of values", FALSE)
      ),
      explain_card("Зачем нужна пространственная карта", "Полевые данные часто содержат позиционный эффект: край ряда, блок, неоднородность почвы или полива. Карта помогает увидеть, где именно концентрируются ошибки и экстремальные значения."),
      card(card_header("Field map by row and column"), plotlyOutput("field_map", height = "560px"))
    )
  ),

  nav_panel(
    nav_label("bug", "Error Report"),
    layout_columns(
      explain_card("Сводка ошибок", "Этот раздел показывает все проблемные записи в одном месте. Его можно использовать как журнал исправлений к статье или как приложение к воспроизводимому workflow."),
      card(card_header("Top problematic traits"), DTOutput("top_errors")),
      col_widths = c(4, 8)
    ),
    card(card_header("Full error log"), DTOutput("full_error_log"))
  ),

  nav_panel(
    nav_label("diagram-3", "Breeding Model"),
    layout_sidebar(
      sidebar = sidebar(
        selectInput("trait_model", "Trait for simple LMM / BLUP", choices = NULL),
        helpText("Model: trait ~ block + (1 | entry). Demonstration only; final publication model may include spatial covariance and genotype data.")
      ),
      explain_card("Для чего это нужно", "После очистки данные становятся пригодными для смешанных линейных моделей. Здесь показан базовый пример ранжирования линий по случайному генетическому эффекту."),
      card(card_header("Entry BLUP ranking"), DTOutput("blup_table")),
      card(card_header("BLUP plot"), plotOutput("blup_plot", height = "420px"))
    )
  ),

  nav_panel(
    nav_label("activity", "DNA + Genomics"),
    layout_columns(
      explain_card("Параллельный геномный модуль", "Для геномной селекции нужен не только чистый фенотип, но и качественный генотип. В проекте CTAB-PK показал более высокий выход и чистоту HMW DNA, чем Qiagen Genomics-Tip."),
      card(card_header("DNA extraction comparison"), DTOutput("dna_table")),
      col_widths = c(4, 8)
    ),
    layout_columns(
      card(card_header("dsDNA yield"), plotOutput("dna_yield_plot", height = "360px")),
      card(card_header("260/230 purity"), plotOutput("dna_purity_plot", height = "360px")),
      col_widths = c(6, 6)
    )
  ),

  nav_panel(
    nav_label("easel", "Presentation Story"),
    card(card_header("Scientific narrative from slides"),
         div(class = "story-grid",
             div(class = "story-card", h4("Problem"), p("Raw field observations + environmental noise can produce inaccurate breeding values.")),
             div(class = "story-card", h4("Barrier"), p("Human entry errors, soil heterogeneity, pests, diseases and arid-climate stress mask genetic potential.")),
             div(class = "story-card", h4("Goal"), p("Build a reproducible R workflow for automated cleaning and model-ready data preparation.")),
             div(class = "story-card", h4("Solution"), p("Digital capture, consolidation, detection, correction and preparation.")),
             div(class = "story-card", h4("Output"), p("Clean phenotype matrix + high-quality HMW DNA for SNP genotyping and GWAS.")),
             div(class = "story-card", h4("Future"), p("Hybrid breeding models, F5-F7 validation and open sesame phenotype database."))
         )
    )
  )
)

server <- function(input, output, session) {
  progress_percent <- reactiveVal(0)
  progress_label <- reactiveVal("Waiting for data")
  run_log_data <- reactiveVal(tibble(step = "Initialized", status = "waiting", percent = 0, time = as.character(Sys.time())))
  custom_limits <- reactiveVal(NULL)
  project_info <- reactiveValues(
    crop = slide_facts$crop,
    location = slide_facts$location,
    scale = slide_facts$scale,
    design = slide_facts$design,
    goal = slide_facts$target
  )
  active_example <- reactiveVal("qc")

  append_log <- function(step, status, percent) {
    old <- run_log_data()
    run_log_data(bind_rows(old, tibble(step = step, status = status, percent = percent, time = as.character(Sys.time()))))
    progress_percent(percent)
    progress_label(paste(step, status))
  }

  input_data <- reactive({
    if (!is.null(input$fieldbook_file)) {
      read_fieldbook(input$fieldbook_file$datapath)
    } else if (isTRUE(input$use_demo) && file.exists(example_path)) {
      read_fieldbook(example_path)
    } else {
      NULL
    }
  })

  observeEvent(input_data(), {
    df <- input_data()
    if (is.null(df)) return()
    cols <- names(df)
    progress_percent(10)
    progress_label("File loaded / columns detected")
    run_log_data(tibble(step = "File upload", status = paste0(nrow(df), " rows, ", ncol(df), " columns"), percent = 10, time = as.character(Sys.time())))
  }, ignoreInit = FALSE)

  output$mapping_ui <- renderUI({
    df <- input_data()
    validate(need(!is.null(df), "Upload a file or use the demo dataset."))
    cols <- names(df)
    tagList(
      selectInput("map_uuid", "UUID / sample ID", choices = choose_vec(cols), selected = detect_col(cols, "uuid")),
      selectInput("map_plant_id", "Plant ID", choices = choose_vec(cols), selected = detect_col(cols, "plant_id")),
      selectInput("map_entry", "Entry / genotype / line", choices = choose_vec(cols), selected = detect_col(cols, "entry")),
      selectInput("map_block", "Block / replicate", choices = choose_vec(cols), selected = detect_col(cols, "block")),
      selectInput("map_row", "Row", choices = choose_vec(cols), selected = detect_col(cols, "row")),
      selectInput("map_column", "Column", choices = choose_vec(cols), selected = detect_col(cols, "column")),
      selectInput("map_trait", "Trait column, for long format", choices = choose_vec(cols), selected = detect_col(cols, "trait")),
      selectInput("map_value", "Value column, for long format", choices = choose_vec(cols), selected = detect_col(cols, "value")),
      selectInput("map_timestamp", "Timestamp", choices = choose_vec(cols), selected = detect_col(cols, "timestamp")),
      selectInput("map_notes", "Notes / comments", choices = choose_vec(cols), selected = detect_col(cols, "notes")),
      selectInput("map_disease", "Disease / pest status", choices = choose_vec(cols), selected = detect_col(cols, "disease"))
    )
  })

  mapping <- reactive({
    list(
      uuid = input$map_uuid, plant_id = input$map_plant_id, entry = input$map_entry, block = input$map_block,
      row = input$map_row, column = input$map_column, trait = input$map_trait, value = input$map_value,
      timestamp = input$map_timestamp, notes = input$map_notes, disease = input$map_disease, person = detect_col(names(input_data()), "person")
    )
  })

  detected_traits <- reactive({
    df <- input_data()
    if (is.null(df)) return(character())
    trait_col <- input$map_trait %||% detect_col(names(df), "trait")
    value_col <- input$map_value %||% detect_col(names(df), "value")
    if (!is_blank(trait_col) && !is_blank(value_col) && trait_col %in% names(df)) {
      sort(unique(janitor::make_clean_names(as.character(df[[trait_col]]))))
    } else {
      meta <- unlist(mapping(), use.names = FALSE)
      sort(setdiff(numericish_cols(df), meta))
    }
  })

  output$trait_selector_ui <- renderUI({
    traits <- detected_traits()
    if (isTRUE(input$auto_select_traits)) return(helpText(paste("Detected", length(traits), "trait(s). All will be analyzed.")))
    checkboxGroupInput("selected_traits", "Select traits", choices = traits, selected = traits)
  })

  selected_traits <- reactive({
    traits <- detected_traits()
    if (isTRUE(input$auto_select_traits)) traits else input$selected_traits
  })

  long_preview <- reactive({
    df <- input_data()
    validate(need(!is.null(df), "Upload data first."))
    prepare_long(df, mapping(), selected_traits())
  })

  observeEvent(long_preview(), {
    long <- long_preview()
    lim <- make_limits_for_traits(long)
    custom_limits(lim)
    updateSelectInput(session, "limit_trait", choices = lim$trait, selected = lim$trait[1])
  }, ignoreInit = FALSE)

  observeEvent(input$limit_trait, {
    lim <- custom_limits()
    req(lim, input$limit_trait)
    row <- lim |> filter(trait == input$limit_trait) |> slice(1)
    updateNumericInput(session, "limit_min", value = row$min_value)
    updateNumericInput(session, "limit_max", value = row$max_value)
    updateTextInput(session, "limit_unit", value = row$unit)
  }, ignoreInit = TRUE)

  observeEvent(input$save_limit, {
    lim <- custom_limits()
    req(lim, input$limit_trait)
    lim <- lim |>
      mutate(
        min_value = if_else(trait == input$limit_trait, as.numeric(input$limit_min), min_value),
        max_value = if_else(trait == input$limit_trait, as.numeric(input$limit_max), max_value),
        unit = if_else(trait == input$limit_trait, input$limit_unit, unit)
      )
    custom_limits(lim)
    showNotification("Range saved for selected trait.", type = "message")
  })

  limits_data <- reactive({
    if (!is.null(input$limits_file)) {
      readr::read_csv(input$limits_file$datapath, show_col_types = FALSE) |> janitor::clean_names()
    } else {
      custom_limits()
    }
  })

  processed <- eventReactive(input$run + input$quick_demo, {
    withProgress(message = "FIELD CLEANER is processing your data", value = 0, {
      append_log("1/7 Upload", "loaded", 15); incProgress(0.15, detail = "Reading CSV and detecting columns")
      Sys.sleep(0.1)
      df <- input_data(); validate(need(!is.null(df), "Upload a CSV file or enable demo data."))

      append_log("2/7 Standardization", "long/wide detection", 30); incProgress(0.15, detail = "Standardizing long/wide structure")
      long <- prepare_long(df, mapping(), selected_traits())

      append_log("3/7 Limits", "biological ranges ready", 45); incProgress(0.15, detail = "Preparing biological limits")
      limits <- limits_data()
      if (is.null(limits)) limits <- make_limits_for_traits(long)

      append_log("4/7 Range QC", "checking plausibility", 60); incProgress(0.15, detail = "Checking missing values and biological ranges")
      cleaned <- clean_long(long, limits, exclude_disease_from_yield = isTRUE(input$exclude_disease))

      append_log("5/7 Summary", "building QC tables", 75); incProgress(0.15, detail = "Building QC summary and error log")
      qc <- summarise_qc(cleaned)
      err_long <- make_error_long(qc)
      error_log <- cleaned |>
        filter(final_status != "kept") |>
        select(any_of(c("uuid", "plant_id", "entry", "block", "row", "column", "trait", "value_raw", "value_num", "final_status", "min_value", "max_value", "range_reason", "iqr_reason", "notes", "disease_status")))

      append_log("6/7 Matrix", "creating model-ready wide table", 90); incProgress(0.15, detail = "Creating clean wide matrix")
      wide <- make_wide(cleaned)
      traits <- sort(unique(cleaned$trait))

      append_log("7/7 Done", "100% complete", 100); incProgress(0.10, detail = "Done")
      list(long = long, cleaned = cleaned, wide = wide, qc = qc, err_long = err_long, error_log = error_log, traits = traits, limits = limits)
    })
  }, ignoreInit = FALSE)

  observe({
    p <- processed()
    traits <- p$traits
    selected <- if ("capsule_count" %in% traits) "capsule_count" else traits[1]
    updateSelectInput(session, "trait_outlier", choices = traits, selected = selected)
    updateSelectInput(session, "trait_spatial", choices = traits, selected = selected)
    updateSelectInput(session, "trait_model", choices = traits, selected = selected)
  })


  observeEvent(input$apply_project, {
    project_info$crop <- input$project_crop
    project_info$location <- input$project_location
    project_info$scale <- input$project_scale
    project_info$design <- input$project_design
    project_info$goal <- input$project_goal
    showNotification("Project card updated.", type = "message")
  })

  output$project_metric_boxes <- renderUI({
    layout_columns(
      metric_box("Crop", project_info$crop, "flower1"),
      metric_box("Location", project_info$location, "geo-alt"),
      metric_box("Scale", project_info$scale, "database"),
      metric_box("Design", project_info$design, "grid-3x3-gap"),
      col_widths = c(3, 3, 3, 3)
    )
  })

  output$project_preview <- renderUI({
    div(class = "preview-card",
        h3("FIELD CLEANER project card"),
        tags$dl(
          tags$dt("Crop"), tags$dd(project_info$crop),
          tags$dt("Location"), tags$dd(project_info$location),
          tags$dt("Scale"), tags$dd(project_info$scale),
          tags$dt("Design"), tags$dd(project_info$design),
          tags$dt("Goal"), tags$dd(project_info$goal)
        ))
  })

  output$data_status_cards <- renderUI({
    df <- input_data()
    if (is.null(df)) {
      return(div(class = "status-grid", div(class = "status-card warn", bsicons::bs_icon("cloud-upload"), strong("No data loaded"), span("Upload CSV or use demo dataset"))))
    }
    cols <- names(df)
    long_ok <- !is_blank(detect_col(cols, "trait")) && !is_blank(detect_col(cols, "value"))
    div(class = "status-grid",
        div(class = "status-card", bsicons::bs_icon("table"), strong(format(nrow(df), big.mark = ",")), span("rows loaded")),
        div(class = "status-card", bsicons::bs_icon("columns-gap"), strong(ncol(df)), span("columns")),
        div(class = "status-card", bsicons::bs_icon(if (long_ok) "check-circle" else "layout-three-columns"), strong(if (long_ok) "Long format" else "Wide / auto"), span("detected structure"))
    )
  })

  output$processing_status <- renderText({
    paste0("Status: ", progress_percent(), "% - ", progress_label())
  })

  observeEvent(input$example_qc, active_example("qc"))
  observeEvent(input$example_trait, active_example("trait"))
  observeEvent(input$example_spatial, active_example("spatial"))
  observeEvent(input$example_errors, active_example("errors"))

  output$example_title <- renderText({
    switch(active_example(), qc = "Example: QC summary", trait = "Example: trait distribution", spatial = "Example: spatial field map", errors = "Example: error log")
  })

  output$example_output <- renderUI({
    switch(active_example(),
      qc = plotOutput("example_qc_plot", height = "430px"),
      trait = plotOutput("example_trait_plot", height = "430px"),
      spatial = plotOutput("example_spatial_plot", height = "430px"),
      errors = DTOutput("example_error_table")
    )
  })

  output$example_qc_plot <- renderPlot({
    processed()$err_long |>
      group_by(trait) |> summarise(flagged = sum(count), .groups = "drop") |>
      arrange(desc(flagged)) |> slice_head(n = 12) |>
      ggplot(aes(x = reorder(trait, flagged), y = flagged)) + geom_col() + coord_flip() +
      labs(x = NULL, y = "Flagged records") + theme_minimal(base_size = 14)
  })

  output$example_trait_plot <- renderPlot({
    tr <- processed()$traits[1]
    if ("plant_height" %in% processed()$traits) tr <- "plant_height"
    processed()$cleaned |> filter(trait == tr) |>
      ggplot(aes(x = final_status, y = value_num)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, alpha = .55) +
      labs(x = "Status", y = tr) + theme_minimal(base_size = 14)
  })

  output$example_spatial_plot <- renderPlot({
    tr <- processed()$traits[1]
    clean <- processed()$cleaned |> filter(trait == tr)
    validate(need(all(c("row", "column") %in% names(clean)), "row and column are required."))
    clean |> group_by(row, column) |> summarise(flagged = sum(final_status != "kept"), .groups = "drop") |>
      ggplot(aes(x = column, y = row, fill = flagged)) + geom_tile(color = "white") + scale_y_reverse() +
      labs(x = "Column", y = "Row", fill = "Flags") + theme_minimal(base_size = 14)
  })

  output$example_error_table <- renderDT({
    datatable(head(processed()$error_log, 30), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10))
  })

  output$progress_bar_ui <- renderUI({
    pct <- progress_percent()
    label <- progress_label()
    div(
      div(class = "progress-label", paste0(pct, "% - ", label)),
      div(class = "progress", div(class = "progress-bar", role = "progressbar", style = paste0("width:", pct, "%"), paste0(pct, "%")))
    )
  })

  output$run_log <- renderDT({
    datatable(run_log_data() |> arrange(desc(time)), rownames = FALSE, options = list(dom = "t", pageLength = 8, scrollX = TRUE))
  })

  output$raw_preview <- renderDT({
    df <- input_data()
    validate(need(!is.null(df), "Upload a CSV file or keep demo data enabled."))
    datatable(head(df, 50), options = list(scrollX = TRUE, pageLength = 8, deferRender = TRUE))
  })

  output$detected_columns <- renderDT({
    df <- input_data(); validate(need(!is.null(df), "Upload data first."))
    fields <- names(synonyms)
    tbl <- tibble(field = fields, detected_column = purrr::map_chr(fields, ~ detect_col(names(df), .x))) |>
      mutate(status = if_else(detected_column == "", "not detected", "detected"))
    datatable(tbl, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  output$trait_preview <- renderDT({
    traits <- detected_traits()
    tibble(trait = traits, selected_for_analysis = traits %in% selected_traits()) |>
      datatable(rownames = FALSE, options = list(pageLength = 15))
  })

  output$limits_table <- renderDT({
    datatable(custom_limits(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15))
  })

  output$n_records <- renderText(nrow(processed()$cleaned))
  output$n_plants <- renderText(length(unique(processed()$cleaned$uuid)))
  output$n_traits <- renderText(length(unique(processed()$cleaned$trait)))
  output$n_flagged <- renderText(sum(processed()$cleaned$final_status != "kept", na.rm = TRUE))

  output$qc_table <- renderDT({
    datatable(processed()$qc, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15, deferRender = TRUE))
  })

  output$error_overview_plot <- renderPlot({
    processed()$err_long |>
      ggplot(aes(x = error_type, y = reorder(trait, count), fill = count)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = if_else(count > 0, as.character(count), "")), size = 3) +
      labs(x = "Error type", y = "Trait", fill = "Count") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())
  })

  output$error_type_plot <- renderPlot({
    processed()$err_long |>
      group_by(error_type) |>
      summarise(count = sum(count), .groups = "drop") |>
      ggplot(aes(x = reorder(error_type, count), y = count)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = "Number of records") +
      theme_minimal(base_size = 13)
  })

  selected_trait_data <- reactive({
    req(input$trait_outlier)
    processed()$cleaned |> filter(trait == input$trait_outlier)
  })

  output$trait_plot_ui <- renderUI({
    if (isTRUE(input$use_plotly_trait)) plotlyOutput("outlier_plot_interactive", height = "480px") else plotOutput("outlier_plot_static", height = "440px")
  })

  make_trait_plot <- reactive({
    df <- selected_trait_data()
    if (!isTRUE(input$show_removed)) df <- df |> filter(final_status == "kept")
    ggplot(df, aes(x = final_status, y = value_num)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(text = paste("Plant:", uuid, "<br>Entry:", entry, "<br>Value:", value_num, "<br>Status:", final_status)), width = 0.15, alpha = 0.65) +
      labs(x = "Record status", y = input$trait_outlier) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
  })

  output$outlier_plot_static <- renderPlot({ make_trait_plot() })
  output$outlier_plot_interactive <- renderPlotly({ ggplotly(make_trait_plot(), tooltip = "text") })

  output$outlier_table <- renderDT({
    tbl <- selected_trait_data() |>
      filter(final_status != "kept") |>
      select(any_of(c("uuid", "plant_id", "entry", "block", "row", "column", "trait", "value_raw", "value_num", "final_status", "min_value", "max_value", "range_reason", "iqr_reason", "notes", "disease_status")))
    datatable(tbl, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10, deferRender = TRUE))
  })

  output$field_map <- renderPlotly({
    req(input$trait_spatial)
    wide <- processed()$wide
    clean <- processed()$cleaned |> filter(trait == input$trait_spatial)
    validate(need(all(c("row", "column") %in% names(wide)), "row and column are required for spatial map."))

    if (isTRUE(input$show_spatial_flags)) {
      map_df <- clean |>
        group_by(uuid, entry, row, column) |>
        summarise(flagged = any(final_status != "kept"), status = paste(unique(final_status), collapse = "; "), .groups = "drop")
      p <- map_df |>
        ggplot(aes(x = column, y = row, fill = flagged, text = paste("Plant:", uuid, "<br>Entry:", entry, "<br>Status:", status))) +
        geom_tile(color = "white") + scale_y_reverse() + labs(x = "Column", y = "Row", fill = "Flagged") + theme_minimal(base_size = 13)
    } else {
      validate(need(input$trait_spatial %in% names(wide), "Selected trait not available in wide matrix."))
      p <- wide |>
        ggplot(aes(x = column, y = row, fill = .data[[input$trait_spatial]], text = paste("Plant:", uuid, "<br>Entry:", entry, "<br>Value:", .data[[input$trait_spatial]]))) +
        geom_tile(color = "white") + scale_y_reverse() + labs(x = "Column", y = "Row", fill = input$trait_spatial) + theme_minimal(base_size = 13)
    }
    ggplotly(p, tooltip = "text")
  })

  output$top_errors <- renderDT({
    processed()$qc |>
      select(trait, records, flagged_total, flagged_pct, range_errors, statistical_outliers, biological_confounders) |>
      arrange(desc(flagged_pct)) |>
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 10))
  })

  output$full_error_log <- renderDT({
    datatable(processed()$error_log, rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 20, deferRender = TRUE))
  })

  blups <- reactive({
    req(input$trait_model)
    fit_simple_lmm(processed()$wide, input$trait_model)
  })

  output$blup_table <- renderDT({ datatable(blups(), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10)) })

  output$blup_plot <- renderPlot({
    b <- blups()
    validate(need(all(c("entry", "blup") %in% names(b)), "Model output is not available for this trait."))
    b |>
      slice_head(n = 30) |>
      ggplot(aes(x = reorder(entry, blup), y = blup)) +
      geom_col() + coord_flip() + labs(x = NULL, y = "Random genetic effect / BLUP") + theme_minimal(base_size = 13)
  })

  output$dna_table <- renderDT({ datatable(dna_table, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5)) })

  output$dna_yield_plot <- renderPlot({
    dna_table |>
      ggplot(aes(x = protocol, y = dsdna_ng_ul)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 3) +
      labs(x = NULL, y = "dsDNA yield, ng/ul") + theme_minimal(base_size = 13)
  })

  output$dna_purity_plot <- renderPlot({
    dna_table |>
      ggplot(aes(x = protocol, y = ratio_260_230)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.12, size = 3) +
      geom_hline(yintercept = 2.0, linetype = "dashed") + labs(x = NULL, y = "260/230 purity ratio") + theme_minimal(base_size = 13)
  })

  output$download_clean <- downloadHandler(
    filename = function() paste0("field_cleaner_clean_matrix_", Sys.Date(), ".csv"),
    content = function(file) readr::write_csv(processed()$wide, file)
  )

  output$download_qc <- downloadHandler(
    filename = function() paste0("field_cleaner_qc_summary_", Sys.Date(), ".csv"),
    content = function(file) readr::write_csv(processed()$qc, file)
  )
}

shinyApp(ui, server)
