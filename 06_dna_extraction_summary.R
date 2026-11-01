# 06_dna_extraction_summary.R
# Summarize HMW DNA extraction comparison

source("R/00_config.R")
load_required_packages()

dna <- tibble::tribble(
  ~protocol, ~sample, ~dsdna_ng_ul, ~ratio_260_280, ~ratio_260_230, ~integrity,
  "CTAB-PK", "Leaf 1", 250.3, 1.83, 2.01, "HMW > 50 kb",
  "CTAB-PK", "Leaf 2", 215.1, 1.79, 2.15, "HMW > 50 kb",
  "Qiagen Genomics-Tip", "Leaf 1", 38.7, 1.45, 1.05, "Medium, partial fragmentation",
  "Qiagen Genomics-Tip", "Leaf 2", 42.4, 1.40, 0.98, "Medium, partial fragmentation"
)

summary <- dna |>
  group_by(protocol) |>
  summarise(
    mean_dsdna_ng_ul = mean(dsdna_ng_ul),
    mean_260_280 = mean(ratio_260_280),
    mean_260_230 = mean(ratio_260_230),
    n = n(),
    .groups = "drop"
  )

readr::write_csv(dna, file.path(project_paths$tables_dir, "dna_extraction_raw.csv"))
readr::write_csv(summary, file.path(project_paths$tables_dir, "dna_extraction_summary.csv"))

p <- dna |>
  ggplot(aes(x = protocol, y = dsdna_ng_ul)) +
  geom_col() +
  labs(x = "Protocol", y = "dsDNA concentration, ng/ul", title = "HMW DNA extraction yield") +
  theme_minimal()

ggsave(file.path(project_paths$figures_dir, "dna_extraction_yield.png"), p, width = 7, height = 5, dpi = 300)
message("Saved DNA extraction summary and figure.")
