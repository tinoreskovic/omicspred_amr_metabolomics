#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(grid)
  library(ggrepel)
})

require_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set environment variable ", name)
  value
}

results <- require_env("OMICSPRED_RESULTS_DIR")
output_dir <- require_env("OMICSPRED_FIGURE_DIR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
primary <- c("Type 2 diabetes", "Ischaemic heart disease", "Chronic kidney disease")
data <- read_csv(file.path(results, "aou_disease_associations.csv"), show_col_types = FALSE) %>%
  mutate(Disease_Name = recode(Disease_Name,
    "Type 2 Diabetes" = "Type 2 diabetes",
    "Ischemic Heart Disease" = "Ischaemic heart disease",
    "Chronic Kidney Disease" = "Chronic kidney disease"
  )) %>%
  filter(model == "Logistic", Disease_Name %in% primary) %>%
  group_by(Disease_Name, score_set) %>% mutate(fdr = p.adjust(p, method = "fdr")) %>% ungroup() %>%
  mutate(score_set = recode(score_set, "INTERVAL-trained" = "INTERVAL", "MCPS-trained" = "MCPS"), beta = log(estimate), low = log(L95), high = log(U95)) %>%
  select(Disease_Name, Biomarker.Name, n_cases, n_controls, score_set, beta, low, high, fdr) %>%
  pivot_wider(names_from = score_set, values_from = c(beta, low, high, fdr))

low_limit <- log(0.85)
high_limit <- log(1.15)
breaks <- c(0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15)
data <- data %>% mutate(
  x = pmin(pmax(beta_INTERVAL, low_limit), high_limit),
  y = pmin(pmax(beta_MCPS, low_limit), high_limit),
  x_low = pmin(pmax(low_INTERVAL, low_limit), high_limit),
  x_high = pmin(pmax(high_INTERVAL, low_limit), high_limit),
  y_low = pmin(pmax(low_MCPS, low_limit), high_limit),
  y_high = pmin(pmax(high_MCPS, low_limit), high_limit),
  category = case_when(
    fdr_INTERVAL < 0.05 & fdr_MCPS < 0.05 ~ "FDR-significant for both",
    fdr_INTERVAL < 0.05 ~ "FDR-significant only for INTERVAL-trained model",
    fdr_MCPS < 0.05 ~ "FDR-significant only for MCPS-trained model",
    TRUE ~ "Not significant"
  ),
  outlier = beta_INTERVAL < low_limit | beta_INTERVAL > high_limit | beta_MCPS < low_limit | beta_MCPS > high_limit
)
colors <- c("FDR-significant for both" = "#E69F00", "FDR-significant only for INTERVAL-trained model" = "#56B4E9", "FDR-significant only for MCPS-trained model" = "#009E73", "Not significant" = "#DFDFDF")

for (disease in primary) {
  panel <- filter(data, Disease_Name == disease)
  p <- ggplot(panel) +
    geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "gray60", linetype = "dashed") +
    geom_smooth(aes(beta_INTERVAL, beta_MCPS), method = "lm", formula = y ~ x, color = "red", se = FALSE, linewidth = 1) +
    geom_segment(aes(x = x_low, xend = x_high, y = y, yend = y, color = category), alpha = 0.4) +
    geom_segment(aes(x = x, xend = x, y = y_low, yend = y_high, color = category), alpha = 0.4) +
    geom_point(aes(x, y, color = category), alpha = 0.8, size = 2.5) +
    geom_text_repel(data = filter(panel, outlier), aes(x, y, label = Biomarker.Name), size = 3.5, fontface = "italic", color = "black", max.overlaps = Inf) +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = log(breaks), labels = breaks) + scale_y_continuous(breaks = log(breaks), labels = breaks) +
    coord_cartesian(xlim = c(low_limit, high_limit), ylim = c(low_limit, high_limit)) +
    theme_minimal(base_size = 20) + theme(legend.position = "none", axis.title = element_text(face = "bold"), panel.grid.minor = element_blank()) +
    labs(title = disease, x = "OR (per s.d. of INTERVAL-based score)", y = "OR (per s.d. of MCPS-based score)")
  stem <- match(disease, primary)
  ggsave(file.path(output_dir, paste0("figure_5", letters[stem], ".pdf")), p, width = 8, height = 8, device = cairo_pdf)
}

counts <- data %>% group_by(Disease_Name) %>% summarise(
  `Either model` = sum(fdr_INTERVAL < 0.05 | fdr_MCPS < 0.05),
  `INTERVAL-trained model` = sum(fdr_INTERVAL < 0.05),
  `MCPS-trained model` = sum(fdr_MCPS < 0.05), .groups = "drop"
) %>% pivot_longer(-Disease_Name, names_to = "category", values_to = "count")

p_counts <- ggplot(counts, aes(count, Disease_Name, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Either model" = "#E69F00", "INTERVAL-trained model" = "#56B4E9", "MCPS-trained model" = "#009E73")) +
  theme_minimal(base_size = 18) + theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "Number of metabolomic traits with FDR-significant associations", y = NULL)
ggsave(file.path(output_dir, "figure_5d.pdf"), p_counts, width = 14, height = 7, device = cairo_pdf)
