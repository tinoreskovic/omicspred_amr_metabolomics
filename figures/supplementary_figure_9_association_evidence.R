#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})

require_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set environment variable ", name)
  value
}

results <- require_env("OMICSPRED_RESULTS_DIR")
output_dir <- require_env("OMICSPRED_FIGURE_DIR")
primary <- c("Type 2 diabetes", "Ischaemic heart disease", "Chronic kidney disease")
data <- read_csv(file.path(results, "aou_disease_associations.csv"), show_col_types = FALSE) %>%
  mutate(Disease_Name = recode(Disease_Name, "Type 2 Diabetes" = "Type 2 diabetes", "Ischemic Heart Disease" = "Ischaemic heart disease", "Chronic Kidney Disease" = "Chronic kidney disease")) %>%
  filter(model == "Logistic", Disease_Name %in% primary) %>%
  mutate(score_set = recode(score_set, "INTERVAL-trained" = "INTERVAL", "MCPS-trained" = "MCPS")) %>%
  select(Disease_Name, Biomarker.Name, score_set, p) %>%
  pivot_wider(names_from = score_set, values_from = p) %>%
  mutate(x = -log10(pmax(INTERVAL, .Machine$double.xmin)), y = -log10(pmax(MCPS, .Machine$double.xmin)), evidence = ifelse(y > x, "MCPS-trained model", ifelse(y < x, "INTERVAL-trained model", "Equal")))
colors <- c("MCPS-trained model" = "#009E73", "INTERVAL-trained model" = "#56B4E9", Equal = "#BDBDBD")

make_panel <- function(disease) {
  panel <- filter(data, Disease_Name == disease)
  limit <- 1.04 * max(c(panel$x, panel$y))
  label <- panel %>% filter(pmax(x, y) >= 10) %>% slice_max(pmax(x, y), n = 1, with_ties = FALSE)
  ggplot(panel, aes(x, y, fill = evidence)) +
    geom_abline(slope = 1, intercept = 0, color = "gray45", linetype = "dashed") +
    geom_point(shape = 21, color = "white", alpha = 0.78, size = 3) +
    geom_text_repel(data = label, aes(label = Biomarker.Name), color = "black", size = 3.2, fontface = "italic", max.overlaps = Inf, show.legend = FALSE) +
    scale_fill_manual(values = colors) + scale_x_sqrt() + scale_y_sqrt() +
    coord_equal(xlim = c(0, limit), ylim = c(0, limit), expand = FALSE) +
    theme_minimal(base_size = 17) + theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank()) +
    labs(title = disease, x = expression(-log[10](italic(P)) ~ "for INTERVAL-trained model"), y = expression(-log[10](italic(P)) ~ "for MCPS-trained model"))
}

combined <- plot_grid(plotlist = lapply(primary, make_panel), labels = c("a)", "b)", "c)"), nrow = 1)
ggsave(file.path(output_dir, "supplementary_figure_9.pdf"), combined, width = 22, height = 8.2, device = cairo_pdf)
