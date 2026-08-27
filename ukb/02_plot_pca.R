#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

require_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set environment variable ", name)
  value
}

dat <- fread(require_env("UKB_ANCESTRY_FILE"))
output_dir <- require_env("OMICSPRED_FIGURE_DIR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
reference <- dat[sampleset == "reference" & !is.na(SuperPop)]
ukb <- dat[sampleset == "UKB-NMR" & RF_P_AMR >= 0.90]
reference[, PC3 := -PC3]
ukb[, PC3 := -PC3]
palette <- c(AFR = "#7C3AED", AMR = "#EF4444", CSA = "#F59E0B", EAS = "#10B981", EUR = "#3B82F6", MID = "#8B5E3C")

p <- ggplot() +
  geom_point(data = ukb, aes(PC3, PC1), color = "grey40", size = 2.5, alpha = 0.4) +
  geom_point(data = reference, aes(PC3, PC1, color = SuperPop), shape = 1, size = 2.1, linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = palette) +
  scale_x_continuous(breaks = seq(-100, 0, 25)) +
  scale_y_continuous(breaks = seq(-40, 40, 20)) +
  coord_cartesian(xlim = c(-100, 25), ylim = c(-60, 40)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey92")
  ) +
  labs(x = "PC3", y = "PC1", color = "Populations")

ggsave(file.path(output_dir, "figure_3d_ukb_pca.pdf"), p, width = 8.5, height = 6, dpi = 600, device = cairo_pdf)
