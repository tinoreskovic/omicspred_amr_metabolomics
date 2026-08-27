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

input_file <- require_env("MCPS_ANCESTRY_FILE")
output_dir <- require_env("OMICSPRED_FIGURE_DIR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(input_file)
analysis <- readRDS(file.path(require_env("OMICSPRED_MCPS_ROOT"), "secure_intermediates", "mcps_interval_validation_dataset.rds"))
dat <- dat[sampleset != "MCPS" | as.character(FID) %chin% as.character(analysis$IID)]
reference <- dat[sampleset == "reference" & !is.na(SuperPop)]
mcps <- dat[sampleset == "MCPS" & RF_P_AMR >= 0.90]
palette <- c(AFR = "#7C3AED", AMR = "#EF4444", CSA = "#F59E0B", EAS = "#10B981", EUR = "#3B82F6", MID = "#8B5E3C")

p <- ggplot() +
  geom_point(data = mcps, aes(PC3, PC1), color = "black", size = 1.5, alpha = 0.4) +
  geom_point(data = reference, aes(PC3, PC1, color = SuperPop), shape = 1, size = 2.1, linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = palette) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3)
  ) +
  labs(x = "PC3", y = "PC1", color = "Populations")

ggsave(file.path(output_dir, "figure_3c_mcps_pca.pdf"), p, width = 8.5, height = 6, dpi = 600, device = cairo_pdf)
