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

results <- require_env("OMICSPRED_RESULTS_DIR")
output_dir <- require_env("OMICSPRED_FIGURE_DIR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
interval <- readRDS(file.path(results, "interval_models_mcps_by_ancestry.rds"))[ancestry %in% c("MCPS_AMR", "MCPS_EUR")]
mcps <- readRDS(file.path(results, "mcps_models_heldout_by_ancestry.rds"))[ancestry %in% c("MCPS_AMR", "MCPS_EUR")]
interval$model <- "INTERVAL-trained models"
mcps$model <- "MCPS-trained models"
plot_data <- rbindlist(list(interval, mcps), fill = TRUE)
plot_data[, ancestry_label := factor(ancestry, c("MCPS_AMR", "MCPS_EUR"), c("MCPS IAM >= 70%", "MCPS EUR >= 70%"))]
plot_data[, group := factor(paste(model, ancestry_label, sep = "
in "), levels = c("INTERVAL-trained models
in MCPS IAM >= 70%", "INTERVAL-trained models
in MCPS EUR >= 70%", "MCPS-trained models
in MCPS IAM >= 70%", "MCPS-trained models
in MCPS EUR >= 70%"))]

p <- ggplot(plot_data, aes(group, R2, fill = ancestry_label)) +
  geom_violin(trim = TRUE, alpha = 0.6, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.1, outlier.shape = 21, outlier.fill = "white", fill = "white", color = "black") +
  scale_fill_manual(values = c("MCPS IAM >= 70%" = "#00BFC4", "MCPS EUR >= 70%" = "#F8766D")) +
  scale_y_continuous(limits = c(0, 0.20), breaks = seq(0, 0.20, 0.05)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 9)) +
  labs(x = "Model training cohort and MCPS evaluation subset", y = expression(R^2))

ggsave(file.path(output_dir, "figure_2.pdf"), p, width = 11, height = 6, device = cairo_pdf)
