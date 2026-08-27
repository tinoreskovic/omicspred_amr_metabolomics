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

scatter <- function(data, x, y, x_label, y_label, output) {
  limits <- c(0, 0.20)
  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.7, size = 2, color = "black") +
    geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
    geom_smooth(method = "lm", formula = y ~ 0 + x, color = "red", se = FALSE, linewidth = 0.7) +
    scale_x_continuous(limits = limits, breaks = seq(0, 0.20, 0.05)) +
    scale_y_continuous(limits = limits, breaks = seq(0, 0.20, 0.05)) +
    coord_equal(expand = FALSE) + theme_bw() + labs(x = x_label, y = y_label)
  ggsave(file.path(output_dir, output), p, width = 8.5, height = 5, device = cairo_pdf)
}

interval_mcps <- as.data.table(readRDS(file.path(results, "interval_models_mcps_heldout_overall.rds")))
mcps_mcps <- as.data.table(readRDS(file.path(results, "mcps_models_heldout_overall.rds")))
panel_a <- merge(interval_mcps[, .(NMR_name, interval_R2 = old_test_R2)], mcps_mcps[, .(NMR_name, mcps_R2 = R2)], by = "NMR_name")
scatter(panel_a, "interval_R2", "mcps_R2", expression("INTERVAL-trained models' " * R^2 * " in MCPS"), expression("MCPS-trained models' " * R^2 * " in MCPS"), "figure_3a.pdf")

ukb <- fread(file.path(results, "ukb_amr_model_performance.tsv"))[subgroup == "All"]
panel_b <- dcast(ukb, trait ~ model, value.var = "R2_pearson")
setnames(panel_b, c("INTERVAL-trained", "MCPS-trained"), c("interval_R2", "mcps_R2"))
scatter(panel_b, "interval_R2", "mcps_R2", expression("INTERVAL-trained models' " * R^2 * " in UKB AMR"), expression("MCPS-trained models' " * R^2 * " in UKB AMR"), "figure_3b.pdf")
