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
heldout <- as.data.table(readRDS(file.path(results, "mcps_models_heldout_overall.rds")))[, .(PRS_base = PRS_name, MCPS_R2 = R2)]
ukb <- fread(file.path(results, "ukb_amr_model_performance.tsv"))[model == "MCPS-trained" & subgroup == "All", .(PRS_base = sub("_MCPS_80_20_fixed$", "", PGS), UKB_R2 = R2_pearson)]
plot_data <- merge(heldout, ukb, by = "PRS_base")

p <- ggplot(plot_data, aes(MCPS_R2, UKB_R2)) +
  geom_point(alpha = 0.7, size = 2, color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dashed") +
  geom_smooth(method = "lm", formula = y ~ 0 + x, color = "red", se = FALSE, linewidth = 0.7) +
  scale_x_continuous(limits = c(0, 0.20), breaks = seq(0, 0.20, 0.05)) +
  scale_y_continuous(limits = c(0, 0.20), breaks = seq(0, 0.20, 0.05)) +
  coord_equal(expand = FALSE) + theme_bw() +
  labs(x = expression("MCPS-trained models' " * R^2 * " in withheld MCPS"), y = expression("MCPS-trained models' " * R^2 * " in UKB AMR"))
ggsave(file.path(output_dir, "supplementary_figure_2.pdf"), p, width = 8.5, height = 5, device = cairo_pdf)
