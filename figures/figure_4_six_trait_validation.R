#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
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
trait_labels <- c(HDL_C = "HDL", LDL_C = "LDL", VLDL_C = "VLDL", Total_TG_serum = "Triglycerides", Creatinine = "Creatinine", Glucose = "Glucose")
colors <- c("INTERVAL-trained" = "#56B4E9", "MCPS-trained" = "#009E73")

make_plot <- function(data, output) {
  order <- data %>% group_by(trait) %>% summarise(max_r2 = max(R2, na.rm = TRUE), .groups = "drop") %>% arrange(max_r2) %>% pull(trait)
  data$trait <- factor(data$trait, levels = order)
  p <- ggplot(data, aes(trait, R2, fill = model)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72, color = "black", linewidth = 0.25) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02), expand = expansion(mult = c(0, 0.03))) +
    theme_minimal(base_size = 15) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = NULL, y = expression(R^2), fill = NULL)
  ggsave(file.path(output_dir, output), p, width = 5, height = 7, device = cairo_pdf)
}

interval <- as.data.table(readRDS(file.path(results, "interval_models_mcps_heldout_overall.rds")))[, .(NMR_name, R2 = old_test_R2, model = "INTERVAL-trained")]
mcps <- as.data.table(readRDS(file.path(results, "mcps_models_heldout_overall.rds")))[, .(NMR_name, R2, model = "MCPS-trained")]
mcps_data <- rbindlist(list(interval, mcps))
mcps_data[, trait := unname(trait_labels[NMR_name])]
make_plot(mcps_data[!is.na(trait)], "figure_4_mcps.pdf")

aou <- fread(file.path(results, "aou_six_trait_model_performance.csv"))[subgroup == "All"]
aou[, model := fifelse(score_type == "OmicsPred", "INTERVAL-trained", "MCPS-trained")]
make_plot(aou[, .(trait, R2 = R2_pearson_RINT, model)], "figure_4_aou.pdf")
