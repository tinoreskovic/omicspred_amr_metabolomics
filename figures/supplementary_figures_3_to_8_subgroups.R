#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
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

make_plot <- function(data, group_column, levels, labels, x_label, output, fills = NULL) {
  plot_data <- data %>% mutate(group = factor(.data[[group_column]], levels = levels, labels = labels)) %>% filter(!is.na(group), !is.na(R2))
  if (is.null(fills)) fills <- setNames(rep("#E5E7EBCC", length(labels)), labels)
  p <- ggplot(plot_data, aes(group, R2, fill = group)) +
    geom_violin(trim = TRUE, alpha = 0.8, color = "black", linewidth = 0.3) +
    geom_point(position = position_jitter(width = 0.06, seed = 1), shape = 21, size = 1.15, alpha = 0.45) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.45) +
    scale_fill_manual(values = fills) + scale_y_continuous(limits = c(0, 0.20), breaks = seq(0, 0.20, 0.05)) +
    theme_bw() + theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = x_label, y = expression(R^2))
  ggsave(file.path(output_dir, output), p, width = ifelse(group_column == "HealthStat", 8, 7), height = 5, device = cairo_pdf)
}

common_specifications <- list(
  list("age", "AgeGroup", c("35~49", "50~64", "65 and older"), c("35–49", "50–64", "≥65"), "Age group"),
  list("health", "HealthStat", c("No diabetes", "Diabetes", "No disease", "Any disease"), c("No diabetes", "Diabetes", "No disease", "Any disease"), "Health status")
)

figure <- 3L
for (model in c("interval", "mcps")) {
  prefix <- if (model == "interval") "interval_models_mcps" else "mcps_models_heldout"
  first_specification <- if (model == "interval") {
    list("sex_ancestry", "subgroup", c("AMR_Female", "EUR_Female", "AMR_Male", "EUR_Male"), c("IAM female", "EUR female", "IAM male", "EUR male"), "Sex and genetic ancestry subset")
  } else {
    list("sex", "Sex", c("Female", "Male"), c("Female", "Male"), "Sex")
  }
  specifications <- c(list(first_specification), common_specifications)
  for (spec in specifications) {
    input <- readRDS(file.path(results, paste0(prefix, "_by_", spec[[1]], ".rds")))
    make_plot(input, spec[[2]], spec[[3]], spec[[4]], spec[[5]], paste0("supplementary_figure_", figure, ".pdf"))
    figure <- figure + 1L
  }
}
