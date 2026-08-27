#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

require_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set environment variable ", name)
  value
}

input <- readRDS(file.path(require_env("OMICSPRED_RESULTS_DIR"), "interval_models_by_cohort.rds"))
output_dir <- require_env("OMICSPRED_FIGURE_DIR")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
labels <- c(MCPS_ALL = "MCPS", MEC_CN = "MEC Chinese", MEC_IN = "MEC Indian", MEC_MA = "MEC Malay", UKB_EUR = "UKB European")

plot_data <- input %>%
  mutate(cohort = unname(labels[str_trim(as.character(cohort_name))])) %>%
  filter(!is.na(cohort))
order <- plot_data %>% group_by(cohort) %>% summarise(value = median(R2, na.rm = TRUE), .groups = "drop") %>% arrange(value) %>% pull(cohort)
plot_data$cohort <- factor(plot_data$cohort, levels = order)
fills <- setNames(rep("#E5E7EBCC", length(order)), order)
fills["MCPS"] <- "#D62828CC"

p <- ggplot(plot_data, aes(cohort, R2, fill = cohort)) +
  geom_violin(alpha = 0.8, linewidth = 0.3, trim = TRUE) +
  geom_boxplot(fill = "white", color = "black", width = 0.1, outlier.shape = 21, outlier.fill = "white", outlier.size = 1) +
  scale_fill_manual(values = fills) +
  scale_y_continuous(limits = c(0, 0.22), breaks = seq(0, 0.20, 0.025)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Cohort used for performance evaluation", y = expression(R^2)) +
  coord_cartesian(ylim = c(0, 0.20))

ggsave(file.path(output_dir, "figure_1.pdf"), p, width = 9, height = 6, device = cairo_pdf)
