#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readxl)
})

require_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set environment variable ", name)
  value
}

results <- require_env("OMICSPRED_RESULTS_DIR")
output_dir <- require_env("OMICSPRED_FIGURE_DIR")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
metadata_file <- normalizePath(file.path(dirname(script_file), "..", "metadata", "omicspred_validation.xlsx"), mustWork = TRUE)
gwas <- fread(require_env("MCPS_GWAS_SUMMARY_FILE")) %>% transmute(trait = sub(" \\(.*$", "", `Phenotype description`), h2 = as.numeric(sub(" .*", "", h2)))
performance <- fread(file.path(results, "mcps_models_heldout_overall.csv"))
mapping <- read_excel(metadata_file, sheet = "Table S2") %>% transmute(PRS_name, trait = Biomarker.Name)
plot_data <- performance %>% inner_join(mapping, by = "PRS_name") %>% inner_join(gwas, by = "trait")
limit <- ceiling(max(c(plot_data$h2, plot_data$R2), na.rm = TRUE) * 10) / 10 + 0.05

p <- ggplot(plot_data, aes(h2, R2)) +
  geom_point(alpha = 0.7, size = 2.5, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, color = "red", linewidth = 0.8) +
  scale_x_continuous(limits = c(0, limit)) + scale_y_continuous(limits = c(0, limit)) +
  coord_equal(expand = FALSE) + theme_bw() +
  labs(x = expression(h^2 * " in MCPS"), y = expression("MCPS-trained models' " * R^2 * " in withheld MCPS"))
ggsave(file.path(output_dir, "supplementary_figure_1.pdf"), p, width = 7, height = 7, device = cairo_pdf)
