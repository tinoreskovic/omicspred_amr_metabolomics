#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

require_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set environment variable ", name)
  value
}

root <- require_env("OMICSPRED_MCPS_ROOT")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
metadata_file <- normalizePath(file.path(dirname(script_file), "..", "..", "metadata", "omicspred_validation.xlsx"), mustWork = TRUE)
input_file <- file.path(root, "secure_intermediates", "mcps_interval_validation_dataset.rds")
output_dir <- file.path(root, "aggregate_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- readRDS(input_file)
trait_map <- as.data.table(read_excel(metadata_file, sheet = "Table S2"))[, .(PRS_name, NMR_name)]

inverse_normal <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - 0.5) / sum(!is.na(r)))
}

evaluate_trait <- function(prs_name, nmr_name) {
  if (!all(c(prs_name, nmr_name) %in% names(dat))) return(NULL)
  true_fit <- lm(
    reformulate(c("AGE", "FEMALE", "COYOACAN", "processing_duration", "donation_month", "donation_hour", paste0("PC", 1:7)), response = nmr_name),
    data = dat, na.action = na.exclude
  )
  pred_fit <- lm(
    reformulate(c("AGE", "FEMALE", "COYOACAN", paste0("PC", 1:7)), response = prs_name),
    data = dat, na.action = na.exclude
  )
  measured <- inverse_normal(residuals(true_fit))
  predicted <- inverse_normal(residuals(pred_fit))
  keep <- is.finite(measured) & is.finite(predicted)
  data.table(
    PRS_name = prs_name,
    NMR_name = nmr_name,
    Rho = cor(measured[keep], predicted[keep], method = "spearman"),
    R2 = cor(measured[keep], predicted[keep], method = "pearson")^2,
    N = sum(keep)
  )
}

results <- rbindlist(Map(evaluate_trait, trait_map$PRS_name, trait_map$NMR_name), fill = TRUE)
saveRDS(results, file.path(output_dir, "interval_models_mcps_overall.rds"))
