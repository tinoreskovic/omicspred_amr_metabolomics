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
output_dir <- file.path(root, "aggregate_results")
mcps <- as.data.table(readRDS(file.path(output_dir, "interval_models_mcps_overall.rds")))
external <- as.data.table(read_excel(metadata_file, sheet = "Table S2"))

cohorts <- list(
  UKB_EUR = c("UKB_R2", "UKB_Rho"),
  MEC_CN = c("MEC_CN_R2", "MEC_CN_Rho"),
  MEC_IN = c("MEC_IN_R2", "MEC_IN_Rho"),
  MEC_MA = c("MEC_MA_R2", "MEC_MA_Rho")
)

combined <- list(mcps[, .(PRS_name, NMR_name, Rho, R2, N, cohort_name = "MCPS_ALL")])
for (cohort in names(cohorts)) {
  columns <- cohorts[[cohort]]
  if (!all(c("PRS_name", columns) %in% names(external))) stop("Validation table is missing fields for ", cohort)
  values <- merge(mcps[, .(PRS_name, NMR_name)], external[, c("PRS_name", columns), with = FALSE], by = "PRS_name")
  setnames(values, columns, c("R2", "Rho"))
  values[, `:=`(N = NA_real_, cohort_name = cohort)]
  combined[[cohort]] <- values[, .(PRS_name, NMR_name, Rho, R2, N, cohort_name)]
}

saveRDS(rbindlist(combined, fill = TRUE), file.path(output_dir, "interval_models_by_cohort.rds"))
