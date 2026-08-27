#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(readxl)
})

require_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set environment variable ", name)
  value
}

root <- require_env("OMICSPRED_MCPS_ROOT")
input_file <- require_env("MCPS_INPUT_RDS")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
metadata_file <- normalizePath(file.path(dirname(script_file), "..", "..", "metadata", "omicspred_validation.xlsx"), mustWork = TRUE)
output_dir <- file.path(root, "secure_intermediates")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260821)
dat <- readRDS(input_file)
trait_map <- as.data.table(read_excel(metadata_file, sheet = "Table S2"))[, .(PRS_name, NMR_name)]
trait_columns <- intersect(trait_map$NMR_name, names(dat))
if (length(trait_columns) == 0L) stop("No mapped metabolite columns were found")

required_covariates <- c(
  "IID", "FID", "FEMALE", "AGE", "COYOACAN",
  paste0("PC", 1:10), "processing_duration", "donation_month", "donation_hour"
)
missing_covariates <- setdiff(required_covariates, names(dat))
if (length(missing_covariates)) stop("Missing covariates: ", paste(missing_covariates, collapse = ", "))

dat <- dat %>%
  filter(rowMeans(is.na(across(all_of(trait_columns)))) <= 0.30) %>%
  filter(if_all(all_of(required_covariates), ~ !is.na(.x)))

replace_zero <- function(x) {
  positive <- x[is.finite(x) & x > 0]
  if (!length(positive)) return(rep(NA_real_, length(x)))
  zero <- which(x == 0)
  if (length(zero)) x[zero] <- runif(length(zero), 0.001, 0.9) * min(positive)
  x
}

dat[trait_columns] <- lapply(dat[trait_columns], replace_zero)
dat[trait_columns] <- lapply(dat[trait_columns], log)
dat[trait_columns] <- lapply(dat[trait_columns], function(x) {
  z <- abs(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  x[z > 10] <- NA_real_
  x
})

saveRDS(dat, file.path(output_dir, "mcps_interval_validation_dataset.rds"))
