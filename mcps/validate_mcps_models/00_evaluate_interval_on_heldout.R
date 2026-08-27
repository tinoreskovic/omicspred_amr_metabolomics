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
baseline_file <- require_env("MCPS_BASELINE_RDS")
test_manifest_file <- file.path(root, "model_training", "heldout_preds_long.tsv.gz")
analysis_file <- file.path(root, "secure_intermediates", "mcps_interval_validation_dataset.rds")
output_dir <- file.path(root, "aggregate_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

test_manifest <- fread(test_manifest_file, select = c("IID", "metabolite"))
test_manifest[, IID := as.character(IID)]
analysis <- as.data.table(readRDS(analysis_file))
baseline <- as.data.table(readRDS(baseline_file))
mapping <- as.data.table(read_excel(metadata_file, sheet = "Table S2"))[, .(PRS_name, NMR_name)]
analysis[, IID := as.character(IID)]
baseline[, IID := as.character(IID)]

baseline_columns <- c(
  "IID", "Score_EUR", "Score_AMR", "BASE_DIABETES", "BASE_HBA1C",
  "BASE_CVD", "BASE_CANCER", "BASE_CKD", "BASE_EMPHYSEMA",
  "BASE_CIRR", "BASE_PEP", "BASE_PAD"
)
missing_baseline <- setdiff(baseline_columns, names(baseline))
if (length(missing_baseline)) stop("Baseline input is missing: ", paste(missing_baseline, collapse = ", "))
replace_columns <- intersect(setdiff(baseline_columns, "IID"), names(analysis))
if (length(replace_columns)) analysis[, (replace_columns) := NULL]
analysis <- merge(analysis, baseline[, ..baseline_columns], by = "IID", all.x = TRUE, sort = FALSE)

analysis[, Sex := fifelse(FEMALE == 1, "Female", "Male")]
analysis[, AgeGroup := cut(AGE, c(35, 50, 65, Inf), c("35~49", "50~64", "65 and older"), right = FALSE)]
analysis[, ancestry := fifelse(Score_EUR >= 0.70, "MCPS_EUR", fifelse(Score_AMR >= 0.70, "MCPS_AMR", NA_character_))]
analysis[, Diabetes := BASE_DIABETES == 1 | BASE_HBA1C > 6.5]
analysis[, No_diabetes := BASE_DIABETES == 0 & BASE_HBA1C <= 6.5]
analysis[, Any_disease := BASE_CVD == 1 | BASE_CANCER == 1 | BASE_CKD == 1 |
  BASE_EMPHYSEMA == 1 | BASE_CIRR == 1 | BASE_PEP == 1 | BASE_PAD == 1 | Diabetes]
analysis[, No_disease := !Any_disease]

inverse_normal <- function(x) {
  keep <- !is.na(x)
  result <- rep(NA_real_, length(x))
  result[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / sum(keep))
  result
}

evaluate <- function(data) {
  result <- lapply(seq_len(nrow(mapping)), function(i) {
    score <- mapping$PRS_name[i]
    trait <- mapping$NMR_name[i]
    if (!all(c(score, trait) %in% names(data))) return(NULL)
    test_ids <- test_manifest[metabolite == trait, IID]
    true_covariates <- c("AGE", "FEMALE", "COYOACAN", "processing_duration", "donation_month", "donation_hour", paste0("PC", 1:7))
    predicted_covariates <- c("AGE", "FEMALE", "COYOACAN", paste0("PC", 1:7))
    required <- unique(c("IID", score, trait, true_covariates, predicted_covariates))
    keep <- data$IID %chin% test_ids & complete.cases(data[, ..required])
    if (sum(keep) < 10L) return(NULL)
    true_fit <- lm(reformulate(true_covariates, response = trait), data = data[keep])
    predicted_fit <- lm(reformulate(predicted_covariates, response = score), data = data[keep])
    measured <- inverse_normal(residuals(true_fit))
    predicted <- inverse_normal(residuals(predicted_fit))
    data.table(
      NMR_name = trait,
      PRS_name = score,
      old_test_R2 = cor(measured, predicted, method = "pearson")^2,
      old_test_Rho = cor(measured, predicted, method = "spearman"),
      N = length(measured)
    )
  })
  rbindlist(result, fill = TRUE)
}

save_grouped <- function(column, labels, output_name) {
  result <- rbindlist(lapply(labels, function(label) {
    subset <- analysis[get(column) == label & !is.na(get(column))]
    values <- evaluate(subset)
    values[, (column) := label]
    values
  }), fill = TRUE)
  saveRDS(result, file.path(output_dir, output_name))
}

overall <- evaluate(analysis)
saveRDS(overall, file.path(output_dir, "interval_models_mcps_heldout_overall.rds"))
save_grouped("Sex", c("Female", "Male"), "interval_models_mcps_heldout_by_sex.rds")
save_grouped("AgeGroup", c("35~49", "50~64", "65 and older"), "interval_models_mcps_heldout_by_age.rds")
save_grouped("ancestry", c("MCPS_AMR", "MCPS_EUR"), "interval_models_mcps_heldout_by_ancestry.rds")

health <- list(
  "No diabetes" = analysis[No_diabetes == TRUE],
  "Diabetes" = analysis[Diabetes == TRUE],
  "No disease" = analysis[No_disease == TRUE],
  "Any disease" = analysis[Any_disease == TRUE]
)
health_results <- rbindlist(lapply(names(health), function(label) {
  values <- evaluate(health[[label]])
  values[, HealthStat := label]
  values
}), fill = TRUE)
saveRDS(health_results, file.path(output_dir, "interval_models_mcps_heldout_by_health.rds"))
