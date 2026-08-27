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
predictions_file <- file.path(root, "model_training", "heldout_preds_long.tsv.gz")
true_file <- file.path(root, "secure_intermediates", "residualized_metabolites.rds")
covariate_file <- file.path(root, "secure_intermediates", "mcps_interval_validation_dataset.rds")
baseline_file <- require_env("MCPS_BASELINE_RDS")
results_dir <- file.path(root, "aggregate_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

inverse_normal <- function(x) {
  keep <- !is.na(x)
  result <- rep(NA_real_, length(x))
  if (any(keep)) result[keep] <- qnorm((rank(x[keep], ties.method = "average") - 0.5) / sum(keep))
  result
}

calculate_stats <- function(data) {
  predicted_covariates <- c("AGE", "FEMALE", "COYOACAN", paste0("PC", 1:7))
  data[!is.na(y_true) & !is.na(y_pred), {
    keep <- complete.cases(.SD[, c("y_pred", predicted_covariates), with = FALSE])
    predicted_residuals <- rep(NA_real_, .N)
    if (any(keep)) {
      fit <- lm(
        y_pred ~ .,
        data = .SD[keep, c("y_pred", predicted_covariates), with = FALSE]
      )
      predicted_residuals[keep] <- residuals(fit)
    }
    predicted <- inverse_normal(predicted_residuals)
    paired <- !is.na(y_true) & !is.na(predicted)
    if (sum(paired) < 3L) {
      .(Rho = NA_real_, R2 = NA_real_)
    } else {
      .(
        Rho = cor(y_true[paired], predicted[paired], method = "spearman"),
        R2 = cor(y_true[paired], predicted[paired], method = "pearson")^2
      )
    }
  }, by = .(NMR_name = metabolite)]
}

add_prs_name <- function(results, mapping) {
  output <- merge(results, mapping, by = "NMR_name", all.x = TRUE, sort = FALSE)
  grouping_columns <- setdiff(names(output), c("PRS_name", "NMR_name", "Rho", "R2"))
  setcolorder(output, c("PRS_name", "NMR_name", grouping_columns, "Rho", "R2"))
  output
}

predictions <- fread(predictions_file, select = c("IID", "metabolite", "y_pred"))
true_values <- as.data.table(readRDS(true_file))
covariates <- as.data.table(readRDS(covariate_file))
baseline <- as.data.table(readRDS(baseline_file))

predictions[, IID := as.character(IID)]
true_values[, IID := as.character(IID)]
covariates[, IID := as.character(IID)]
baseline[, IID := as.character(IID)]

mapping <- as.data.table(read_excel(metadata_file, sheet = "Table S2"))[, .(NMR_name, PRS_name)]
true_long <- melt(
  true_values,
  id.vars = "IID",
  measure.vars = setdiff(names(true_values), "IID"),
  variable.name = "metabolite",
  value.name = "y_true"
)
data <- merge(predictions, true_long, by = c("IID", "metabolite"), all.x = TRUE, sort = FALSE)

covariate_columns <- c("IID", "AGE", "FEMALE", "COYOACAN", paste0("PC", 1:7))
baseline_columns <- c(
  "IID", "Score_EUR", "Score_AMR", "BASE_DIABETES", "BASE_HBA1C",
  "BASE_CVD", "BASE_CANCER", "BASE_CKD", "BASE_EMPHYSEMA",
  "BASE_CIRR", "BASE_PEP", "BASE_PAD"
)
missing_covariates <- setdiff(covariate_columns, names(covariates))
missing_baseline <- setdiff(baseline_columns, names(baseline))
if (length(missing_covariates)) stop("Covariate input is missing: ", paste(missing_covariates, collapse = ", "))
if (length(missing_baseline)) stop("Baseline input is missing: ", paste(missing_baseline, collapse = ", "))

data <- merge(data, covariates[, ..covariate_columns], by = "IID", all.x = TRUE, sort = FALSE)
data <- merge(data, baseline[, ..baseline_columns], by = "IID", all.x = TRUE, sort = FALSE)
data[, Diabetes := BASE_DIABETES == 1 | BASE_HBA1C > 6.5]
data[, No_diabetes := BASE_DIABETES == 0 & BASE_HBA1C <= 6.5]
data[, Any_disease := BASE_CVD == 1 | BASE_CANCER == 1 | BASE_CKD == 1 |
  BASE_EMPHYSEMA == 1 | BASE_CIRR == 1 | BASE_PEP == 1 |
  BASE_PAD == 1 | Diabetes]
data[, No_disease := !Any_disease]
data[, Sex := fifelse(FEMALE == 1, "Female", "Male")]
data[, AgeGroup := cut(AGE, c(35, 50, 65, Inf), c("35~49", "50~64", "65 and older"), right = FALSE)]
data[, ancestry := fifelse(Score_EUR >= 0.70, "MCPS_EUR", fifelse(Score_AMR >= 0.70, "MCPS_AMR", NA_character_))]

overall <- add_prs_name(calculate_stats(data), mapping)
saveRDS(overall, file.path(results_dir, "mcps_models_heldout_overall.rds"))
fwrite(overall, file.path(results_dir, "mcps_models_heldout_overall.csv"))

sex <- add_prs_name(data[, calculate_stats(.SD), by = Sex], mapping)
age <- add_prs_name(data[!is.na(AgeGroup), calculate_stats(.SD), by = AgeGroup], mapping)
ancestry <- add_prs_name(data[!is.na(ancestry), calculate_stats(.SD), by = ancestry], mapping)
saveRDS(sex, file.path(results_dir, "mcps_models_heldout_by_sex.rds"))
saveRDS(age, file.path(results_dir, "mcps_models_heldout_by_age.rds"))
saveRDS(ancestry, file.path(results_dir, "mcps_models_heldout_by_ancestry.rds"))

health_groups <- list(
  "No diabetes" = data[No_diabetes == TRUE],
  "Diabetes" = data[Diabetes == TRUE],
  "No disease" = data[No_disease == TRUE],
  "Any disease" = data[Any_disease == TRUE]
)
health <- rbindlist(lapply(names(health_groups), function(label) {
  calculate_stats(health_groups[[label]])[, HealthStat := label]
}), fill = TRUE)
health <- add_prs_name(health, mapping)
saveRDS(health, file.path(results_dir, "mcps_models_heldout_by_health.rds"))
