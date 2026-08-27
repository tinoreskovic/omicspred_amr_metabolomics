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
dat <- readRDS(file.path(root, "secure_intermediates", "mcps_interval_validation_dataset.rds"))
baseline <- as.data.table(readRDS(require_env("MCPS_BASELINE_RDS")))
trait_map <- as.data.table(read_excel(metadata_file, sheet = "Table S2"))[, .(PRS_name, NMR_name)]
output_dir <- file.path(root, "aggregate_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

baseline_columns <- c(
  "IID", "Score_EUR", "Score_AMR", "BASE_DIABETES", "BASE_HBA1C",
  "BASE_CVD", "BASE_CANCER", "BASE_CKD", "BASE_EMPHYSEMA",
  "BASE_CIRR", "BASE_PEP", "BASE_PAD"
)
missing_baseline <- setdiff(baseline_columns, names(baseline))
if (length(missing_baseline)) stop("Baseline input is missing: ", paste(missing_baseline, collapse = ", "))
dat <- as.data.table(dat)
dat[, IID := as.character(IID)]
baseline[, IID := as.character(IID)]
replace_columns <- intersect(setdiff(baseline_columns, "IID"), names(dat))
if (length(replace_columns)) dat[, (replace_columns) := NULL]
dat <- merge(dat, baseline[, ..baseline_columns], by = "IID", all.x = TRUE, sort = FALSE)

inverse_normal <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - 0.5) / sum(!is.na(r)))
}

evaluate_group <- function(group_data, label, group_name) {
  rbindlist(lapply(seq_len(nrow(trait_map)), function(i) {
    prs <- trait_map$PRS_name[i]
    trait <- trait_map$NMR_name[i]
    if (!all(c(prs, trait) %in% names(group_data))) return(NULL)
    true_formula <- reformulate(c("AGE", "FEMALE", "COYOACAN", "processing_duration", "donation_month", "donation_hour", paste0("PC", 1:7)), trait)
    pred_formula <- reformulate(c("AGE", "FEMALE", "COYOACAN", paste0("PC", 1:7)), prs)
    measured <- inverse_normal(residuals(lm(true_formula, group_data, na.action = na.exclude)))
    predicted <- inverse_normal(residuals(lm(pred_formula, group_data, na.action = na.exclude)))
    keep <- is.finite(measured) & is.finite(predicted)
    out <- data.table(PRS_name = prs, NMR_name = trait, Rho = NA_real_, R2 = NA_real_, N = sum(keep))
    if (sum(keep) >= 10L) {
      out[, `:=`(
        Rho = cor(measured[keep], predicted[keep], method = "spearman"),
        R2 = cor(measured[keep], predicted[keep], method = "pearson")^2
      )]
    }
    out[, (group_name) := label]
    out
  }), fill = TRUE)
}

groups <- list(
  sex = list(column = "Sex", values = list(Female = dat$FEMALE == 1, Male = dat$FEMALE == 0)),
  age = list(column = "AgeGroup", values = list("35~49" = dat$AGE >= 35 & dat$AGE < 50, "50~64" = dat$AGE >= 50 & dat$AGE < 65, "65 and older" = dat$AGE >= 65)),
  ancestry = list(column = "ancestry", values = list(MCPS_AMR = dat$Score_AMR >= 0.70, MCPS_EUR = dat$Score_EUR >= 0.70)),
  health = list(column = "HealthStat", values = list(
    "No diabetes" = dat$BASE_DIABETES == 0 & dat$BASE_HBA1C <= 6.5,
    "Diabetes" = dat$BASE_DIABETES == 1 | dat$BASE_HBA1C > 6.5,
    "No disease" = !(dat$BASE_CVD == 1 | dat$BASE_CANCER == 1 | dat$BASE_CKD == 1 |
      dat$BASE_EMPHYSEMA == 1 | dat$BASE_CIRR == 1 | dat$BASE_PEP == 1 |
      dat$BASE_PAD == 1 | dat$BASE_DIABETES == 1 | dat$BASE_HBA1C > 6.5),
    "Any disease" = dat$BASE_CVD == 1 | dat$BASE_CANCER == 1 | dat$BASE_CKD == 1 |
      dat$BASE_EMPHYSEMA == 1 | dat$BASE_CIRR == 1 | dat$BASE_PEP == 1 |
      dat$BASE_PAD == 1 | dat$BASE_DIABETES == 1 | dat$BASE_HBA1C > 6.5
  ))
)

for (name in names(groups)) {
  specification <- groups[[name]]
  result <- rbindlist(lapply(names(specification$values), function(label) {
    keep <- specification$values[[label]]
    evaluate_group(dat[which(!is.na(keep) & keep), ], label, specification$column)
  }), fill = TRUE)
  saveRDS(result, file.path(output_dir, paste0("interval_models_mcps_by_", name, ".rds")))
}

sex_ancestry <- list(
  AMR_Female = dat$FEMALE == 1 & dat$Score_AMR >= 0.70,
  EUR_Female = dat$FEMALE == 1 & dat$Score_EUR >= 0.70,
  AMR_Male = dat$FEMALE == 0 & dat$Score_AMR >= 0.70,
  EUR_Male = dat$FEMALE == 0 & dat$Score_EUR >= 0.70
)
result <- rbindlist(lapply(names(sex_ancestry), function(label) {
  keep <- sex_ancestry[[label]]
  evaluate_group(dat[which(!is.na(keep) & keep), ], label, "subgroup")
}), fill = TRUE)
saveRDS(result, file.path(output_dir, "interval_models_mcps_by_sex_ancestry.rds"))
