#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lubridate)
})

require_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop("Set environment variable ", name)
  value
}

nmr_file <- require_env("UKB_NMR_FILE")
interval_file <- require_env("UKB_INTERVAL_SCORES")
mcps_file <- require_env("UKB_MCPS_SCORES")
ancestry_file <- require_env("UKB_ANCESTRY_FILE")
output_dir <- file.path(require_env("OMICSPRED_UKB_ROOT"), "aggregate_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

interval <- fread(interval_file)
mcps <- fread(mcps_file)
ancestry <- fread(ancestry_file)

score_columns <- c("IID", "PGS", "SUM", "nmr_field", "trait")
stopifnot(all(score_columns %in% names(interval)))
stopifnot(all(score_columns %in% names(mcps)))
stopifnot(all(c("IID", "RF_P_AMR") %in% names(ancestry)))

amr_ids <- unique(as.character(ancestry[RF_P_AMR >= 0.90, IID]))
interval[, IID := as.character(IID)]
mcps[, IID := as.character(IID)]
interval <- interval[IID %chin% amr_ids]
mcps <- mcps[IID %chin% amr_ids]

needed_nmr_fields <- sort(unique(c(interval$nmr_field, mcps$nmr_field)))
needed_nmr_fields <- needed_nmr_fields[!is.na(needed_nmr_fields) & nzchar(needed_nmr_fields)]
raw_covariates <- c(
  "eid", "p21003_i0", "p22001", "p54_i0",
  paste0("p22009_a", 1:10),
  "p23658_i0", "p23659_i0",
  paste0("p3166_i0_a", 0:5)
)

nmr_header <- names(fread(nmr_file, nrows = 0))
missing_columns <- setdiff(c(raw_covariates, needed_nmr_fields), nmr_header)
if (length(missing_columns)) {
  stop("UKB NMR input is missing: ", paste(missing_columns, collapse = ", "))
}
nmr <- fread(nmr_file, select = c(raw_covariates, needed_nmr_fields))
nmr[, eid := as.character(eid)]
nmr <- nmr[eid %chin% amr_ids]

nmr[, AGE := as.numeric(p21003_i0)]
nmr[, FEMALE := fifelse(p22001 == 0, 1L, fifelse(p22001 == 1, 0L, NA_integer_))]
nmr[, CENTRE := factor(p54_i0)]
for (index in 1:10) {
  setnames(nmr, paste0("p22009_a", index), paste0("PC", index))
  nmr[, (paste0("PC", index)) := as.numeric(get(paste0("PC", index)))]
}

parse_datetime <- function(x) {
  x <- as.character(x)
  value <- suppressWarnings(ymd_hms(x, tz = "UTC"))
  missing <- is.na(value)
  if (any(missing)) value[missing] <- suppressWarnings(dmy_hms(x[missing], tz = "UTC"))
  missing <- is.na(value)
  if (any(missing)) value[missing] <- suppressWarnings(ymd(x[missing], tz = "UTC"))
  missing <- is.na(value)
  if (any(missing)) value[missing] <- suppressWarnings(dmy(x[missing], tz = "UTC"))
  value
}

measured_time <- parse_datetime(nmr$p23658_i0)
prepared_time <- parse_datetime(nmr$p23659_i0)
nmr[, processing_duration := as.numeric(difftime(measured_time, prepared_time, units = "hours"))]
nmr[!is.finite(processing_duration) | processing_duration < -1e-6, processing_duration := NA_real_]

donation_columns <- paste0("p3166_i0_a", 0:5)
donation_time <- do.call(dplyr::coalesce, lapply(donation_columns, function(column) as.character(nmr[[column]])))
donation_time[grepl("^\\s*1900-01-01", donation_time)] <- NA_character_

donation_hour <- rep(NA_integer_, length(donation_time))
time_only <- grepl("^\\s*\\d{1,2}:\\d{2}(:\\d{2})?\\s*$", donation_time)
donation_hour[time_only] <- hour(hms(paste0(
  trimws(donation_time[time_only]),
  ifelse(grepl(":\\d{2}:\\d{2}$", donation_time[time_only]), "", ":00")
)))
remaining <- is.na(donation_hour) & !is.na(donation_time)
if (any(remaining)) donation_hour[remaining] <- hour(parse_datetime(donation_time[remaining]))
donation_datetime <- parse_datetime(donation_time)
donation_month <- rep(NA_character_, length(donation_time))
donation_month[!is.na(donation_datetime)] <- format(donation_datetime[!is.na(donation_datetime)], "%m/%Y")

nmr[, donation_hour := factor(donation_hour)]
nmr[, donation_month := factor(donation_month)]
nmr[, missing_donation_time := ifelse(is.na(donation_hour), "Yes", "No")]
setnames(nmr, "eid", "IID")

inverse_normal <- function(x) {
  ranks <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((ranks - 0.5) / sum(!is.na(ranks)))
}

pvalue_from_squared_correlation <- function(r2, n) {
  if (!is.finite(r2) || n < 3L || r2 < 0) {
    return(list(p = NA_character_, log10_p = NA_real_))
  }
  log_p <- if (r2 >= 1) -Inf else pf(
    r2 * (n - 2) / (1 - r2), 1, n - 2,
    lower.tail = FALSE, log.p = TRUE
  )
  log10_p <- log_p / log(10)
  p_text <- if (is.infinite(log_p) && log_p < 0) {
    "0"
  } else {
    exponent <- floor(log10_p)
    sprintf("%.15fe%d", 10^(log10_p - exponent), as.integer(exponent))
  }
  list(p = p_text, log10_p = log10_p)
}

prepare_measured <- function(x) {
  x <- as.numeric(x)
  positive <- x[is.finite(x) & x > 0]
  if (!length(positive)) return(rep(NA_real_, length(x)))
  zero <- which(x == 0)
  if (length(zero)) x[zero] <- runif(length(zero), 0.001, 0.9) * min(positive)
  log(x)
}

safe_residuals <- function(y, data, covariates) {
  model_data <- data.table::copy(data[, ..covariates])
  model_data[, y := y]
  keep <- complete.cases(model_data)
  result <- rep(NA_real_, nrow(model_data))
  if (!any(keep)) return(result)

  analysis_data <- model_data[keep]
  usable <- covariates[vapply(covariates, function(column) {
    length(unique(analysis_data[[column]])) > 1L
  }, logical(1))]
  formula <- if (length(usable)) reformulate(usable, response = "y") else y ~ 1
  result[keep] <- residuals(lm(formula, data = analysis_data))
  result
}

true_covariates_overall <- c(
  "AGE", "FEMALE", "CENTRE", "processing_duration",
  "donation_month", "donation_hour", paste0("PC", 1:10)
)
true_covariates_subgroups <- c(
  "AGE", "FEMALE", "CENTRE", "processing_duration",
  "donation_month", "donation_hour", "missing_donation_time", paste0("PC", 1:10)
)
predicted_covariates <- c("AGE", "FEMALE", "CENTRE", paste0("PC", 1:10))

measured_residuals <- function(covariates) {
  set.seed(20260821)
  setNames(lapply(needed_nmr_fields, function(field) {
    safe_residuals(prepare_measured(nmr[[field]]), nmr, covariates)
  }), needed_nmr_fields)
}

predicted_residuals <- function(scores) {
  pgs_names <- unique(scores$PGS)
  result <- setNames(vector("list", length(pgs_names)), pgs_names)
  covariate_data <- nmr[, c("IID", predicted_covariates), with = FALSE]
  for (pgs in pgs_names) {
    score <- scores[PGS == pgs, .(IID, SUM = as.numeric(SUM))]
    score <- merge(score, covariate_data, by = "IID", all.x = TRUE, sort = FALSE)
    residual <- safe_residuals(score$SUM, score, predicted_covariates)
    aligned <- rep(NA_real_, nrow(nmr))
    position <- match(score$IID, nmr$IID)
    aligned[position[!is.na(position)]] <- residual[!is.na(position)]
    result[[pgs]] <- aligned
  }
  result
}

evaluate <- function(scores, model_label, true_residuals, pred_residuals, index, subgroup) {
  map <- unique(scores[, .(PGS, nmr_field, trait)])
  if (any(map[, .N, by = PGS]$N != 1L)) stop("Each PGS must map to one Nightingale NMR trait")

  rbindlist(lapply(seq_len(nrow(map)), function(i) {
    measured <- true_residuals[[map$nmr_field[i]]][index]
    predicted <- pred_residuals[[map$PGS[i]]][index]
    keep <- is.finite(measured) & is.finite(predicted)
    n_overlap <- sum(keep)
    rho <- r2 <- NA_real_
    rho_p <- r2_p <- list(p = NA_character_, log10_p = NA_real_)
    if (n_overlap >= 10L) {
      measured <- inverse_normal(measured[keep])
      predicted <- inverse_normal(predicted[keep])
      rho <- cor(measured, predicted, method = "spearman")
      r2 <- cor(measured, predicted, method = "pearson")^2
      rho_p <- pvalue_from_squared_correlation(rho^2, n_overlap)
      r2_p <- pvalue_from_squared_correlation(r2, n_overlap)
    }
    data.table(
      model = model_label,
      subgroup = subgroup,
      PGS = map$PGS[i],
      trait = map$trait[i],
      n_overlap = n_overlap,
      rho_spearman = rho,
      p_spearman = rho_p$p,
      log10_p_spearman = rho_p$log10_p,
      R2_pearson = r2,
      p_R2 = r2_p$p,
      log10_p_R2 = r2_p$log10_p
    )
  }))
}

overall_residuals <- measured_residuals(true_covariates_overall)
subgroup_residuals <- measured_residuals(true_covariates_subgroups)
interval_predicted <- predicted_residuals(interval)
mcps_predicted <- predicted_residuals(mcps)

groups <- list(
  All = seq_len(nrow(nmr)),
  Female = which(nmr$FEMALE == 1L),
  Male = which(nmr$FEMALE == 0L),
  `35~49` = which(nmr$AGE >= 35 & nmr$AGE < 50),
  `50~64` = which(nmr$AGE >= 50 & nmr$AGE < 65),
  `65 and older` = which(nmr$AGE >= 65)
)

results <- list()
for (subgroup in names(groups)) {
  true_residuals <- if (subgroup == "All") overall_residuals else subgroup_residuals
  results <- c(results, list(
    evaluate(interval, "INTERVAL-trained", true_residuals, interval_predicted, groups[[subgroup]], subgroup),
    evaluate(mcps, "MCPS-trained", true_residuals, mcps_predicted, groups[[subgroup]], subgroup)
  ))
}

fwrite(
  rbindlist(results),
  file.path(output_dir, "ukb_amr_model_performance.tsv"),
  sep = "\t",
  quote = TRUE
)
