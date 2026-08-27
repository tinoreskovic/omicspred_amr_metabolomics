#!/usr/bin/env Rscript

library(data.table)

wd <- Sys.getenv("OMICSPRED_MCPS_ROOT")
if (!nzchar(wd)) stop("Set OMICSPRED_MCPS_ROOT")
dir.create(file.path(wd, "training_inputs"), recursive = TRUE, showWarnings = FALSE)

infile      <- file.path(wd, "secure_intermediates", "mcps_interval_validation_dataset.rds")

outfile     <- file.path(wd, "secure_intermediates", "residualized_metabolites.rds")

prs <- readRDS(infile)

meta_cols <- 161:301
stopifnot(length(meta_cols) == 141)

covars <- c("AGE", "FEMALE", "COYOACAN",
            "processing_duration", "donation_month", "donation_hour",
             paste0("PC", 1:7))

irnt <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - 0.5) / sum(!is.na(r)))
}

N          <- nrow(prs)
res_matrix <- matrix(NA_real_, N, length(meta_cols))

for (i in seq_along(meta_cols)) {
  y <- prs[[meta_cols[i]]]
  lm_fit <- lm(y ~ ., data = prs[, c(covars)], na.action = na.exclude)
  res_matrix[, i] <- irnt(residuals(lm_fit))
}

out_df <- data.frame(
  IID = prs$IID,
  setNames(as.data.frame(res_matrix), names(prs)[meta_cols])
)

saveRDS(out_df, outfile)
message("Written: ", outfile)

fwrite(out_df, file.path(wd, "training_inputs", "metabolite_phenos.tsv.gz"), sep="\t")
