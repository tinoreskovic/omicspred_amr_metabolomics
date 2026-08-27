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

interval_archive <- require_env("UKB_INTERVAL_SCORE_ARCHIVE")
mcps_archive <- require_env("UKB_MCPS_SCORE_ARCHIVE")
ancestry_archive <- require_env("UKB_ANCESTRY_ARCHIVE")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
metadata_file <- normalizePath(file.path(dirname(script_file), "..", "metadata", "omicspred_validation.xlsx"), mustWork = TRUE)
root <- require_env("OMICSPRED_UKB_ROOT")
output_dir <- file.path(root, "secure_intermediates")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

extract_member <- function(archive, patterns) {
  members <- utils::untar(archive, list = TRUE)
  for (pattern in patterns) {
    hits <- members[grepl(pattern, members, perl = TRUE)]
    if (length(hits)) {
      directory <- tempfile("ukb_archive_")
      dir.create(directory)
      utils::untar(archive, files = hits[1], exdir = directory)
      return(file.path(directory, hits[1]))
    }
  }
  stop("Required member was not found in ", basename(archive))
}

read_scores <- function(archive) {
  fread(extract_member(archive, c("(^|/)score/aggregated_scores\\.txt\\.gz$", "aggregated_scores\\.txt\\.gz$")))
}

read_ancestry <- function(archive) {
  fread(extract_member(archive, c("popsimilarity\\.txt\\.gz$", "ancestry\\.tsv(\\.gz)?$")))
}

interval <- read_scores(interval_archive)
mcps <- read_scores(mcps_archive)
ancestry <- read_ancestry(ancestry_archive)
stopifnot(all(c("IID", "PGS", "SUM") %in% names(interval)))
stopifnot(all(c("IID", "PGS", "SUM") %in% names(mcps)))
stopifnot(all(c("IID", "RF_P_AMR") %in% names(ancestry)))

amr_ids <- unique(as.character(ancestry[RF_P_AMR >= 0.90, IID]))
interval[, IID := as.character(IID)]
mcps[, IID := as.character(IID)]
interval <- interval[IID %chin% amr_ids]
mcps <- mcps[IID %chin% amr_ids]

metadata <- as.data.table(read_excel(metadata_file, sheet = "Table S2"))
stopifnot(all(c("PRS_name", "Biomarker.Name", "UKB_field_id") %in% names(metadata)))
mapping <- metadata[, .(
  PRS_name = trimws(as.character(PRS_name)),
  trait = trimws(as.character(Biomarker.Name)),
  nmr_field = paste0("p", as.integer(UKB_field_id), "_i0")
)]
interval_map <- mapping[, .(PGS = PRS_name, trait, nmr_field)]
mcps_map <- mapping[, .(PGS = paste0(PRS_name, "_MCPS_80_20_fixed"), trait, nmr_field)]
interval <- merge(interval, interval_map, by = "PGS")
mcps <- merge(mcps, mcps_map, by = "PGS")

fwrite(interval, file.path(output_dir, "ukb_interval_scores_annotated.tsv.gz"), sep = "	")
fwrite(mcps, file.path(output_dir, "ukb_mcps_scores_annotated.tsv.gz"), sep = "	")
