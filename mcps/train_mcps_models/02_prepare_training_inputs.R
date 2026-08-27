#!/usr/bin/env Rscript
library(data.table)
library(readxl)

wk <- Sys.getenv("OMICSPRED_MCPS_ROOT")
if (!nzchar(wk)) stop("Set OMICSPRED_MCPS_ROOT")
gwas_root <- Sys.getenv("MCPS_GWAS_RESULTS_DIR")
gwas_sum <- Sys.getenv("MCPS_GWAS_SUMMARY_FILE")
sample_csv <- Sys.getenv("MCPS_SAMPLE_SHEET")
if (!nzchar(gwas_root)) stop("Set MCPS_GWAS_RESULTS_DIR")
if (!nzchar(gwas_sum)) stop("Set MCPS_GWAS_SUMMARY_FILE")
if (!nzchar(sample_csv)) stop("Set MCPS_SAMPLE_SHEET")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
metadata_file <- normalizePath(file.path(dirname(script_file), "..", "..", "metadata", "omicspred_validation.xlsx"), mustWork = TRUE)

out_dir <- file.path(wk, "training_inputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ss <- fread(sample_csv)
if (!("path_prefix" %in% names(ss))) stop("column 'path_prefix' not found in sample sheet")

flip_strand <- function(x) {
  x <- gsub("A", "V", x)
  x <- gsub("T", "X", x)
  x <- gsub("C", "Y", x)
  x <- gsub("G", "Z", x)
  x <- gsub("V", "T", x)
  x <- gsub("X", "A", x)
  x <- gsub("Y", "G", x)
  x <- gsub("Z", "C", x)
  return(x)
}

t2 <- as.data.table(read_excel(
  metadata_file,
  sheet = "Table S2"
))
need_t2 <- c("PRS_name", "NMR_name", "Biomarker.Name")
if (!all(need_t2 %in% names(t2))) {
  stop("Table S2 is missing one or more required columns: PRS_name, NMR_name, Biomarker.Name")
}

target_traits <- unique(t2$Biomarker.Name)
target_traits <- target_traits[!is.na(target_traits) & nzchar(target_traits)]
cat("Target traits (Biomarker.Name):", length(target_traits), "\n")

sum_dt <- fread(gwas_sum)
if (!all(c("Phenotype name","Phenotype description") %in% names(sum_dt))) {
  stop("GWAS summary missing required columns: Phenotype name, Phenotype description")
}

clean_desc <- function(x) sub("\\s*\\(.*$", "", x)
sum_dt[, desc_clean := clean_desc(`Phenotype description`)]

map_dt <- sum_dt[desc_clean %in% target_traits,
                 .(pheno_code = `Phenotype name`,
                   Phenotype.description = `Phenotype description`,
                   Biomarker.Name = desc_clean)]
cat("Matched traits to GWAS folders:", nrow(map_dt), "\n")
if (nrow(map_dt) == 0) stop("No traits matched phenotype_GWAS_summary.txt cleaning rule")

pvar_list <- vector("list", 22)
for (k in 1:22) {
  pv_file <- paste0(ss$path_prefix[k], ".pvar")
  pvar <- fread(pv_file, select = c("#CHROM","POS","ID","REF","ALT"))
  setnames(pvar, c("chr","pos","ID","REF","ALT"))
  pvar[, chr := as.integer(chr)]
  pvar[, pos := as.integer(pos)]
  pvar_list[[k]] <- pvar
}
pvar_all_raw <- rbindlist(pvar_list, use.names = TRUE)

pvar_all_filt <- copy(pvar_all_raw)

multi <- pvar_all_filt[grepl("^rs", ID), .N, by=ID][N > 1]
if (nrow(multi) > 0) pvar_all_filt <- pvar_all_filt[!multi, on = .(ID)]

multi <- pvar_all_filt[!grepl("^rs", ID), .N, by=.(chr, pos)][N > 1]
if (nrow(multi) > 0) pvar_all_filt <- pvar_all_filt[!multi, on = .(chr, pos)]

pvar_all_filt <- pvar_all_filt[REF != flip_strand(ALT)]

pvar_all_filt <- pvar_all_filt[nchar(REF) == 1 & nchar(ALT) == 1]

setkey(pvar_all_raw,  chr, pos)
setkey(pvar_all_filt, chr, pos)

variant_id_map <- unique(pvar_all_filt[, .(ID, chr, pos, REF, ALT)])
fwrite(
  variant_id_map,
  file = file.path(out_dir, "variant_id_map.tsv"),
  sep = "\t",
  quote = FALSE
)
cat("Wrote variant_id_map.tsv (rows=", nrow(variant_id_map), ")\n")

metab   <- readRDS(file.path(wk, "secure_intermediates", "residualized_metabolites.rds"))
keepids <- data.table(FID = metab$IID, IID = metab$IID)
fwrite(keepids, file = file.path(out_dir, "samples.keep"), sep = "\t", quote = FALSE)
cat("Wrote samples.keep\n")

t2n <- unique(t2[, .(PRS_name, Biomarker.Name, NMR_name)])

missing_nmr <- t2n[is.na(NMR_name) | !nzchar(NMR_name), .N]
if (missing_nmr > 0) cat("WARN: Table S2 rows with no NMR_name:", missing_nmr, "\n")

trait_map <- merge(
  map_dt[, .(Biomarker.Name, pheno_code)],
  t2n[, .(Biomarker.Name, PRS_name, NMR_name)],
  by="Biomarker.Name",
  all.x=TRUE
)

trait_map[, safe_trait := gsub("[^A-Za-z0-9_]+", "_", Biomarker.Name)]

trait_map[, base := paste0(safe_trait, "__", pheno_code)]

dup_nmr <- trait_map[!is.na(NMR_name) & nzchar(NMR_name), .N, by=NMR_name][N > 1]
if (nrow(dup_nmr) > 0) {
  cat("WARN: duplicated NMR_name entries in trait_map.tsv\n")
}

fwrite(trait_map,
       file=file.path(out_dir, "trait_map.tsv"),
       sep="\t", quote=FALSE)

cat("Wrote trait_map.tsv (rows=", nrow(trait_map), ")\n")
cat("  trait_map non-missing NMR_name:", trait_map[!is.na(NMR_name) & nzchar(NMR_name), .N], "\n")

hits_out_dir <- file.path(out_dir, "per_trait")
dir.create(hits_out_dir, showWarnings = FALSE)

all_hits_mapped <- list()
count_rows <- list()

for (i in seq_len(nrow(map_dt))) {
  pheno_code <- map_dt$pheno_code[i]
  biomarker  <- map_dt$Biomarker.Name[i]
  gws_file   <- file.path(gwas_root, pheno_code, "All/output_files/genome-wide-significant.txt.gz")

  if (!file.exists(gws_file)) {
    message("Missing: ", gws_file, "  (skip)")
    count_rows[[length(count_rows)+1]] <- data.table(
      Biomarker.Name = biomarker,
      pheno_code = pheno_code,
      n_gws_raw = NA_integer_,
      n_gws_raw_uniqpos = NA_integer_,
      n_mapped_prefilter = 0L,
      n_mapped_postfilter = 0L,
      n_ids_prefilter = 0L,
      n_ids_postfilter = 0L
    )
    next
  }

  gws <- fread(gws_file, select = c("CHROM","GENPOS","ALLELE1","ALLELE0","BETA","SE","P","N","INFO","ID"))
  setnames(gws, c("CHROM","GENPOS"), c("chr","pos"))
  gws[, chr := as.integer(chr)]
  gws[, pos := as.integer(pos)]

  n_raw <- nrow(gws)
  n_raw_uniqpos <- uniqueN(gws, by=c("chr","pos"))

  mapped_pre <- merge(gws, pvar_all_raw,  by=c("chr","pos"), all=FALSE)
  mapped     <- merge(gws, pvar_all_filt, by=c("chr","pos"), all=FALSE)

  n_pre  <- nrow(mapped_pre)
  n_post <- nrow(mapped)
  n_ids_pre  <- uniqueN(mapped_pre$ID.y)
  n_ids_post <- uniqueN(mapped$ID.y)

  count_rows[[length(count_rows)+1]] <- data.table(
    Biomarker.Name = biomarker,
    pheno_code = pheno_code,
    n_gws_raw = n_raw,
    n_gws_raw_uniqpos = n_raw_uniqpos,
    n_mapped_prefilter = n_pre,
    n_mapped_postfilter = n_post,
    n_ids_prefilter = n_ids_pre,
    n_ids_postfilter = n_ids_post
  )

  if (n_post == 0) {
    message("No MCPS pvar matches after filtering for ", biomarker, " / ", pheno_code)
    next
  }

  mapped[, trait := biomarker]
  mapped[, pheno_code := pheno_code]

  all_hits_mapped[[length(all_hits_mapped) + 1]] <- mapped

  safe_name <- gsub("[^A-Za-z0-9_]+", "_", biomarker)

  fwrite(
    mapped[, .(trait, pheno_code, ID = ID.y, chr, pos, REF, ALT, ALLELE1, ALLELE0, BETA, SE, P, N, INFO)],
    file = file.path(hits_out_dir, paste0(safe_name, "__", pheno_code, "__hits_mapped.tsv.gz")),
    sep  = "\t"
  )

  fwrite(
    unique(mapped[, .(trait, pheno_code, ID = ID.y, chr, pos)]),
    file = file.path(hits_out_dir, paste0(safe_name, "__", pheno_code, "__variants_unthinned.tsv.gz")),
    sep  = "\t"
  )

  message("OK: ", biomarker,
          " | raw=", n_raw,
          " | postfilter IDs=", n_ids_post)
}

all_dt <- rbindlist(all_hits_mapped, use.names = TRUE, fill = TRUE)
fwrite(all_dt, file = file.path(out_dir, "gws_all_traits_hits_mapped.tsv.gz"), sep = "\t")
cat("Wrote: gws_all_traits_hits_mapped.tsv.gz (rows=", nrow(all_dt), ")\n")

counts_dt <- rbindlist(count_rows, use.names=TRUE, fill=TRUE)
fwrite(counts_dt, file=file.path(out_dir, "gws_counts_raw_filter_map.tsv"), sep="\t")
cat("Wrote counts: ", file.path(out_dir, "gws_counts_raw_filter_map.tsv"), "\n")
