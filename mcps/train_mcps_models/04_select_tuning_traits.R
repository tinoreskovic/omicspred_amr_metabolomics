#!/usr/bin/env Rscript
library(data.table)

set.seed(1979)

wk <- Sys.getenv("OMICSPRED_MCPS_ROOT")
if (!nzchar(wk)) stop("Set OMICSPRED_MCPS_ROOT")
outroot <- file.path(wk, "training_inputs")

trait_map <- fread(file.path(outroot, "trait_map.tsv"))
if (!("NMR_name" %in% names(trait_map))) stop("trait_map.tsv missing NMR_name")

ph <- fread(file.path(outroot, "metabolite_phenos.tsv.gz"), nrows = 1)
ph_cols <- setdiff(names(ph), "IID")

cand <- unique(trait_map[!is.na(NMR_name) & nzchar(NMR_name), NMR_name])
cand <- intersect(cand, ph_cols)

if (length(cand) == 0) stop("No candidates after intersecting trait_map NMR_name with phenotype columns")

tune <- sample(cand, size = min(10, length(cand)))

fwrite(data.table(tune),
       file = file.path(outroot, "tuning_list.txt"),
       sep = "\n", col.names = FALSE, quote = FALSE)

cat("Wrote tuning_list.txt with n=", length(tune), "\n")
