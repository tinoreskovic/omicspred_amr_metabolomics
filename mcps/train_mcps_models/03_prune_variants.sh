#!/usr/bin/env bash

set -euo pipefail

WK="${OMICSPRED_MCPS_ROOT:?Set OMICSPRED_MCPS_ROOT}"
OUT="$WK/training_inputs"
PER="$OUT/per_trait"
PRUNED="$OUT/pruned"

GENO="${MCPS_GENOTYPE_DIR:?Set MCPS_GENOTYPE_DIR}"
PFILE_PREF="$GENO/mcps-freeze150k_qcd_chr"

mkdir -p "$PRUNED"

shopt -s nullglob
files=("$PER"/*__variants_unthinned.tsv.gz)
index="${TASK_INDEX:?Set TASK_INDEX to a zero-based trait index}"
if (( index < 0 || index >= ${#files[@]} )); then
    echo "ERROR: TASK_INDEX is outside the available trait range"
    exit 1
fi
f="${files[$index]}"
base=$(basename "$f" "__variants_unthinned.tsv.gz")
outpref="$PRUNED/$base"

echo "[$(date)] LD pruning trait: $base"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

tmp_keep_all="$tmpdir/all_pruned.ids"
> "$tmp_keep_all"

for chr in {1..22}; do
    zcat "$f" | awk -v c="$chr" 'BEGIN{FS="\t"} $4==c {print $3}' > "$tmpdir/chr${chr}.ids"

    [[ -s "$tmpdir/chr${chr}.ids" ]] || continue

    v_count=$(wc -l < "$tmpdir/chr${chr}.ids")

    if [ "$v_count" -eq 1 ]; then
        cat "$tmpdir/chr${chr}.ids" >> "$tmp_keep_all"
        continue
    fi

    plink2 \
        --pfile "${PFILE_PREF}${chr}" \
        --keep "$OUT/samples.keep" \
        --extract "$tmpdir/chr${chr}.ids" \
        --maf 0.005 \
        --indep-pairwise 1000kb 0.8 \
        --threads "${N_THREADS:-1}" \
        --out "$tmpdir/chr${chr}" \
        --silent || true

    if [[ -f "$tmpdir/chr${chr}.prune.in" ]]; then
        cat "$tmpdir/chr${chr}.prune.in" >> "$tmp_keep_all"
    else
        cat "$tmpdir/chr${chr}.ids" >> "$tmp_keep_all"
    fi
done

if [[ -s "$tmp_keep_all" ]]; then
    sort -u "$tmp_keep_all" > "${outpref}__variants_pruned.keep"
    n_keep=$(wc -l < "${outpref}__variants_pruned.keep")
    echo "[$(date)] Final pruned variants: $n_keep"
    echo "[$(date)] Done: ${outpref}__variants_pruned.keep"
else
    echo "ERROR: zero variants found for $base"
    exit 1
fi
