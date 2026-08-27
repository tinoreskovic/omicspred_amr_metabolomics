#!/usr/bin/env python3
"""
Aggregate per-metabolite training metadata.
"""

import os
from pathlib import Path
import pandas as pd


WK = Path(os.environ["OMICSPRED_MCPS_ROOT"]).expanduser()
OUT_ROOT = WK / "model_training"
STAT_DIR = OUT_ROOT / "stats"
OUT_FILE = OUT_ROOT / "training_summary_all_metabolites.tsv"


print(f"Searching for stat files in {STAT_DIR}...")
files = sorted(STAT_DIR.glob("*_fixed_summary.tsv"))

if not files:
    print("Error: No stats files found. Check that model fitting finished.")
    exit(1)

dfs = []
for f in files:
    try:
        dfs.append(pd.read_csv(f, sep="\t"))
    except Exception as e:
        print(f"Skipping {f.name} due to error: {e}")

combined_df = pd.concat(dfs, ignore_index=True)


combined_df.to_csv(OUT_FILE, sep="\t", index=False)
print(f"Success! Wrote summary for {len(files)} metabolites to:")
print(f"  {OUT_FILE}")
