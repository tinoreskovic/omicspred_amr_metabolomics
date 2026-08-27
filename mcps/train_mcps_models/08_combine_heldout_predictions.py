#!/usr/bin/env python3
"""
Combine the participant-indexed held-out predictions inside the secure workspace.
"""

import gzip
import os
from pathlib import Path

import pandas as pd


WK = Path(os.environ["OMICSPRED_MCPS_ROOT"]).expanduser()

OUT_ROOT = WK / "model_training"
PRED_DIR = OUT_ROOT / "predictions"
OUT_FILE = OUT_ROOT / "heldout_preds_long.tsv.gz"


files = sorted(PRED_DIR.glob("*_heldout.tsv.gz"))
if not files:
    raise SystemExit(f"No _heldout.tsv.gz files found in {PRED_DIR}")

dfs = []
for f in files:
    dfs.append(pd.read_csv(gzip.open(f, "rt"), sep="\t"))

long_df = pd.concat(dfs, ignore_index=True)


long_df.to_csv(OUT_FILE, sep="\t", index=False, compression="gzip")
print(f"Wrote {OUT_FILE}  (n = {len(long_df):,} rows; "
      f"{len(files)} metabolites)")
