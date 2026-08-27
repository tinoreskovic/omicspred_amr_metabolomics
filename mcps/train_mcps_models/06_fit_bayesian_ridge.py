#!/usr/bin/env python3
"""
Train the final MCPS genetic model for one Nightingale NMR metabolomic trait
using Bayesian ridge regression.

The script uses the pruned per-trait variant list, fixed non-informative
priors, and the first of five reproducible KFold splits as the held-out 20%
subset. It writes the fitted model, summary statistics, held-out predictions,
and an OmicsPred-format score file.
"""

import os, sys, gc, shutil, subprocess, tempfile, datetime, warnings, re
from pathlib import Path
from typing import Optional, Tuple, List, Dict

import joblib
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import BayesianRidge
from sklearn.model_selection import KFold
from sklearn.metrics import explained_variance_score

warnings.filterwarnings("ignore", message="overflow encountered", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="divide by zero encountered", category=RuntimeWarning)


WK = Path(os.environ["OMICSPRED_MCPS_ROOT"]).expanduser()
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]

OUTROOT = WK / "training_inputs"
PRUNED    = OUTROOT / "pruned"
KEEP_SAMP = OUTROOT / "samples.keep"
VARMAP    = OUTROOT / "variant_id_map.tsv"
PHENO     = OUTROOT / "metabolite_phenos.tsv.gz"
TRAITMAP  = OUTROOT / "trait_map.tsv"
NAME_MAP  = REPOSITORY_ROOT / "metadata" / "omicspred_validation.xlsx"

GENO_DIR = Path(os.environ["MCPS_GENOTYPE_DIR"]).expanduser()
PFILE_PREF = GENO_DIR / "mcps-freeze150k_qcd_chr{chr}"


OUT_ROOT = WK / "model_training"
MODEL_DIR = OUT_ROOT / "models"
STAT_DIR  = OUT_ROOT / "stats"
SCORE_DIR = OUT_ROOT / "scores"
PRED_DIR  = OUT_ROOT / "predictions"
for d in (MODEL_DIR, STAT_DIR, SCORE_DIR, PRED_DIR):
    d.mkdir(parents=True, exist_ok=True)


PREFER_PGEN = True


PRIOR = dict(alpha_1=1e-5, alpha_2=1e-5, lambda_1=1e-5, lambda_2=1e-5)
kf_outer = KFold(n_splits=5, shuffle=True, random_state=21)


COLS = ["chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight"]
TAB  = "\t"

def log(msg: str) -> None:
    print(f"[{datetime.datetime.now():%F %T}] {msg}", flush=True)


def _plink_common_args() -> list[str]:
    threads = os.environ.get("N_THREADS")
    mem_mb  = os.environ.get("MEMORY_MB")
    extra: list[str] = []
    if threads:
        extra += ["--threads", str(int(threads))]
    if mem_mb:
        try:
            extra += ["--memory", str(int(mem_mb))]
        except Exception:
            pass
    return extra

def run_plink_make_pgen(chr_num: int, ids_file: Path, outpref: Path) -> None:
    prefix = str(PFILE_PREF).format(chr=chr_num)
    cmd = [
        "plink2",
        "--pfile", prefix,
        "--keep", str(KEEP_SAMP),
        "--extract", str(ids_file),
        "--make-pgen",
        "--out", str(outpref),
        "--silent",
    ] + _plink_common_args()
    subprocess.run(cmd, check=True)

def run_plink_export_raw(chr_num: int, ids_file: Path, outpref: Path) -> None:
    prefix = str(PFILE_PREF).format(chr=chr_num)
    cmd = [
        "plink2",
        "--pfile", prefix,
        "--keep", str(KEEP_SAMP),
        "--extract", str(ids_file),
        "--export", "A",
        "--out", str(outpref),
        "--silent",
    ] + _plink_common_args()
    subprocess.run(cmd, check=True)


def _find_header_line(path: Path, required: Tuple[str, ...], max_lines: int = 200) -> Tuple[int, List[str]]:
    """Find the PLINK header row."""
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for i in range(max_lines):
            line = f.readline()
            if not line:
                break
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            toks = line.split("\t")
            norm = [t.lstrip("#") for t in toks]
            if all(r in norm for r in required):

                return i, norm
    raise RuntimeError(f"Could not find header with {required} within first {max_lines} lines of {path}")

def read_psam_iids(psam_path: Path) -> np.ndarray:
    skip, _ = _find_header_line(psam_path, required=("FID", "IID"))
    df = pd.read_csv(
        psam_path,
        sep="\t",
        header=0,
        skiprows=skip,
        dtype=str,
        engine="python",
        comment=None,
    )
    df.columns = [c.lstrip("#") for c in df.columns]
    if "IID" not in df.columns:
        raise RuntimeError(f"{psam_path} missing IID column; columns={df.columns.tolist()[:10]}")
    return df["IID"].astype(str).to_numpy()

def read_pvar_table(pvar_path: Path) -> pd.DataFrame:
    """Read CHROM, POS, ID, REF and ALT in file order."""
    skip, _ = _find_header_line(pvar_path, required=("CHROM", "POS", "REF", "ALT"))
    df = pd.read_csv(
        pvar_path,
        sep="\t",
        header=0,
        skiprows=skip,
        dtype=str,
        engine="python",
        comment=None,
    )
    df.columns = [c.lstrip("#") for c in df.columns]



    if "CHROM" not in df.columns or "POS" not in df.columns:
        raise RuntimeError(f"{pvar_path}: could not find CHROM/POS after normalization; columns={df.columns.tolist()[:12]}")
    if "REF" not in df.columns or "ALT" not in df.columns:
        raise RuntimeError(f"{pvar_path}: could not find REF/ALT after normalization; columns={df.columns.tolist()[:12]}")

    if "ID" not in df.columns:

        if df.shape[1] >= 3:
            df["ID"] = df.iloc[:, 2].astype(str)
        else:
            raise RuntimeError(f"{pvar_path} missing ID column and too few columns to recover it.")

    out = df[["CHROM", "POS", "ID", "REF", "ALT"]].copy()
    return out


def try_import_pgenlib():
    try:
        import pgenlib
        return pgenlib
    except Exception:
        return None

def read_pgen_as_dosage(pgenlib_mod, pgen_path: Path, n_samples: int) -> np.ndarray:
    """Read dosages or hard calls into a sample-by-variant matrix."""
    reader = pgenlib_mod.PgenReader(bytes(str(pgen_path), "utf-8"))
    v_ct = reader.get_variant_ct()
    s_ct = reader.get_raw_sample_ct()
    if s_ct != n_samples:
        raise RuntimeError(f"pgen sample_ct mismatch: pgen has {s_ct}, expected {n_samples}")


    dosage_methods = [m for m in ["read_dosages", "read_dosage", "read_dosage16"] if hasattr(reader, m)]
    hard_methods   = [m for m in ["read", "read_hardcalls", "read_alleles"] if hasattr(reader, m)]

    X = np.empty((n_samples, v_ct), dtype=np.float32)

    if dosage_methods:
        mname = dosage_methods[0]
        meth = getattr(reader, mname)
        buf = np.empty(n_samples, dtype=np.float32)
        for j in range(v_ct):
            ok = False

            try:
                meth(buf, j)
                ok = True
            except Exception:
                pass
            if not ok:

                meth(j, buf)
            X[:, j] = buf
        return X

    if hard_methods:
        mname = hard_methods[0]
        meth = getattr(reader, mname)
        buf_i8 = np.empty(n_samples, dtype=np.int8)
        for j in range(v_ct):
            ok = False
            try:
                meth(buf_i8, j)
                ok = True
            except Exception:
                pass
            if not ok:
                meth(j, buf_i8)
            col = buf_i8.astype(np.float32)
            col[col < 0] = np.nan
            X[:, j] = col
        return X

    raise RuntimeError("pgenlib: no usable read method found on this build")


_ID_RE = re.compile(r"^(?:chr)?(?P<chr>\d+):(?P<pos>\d+):(?P<ref>[ACGT]+):(?P<alt>[ACGT]+)(?:_(?P<a1>[ACGT]+))?$")

def parse_raw_to_numpy(raw_path: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], List[str], Optional[pd.DataFrame]]:
    """Read PLINK raw output and variant metadata."""
    if not raw_path.exists():
        return None, None, [], None


    with raw_path.open("r") as f:
        header = f.readline().strip().split()
    if len(header) <= 6:
        return None, None, [], None
    ids = [str(x) for x in header[6:]]

    df = pd.read_csv(raw_path, sep=r"\s+", engine="python")
    if df.shape[1] <= 6:
        return None, None, [], None

    iid = df["IID"].astype(str).to_numpy()
    X = df.iloc[:, 6:].to_numpy(dtype=np.float32, copy=False)
    del df

    rows = []
    for vid in ids:
        m = _ID_RE.match(vid)
        if not m:
            rows = []
            break
        chr_ = int(m.group("chr"))
        pos  = int(m.group("pos"))
        ref  = m.group("ref")
        alt  = m.group("alt")
        a1   = m.group("a1") or alt
        other = alt if a1 == ref else ref
        rows.append((chr_, pos, ref, alt, a1, other))

    meta_df = None
    if rows:
        meta_df = pd.DataFrame(rows, columns=["chr", "pos", "ref", "alt", "a1", "other"])
    return iid, X, ids, meta_df


def load_all141_nmr() -> list[str]:
    nm = pd.read_excel(NAME_MAP, sheet_name="Table S2", dtype=str)
    if not {"PRS_name", "NMR_name"}.issubset(nm.columns):
        raise RuntimeError("The metadata workbook must contain columns: PRS_name, NMR_name")
    metas = nm["NMR_name"].dropna().astype(str).tolist()

    seen = set()
    out = []
    for m in metas:
        if m not in seen:
            out.append(m); seen.add(m)
    return out

def pick_metab() -> str:
    metas = load_all141_nmr()
    if len(sys.argv) > 1:
        idx = int(sys.argv[1])
    else:
        idx = int(os.environ.get("TASK_INDEX", "0"))
    if idx < 0 or idx >= len(metas):
        raise SystemExit(f"Index {idx} out of range (0..{len(metas)-1})")
    return metas[idx]

def load_prs_name_map() -> dict:
    nm = pd.read_excel(NAME_MAP, sheet_name="Table S2", dtype=str)
    return nm.set_index("NMR_name")["PRS_name"].to_dict()

def load_traitmap_base(metab: str) -> str:
    tm = pd.read_csv(TRAITMAP, sep="\t", dtype=str)
    if not {"NMR_name", "base"}.issubset(tm.columns):
        raise RuntimeError("trait_map.tsv must contain columns: NMR_name, base")
    hit = tm.loc[tm["NMR_name"] == metab]
    if hit.shape[0] == 0:
        raise FileNotFoundError(f"{metab} not found in trait_map.tsv")
    return str(hit.iloc[0]["base"])

def build_id_to_chr() -> dict:
    vm = pd.read_csv(VARMAP, sep="\t", dtype={"ID": str})
    if not {"ID", "chr"}.issubset(vm.columns):
        raise RuntimeError("variant_id_map.tsv must contain columns ID and chr")
    vm["chr"] = vm["chr"].astype(int)
    return dict(zip(vm["ID"].astype(str), vm["chr"].astype(int)))

def read_ids(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]

def load_pheno_arrays(metab: str) -> Tuple[np.ndarray, np.ndarray]:
    ph = pd.read_csv(PHENO, sep="\t", usecols=["IID", metab])
    iid = ph["IID"].astype(str).to_numpy()
    y = ph[metab].to_numpy(dtype=np.float32)
    del ph
    return iid, y

def align_y_to_iids(iid_ref: np.ndarray, iid_ph: np.ndarray, y_ph: np.ndarray) -> np.ndarray:
    order = np.argsort(iid_ph)
    iid_ph_sorted = iid_ph[order]
    y_sorted = y_ph[order]
    idx = np.searchsorted(iid_ph_sorted, iid_ref)
    y = np.full((len(iid_ref),), np.nan, dtype=np.float32)
    candidates = np.flatnonzero(idx < len(iid_ph_sorted))
    matched = candidates[iid_ph_sorted[idx[candidates]] == iid_ref[candidates]]
    y[matched] = y_sorted[idx[matched]]
    return y

def first_fold_indices(n: int):
    return next(iter(kf_outer.split(np.empty(n))))


def make_score_table_from_pvar(pvar_df: pd.DataFrame, coefs: np.ndarray) -> pd.DataFrame:
    """Build the score table from PVAR allele metadata."""
    if len(pvar_df) != len(coefs):
        raise RuntimeError(f"PVAR length {len(pvar_df)} != number of coefs {len(coefs)}")
    chr_name = pvar_df["CHROM"].astype(str).str.replace("^chr", "", regex=True).astype("int64")
    pos = pvar_df["POS"].astype("int64")
    alt = pvar_df["ALT"].astype(str).str.strip()
    ref = pvar_df["REF"].astype(str).str.strip()
    df = pd.DataFrame({
        "chr_name": chr_name.to_numpy(),
        "chr_position": pos.to_numpy(),
        "effect_allele": alt.to_numpy(),
        "other_allele": ref.to_numpy(),
        "effect_weight": np.asarray(coefs, dtype="float64"),
    })[COLS]
    return df

def make_score_table_from_raw_meta(meta_df: pd.DataFrame, coefs: np.ndarray) -> pd.DataFrame:
    """Build the score table from raw variant metadata."""
    if meta_df is None or meta_df.shape[0] == 0:
        raise RuntimeError("Cannot build scorefile from .raw because variant IDs were not parseable into chr/pos/alleles.")
    if len(meta_df) != len(coefs):
        raise RuntimeError(f"RAW meta length {len(meta_df)} != number of coefs {len(coefs)}")
    df = pd.DataFrame({
        "chr_name": meta_df["chr"].astype("int64").to_numpy(),
        "chr_position": meta_df["pos"].astype("int64").to_numpy(),
        "effect_allele": meta_df["a1"].astype(str).to_numpy(),
        "other_allele": meta_df["other"].astype(str).to_numpy(),
        "effect_weight": np.asarray(coefs, dtype="float64"),
    })[COLS]
    return df

def write_scorefile(score_path: Path, header_lines: list[str], table: pd.DataFrame) -> None:
    with open(score_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(header_lines) + "\n")
        fh.write(TAB.join(COLS) + "\n")
        table.to_csv(
            fh,
            sep=TAB,
            index=False,
            header=False,
            float_format="%.10g",
            lineterminator="\n",
        )


def main():
    metab = pick_metab()
    log(f"Metabolite (NMR_name / phenotype column): {metab}")

    prs_map = load_prs_name_map()
    prs_name = prs_map.get(metab, metab)
    pgs_id = f"{prs_name}_MCPS_80_20_fixed"

    base = load_traitmap_base(metab)
    pruned_keep = PRUNED / f"{base}__variants_pruned.keep"
    if not pruned_keep.exists():
        log(f"Missing keepfile for {metab}: {pruned_keep} -> skipping")
        return

    ids_all = read_ids(pruned_keep)
    if not ids_all:
        log("No SNPs in keepfile -> skipping")
        return


    id_to_chr_lookup = build_id_to_chr()
    by_chr: Dict[int, List[str]] = {}
    miss = 0
    for vid in ids_all:
        c = id_to_chr_lookup.get(str(vid))
        if c is None:
            miss += 1
            continue
        by_chr.setdefault(int(c), []).append(str(vid))
    if miss:
        log(f"WARN: {miss} IDs not found in variant_id_map.tsv (ignored)")
    if not by_chr:
        log("No SNPs mapped -> skipping")
        return

    iid_ph, y_ph = load_pheno_arrays(metab)

    pgenlib_mod = try_import_pgenlib() if (PREFER_PGEN) else None
    if pgenlib_mod is None:
        log("pgenlib not available (or disabled) -> will use text .raw export where needed")

    tmpdir = Path(tempfile.mkdtemp(prefix=f"gws_fit_{base[:60]}_"))
    try:
        iid_ref: Optional[np.ndarray] = None
        mask = None
        tr_idx = te_idx = None

        Xtr_blocks: List[np.ndarray] = []
        Xte_blocks: List[np.ndarray] = []


        score_tables_parts: List[pd.DataFrame] = []
        raw_meta_parts: List[pd.DataFrame] = []

        for c in sorted(by_chr.keys()):
            ids_chr = by_chr[c]
            if not ids_chr:
                continue

            ids_file = tmpdir / f"chr{c}.ids"
            ids_file.write_text("\n".join(ids_chr) + "\n")

            used_pvar_df = None
            used_raw_meta = None

            if pgenlib_mod is not None:
                outpref = tmpdir / f"chr{c}_subset"
                log(f"PLINK make-pgen chr{c}: n_ids={len(ids_chr)}")
                run_plink_make_pgen(c, ids_file, outpref)

                pgen = outpref.with_suffix(".pgen")
                psam = outpref.with_suffix(".psam")
                pvar = outpref.with_suffix(".pvar")

                if not (pgen.exists() and psam.exists() and pvar.exists()):
                    log(f"WARN: missing pgen trio chr{c} -> fallback to .raw for this chr")
                    outpref2 = tmpdir / f"chr{c}_dos"
                    run_plink_export_raw(c, ids_file, outpref2)
                    raw = outpref2.with_suffix(".raw")
                    iid_chr, X_chr, _ids_order, meta_df = parse_raw_to_numpy(raw)
                    used_raw_meta = meta_df
                else:
                    iid_chr = read_psam_iids(psam)
                    used_pvar_df = read_pvar_table(pvar)
                    try:
                        X_chr = read_pgen_as_dosage(pgenlib_mod, pgen, n_samples=len(iid_chr))
                    except Exception as e:
                        log(f"WARN: pgen read failed ({e}) -> fallback to .raw for chr{c}")
                        outpref2 = tmpdir / f"chr{c}_dos"
                        run_plink_export_raw(c, ids_file, outpref2)
                        raw = outpref2.with_suffix(".raw")
                        iid_chr, X_chr, _ids_order, meta_df = parse_raw_to_numpy(raw)
                        used_raw_meta = meta_df

            else:
                outpref = tmpdir / f"chr{c}_dos"
                log(f"PLINK export .raw chr{c}: n_ids={len(ids_chr)}")
                run_plink_export_raw(c, ids_file, outpref)
                raw = outpref.with_suffix(".raw")
                iid_chr, X_chr, _ids_order, meta_df = parse_raw_to_numpy(raw)
                used_raw_meta = meta_df

            if iid_chr is None or X_chr is None or X_chr.shape[1] == 0:
                continue

            iid_chr = np.asarray(iid_chr, dtype=str)

            if iid_ref is None:
                iid_ref = iid_chr
                y_full = align_y_to_iids(iid_ref=iid_ref, iid_ph=iid_ph, y_ph=y_ph)
                mask = ~np.isnan(y_full)
                y = y_full[mask]
                if len(y) < 10:
                    log("Too few non-missing phenotypes -> skipping")
                    return
                tr_idx, te_idx = first_fold_indices(len(y))
                log(f"After NA drop: n={len(y):,} (train={len(tr_idx):,}, test={len(te_idx):,})")
            else:
                if len(iid_chr) != len(iid_ref) or not np.all(iid_chr == iid_ref):
                    raise RuntimeError("IID mismatch across chromosome exports (unexpected).")


            X_chr = X_chr[mask, :].astype(np.float32, copy=False)
            Xtr_blocks.append(X_chr[tr_idx, :])
            Xte_blocks.append(X_chr[te_idx, :])


            if used_pvar_df is not None:
                score_tables_parts.append(used_pvar_df)
            else:
                if used_raw_meta is None or used_raw_meta.shape[0] == 0:
                    raise RuntimeError(
                        f"chr{c}: had to fallback to .raw but could not parse variant IDs into alleles/pos; "
                        "cannot write a correct scorefile."
                    )
                raw_meta_parts.append(used_raw_meta)

            del iid_chr, X_chr
            gc.collect()

        if iid_ref is None or not Xtr_blocks:
            log("No genotype blocks produced -> skipping")
            return

        Xtr = np.concatenate(Xtr_blocks, axis=1).astype(np.float32, copy=False)
        Xte = np.concatenate(Xte_blocks, axis=1).astype(np.float32, copy=False)
        del Xtr_blocks, Xte_blocks
        gc.collect()


        y_full = align_y_to_iids(iid_ref=iid_ref, iid_ph=iid_ph, y_ph=y_ph)
        y = y_full[mask]
        ytr = y[tr_idx]
        yte = y[te_idx]

        log(f"Final design: Xtr={Xtr.shape} Xte={Xte.shape}")

        br = BayesianRidge(**PRIOR).fit(Xtr, ytr)

        pred = br.predict(Xte)
        r,  p_r    = pearsonr(yte, pred)
        rho, p_rho = spearmanr(yte, pred)
        evs        = explained_variance_score(yte, pred)
        r2         = float(r ** 2)







        parts = []
        coef_offset = 0


        pvar_i = 0
        raw_i = 0
        for c in sorted(by_chr.keys()):
            ids_chr = by_chr[c]
            if not ids_chr:
                continue




            use_pvar = None
            if pvar_i < len(score_tables_parts):
                cand = score_tables_parts[pvar_i]

                try:
                    cand_chr = str(cand.iloc[0]["CHROM"]).replace("chr", "")
                    if int(cand_chr) == int(c):
                        use_pvar = cand
                        pvar_i += 1
                except Exception:
                    pass

            if use_pvar is not None:
                v_ct = len(use_pvar)
                coefs_chr = br.coef_[coef_offset:coef_offset + v_ct]
                coef_offset += v_ct
                parts.append(make_score_table_from_pvar(use_pvar, coefs_chr))
            else:
                if raw_i >= len(raw_meta_parts):
                    raise RuntimeError("Internal mismatch: expected raw_meta part but none left.")
                meta = raw_meta_parts[raw_i]
                raw_i += 1
                v_ct = len(meta)
                coefs_chr = br.coef_[coef_offset:coef_offset + v_ct]
                coef_offset += v_ct
                parts.append(make_score_table_from_raw_meta(meta, coefs_chr))

        var_table = pd.concat(parts, axis=0, ignore_index=True)
        if var_table.shape[0] != len(br.coef_):
            raise RuntimeError(f"Score table rows {var_table.shape[0]} != number of coefficients {len(br.coef_)}")

        score_path = SCORE_DIR / f"{pgs_id}.txt"
        header = [
            "##OMICSPRED SCORE INFORMATION",
            f"#omicspred_id={pgs_id}",
            f"#pgs_id={pgs_id}",
            f"#pgs_name={prs_name}",
            "#trait_type=metabolomics",
            "#measurement_tissue=serum",
            "#measurement_platform=Nightingale",
            f"#trait_reported={metab}",
            "#genome_build=GRCh38",
            f"#variants_number={len(var_table)}",
            "##SOURCE INFORMATION",
            "#citation=This study (MCPS GWAS hits; BR fixed priors)",
            "#license=CC BY",
        ]
        write_scorefile(score_path, header, var_table)


        model_path = MODEL_DIR / f"{pgs_id}_BR_fixed.pkl"
        stats_path = STAT_DIR  / f"{pgs_id}_fixed_summary.tsv"
        pred_path  = PRED_DIR  / f"{pgs_id}_heldout.tsv.gz"

        joblib.dump(br, model_path)

        pd.DataFrame([{
            "opgs_id": pgs_id,
            "pgs_name": prs_name,
            "metabolite": metab,
            "base": base,
            "n_variants": int(Xtr.shape[1]),
            "n_train": int(len(tr_idx)),
            "n_test": int(len(te_idx)),
            **PRIOR,
            "outer_r": float(r),
            "outer_p_r": float(p_r),
            "outer_r2": float(r2),
            "outer_rho": float(rho),
            "outer_p_rho": float(p_rho),
            "outer_evs": float(evs),
        }]).to_csv(stats_path, sep="\t", index=False)

        heldout_iids = iid_ref[mask][te_idx]
        pd.DataFrame({
            "IID": heldout_iids.astype(str),
            "opgs_id": pgs_id,
            "pgs_name": prs_name,
            "metabolite": metab,
            "y_true": yte.astype(np.float32),
            "y_pred": pred.astype(np.float32),
        }).to_csv(pred_path, sep="\t", index=False, compression="gzip")

        log(f"Done: R²={r2:.4f}  ρ={float(rho):.4f}  | wrote score/model/stats/preds")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
        gc.collect()

if __name__ == "__main__":
    main()
