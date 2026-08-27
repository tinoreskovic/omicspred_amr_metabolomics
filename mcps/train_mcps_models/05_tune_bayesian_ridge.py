#!/usr/bin/env python3
import os, sys, gc, shutil, subprocess, tempfile, datetime, warnings
from pathlib import Path
from typing import Optional, Tuple, List, Dict

import numpy as np
import pandas as pd
from sklearn.linear_model import BayesianRidge
from sklearn.model_selection import KFold, RandomizedSearchCV, cross_val_score
from sklearn.metrics import r2_score
from scipy.stats import spearmanr, loguniform

warnings.filterwarnings("ignore", message="overflow encountered", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="divide by zero encountered", category=RuntimeWarning)


WK = Path(os.environ["OMICSPRED_MCPS_ROOT"]).expanduser()
OUTROOT = WK / "training_inputs"
PRUNED    = OUTROOT / "pruned"
KEEP_SAMP = OUTROOT / "samples.keep"
VARMAP    = OUTROOT / "variant_id_map.tsv"
PHENO     = OUTROOT / "metabolite_phenos.tsv.gz"
LIST_FILE = OUTROOT / "tuning_list.txt"
TRAITMAP  = OUTROOT / "trait_map.tsv"

OUT_DIR   = OUTROOT / "tuning_results"
OUT_DIR.mkdir(exist_ok=True, parents=True)

GENO_DIR = Path(os.environ["MCPS_GENOTYPE_DIR"]).expanduser()
PFILE_PREF = GENO_DIR / "mcps-freeze150k_qcd_chr{chr}"

PREFER_PGEN = True

THREADS = int(os.environ.get("N_THREADS", 1))

def log(msg: str):
    print(f"[{datetime.datetime.now():%F %T}] {msg}", flush=True)


def pick_metab() -> str:
    metas = [ln.strip() for ln in LIST_FILE.read_text().splitlines() if ln.strip()]
    if not metas:
        raise SystemExit("tuning_list.txt empty")

    if len(sys.argv) > 1:
        arg = sys.argv[1]
        try:
            idx = int(arg)
            return metas[idx]
        except ValueError:
            return arg.strip()

    if "TASK_INDEX" in os.environ:
        return metas[int(os.environ["TASK_INDEX"])]
    return metas[0]

def load_traitmap_base(metab: str) -> str:
    tm = pd.read_csv(TRAITMAP, sep="\t", dtype=str)
    hit = tm.loc[tm["NMR_name"] == metab]
    if hit.shape[0] == 0:
        raise FileNotFoundError(f"{metab} not found in trait_map.tsv")
    return str(hit.iloc[0]["base"])

def build_id_to_chr() -> dict:
    vm = pd.read_csv(VARMAP, sep="\t", dtype={"ID": str})
    return dict(zip(vm["ID"].astype(str), vm["chr"].astype(int)))

def read_ids(path: Path):
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]

def load_pheno_arrays(metab: str) -> Tuple[np.ndarray, np.ndarray]:
    ph = pd.read_csv(PHENO, sep="\t", usecols=["IID", metab])
    iid = ph["IID"].astype(str).to_numpy()
    y = ph[metab].to_numpy(dtype=np.float32)
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


def _plink_common_args() -> list[str]:
    extra = []
    if "N_THREADS" in os.environ:
        extra += ["--threads", os.environ["N_THREADS"]]
    if "MEMORY_MB" in os.environ:
        extra += ["--memory", os.environ["MEMORY_MB"]]
    return extra

def run_plink_make_pgen(chr_num: int, ids_file: Path, outpref: Path):
    cmd = ["plink2", "--pfile", str(PFILE_PREF).format(chr=chr_num),
           "--keep", str(KEEP_SAMP), "--extract", str(ids_file),
           "--make-pgen", "--out", str(outpref), "--silent"] + _plink_common_args()
    subprocess.run(cmd, check=True)

def run_plink_export_raw(chr_num: int, ids_file: Path, outpref: Path):
    cmd = ["plink2", "--pfile", str(PFILE_PREF).format(chr=chr_num),
           "--keep", str(KEEP_SAMP), "--extract", str(ids_file),
           "--export", "A", "--out", str(outpref), "--silent"] + _plink_common_args()
    subprocess.run(cmd, check=True)

def _find_header_line(path: Path, required: Tuple[str, ...], max_lines: int = 200):
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for i in range(max_lines):
            line = f.readline()
            if not line: break
            line = line.rstrip("\n")
            if not line or line.startswith("##"): continue
            toks = line.split("\t")
            norm = [t.lstrip("#") for t in toks]
            if all(r in norm for r in required): return i, norm
    raise RuntimeError("Header not found")

def read_psam_iids(psam_path: Path) -> np.ndarray:
    skip, _ = _find_header_line(psam_path, required=("FID", "IID"))
    df = pd.read_csv(psam_path, sep="\t", header=0, skiprows=skip, dtype=str, engine="python", comment=None)
    df.columns = [c.lstrip("#") for c in df.columns]
    return df["IID"].astype(str).to_numpy()

def try_import_pgenlib():
    try:
        import pgenlib
        return pgenlib
    except Exception:
        return None

def read_pgen_as_dosage(pgenlib_mod, pgen_path: Path, n_samples: int) -> np.ndarray:
    reader = pgenlib_mod.PgenReader(bytes(str(pgen_path), "utf-8"))
    v_ct = reader.get_variant_ct()
    X = np.empty((n_samples, v_ct), dtype=np.float32)

    dosage_methods = [m for m in ["read_dosages", "read_dosage", "read_dosage16"] if hasattr(reader, m)]
    if dosage_methods:
        meth = getattr(reader, dosage_methods[0])
        buf = np.empty(n_samples, dtype=np.float32)
        for j in range(v_ct):
            try: meth(buf, j)
            except: meth(j, buf)
            X[:, j] = buf
        return X

    hard_methods = [m for m in ["read", "read_hardcalls", "read_alleles"] if hasattr(reader, m)]
    if hard_methods:
        meth = getattr(reader, hard_methods[0])
        buf_i8 = np.empty(n_samples, dtype=np.int8)
        for j in range(v_ct):
            try: meth(buf_i8, j)
            except: meth(j, buf_i8)
            col = buf_i8.astype(np.float32)
            col[col < 0] = np.nan
            X[:, j] = col
        return X
    raise RuntimeError("No usable pgenlib read method found")

def parse_raw_to_numpy(raw_path: Path):
    if not raw_path.exists(): return None, None
    df = pd.read_csv(raw_path, sep=r"\s+", engine="python")
    if df.shape[1] <= 6: return None, None
    return df["IID"].astype(str).to_numpy(), df.iloc[:, 6:].to_numpy(dtype=np.float32, copy=False)


def main():
    metab = pick_metab()
    log(f"Metabolite (NMR_name): {metab}")

    base = load_traitmap_base(metab)
    pruned_keep = PRUNED / f"{base}__variants_pruned.keep"
    if not pruned_keep.exists():
        log(f"Missing keepfile: {pruned_keep} -> skipping")
        return

    ids_all = read_ids(pruned_keep)
    id_to_chr_lookup = build_id_to_chr()
    by_chr: Dict[int, List[str]] = {}
    for vid in ids_all:
        c = id_to_chr_lookup.get(str(vid))
        if c is not None:
            by_chr.setdefault(int(c), []).append(str(vid))

    iid_ph, y_ph = load_pheno_arrays(metab)
    pgenlib_mod = try_import_pgenlib() if PREFER_PGEN else None

    tmpdir = Path(tempfile.mkdtemp(prefix=f"gws_tune_{base[:60]}_"))
    try:
        iid_ref: Optional[np.ndarray] = None
        mask = tr_idx = te_idx = None
        Xtr_blocks, Xte_blocks = [], []

        for c in sorted(by_chr.keys()):
            ids_chr = by_chr[c]
            if not ids_chr: continue

            ids_file = tmpdir / f"chr{c}.ids"
            ids_file.write_text("\n".join(ids_chr) + "\n")

            if pgenlib_mod is not None:
                outpref = tmpdir / f"chr{c}_subset"
                run_plink_make_pgen(c, ids_file, outpref)
                pgen, psam, pvar = outpref.with_suffix(".pgen"), outpref.with_suffix(".psam"), outpref.with_suffix(".pvar")

                if pgen.exists() and psam.exists() and pvar.exists():
                    iid_chr = read_psam_iids(psam)
                    try:
                        X_chr = read_pgen_as_dosage(pgenlib_mod, pgen, len(iid_chr))
                    except Exception as e:
                        log(f"Fallback .raw chr{c} ({e})")
                        run_plink_export_raw(c, ids_file, tmpdir / f"chr{c}_dos")
                        iid_chr, X_chr = parse_raw_to_numpy(tmpdir / f"chr{c}_dos.raw")
                else:
                    run_plink_export_raw(c, ids_file, tmpdir / f"chr{c}_dos")
                    iid_chr, X_chr = parse_raw_to_numpy(tmpdir / f"chr{c}_dos.raw")
            else:
                run_plink_export_raw(c, ids_file, tmpdir / f"chr{c}_dos")
                iid_chr, X_chr = parse_raw_to_numpy(tmpdir / f"chr{c}_dos.raw")

            if iid_chr is None or X_chr is None or X_chr.shape[1] == 0: continue
            iid_chr = np.asarray(iid_chr, dtype=str)

            if iid_ref is None:
                iid_ref = iid_chr
                y_full = align_y_to_iids(iid_ref=iid_ref, iid_ph=iid_ph, y_ph=y_ph)
                mask = ~np.isnan(y_full)
                y = y_full[mask]

                kf_outer = KFold(n_splits=5, shuffle=True, random_state=21)
                tr_idx, te_idx = next(iter(kf_outer.split(y)))

            X_chr = X_chr[mask, :].astype(np.float32, copy=False)
            Xtr_blocks.append(X_chr[tr_idx, :])
            Xte_blocks.append(X_chr[te_idx, :])
            del iid_chr, X_chr
            gc.collect()

        if iid_ref is None or not Xtr_blocks: return

        Xtr = np.concatenate(Xtr_blocks, axis=1).astype(np.float32, copy=False)
        Xte = np.concatenate(Xte_blocks, axis=1).astype(np.float32, copy=False)
        del Xtr_blocks, Xte_blocks
        gc.collect()

        y_full = align_y_to_iids(iid_ref=iid_ref, iid_ph=iid_ph, y_ph=y_ph)
        y = y_full[mask]
        ytr, yte = y[tr_idx], y[te_idx]

        log(f"Matrix ready: Xtr={Xtr.shape} Xte={Xte.shape}")


        param_dist = {
            'alpha_1': loguniform(1e-10, 1e10),
            'alpha_2': loguniform(1e-10, 1e10),
            'lambda_1': loguniform(1e-10, 1e10),
            'lambda_2': loguniform(1e-10, 1e10)
        }

        inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)

        log(f"Starting RandomizedSearchCV (200 iters) across {THREADS} CPUs...")
        search = RandomizedSearchCV(
            estimator=BayesianRidge(),
            param_distributions=param_dist,
            n_iter=200,
            cv=inner_cv,
            scoring='r2',
            n_jobs=THREADS,
            pre_dispatch='1.5*n_jobs',
            verbose=1,
            random_state=42,
            return_train_score=False
        )
        search.fit(Xtr, ytr)



        res_df = pd.DataFrame(search.cv_results_)
        rename_map = {
            'param_alpha_1': 'alpha_1', 'param_alpha_2': 'alpha_2',
            'param_lambda_1': 'lambda_1', 'param_lambda_2': 'lambda_2',
            'mean_test_score': 'mean_inner_R2'
        }
        tune_df = res_df[list(rename_map.keys())].rename(columns=rename_map)


        pred_best = search.predict(Xte)
        outer_r2_best = r2_score(yte, pred_best)
        outer_sp_best = float(spearmanr(yte, pred_best).statistic)
        best_p = search.best_params_


        log("Evaluating fixed 1e-5 baseline...")
        model_fixed = BayesianRidge(alpha_1=1e-5, alpha_2=1e-5, lambda_1=1e-5, lambda_2=1e-5)


        scores_fixed = cross_val_score(model_fixed, Xtr, ytr, cv=inner_cv, scoring='r2', n_jobs=THREADS)
        inner_r2_fixed = float(np.mean(scores_fixed))


        model_fixed.fit(Xtr, ytr)
        pred_fixed = model_fixed.predict(Xte)
        outer_r2_fixed = r2_score(yte, pred_fixed)
        outer_sp_fixed = float(spearmanr(yte, pred_fixed).statistic)


        fixed_row = pd.DataFrame([{
            'alpha_1': 1e-5, 'alpha_2': 1e-5, 'lambda_1': 1e-5, 'lambda_2': 1e-5,
            'mean_inner_R2': inner_r2_fixed
        }])
        tune_df = pd.concat([tune_df, fixed_row], ignore_index=True)

        tune_out = OUT_DIR / f"{metab}_tune.tsv"
        tune_df.to_csv(tune_out, sep="\t", index=False)


        summ = pd.DataFrame([
            {
                "metabolite": metab, "type": "best", "base": base,
                "alpha_1": best_p['alpha_1'], "alpha_2": best_p['alpha_2'],
                "lambda_1": best_p['lambda_1'], "lambda_2": best_p['lambda_2'],
                "inner_R2": float(search.best_score_),
                "outer_R2": float(outer_r2_best),
                "outer_spearman": float(outer_sp_best),
                "n_samples": int(len(y)), "n_snps": int(Xtr.shape[1]),
                "n_train": int(len(tr_idx)), "n_test": int(len(te_idx))
            },
            {
                "metabolite": metab, "type": "fixed", "base": base,
                "alpha_1": 1e-5, "alpha_2": 1e-5,
                "lambda_1": 1e-5, "lambda_2": 1e-5,
                "inner_R2": float(inner_r2_fixed),
                "outer_R2": float(outer_r2_fixed),
                "outer_spearman": float(outer_sp_fixed),
                "n_samples": int(len(y)), "n_snps": int(Xtr.shape[1]),
                "n_train": int(len(tr_idx)), "n_test": int(len(te_idx))
            }
        ])
        summ_out = OUT_DIR / f"{metab}_tune_summary.tsv"
        summ.to_csv(summ_out, sep="\t", index=False)

        log(f"Done! Best R2_inner={search.best_score_:.4f}, Fixed R2_inner={inner_r2_fixed:.4f}")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
        gc.collect()

if __name__ == "__main__":
    main()
