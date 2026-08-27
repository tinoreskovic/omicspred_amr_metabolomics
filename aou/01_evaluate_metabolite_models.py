#!/usr/bin/env python3
"""Evaluate the six manuscript traits inside the All of Us Workbench.

The script applies age, sex, and PC adjustment and reports R-squared and
Spearman rho with their p-values.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import betaln, hyp2f1
from scipy.stats import f as f_dist
from scipy.stats import norm, rankdata, spearmanr



AOU_ZERO_JITTER_SEED = 20260821

TRAIT_TO_OPGS = {
    "HDL_cholesterol": "OPGS003520",
    "LDL_cholesterol": "OPGS003423",
    "VLDL_cholesterol": "OPGS003488",
    "Triglycerides": "OPGS003486",
    "Creatinine": "OPGS003552",
    "Glucose": "OPGS003557",
}

SIMPLE_NAMES = {
    "HDL_cholesterol": "HDL",
    "LDL_cholesterol": "LDL",
    "VLDL_cholesterol": "VLDL",
    "Triglycerides": "Triglycerides",
    "Creatinine": "Creatinine",
    "Glucose": "Glucose",
}

SUBGROUPS = ("All", "Male", "Female", "Age 35-49", "Age 50-64", "Age >= 65")
COVARIATE_COLUMNS = ["sex_at_birth", "age_at_measurement"] + [f"pc{i}" for i in range(1, 11)]

def pvalue_from_squared_correlation(r2: float, n: int) -> tuple[str | None, float]:
    """Return scientific-notation p text and its uncapped log10 value."""
    if not np.isfinite(r2) or n < 3 or r2 < 0:
        return None, np.nan
    if r2 >= 1.0:
        log_p = -np.inf
    else:
        f_stat = r2 * (n - 2) / (1.0 - r2)
        p_direct = float(f_dist.sf(f_stat, 1, n - 2))
        if np.isfinite(p_direct) and p_direct > 0.0:
            log_p = float(np.log(p_direct))
        else:


            x_beta = 1.0 - r2
            a_beta = (n - 2) / 2.0
            hypergeom = hyp2f1(a_beta, 0.5, a_beta + 1.0, x_beta)
            log_p = float(
                a_beta * np.log(x_beta)
                - np.log(a_beta)
                - betaln(a_beta, 0.5)
                + np.log(hypergeom)
            )
    log10_p = float(log_p / np.log(10.0))
    if np.isneginf(log_p):
        p_text = "0"
    else:
        exponent = int(np.floor(log10_p))
        mantissa = 10.0 ** (log10_p - exponent)
        p_text = f"{mantissa:.15f}e{exponent}"
    return p_text, log10_p

def find_omics_root() -> Path:
    configured = os.environ.get("AOU_OMICS_PRED_ROOT", "")
    if not configured:
        raise FileNotFoundError("Set AOU_OMICS_PRED_ROOT")
    path = Path(configured).expanduser()
    if not path.is_dir():
        raise FileNotFoundError(f"AOU_OMICS_PRED_ROOT does not exist: {path}")
    return path

def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    return path

def inverse_normal(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    output = pd.Series(np.nan, index=series.index, dtype=float)
    mask = values.notna() & np.isfinite(values)
    if not mask.any():
        return output
    ranks = rankdata(values.loc[mask], method="average")
    output.loc[mask] = norm.ppf((ranks - 0.5) / mask.sum())
    return output

def log_transform_positive(series: pd.Series, rng: np.random.Generator) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").copy()
    positive = values[(values > 0) & np.isfinite(values)]
    if positive.empty:
        return pd.Series(np.nan, index=series.index, dtype=float)
    zeros = values == 0
    if zeros.any():
        values.loc[zeros] = rng.uniform(0.001, 0.9, zeros.sum()) * positive.min()
    values.loc[~np.isfinite(values) | (values <= 0)] = np.nan
    return np.log(values)

def covariates(df: pd.DataFrame) -> pd.DataFrame:
    missing = [column for column in COVARIATE_COLUMNS if column not in df.columns]
    if missing:
        raise KeyError(f"Measured-trait input lacks covariates: {', '.join(missing)}")
    out = pd.DataFrame(index=df.index)
    out["FEMALE"] = pd.to_numeric(df["sex_at_birth"], errors="coerce")
    out["AGE"] = pd.to_numeric(df["age_at_measurement"], errors="coerce")
    for index in range(1, 11):
        out[f"pc{index}"] = pd.to_numeric(df[f"pc{index}"], errors="coerce")
    return out

def residualize(y: pd.Series, x: pd.DataFrame) -> pd.Series:
    numeric_y = pd.to_numeric(y, errors="coerce")
    numeric_x = x.apply(pd.to_numeric, errors="coerce")
    model_data = pd.concat([numeric_y.rename("y"), numeric_x], axis=1).dropna()
    output = pd.Series(np.nan, index=y.index, dtype=float)
    if model_data.shape[0] <= numeric_x.shape[1] + 1:
        return output
    usable = [column for column in numeric_x if model_data[column].nunique() > 1]
    design = sm.add_constant(model_data[usable], has_constant="add")
    output.loc[model_data.index] = sm.OLS(model_data["y"], design).fit().resid
    return output

def subgroup_mask(label: str, age: pd.Series, sex: pd.Series) -> pd.Series:
    if label == "All":
        return pd.Series(True, index=age.index)
    if label == "Male":
        return sex == 0
    if label == "Female":
        return sex == 1
    if label == "Age 35-49":
        return (age >= 35) & (age < 50)
    if label == "Age 50-64":
        return (age >= 50) & (age < 65)
    if label == "Age >= 65":
        return age >= 65
    raise ValueError(f"Unknown subgroup: {label}")

def load_scores(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_dir = root / "filtered_scores"
    mcps = pd.read_csv(
        require_file(score_dir / "MCPS_AMR_scores.txt.gz"),
        sep="\t",
        compression="gzip",
        dtype=str,
    )
    interval = pd.read_csv(
        require_file(score_dir / "OmicsPred_AMR_scores.txt.gz"),
        sep="\t",
        compression="gzip",
        dtype=str,
    )
    required = {"IID", "PGS", "SUM"}
    for label, frame in (("MCPS", mcps), ("INTERVAL", interval)):
        missing = required.difference(frame.columns)
        if missing:
            raise KeyError(f"{label} score file lacks columns: {sorted(missing)}")
    return mcps, interval

def find_mcps_pgs(mcps_scores: pd.DataFrame, opgs_id: str) -> str:
    available = set(mcps_scores["PGS"].dropna().astype(str))
    pgs_id = f"{opgs_id}_MCPS_80_20_fixed"
    if pgs_id not in available:
        raise KeyError(f"No MCPS score was found for {pgs_id}")
    return pgs_id

def prepare_trait(
    root: Path,
    trait: str,
    opgs_id: str,
    mcps_scores: pd.DataFrame,
    interval_scores: pd.DataFrame,
) -> tuple[pd.DataFrame, str, str]:
    measured = pd.read_csv(
        require_file(root / "measured_traits" / f"{trait}.csv.gz"),
        compression="gzip",
        dtype={"person_id": str},
    )
    measured["person_id"] = measured["person_id"].astype(str)
    mcps_pgs = find_mcps_pgs(mcps_scores, opgs_id)
    interval_column = f"{opgs_id}_SUM_INTERVAL"
    mcps_column = f"{mcps_pgs}_SUM_MCPS"

    interval_part = interval_scores.loc[
        interval_scores["PGS"] == opgs_id, ["IID", "SUM"]
    ].rename(columns={"IID": "person_id", "SUM": interval_column})
    mcps_part = mcps_scores.loc[
        mcps_scores["PGS"] == mcps_pgs, ["IID", "SUM"]
    ].rename(columns={"IID": "person_id", "SUM": mcps_column})
    interval_part["person_id"] = interval_part["person_id"].astype(str)
    mcps_part["person_id"] = mcps_part["person_id"].astype(str)

    merged = measured.merge(interval_part, on="person_id", how="left").merge(
        mcps_part, on="person_id", how="left"
    )
    merged = merged.dropna(subset=[interval_column, mcps_column])
    return merged, interval_column, mcps_column

def compute_results(root: Path, rng: np.random.Generator) -> pd.DataFrame:
    mcps_scores, interval_scores = load_scores(root)
    records = []

    for trait, opgs_id in TRAIT_TO_OPGS.items():
        df, interval_column, mcps_column = prepare_trait(
            root, trait, opgs_id, mcps_scores, interval_scores
        )
        age = pd.to_numeric(df["age_at_measurement"], errors="coerce")
        df = df.loc[age >= 35].copy()
        if df.empty:
            continue
        age = pd.to_numeric(df["age_at_measurement"], errors="coerce")
        sex = pd.to_numeric(df["sex_at_birth"], errors="coerce")
        x_covariates = covariates(df)
        measured_residual = residualize(
            log_transform_positive(df[trait], rng), x_covariates
        )
        predicted_residuals = {
            "OmicsPred": residualize(df[interval_column], x_covariates),
            "MCPS": residualize(df[mcps_column], x_covariates),
        }

        for subgroup in SUBGROUPS:
            subset = subgroup_mask(subgroup, age, sex)
            measured_rint = inverse_normal(measured_residual.loc[subset])
            for score_type, predicted_residual in predicted_residuals.items():
                predicted_rint = inverse_normal(predicted_residual.loc[subset])
                paired = measured_rint.notna() & predicted_rint.notna()
                n_pairs = int(paired.sum())
                if n_pairs < 3:
                    continue
                x = measured_rint.loc[paired].to_numpy()
                y = predicted_rint.loc[paired].to_numpy()

                r2 = float(np.corrcoef(x, y)[0, 1] ** 2)
                p_r2, log10_p_r2 = pvalue_from_squared_correlation(r2, n_pairs)
                spearman = spearmanr(x, y)
                rho = float(spearman.statistic)
                p_spearman, log10_p_spearman = pvalue_from_squared_correlation(
                    rho**2, n_pairs
                )
                records.append(
                    {
                        "trait": SIMPLE_NAMES[trait],
                        "subgroup": subgroup,
                        "OPGS_ID": opgs_id,
                        "score_type": score_type,
                        "p_R2_RINT": p_r2,
                        "log10_p_R2_RINT": log10_p_r2,
                        "rho_spearman_RINT": rho,
                        "p_spearman_RINT": p_spearman,
                        "log10_p_spearman_RINT": log10_p_spearman,
                        "R2_pearson_RINT": r2,
                        "n_pairs": n_pairs,
                    }
                )
    result = pd.DataFrame.from_records(records)
    if result.empty:
        raise RuntimeError("No AoU trait/model results were produced.")
    return result

def save_outputs(results: pd.DataFrame, root: Path) -> None:
    result_dir = root / "aggregate_results"
    result_dir.mkdir(parents=True, exist_ok=True)
    overall = results.loc[results["subgroup"] == "All"].drop(columns="subgroup")
    overall.to_csv(result_dir / "aou_six_trait_model_performance.csv", index=False)
    results.to_csv(result_dir / "aou_six_trait_model_performance_subgroups.csv", index=False)

def main() -> None:
    root = find_omics_root()
    seed_text = os.environ.get(
        "AOU_ZERO_JITTER_SEED", str(AOU_ZERO_JITTER_SEED)
    )
    seed = int(seed_text)
    rng = np.random.default_rng(seed)
    save_outputs(compute_results(root, rng), root)

if __name__ == "__main__":
    main()
