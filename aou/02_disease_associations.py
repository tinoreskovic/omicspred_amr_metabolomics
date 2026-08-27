#!/usr/bin/env python3
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from google.cloud import bigquery
from statsmodels.stats.multitest import multipletests

from phetk.cohort import Cohort
from phetk.phecode import Phecode


TARGET_PHECODES = {
    "CV_404": "Ischemic Heart Disease",
    "EM_202.2": "Type 2 Diabetes",
    "GU_582.2": "Chronic Kidney Disease",
}

METADATA_FILE = Path(__file__).resolve().parents[1] / "metadata" / "omicspred_validation.xlsx"


def require_env(name: str) -> str:
    value = os.environ.get(name, "")
    if not value:
        raise RuntimeError(f"Set environment variable {name}")
    return value


def run_model(
    cohort: pd.DataFrame,
    cases: set[str],
    scores: pd.DataFrame,
    pgs_id: str,
    score_label: str,
    phecode: str,
) -> dict | None:
    score = scores.loc[scores["PGS"] == pgs_id, ["IID", "SUM"]].rename(
        columns={"IID": "person_id", "SUM": "score"}
    )
    score["person_id"] = score["person_id"].astype(str)
    data = cohort.merge(score, on="person_id", how="inner")
    data["status"] = data["person_id"].isin(cases).astype(int)
    data["sex_cov"] = pd.to_numeric(data["sex_at_birth"], errors="coerce").fillna(0).astype(int)
    pc_columns = [
        column
        for column in data
        if column.lower().startswith("principal_component")
        or column.lower().startswith("pc")
    ]
    columns = ["score", "age_at_last_ehr_event", "sex_cov", *pc_columns]
    model_data = data[["status", *columns]].apply(pd.to_numeric, errors="coerce").dropna()
    if model_data["status"].sum() < 20:
        return None
    model_data["score"] = (model_data["score"] - model_data["score"].mean()) / model_data["score"].std()
    for column in ["age_at_last_ehr_event", *pc_columns]:
        if model_data[column].std() > 0:
            model_data[column] = (model_data[column] - model_data[column].mean()) / model_data[column].std()
    usable = [column for column in columns if model_data[column].nunique() > 1]
    fit = sm.Logit(
        model_data["status"], sm.add_constant(model_data[usable])
    ).fit(disp=False, maxiter=100)
    beta = fit.params["score"]
    standard_error = fit.bse["score"]
    return {
        "Disease": phecode,
        "PGS": pgs_id,
        "score_set": score_label,
        "model": "Logistic",
        "n_cases": int(model_data["status"].sum()),
        "n_controls": int((model_data["status"] == 0).sum()),
        "z": float(fit.tvalues["score"]),
        "p": float(fit.pvalues["score"]),
        "estimate": float(np.exp(beta)),
        "L95": float(np.exp(beta - 1.96 * standard_error)),
        "U95": float(np.exp(beta + 1.96 * standard_error)),
    }


def main() -> None:
    root = Path(require_env("AOU_OMICS_PRED_ROOT"))
    output_dir = Path(require_env("AOU_DISEASE_RESULTS_DIR"))
    output_dir.mkdir(parents=True, exist_ok=True)
    project = require_env("GOOGLE_PROJECT")
    cdr = require_env("WORKSPACE_CDR")
    flagged_samples_file = require_env("AOU_FLAGGED_SAMPLES_FILE")

    ancestry = pd.read_csv(root / "amr_ids.tsv", sep="	", dtype=str)
    flagged = pd.read_csv(flagged_samples_file, sep="	", dtype=str)
    cohort_ids = ancestry.loc[
        ~ancestry["research_id"].astype(str).isin(flagged["s"].astype(str)),
        "research_id",
    ].astype(str).drop_duplicates()

    with tempfile.TemporaryDirectory(prefix="aou_phetk_") as temporary:
        temporary_path = Path(temporary)
        cohort_path = temporary_path / "cohort.tsv"
        covariate_path = temporary_path / "cohort_covariates.tsv"
        pd.DataFrame({"person_id": cohort_ids, "ancestry": "AMR"}).to_csv(cohort_path, sep="	", index=False)

        Cohort(platform="aou").add_covariates(
            cohort_file_path=str(cohort_path),
            date_of_birth=True,
            current_age=True,
            sex_at_birth=True,
            age_at_last_ehr_event=True,
            first_n_pcs=10,
            output_file_path=str(covariate_path),
        )
        cohort = pd.read_csv(covariate_path, sep="	", dtype={"person_id": str})
        cohort["age_at_last_ehr_event"] = pd.to_numeric(cohort["age_at_last_ehr_event"], errors="coerce")
        cohort = cohort.loc[cohort["age_at_last_ehr_event"] >= 35].copy()

        query = f"""
        SELECT CAST(co.person_id AS STRING) AS person_id,
               co.condition_start_date AS date,
               co.condition_source_value AS ICD,
               concept.vocabulary_id
        FROM `{cdr}.condition_occurrence` AS co
        JOIN `{cdr}.concept` AS concept
          ON co.condition_source_concept_id = concept.concept_id
        JOIN UNNEST(@cohort_ids) AS cohort_id
          ON CAST(co.person_id AS STRING) = cohort_id
        WHERE concept.vocabulary_id IN ('ICD9CM', 'ICD10CM')
        """
        job_config = bigquery.QueryJobConfig(
            query_parameters=[bigquery.ArrayQueryParameter("cohort_ids", "STRING", cohort["person_id"].tolist())]
        )
        icd = bigquery.Client(project=project).query(query, job_config=job_config).to_dataframe()
        icd_path = temporary_path / "icd.tsv"
        icd.to_csv(icd_path, sep="	", index=False)
        counts_path = temporary_path / "phecode_counts.tsv"
        module = Phecode(platform="custom", icd_file_path=str(icd_path))
        module.cdr = cdr
        module.project = project
        module.count_phecode(phecode_version="X", output_file_path=str(counts_path))
        counts = pd.read_csv(counts_path, sep="	", dtype={"person_id": str})

    cases = counts.groupby("phecode")["person_id"].apply(set).to_dict()
    interval_scores = pd.read_csv(root / "filtered_scores" / "OmicsPred_AMR_scores.txt.gz", sep="	", dtype=str)
    mcps_scores = pd.read_csv(root / "filtered_scores" / "MCPS_AMR_scores.txt.gz", sep="	", dtype=str)
    metadata = pd.read_excel(METADATA_FILE, sheet_name="Table S2")[[
        "PRS_name",
        "Biomarker.Name",
        "Group",
        "Subgroup",
    ]]
    expected_ids = metadata["PRS_name"].dropna().astype(str).tolist()
    interval_ids = set(interval_scores["PGS"].dropna())
    mcps_ids = set(mcps_scores["PGS"].dropna())
    missing_interval = [pgs for pgs in expected_ids if pgs not in interval_ids]
    missing_mcps = [
        f"{pgs}_MCPS_80_20_fixed"
        for pgs in expected_ids
        if f"{pgs}_MCPS_80_20_fixed" not in mcps_ids
    ]
    if missing_interval or missing_mcps:
        raise RuntimeError(
            f"Missing {len(missing_interval)} INTERVAL-trained and {len(missing_mcps)} MCPS-trained scores"
        )
    pairs = [(pgs, f"{pgs}_MCPS_80_20_fixed") for pgs in expected_ids]

    records = []
    for phecode in TARGET_PHECODES:
        for interval_id, mcps_id in pairs:
            for result in (
                run_model(cohort, cases.get(phecode, set()), interval_scores, interval_id, "INTERVAL-trained", phecode),
                run_model(cohort, cases.get(phecode, set()), mcps_scores, mcps_id, "MCPS-trained", phecode),
            ):
                if result is not None:
                    records.append(result)

    results = pd.DataFrame.from_records(records)
    results["Disease_Name"] = results["Disease"].map(TARGET_PHECODES)
    results["fdr"] = np.nan
    for _, index in results.groupby(["Disease_Name", "score_set"]).groups.items():
        results.loc[index, "fdr"] = multipletests(results.loc[index, "p"], method="fdr_bh")[1]
    results["PRS_base"] = results["PGS"].str.replace(r"_MCPS_80_20_fixed$", "", regex=True)
    results = results.merge(metadata, left_on="PRS_base", right_on="PRS_name", how="left").drop(columns="PRS_name")
    results.to_csv(output_dir / "aou_disease_associations.csv", index=False)


if __name__ == "__main__":
    main()
