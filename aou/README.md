# All of Us analyses

These scripts are run in the All of Us Researcher Workbench and use the files beneath `AOU_OMICS_PRED_ROOT`.

- `01_evaluate_metabolite_models.py` externally validates the INTERVAL-trained and MCPS-trained genetic models for the six Nightingale NMR metabolomic traits used in the manuscript and records both R² and Spearman rho.
- `02_disease_associations.py` tests associations of the predicted metabolomic traits with ischemic heart disease, type 2 diabetes and chronic kidney disease, as used in Figure 5 and Supplementary Figures 9–12.

The disease-association script reads model and Nightingale trait names from the metadata workbook included with the repository. Both scripts use precomputed genetic-score files. Genetic score imputation used `pgsc_calc` 2.0.0 (pipeline revision `a10e219a42`) with Nextflow 24.10.2.

`02_disease_associations.py` uses the `GOOGLE_PROJECT` and `WORKSPACE_CDR` variables and requires PheTK in the selected Workbench environment. The AMR ID file, filtered genetic-score files and measured-trait files remain beneath `AOU_OMICS_PRED_ROOT`; the platform QC file is set with `AOU_FLAGGED_SAMPLES_FILE`, and disease-association results are written to `AOU_DISEASE_RESULTS_DIR`. See [`../config/env.example`](../config/env.example) for the path settings.

PheTK 0.3.6 defaults to v9 genetic PCs, and `02_disease_associations.py` calls `Cohort(platform="aou")`. To select v8 PCs, change that call to `Cohort(platform="aou", aou_db_version=8)` in the covariate-preparation step. `WORKSPACE_CDR` selects the clinical dataset separately.

## Environment and installation

Python 3.12.14 on Ubuntu 24.04.5 LTS, x86_64. The Python package versions are specified in [`../environment.yml`](../environment.yml). Run the following commands in a Researcher Workbench terminal with Python 3.12 available:

```bash
python3.12 -m venv "$HOME/.venvs/omicspred-aou"
source "$HOME/.venvs/omicspred-aou/bin/activate"
python -m pip install \
  "numpy==1.26.4" \
  "pandas==2.3.3" \
  "scipy==1.17.1" \
  "statsmodels==0.15.0" \
  "matplotlib==3.11.2" \
  "seaborn==0.13.2" \
  "phetk==0.3.6" \
  "openpyxl==3.1.5" \
  "pyarrow==24.0.0" \
  "db-dtypes==1.4.3" \
  "google-cloud-bigquery==3.38.0" \
  "google-cloud-bigquery-storage==2.34.0" \
  "google-cloud-storage==3.4.1" \
  "gcsfs==2025.9.0" \
  "google-api-core==2.25.1" \
  "google-auth==2.40.3" \
  "protobuf==6.32.0" \
  "dsub==0.5.6" \
  "ipython==9.17.1" \
  "ipykernel==7.4.0"
python -m pip check
```

To use this environment in Jupyter, register its kernel:

```bash
python -m ipykernel install --user --name omicspred-aou --display-name "Python (omicspred-aou)"
```

Select **Python (omicspred-aou)** in Jupyter. The same pinned Python packages are provided in `environment.yml` for a Conda-based setup.

## Running the scripts

Prepare the required input files using the study's cohort, data release and phenotype mappings. Export the path settings above, then run these commands from the repository root with the project environment activated:

```bash
python aou/01_evaluate_metabolite_models.py
python aou/02_disease_associations.py
```

The evaluation script writes its results beneath `AOU_OMICS_PRED_ROOT/aggregate_results`; the disease-association script writes to `AOU_DISEASE_RESULTS_DIR`.
