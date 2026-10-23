# MCPS analyses

The MCPS code is organised around the three parts of the analysis:

1. `validate_interval/` prepares the Nightingale NMR metabolomic traits and validates the INTERVAL Study–trained genetic models in MCPS, overall and in the subgroups referred to in the manuscript.
2. `train_mcps_models/` prepares the training inputs, selects and prunes genome-wide significant variants, tunes the Bayesian ridge prior on a reproducible subset of traits, and trains the new genetic models for prediction/imputation of the 141 Nightingale NMR metabolomic traits.
3. `validate_mcps_models/` evaluates the MCPS-trained and INTERVAL-trained models in the same withheld MCPS subset.

The model-to-trait mapping is read from `../metadata/omicspred_validation.xlsx`. Working inputs and the intermediate files passed between scripts are written beneath `OMICSPRED_MCPS_ROOT`.

Run `validate_interval/01` through `03` first. `validate_interval/05_combine_cohort_results.R` combines the MCPS results with the published external-validation results used for Figure 1. Run `train_mcps_models/01` through `08` in order, distributing steps `05` and `06` over `TASK_INDEX` where appropriate. Finally, run `validate_mcps_models/00` and `01` to compare both sets of genetic models in the withheld subset.

## Software and setup

Rocky Linux 8.10, x86_64; Python 3.11.3 and R 4.3.2 in separate environments. Versions were recorded on 23 October 2026 without rerunning the analyses.

The genetic scoring used `pgsc_calc` 2.0.0-alpha.4 with Nextflow 23.10.1.

### Python training and summaries

The study Python environment uses module `Python/3.11.3-GCCcore-12.3.0`.

| Software | Version |
| --- | --- |
| numpy | 2.4.0 |
| pandas | 2.3.3 |
| scipy | 1.16.3 |
| scikit-learn | 1.8.0 |
| joblib | 1.5.3 |
| pgenlib (optional) | 0.93.0 |
| openpyxl | 3.1.5 |
| et-xmlfile | 2.0.0 |

The workbook-reading steps require openpyxl. After activating the study Python environment, install it with:

```bash
python -m pip install "openpyxl==3.1.5" "et-xmlfile==2.0.0"
```

The workbook-reader imports passed in the study Python environment with the NumPy and pandas versions above.

PLINK `2.00a2.3LM` (24 January 2020) is provided by module `PLINK/2.00a2.3_x86_64`. The tuning and fitting code can use pgenlib for genotype access and otherwise falls back to PLINK text export.

### R preparation, validation and cohort plots

R 4.3.2 is provided by `RStudio-Server/2023.09.1+494-foss-2023a-Java-11-R-4.3.2` and `R-bundle-CRAN/2023.12-foss-2023a`, together with the study R library.

| Software | Version |
| --- | --- |
| data.table | 1.17.0 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| readr | 2.1.5 |
| readxl | 1.4.5 |
| ggplot2 | 3.5.2 |
| lubridate | 1.9.4 |
| stringr | 1.5.1 |
| cowplot | 1.1.3 |
| ggrepel | 0.9.4 |
| circlize | 0.4.15 |

Base R packages and plotting devices are supplied with R. See [`../figures/README.md`](../figures/README.md) for the local environment used with aggregate figure inputs.
