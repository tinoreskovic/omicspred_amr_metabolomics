# MCPS analyses

The MCPS code is organised around the three parts of the analysis:

1. `validate_interval/` prepares the Nightingale NMR metabolomic traits and validates the INTERVAL Study–trained genetic models in MCPS, overall and in the subgroups referred to in the manuscript.
2. `train_mcps_models/` prepares the training inputs, selects and prunes genome-wide significant variants, tunes the Bayesian ridge prior on a reproducible subset of traits, and trains the new genetic models for prediction/imputation of the 141 Nightingale NMR metabolomic traits.
3. `validate_mcps_models/` evaluates the MCPS-trained and INTERVAL-trained models in the same withheld MCPS subset.

The model-to-trait mapping is read from `../metadata/omicspred_validation.xlsx`. Working inputs and the intermediate files passed between scripts are written beneath `OMICSPRED_MCPS_ROOT`.

Run `validate_interval/01` through `03` first. `validate_interval/05_combine_cohort_results.R` combines the MCPS results with the published external-validation results used for Figure 1. Run `train_mcps_models/01` through `08` in order, distributing steps `05` and `06` over `TASK_INDEX` where appropriate. Finally, run `validate_mcps_models/00` and `01` to compare both sets of genetic models in the withheld subset.
