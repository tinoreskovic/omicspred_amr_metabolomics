# UK Biobank analyses

These scripts externally validate the INTERVAL-trained and MCPS-trained genetic models among UK Biobank participants with admixed American genetic ancestries. They are run within the UK Biobank Research Analysis Platform after the genetic scores, Nightingale NMR traits and ancestry results have been made available to the analysis job.

- `00_prepare_analysis_inputs.R` matches the genetic scores to the Nightingale NMR traits using the metadata workbook included with the repository and prepares the inputs for evaluation.
- `01_evaluate_models.R` evaluates both sets of genetic models and records both R² and Spearman rho. This is the analysis without adjustment for BMI.
- `02_plot_pca.R` makes the UK Biobank PCA panel used in Figure 3d.

The UK Biobank file locations are set in `config/env.example`.
