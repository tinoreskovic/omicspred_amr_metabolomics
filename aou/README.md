# All of Us analyses

These scripts are run in the All of Us Researcher Workbench and use the files beneath `AOU_OMICS_PRED_ROOT`.

- `01_evaluate_metabolite_models.py` externally validates the INTERVAL-trained and MCPS-trained genetic models for the six Nightingale NMR metabolomic traits used in the manuscript and records both R² and Spearman rho.
- `02_disease_associations.py` tests associations of the predicted metabolomic traits with ischemic heart disease, type 2 diabetes and chronic kidney disease, as used in Figure 5 and Supplementary Figures 9–12.

The model and Nightingale trait names are read from the metadata workbook included with the repository.

`02_disease_associations.py` uses the `GOOGLE_PROJECT` and `WORKSPACE_CDR` variables and PheTK installation available in the Workbench. The AMR ID file, filtered genetic-score files and measured-trait files remain beneath `AOU_OMICS_PRED_ROOT`; the platform QC file is set with `AOU_FLAGGED_SAMPLES_FILE`.
