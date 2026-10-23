# Manuscript figures

This folder contains the R code for the main and supplementary figures generated from the aggregate results of the cohort-specific analyses.

- `figure_1_interval_external_validation.R`: Figure 1.
- `figure_2_mcps_ancestry_comparison.R`: Figure 2.
- `figure_3_model_comparison.R`: Figure 3a and 3b. The MCPS and UKB PCA panels in Figure 3c and 3d are made by the cohort-specific PCA scripts.
- `figure_4_six_trait_validation.R`: Figure 4.
- `figure_5_disease_associations.R`: Figure 5.
- `supplementary_figure_1_heritability.R` and `supplementary_figure_2_external_validation.R`: Supplementary Figures 1 and 2.
- `supplementary_figures_3_to_8_subgroups.R`: Supplementary Figures 3–8.
- `supplementary_figure_9_association_evidence.R`: Supplementary Figure 9.
- `supplementary_figures_10_to_12_wheelplots.R`: Supplementary Figures 10–12.

The scripts read the named summary-level outputs from `OMICSPRED_RESULTS_DIR` and write the figure panels to `OMICSPRED_FIGURE_DIR`. Extended Data Figure 1 was assembled in PowerPoint and is not generated here.

## Local plotting environment

These scripts use aggregate cohort results and can run locally; participant-level preparation, validation and PCA use the corresponding cohort platforms.

R 4.5.3 on macOS Tahoe 26.6.2, arm64 (Darwin 25.6.0). 

| Software | Version |
| --- | --- |
| data.table | 1.18.2.1 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| readr | 2.2.0 |
| readxl | 1.4.5 |
| ggplot2 | 4.0.3 |
| stringr | 1.6.0 |
| cowplot | 1.2.0 |
| ggrepel | 0.9.8 |
| circlize | 0.4.18 |

