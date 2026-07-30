# GC_ICI_MR_Niche_REPRO

## Public repository

The canonical public repository for this project is:

`https://github.com/jianfeng-lisysu/GC_ICI_MR_Niche_REPRO`

Before manuscript submission, confirm that this URL opens in a signed-out or private browser window and that the repository visibility is **Public**.

Reproducibility repository for the manuscript on FCN1/TNFSF12-associated stromal-myeloid patterns in gastric cancer.

## Repository status

This cumulative cleaned package contains the deposited MR, single-cell and colocalization scripts; the corrected nine-patient GSE251950 spatial-inference workflow; locked PRJEB25780/PRJEB40416 signature analyses; expression-derived composition sensitivity; PRJEB40416 Firth response models; ROC/DeLong/nested-Firth/LOOCV/bootstrap incremental-value analyses; CPTAC-STAD FCN1 proteomic analyses; verified historical PRJEB25780 DESeq2/GSEA outputs; and the locked Supplementary Figure S7 heatmap workflow.

The core statistical analyses described in the manuscript are now represented. Some Word-table formatting and final composite-main-figure assembly remain presentation-layer tasks rather than missing statistical analyses.

## Running the code

Run scripts from the repository root and define the complete local analysis project directory before starting R:

```r
Sys.setenv(GC_ICI_MR_NICHE_PROJECT = "E:/Bioinformatics analysis/GC_ICI_MR_Niche_REPRO_20260617")
```

The scripts read `config/project_config.R`. Public raw data and locally generated intermediate files remain outside this repository.

## Current script order

1. `scripts/01_format_finngen_outcomes.R`
2. `scripts/02_opengwas_protein_mr_discovery.R`
3. `scripts/03_candidate_cis_only_mr.R`
4. `scripts/04_GSE167297_singlecell_expression.R`
5. `scripts/05_GSE251950_primary9_spatial_inference.R`
6. `scripts/06_FCN1_TNFSF12_colocalization.R`
7. `scripts/07_anti_PD1_locked_signature_analysis.R`
8. `scripts/08_anti_PD1_supplementary_figures_S3_S5.R`
9. `scripts/09_anti_PD1_supplementary_figure_S4_response_groups.R`
10. `scripts/10_PRJEB25780_negative_controls_Figure_S8.R`
11. `scripts/11_tumor_composition_sensitivity_analysis.R`
12. `scripts/12_PRJEB40416_Firth_response_sensitivity.R`
13. `scripts/13_incremental_value_ROC_Firth_bootstrap.R`
14. `scripts/14_CPTAC_STAD_FCN1_proteomic_analysis.R`
15. `scripts/15_PRJEB25780_DESeq2_GSEA.R`
16. `scripts/16_PRJEB25780_supplementary_figure_S7_heatmap.R`
17. `scripts/17_audit_figure_resolution.R`
18. `scripts/99_export_software_versions.R`

Script 10 reconstructs the final negative-control graphic from locked summary statistics because the patient-level negative-control vectors were not preserved in the repository audit bundle. Script 15 validates and packages the verified historical PRJEB25780 DESeq2/GSEA outputs; the reconstructed count matrix is used only for the VST heatmap. Script 16 regenerates Supplementary Figure S7 only. The already locked formal Figure S6 is intentionally not overwritten because its exact historical volcano-plot layer is not recoverable from the archived all-gene table. Script 17 audits the source raster files in `05_figures` at a conservative 600-dpi gate without modifying them.

## Resampling and seeds

Script 13 uses 2,000 bootstrap replicates for AUC confidence intervals and 1,000 replicates for optimism correction with base seed `20260721`. The deposited Script 15 loads verified historical GSEA results and records the historical seed rather than rerunning fgsea. Script 16 is deterministic. Most other analyses are deterministic.

## Authentication

OpenGWAS authentication is read from the user-level `OPENGWAS_JWT` environment variable. Tokens and `.Renviron` files must never be committed.

## Data

Raw sequencing data, individual-level clinical data, large GWAS files and intermediate matrices are not distributed. Public datasets must be obtained from their original repositories. See `docs/DATA_SOURCE_NOTES.md` and `docs/REQUIRED_LOCAL_INPUTS.md`.

## Software

The analysis environment used R 4.6.0. Recorded package versions are in `docs/SOFTWARE_VERSIONS.tsv`. Historical versions of TwoSampleMR and pROC could not be recovered and are therefore not guessed.
