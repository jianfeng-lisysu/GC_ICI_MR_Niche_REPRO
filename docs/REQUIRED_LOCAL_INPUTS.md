# Required local inputs

The GitHub repository intentionally excludes raw sequencing data, individual-level clinical files, large GWAS files and large intermediate matrices. Set `GC_ICI_MR_NICHE_PROJECT` to the complete local analysis project.

## Scripts 7-13

Required locked object:

- `04_results/PRJEB40416_locked_signature_reanalysis_results.rds`

Script 11 additionally uses the ESTIMATE-derived columns embedded in the locked score tables. Its leave-one-Stroma-gene-out component also uses:

- `04_results/PRJEB25780_gene_TPM.csv`
- `04_results/PRJEB40416_gene_TPM.csv`

Script 12 requires poor-differentiation status and tumor location. It reads them from the locked score table or from:

- `04_results/GitHub_reproduction/PRJEB40416_clinical_covariates.tsv`

When missing, the script writes a template rather than guessing clinical coding.

## Script 14

Preferred processed input:

- `04_results/PDC000614_FCN1_tumor_analysis_dataset.rds`

For the target-availability audit, sample accounting and paired tumor-normal analysis, retain the PDC000614 unshared-log-ratio GCT matrix and `CPTAC4_Gastric_Cancer_JHU_Proteome.sample.txt` under:

- `01_raw_data/08_PDC000614_CPTAC_STAD_proteome/`

The script finds these files recursively. TNFSF12-TNFSF13 is audit-only and is never substituted for canonical TNFSF12.

## Scripts 15-16

Script 15 validates and packages the verified historical PRJEB25780 outputs. It requires:

- `04_results/DESeq2_PRJEB25780_clean_all.csv`
- `04_results/DESeq2_PRJEB25780_clean_DEGs.csv`
- `04_results/GSEA_GO_BP_PRJEB25780_Responder_vs_Non_responder.csv`
- `04_results/PRJEB40416_locked_signature_reanalysis_results.rds`

To reconstruct the VST expression matrix used for Supplementary Figure S7, Script 15 also searches for the gene-level PRJEB25780 count matrix under the documented raw-data and results paths. The count matrix is not used to replace the locked historical differential-expression statistics.

Script 16 reads the packaged object from either:

- `04_results/GitHub_reproduction/PRJEB25780_DESeq2_GSEA_results.rds`, or
- `results_summary/PRJEB25780_DESeq2_GSEA_results.rds` inside this repository.

Script 16 overwrites only the formal Supplementary Figure S7 PNG/TIFF at 600 dpi. It intentionally does not regenerate or overwrite Figure S6.
