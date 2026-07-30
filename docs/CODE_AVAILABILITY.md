# Code availability

The repository contains analysis scripts for MR, colocalization, single-cell localization, nine-patient spatial inference, anti-PD-1 signature analyses, expression-derived composition sensitivity, Firth response modeling, ROC/DeLong/LOOCV/bootstrap incremental-value analyses, CPTAC-STAD FCN1 proteomic analyses, validation of the verified historical PRJEB25780 DESeq2/GSEA outputs and the final Supplementary Figure S7 heatmap workflow.

Raw sequencing data, individual-level clinical data, large GWAS summary-statistic files and large intermediate matrices are not redistributed. The scripts require the public source datasets and locally generated locked intermediate objects listed in `REQUIRED_LOCAL_INPUTS.md`.

AUC confidence intervals use 2,000 bootstrap replicates; optimism correction uses 1,000 replicates with base seed 20260721. The historical GO BP GSEA output records seed 20260721; the deposited workflow validates and packages that result rather than rerunning fgsea. Other stochastic operations must declare their seeds in the corresponding script.

The negative-control S8 graphics script uses locked summary statistics because the underlying patient-level negative-control vectors were not preserved in the audit bundle. It is therefore a transparent figure reconstruction, not a patient-level recomputation.

The exact historical volcano-plot layer used in Supplementary Figure S6 could not be reconstructed from the archived all-gene table without changing the locked presentation. The approved formal S6 is therefore retained unchanged and is not overwritten by the deposited figure script.
