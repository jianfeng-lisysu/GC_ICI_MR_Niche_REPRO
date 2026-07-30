# Code execution order

1. 01_format_finngen_outcomes.R
2. 02_opengwas_protein_mr_discovery.R
3. 03_candidate_cis_only_mr.R
4. 04_GSE167297_singlecell_expression.R
5. 05_GSE251950_primary9_spatial_inference.R
6. 06_FCN1_TNFSF12_colocalization.R
7. 07_anti_PD1_locked_signature_analysis.R
8. 08_anti_PD1_supplementary_figures_S3_S5.R
9. 09_anti_PD1_supplementary_figure_S4_response_groups.R
10. 10_PRJEB25780_negative_controls_Figure_S8.R
11. 11_tumor_composition_sensitivity_analysis.R
12. 12_PRJEB40416_Firth_response_sensitivity.R
13. 13_incremental_value_ROC_Firth_bootstrap.R
14. 14_CPTAC_STAD_FCN1_proteomic_analysis.R
15. 15_PRJEB25780_DESeq2_GSEA.R
16. 16_PRJEB25780_supplementary_figure_S7_heatmap.R
17. 17_audit_figure_resolution.R
18. 99_export_software_versions.R

Scripts 7-13 require locally generated locked score/intermediate objects. Script 13 is stochastic only through explicitly seeded bootstrap resampling. Script 14 uses the canonical FCN1 protein row and excludes the TNFSF12-TNFSF13 readthrough from canonical TNFSF12 validation. Script 15 validates and packages the verified historical PRJEB25780 DESeq2/GSEA outputs and reconstructs the VST matrix used for the heatmap. Script 16 generates the final 600-dpi Supplementary Figure S7 only; the locked formal Figure S6 is retained unchanged. Script 17 audits the local source PNG/TIFF/JPEG files and writes a resolution report without resaving any image.
