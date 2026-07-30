# PRJEB25780 DESeq2/GSEA reproducibility boundary

The manuscript and supplementary material use the verified historical files:

- `DESeq2_PRJEB25780_clean_all.csv`: 36,921 rows;
- `DESeq2_PRJEB25780_clean_DEGs.csv`: 369 cleaned unique genes at adjusted P < 0.1;
- `GSEA_GO_BP_PRJEB25780_Responder_vs_Non_responder.csv`: 5,880 GO Biological Process pathways.

The recovered count matrix contained 49,067 genes with at least one non-zero count and 36,921 genes after the historical row-sum >= 10 gate. However, rerunning DESeq2 on that reconstructed matrix yielded 696 rows at adjusted P < 0.1 rather than the locked historical result. The difference indicates that the exact historical DESeq2 object or an upstream processing detail was not fully preserved.

To avoid silently replacing the published analysis, Script 15 validates and packages the verified historical outputs. It uses the reconstructed count matrix only to calculate the VST expression matrix needed for Supplementary Figure S7. The deposited audit records FBLN1 and the manuscript-highlighted GSEA pathways.

The newly reconstructed Figure S7 matches the locked response groups and top-gene directions. The exact historical volcano layer in Figure S6 is retained unchanged because regenerating it from the archived all-gene table altered the visual scale and significant-point mapping.
