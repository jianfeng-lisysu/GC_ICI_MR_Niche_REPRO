# Build history

This repository is cumulative. Each later working build includes all content from the earlier build plus newly recovered or reconstructed analysis code.

- Stage 1: cleaned the original repository; removed obsolete files and absolute paths; replaced the old ten-section spatial script with the nine-primary-tumour inference workflow; added software, security, and reproducibility documentation.
- Stage 2: retained all Stage 1 content and added the locked anti-PD-1 signature analysis and the final scripts for Supplementary Figures S3, S4, S5, and S8.

Only the latest cumulative repository should be retained or pushed. The numbered stage archives are audit snapshots and are not independent dependencies.

## Cumulative update: composition, Firth, incremental value and CPTAC

Added scripts 11-14, updated execution order, corrected bootstrap/seed documentation, documented required local inputs and retained explicit validation against locked manuscript results.


## Cumulative update: PRJEB25780 DESeq2, GSEA and S7

Script 15 now validates and packages the verified historical 36,921-row DESeq2 output, the authoritative 369-gene cleaned DEG file and the 5,880-pathway GO BP GSEA output. The reconstructed count matrix is used only to generate the VST matrix required for the heatmap. The exact locked Figure S6 volcano layer was not reproduced from the archived all-gene table and is intentionally retained unchanged. Script 16 now generates only the final 600-dpi Supplementary Figure S7 top-DEG heatmap. The verified RDS and reproduction audit are deposited under `results_summary/`. This cumulative archive supersedes the earlier stage files.
