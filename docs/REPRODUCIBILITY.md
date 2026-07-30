# Reproducibility notes

Most deposited analyses are deterministic. Explicitly seeded operations are:

- bootstrap AUC confidence intervals: 2,000 replicates;
- bootstrap optimism correction: 1,000 replicates;
- historical GO Biological Process GSEA, whose archived seed record is `20260721`.

The deposited scripts do not rerun the historical PRJEB25780 GSEA or the locked Figure S6 volcano layer. The prespecified bootstrap base seed is `20260721`.

The repository never infers or fabricates missing package versions. TwoSampleMR and pROC versions were not recoverable from the historical environment and are reported without version numbers.

Scripts validate cohort sizes and manuscript-locked numerical results before writing formal outputs. A failed validation stops execution rather than silently replacing locked values. Script 15 checks the locked PRJEB25780 sample composition, the 36,921-row historical all-gene output, the authoritative 369-gene cleaned DEG file, FBLN1, the 5,880-pathway GSEA file and the directions/approximate NES values of three manuscript-highlighted GO BP pathways. The reconstructed count matrix is used only for the VST heatmap because rerunning DESeq2 on the recovered matrix yielded a different significant-row count.

CPTAC rules are fixed: canonical FCN1 is analyzed; canonical TNFSF12 and TNFRSF12A were absent from the released matrix; the TNFSF12-TNFSF13 readthrough is audit-only.
