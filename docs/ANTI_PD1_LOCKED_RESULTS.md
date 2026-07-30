# Locked anti-PD-1 cohort results

The deposited anti-PD-1 scripts read the final locked score object:

`04_results/PRJEB40416_locked_signature_reanalysis_results.rds`

The score object contains baseline patient-level signature scores for PRJEB25780 and PRJEB40416. Signature genes were standardized across samples after log2(TPM + 1) transformation and averaged with equal weights.

## Signature definitions

- FCN1/TNFSF12: FCN1, TNFSF12
- Stroma-ECM: COL1A1, COL1A2, COL3A1, DCN, LUM, FBLN1, FN1
- CD8 T: CD8A, CD8B, GZMK, CXCL13
- Cytotoxic: GZMB, GZMA, PRF1, NKG7, GNLY

## Locked sample composition

- PRJEB25780: 23 baseline and response-evaluable patients; 7 responders and 16 non-responders.
- PRJEB40416: 16 baseline-RNA patients; 15 response-evaluable patients, comprising 7 responders and 8 non-responders; one patient was not evaluable.

## Locked correlations

| Cohort | Comparison | n | Spearman rho | P value |
|---|---|---:|---:|---:|
| PRJEB25780 | FCN1/TNFSF12 vs Stroma-ECM | 23 | 0.5691700 | 0.004589741 |
| PRJEB25780 | FCN1/TNFSF12 vs CD8 T | 23 | 0.3893281 | 0.066320621 |
| PRJEB25780 | FCN1/TNFSF12 vs Cytotoxic | 23 | 0.3675889 | 0.084420641 |
| PRJEB40416 | FCN1/TNFSF12 vs Stroma-ECM | 16 | 0.3911765 | 0.134071286 |
| PRJEB40416 | FCN1/TNFSF12 vs CD8 T | 16 | 0.1294118 | 0.632882332 |
| PRJEB40416 | FCN1/TNFSF12 vs Cytotoxic | 16 | 0.2323529 | 0.386511350 |

## Locked response-group P values

| Cohort | Score | n | Wilcoxon P |
|---|---|---:|---:|
| PRJEB25780 | FCN1/TNFSF12 | 23 | 0.367053313 |
| PRJEB25780 | Stroma-ECM | 23 | 0.014737800 |
| PRJEB25780 | CD8 T | 23 | 0.014737800 |
| PRJEB25780 | Cytotoxic | 23 | 0.048717967 |
| PRJEB40416 | FCN1/TNFSF12 | 15 | 0.056197112 |
| PRJEB40416 | Stroma-ECM | 15 | 0.524449722 |
| PRJEB40416 | CD8 T | 15 | 0.602524352 |
| PRJEB40416 | Cytotoxic | 15 | 0.953857153 |

## Negative-control Figure S8

The final Figure S8 was redrawn from locked Supplementary Table S15 summary statistics. Patient-level negative-control vectors were not preserved in the repository audit bundle, so `10_PRJEB25780_negative_controls_Figure_S8.R` transparently reconstructs the graphic and does not claim to recompute those correlations.
