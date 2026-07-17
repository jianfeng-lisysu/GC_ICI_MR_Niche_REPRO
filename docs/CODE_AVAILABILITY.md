# Code availability

The analysis code is available at: https://github.com/jianfeng-lisysu/GC_ICI_MR_Niche_REPRO

The repository contains reproducible R scripts for:

- FinnGen R13 outcome formatting
- Proteome-wide Mendelian randomization
- Candidate cis-pQTL Mendelian randomization
- GSE167297 single-cell expression analysis
- GSE251950 spatial transcriptomic analysis

Raw sequencing data, individual-level clinical data, large GWAS summary-statistic files,
intermediate expression matrices, and authentication credentials are not included.

Public datasets should be downloaded from their original repositories.

OpenGWAS authentication is read from the user-level `OPENGWAS_JWT` environment variable.
The JWT token itself is never stored in the repository.

All currently deposited analysis scripts are deterministic and do not invoke stochastic
sampling, permutation, bootstrap resampling, random initialization, or clustering procedures.
Therefore, no random seed was required for these scripts.

The Bayesian colocalization analysis performed with `coloc.abf` is also deterministic
and does not require a random seed.
