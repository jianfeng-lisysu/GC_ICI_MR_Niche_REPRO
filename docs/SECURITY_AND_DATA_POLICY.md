# Security and data policy

This repository contains analysis code only.

The following items must not be committed:

- OpenGWAS JWT tokens or other authentication credentials
- `.Renviron` or environment-variable files
- Individual-level clinical or patient data
- FASTQ, BAM, VCF, RDS, RData, H5AD, or large GWAS files
- Raw or intermediate expression matrices

The OpenGWAS token is read from the user-level environment variable `OPENGWAS_JWT`.
The token itself is never stored in this repository.

Publicly available datasets must be downloaded from their original repositories.
File paths in the scripts describe the expected local directory structure only.
