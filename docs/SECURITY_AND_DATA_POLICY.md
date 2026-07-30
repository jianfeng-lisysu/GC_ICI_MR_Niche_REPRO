# Security and data policy

Do not commit:

- OpenGWAS JWT tokens, passwords or authentication files
- `.Renviron`, private keys or access credentials
- patient-level confidential data
- FASTQ, BAM, VCF, RDS/RData, H5/H5AD or large expression matrices
- proprietary or access-controlled raw files

The `.gitignore` supplied with this repository blocks common credential, raw-data and intermediate-file types. Public data should be downloaded from their original repositories.
