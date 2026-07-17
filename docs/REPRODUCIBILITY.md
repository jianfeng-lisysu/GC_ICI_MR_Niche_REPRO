# Reproducibility and random-seed declaration

## Repository

https://github.com/jianfeng-lisysu/GC_ICI_MR_Niche_REPRO

## Random-seed declaration

The deposited R scripts were reviewed for stochastic procedures.
No calls to random-number generation, random sampling, bootstrap resampling,
permutation procedures, stochastic dimensionality reduction, or random clustering
were identified.

Accordingly, the deposited analyses are deterministic and no random seed was required.

The Bayesian colocalization analysis using `coloc.abf` is a deterministic calculation
and likewise does not require a random seed.

A fixed random seed must be declared before execution if stochastic analyses are added
to the repository in the future.

## Authentication

The OpenGWAS JWT is supplied through the user-level `OPENGWAS_JWT` environment variable.
Authentication tokens and `.Renviron` files must never be committed.

## Data availability

The repository distributes analysis code only.
Raw data and large summary-statistic files must be obtained from their original sources.
