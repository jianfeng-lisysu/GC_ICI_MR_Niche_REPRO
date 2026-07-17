# Reproduction Code Order

Run scripts in this order:

1. `01_format_finngen_outcomes.R`
2. `02_opengwas_protein_mr_discovery.R`
3. `03_candidate_cis_only_mr.R`
4. `04_GSE167297_singlecell_expression.R`

Before running script 2, make sure `OPENGWAS_JWT` is stored in `~/.Renviron`.

Do not store the JWT token in the project folder.
