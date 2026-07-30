# 03_candidate_cis_only_mr.R
# Candidate cis-only MR recalculation.

library(data.table)

config_file <- file.path("config", "project_config.R")
if (!file.exists(config_file)) {
  stop("Run this script from the repository/project root, or copy config/project_config.R into the project root.")
}
source(config_file)
harm_dir <- file.path(project_dir, "03_intermediate/02_harmonised_data")
cand_dir <- file.path(project_dir, "03_intermediate/03_candidate_tables")
dir.create(cand_dir, recursive = TRUE, showWarnings = FALSE)

adeno_dat <- fread(file.path(harm_dir, "prot_a_ADENO_harmonised.tsv"))
wide_dat <- fread(file.path(harm_dir, "prot_a_WIDE_harmonised.tsv"))

candidate_gene_map <- data.table(
  exposure_id = c("prot-a-431", "prot-a-1077", "prot-a-2642", "prot-a-3055", "prot-a-421", "prot-a-830", "prot-a-1882", "prot-a-2729", "prot-a-426"),
  gene_symbol = c("CD274", "FCN1", "SCARF2", "TNFSF12", "CD1D", "DLL1", "MET", "SIGLEC7", "CD209"),
  expected_chr = c("9", "9", "22", "17", "1", "6", "7", "19", "19")
)

run_cis_mr <- function(dat, label) {
  d <- merge(candidate_gene_map, dat, by = "exposure_id")
  d <- d[as.character(chr.exposure) == expected_chr]
  d[, ratio := beta.outcome.aligned / beta.exposure]
  d[, ratio_se := sqrt((se.outcome^2 / beta.exposure^2) + ((beta.outcome.aligned^2 * se.exposure^2) / beta.exposure^4))]
  mr_one <- function(x) {
    if (nrow(x) == 1) {
      beta_mr <- x$ratio
      se_mr <- x$ratio_se
      method <- "cis Wald ratio"
    } else {
      w <- 1 / x$ratio_se^2
      beta_mr <- sum(w * x$ratio) / sum(w)
      se_mr <- sqrt(1 / sum(w))
      method <- "cis IVW fixed"
    }
    data.table(
      outcome = label,
      exposure_id = x$exposure_id[1],
      gene_symbol = x$gene_symbol[1],
      exposure = x$exposure[1],
      nsnp_cis = nrow(x),
      method = method,
      beta = beta_mr,
      se = se_mr,
      pval = 2 * pnorm(abs(beta_mr / se_mr), lower.tail = FALSE),
      OR = exp(beta_mr),
      OR_lci95 = exp(beta_mr - 1.96 * se_mr),
      OR_uci95 = exp(beta_mr + 1.96 * se_mr),
      SNPs = paste(x$SNP, collapse = ";")
    )
  }
  rbindlist(lapply(split(d, d$exposure_id), mr_one), fill = TRUE)
}

cis_adeno <- run_cis_mr(adeno_dat, "ADENO")
cis_wide <- run_cis_mr(wide_dat, "WIDE")

cis_compare <- merge(
  cis_adeno[, .(exposure_id, gene_symbol, exposure, nsnp_cis_ADENO = nsnp_cis, beta_ADENO = beta, p_ADENO = pval, OR_ADENO = OR, SNPs_ADENO = SNPs)],
  cis_wide[, .(exposure_id, nsnp_cis_WIDE = nsnp_cis, beta_WIDE = beta, p_WIDE = pval, OR_WIDE = OR, SNPs_WIDE = SNPs)],
  by = "exposure_id"
)
cis_compare[, same_direction := sign(beta_ADENO) == sign(beta_WIDE)]
cis_compare <- cis_compare[order(p_ADENO)]
fwrite(cis_compare, file.path(cand_dir, "candidate_cis_only_MR_compare.tsv"), sep = "\t")
print(cis_compare)
