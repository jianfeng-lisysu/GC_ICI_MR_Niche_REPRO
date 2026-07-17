# 01_format_finngen_outcomes.R
# Standardize FinnGen R13 gastric cancer GWAS summary statistics for MR.

library(data.table)

project_dir <- path.expand("~/Desktop/Scientific Research/Bioinformatics analysis/GC_ICI_MR_Niche_REPRO_20260617")
raw_dir <- file.path(project_dir, "01_raw_data/01_FinnGen_R13_outcome/raw")
fmt_dir <- file.path(project_dir, "01_raw_data/01_FinnGen_R13_outcome/formatted")
dir.create(fmt_dir, recursive = TRUE, showWarnings = FALSE)

format_finngen <- function(infile, outfile) {
  fg <- fread(infile)
  fg2 <- fg[
    !is.na(rsids) & rsids != "." & !is.na(beta) & !is.na(sebeta) & !is.na(pval),
    .(
      SNP = rsids,
      chr = `#chrom`,
      pos = pos,
      effect_allele = alt,
      other_allele = ref,
      beta = beta,
      se = sebeta,
      pval = pval,
      eaf = af_alt,
      nearest_genes = nearest_genes
    )
  ]
  fg2 <- fg2[grepl("^rs", SNP)]
  fg2 <- unique(fg2, by = "SNP")
  fwrite(fg2, outfile, sep = "\t")
  fg2
}

adeno <- format_finngen(
  infile = file.path(raw_dir, "finngen_R13_C3_STOMACH_ADENO.gz"),
  outfile = file.path(fmt_dir, "finngen_R13_C3_STOMACH_ADENO_formatted.tsv.gz")
)

wide <- format_finngen(
  infile = file.path(raw_dir, "finngen_R13_C3_STOMACH_WIDE.gz"),
  outfile = file.path(fmt_dir, "finngen_R13_C3_STOMACH_WIDE_formatted.tsv.gz")
)

cat("ADENO formatted:", nrow(adeno), "SNPs\n")
cat("WIDE formatted:", nrow(wide), "SNPs\n")
