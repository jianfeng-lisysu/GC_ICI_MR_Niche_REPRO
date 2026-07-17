# 02_opengwas_protein_mr_discovery.R
# OpenGWAS protein-QTL instrument extraction and MR discovery.
# Requires OPENGWAS_JWT in ~/.Renviron.

library(data.table)
library(ieugwasr)

project_dir <- path.expand("~/Desktop/Scientific Research/Bioinformatics analysis/GC_ICI_MR_Niche_REPRO_20260617")
meta_dir <- file.path(project_dir, "01_raw_data/02_OpenGWAS_metadata")
inst_dir <- file.path(project_dir, "03_intermediate/01_OpenGWAS_instruments")
harm_dir <- file.path(project_dir, "03_intermediate/02_harmonised_data")
mr_dir <- file.path(project_dir, "04_results/01_MR_results")
fmt_dir <- file.path(project_dir, "01_raw_data/01_FinnGen_R13_outcome/formatted")

dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(inst_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(harm_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mr_dir, recursive = TRUE, showWarnings = FALSE)

if (Sys.getenv("OPENGWAS_JWT") == "") {
  stop("OPENGWAS_JWT is empty. Add OPENGWAS_JWT=<token> to ~/.Renviron and restart R.")
}

if (!file.exists(file.path(meta_dir, "opengwas_prot_a_exposures.tsv"))) {
  gw <- as.data.table(ieugwasr::gwasinfo())
  fwrite(gw, file.path(meta_dir, "opengwas_all_datasets.tsv"), sep = "\t")
  prot_a <- gw[grepl("^prot-a-", id)]
  prot_a <- prot_a[, .(id, trait, author, sample_size, population, category, subcategory, year, pmid)]
  fwrite(prot_a, file.path(meta_dir, "opengwas_prot_a_exposures.tsv"), sep = "\t")
} else {
  prot_a <- fread(file.path(meta_dir, "opengwas_prot_a_exposures.tsv"))
}

inst_file <- file.path(inst_dir, "prot_a_all_instruments_checkpoint.tsv")
fail_file <- file.path(inst_dir, "prot_a_failed_ids.tsv")

if (file.exists(inst_file)) {
  all_inst <- fread(inst_file)
  done_ids <- unique(all_inst$exposure_id)
} else {
  all_inst <- data.table()
  done_ids <- character()
}

if (file.exists(fail_file)) {
  failed <- fread(fail_file)
  failed_ids <- unique(failed$exposure_id)
} else {
  failed <- data.table(exposure_id = character(), error = character())
  failed_ids <- character()
}

ids_to_run <- setdiff(prot_a$id, union(done_ids, failed_ids))
message("Total proteins: ", nrow(prot_a))
message("Done: ", length(done_ids))
message("Failed: ", length(failed_ids))
message("Remaining: ", length(ids_to_run))

get_one_instrument <- function(id) {
  message("Getting instruments for: ", id)
  x <- tryCatch(
    ieugwasr::tophits(id = id, pval = 5e-8, clump = 1),
    error = function(e) list(error = e$message)
  )
  if (is.list(x) && !is.null(x$error)) {
    fwrite(data.table(exposure_id = id, error = x$error), fail_file, sep = "\t", append = file.exists(fail_file))
    return(NULL)
  }
  if (is.null(x) || nrow(x) == 0) {
    fwrite(data.table(exposure_id = id, error = "No genome-wide significant instrument"), fail_file, sep = "\t", append = file.exists(fail_file))
    return(NULL)
  }
  x <- as.data.table(x)
  x[, exposure_id := id]
  fwrite(x, inst_file, sep = "\t", append = file.exists(inst_file))
  Sys.sleep(0.2)
  invisible(x)
}

for (i in seq_along(ids_to_run)) {
  message("Progress: ", i, "/", length(ids_to_run))
  get_one_instrument(ids_to_run[i])
}

run_mr_for_outcome <- function(outcome_file, label) {
  inst <- fread(inst_file)
  outcome <- fread(outcome_file)
  exp_dat <- inst[, .(
    exposure_id = exposure_id,
    exposure = trait,
    SNP = rsid,
    chr.exposure = chr,
    pos.exposure = position,
    effect_allele.exposure = ea,
    other_allele.exposure = nea,
    eaf.exposure = eaf,
    beta.exposure = beta,
    se.exposure = se,
    pval.exposure = p,
    samplesize.exposure = n
  )]
  out_dat <- outcome[SNP %in% exp_dat$SNP, .(
    SNP,
    effect_allele.outcome = effect_allele,
    other_allele.outcome = other_allele,
    eaf.outcome = eaf,
    beta.outcome = beta,
    se.outcome = se,
    pval.outcome = pval
  )]
  dat <- merge(exp_dat, out_dat, by = "SNP")
  dat <- unique(dat, by = c("exposure_id", "SNP"))
  aligned <- dat$effect_allele.exposure == dat$effect_allele.outcome &
    dat$other_allele.exposure == dat$other_allele.outcome
  flipped <- dat$effect_allele.exposure == dat$other_allele.outcome &
    dat$other_allele.exposure == dat$effect_allele.outcome
  dat[, mr_keep := aligned | flipped]
  dat[, beta.outcome.aligned := beta.outcome]
  dat[flipped == TRUE, beta.outcome.aligned := -beta.outcome]
  dat[, palindromic := paste0(effect_allele.exposure, other_allele.exposure) %in% c("AT", "TA", "CG", "GC")]
  dat <- dat[mr_keep == TRUE]
  mr_one <- function(d) {
    d <- d[!is.na(beta.exposure) & !is.na(se.exposure) & !is.na(beta.outcome.aligned) & !is.na(se.outcome)]
    if (nrow(d) == 0) return(NULL)
    d[, ratio := beta.outcome.aligned / beta.exposure]
    d[, ratio_se := sqrt((se.outcome^2 / beta.exposure^2) + ((beta.outcome.aligned^2 * se.exposure^2) / beta.exposure^4))]
    if (nrow(d) == 1) {
      beta_mr <- d$ratio[1]
      se_mr <- d$ratio_se[1]
      method <- "Wald ratio"
    } else {
      w <- 1 / d$ratio_se^2
      beta_mr <- sum(w * d$ratio) / sum(w)
      se_mr <- sqrt(1 / sum(w))
      method <- "IVW fixed"
    }
    data.table(
      exposure_id = d$exposure_id[1],
      exposure = d$exposure[1],
      outcome = label,
      nsnp = nrow(d),
      method = method,
      beta = beta_mr,
      se = se_mr,
      pval = 2 * pnorm(abs(beta_mr / se_mr), lower.tail = FALSE),
      OR = exp(beta_mr),
      OR_lci95 = exp(beta_mr - 1.96 * se_mr),
      OR_uci95 = exp(beta_mr + 1.96 * se_mr),
      min_pval_exposure = min(d$pval.exposure, na.rm = TRUE),
      min_pval_outcome = min(d$pval.outcome, na.rm = TRUE),
      n_palindromic = sum(d$palindromic, na.rm = TRUE)
    )
  }
  mr_res <- rbindlist(lapply(split(dat, dat$exposure_id), mr_one), fill = TRUE)
  mr_res[, FDR := p.adjust(pval, method = "BH")]
  mr_res <- mr_res[order(pval)]
  fwrite(dat, file.path(harm_dir, paste0("prot_a_", label, "_harmonised.tsv")), sep = "\t")
  fwrite(mr_res, file.path(mr_dir, paste0("prot_a_", label, "_mr_results.tsv")), sep = "\t")
  mr_res
}

mr_adeno <- run_mr_for_outcome(file.path(fmt_dir, "finngen_R13_C3_STOMACH_ADENO_formatted.tsv.gz"), "ADENO")
mr_wide <- run_mr_for_outcome(file.path(fmt_dir, "finngen_R13_C3_STOMACH_WIDE_formatted.tsv.gz"), "WIDE")

compare <- merge(
  mr_adeno[, .(exposure_id, exposure, beta_ADENO = beta, p_ADENO = pval, OR_ADENO = OR, FDR_ADENO = FDR)],
  mr_wide[, .(exposure_id, beta_WIDE = beta, p_WIDE = pval, OR_WIDE = OR, FDR_WIDE = FDR)],
  by = "exposure_id"
)
compare[, same_direction := sign(beta_ADENO) == sign(beta_WIDE)]
compare <- compare[order(p_ADENO)]
fwrite(compare, file.path(mr_dir, "prot_a_ADENO_WIDE_MR_compare.tsv"), sep = "\t")
print(head(compare, 30))
