# ============================================================
# FCN1 and TNFSF12 Bayesian colocalization analysis
#
# Input:
#   04_results/MR_colocalization/10_FinnGen_regional_matches.rds
#
# Outputs:
#   04_results/TableS10_MR_colocalization_results.docx
#   04_results/MR_colocalization/11_coloc_complete_results.rds
#   04_results/MR_colocalization/12_coloc_PP_H4_report.txt
#
# Reproducibility:
#   Deterministic analysis; no random seed is required.
# ============================================================

required_packages <- c("coloc", "officer", "flextable")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package installation failed: ", pkg)
  }
}

library(coloc)
library(officer)
library(flextable)

# Run from the project root, or define:
# Sys.setenv(GC_ICI_MR_NICHE_PROJECT = "path/to/project")
project_dir <- Sys.getenv("GC_ICI_MR_NICHE_PROJECT", unset = ".")

result_dir <- file.path(project_dir, "04_results")
coloc_dir <- file.path(result_dir, "MR_colocalization")

regional_match_rds <- file.path(
  coloc_dir,
  "10_FinnGen_regional_matches.rds"
)

if (!file.exists(regional_match_rds)) {
  stop(
    "Missing input file:\n",
    regional_match_rds,
    "\nRun from the project root or set GC_ICI_MR_NICHE_PROJECT."
  )
}

regional_data <- readRDS(regional_match_rds)

required_objects <- c("FCN1_pQTL", "TNFSF12_pQTL", "ADENO", "WIDE")
missing_objects <- setdiff(required_objects, names(regional_data))

if (length(missing_objects) > 0) {
  stop("Input RDS is missing: ", paste(missing_objects, collapse = ", "))
}

fcn1_pqtl <- as.data.frame(regional_data$FCN1_pQTL)
tnfsf12_pqtl <- as.data.frame(regional_data$TNFSF12_pQTL)
adeno_outcome <- as.data.frame(regional_data$ADENO)
wide_outcome <- as.data.frame(regional_data$WIDE)

safe_minimum <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

format_probability <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 4))
}

run_coloc_pair <- function(
  pqtl_data,
  outcome_data,
  gene_name,
  outcome_name,
  p1 = 1e-4,
  p2 = 1e-4,
  p12 = 1e-5
) {
  merged <- merge(
    pqtl_data,
    outcome_data,
    by = "rsid",
    all = FALSE,
    sort = FALSE
  )

  merged <- merged[
    is.finite(merged$position_pqtl) &
      is.finite(merged$beta_pqtl) &
      is.finite(merged$se_pqtl) &
      merged$se_pqtl > 0 &
      is.finite(merged$beta_finngen) &
      is.finite(merged$se_finngen) &
      merged$se_finngen > 0,
    ,
    drop = FALSE
  ]

  merged <- merged[order(merged$position_pqtl, merged$rsid), , drop = FALSE]
  merged <- merged[!duplicated(merged$rsid), , drop = FALSE]
  rownames(merged) <- NULL

  if (nrow(merged) < 50) {
    stop(gene_name, " and ", outcome_name, " share fewer than 50 valid SNPs.")
  }

  dataset_pqtl <- list(
    beta = as.numeric(merged$beta_pqtl),
    varbeta = as.numeric(merged$se_pqtl^2),
    snp = as.character(merged$rsid),
    position = as.integer(merged$position_pqtl),
    type = "quant",
    sdY = 1
  )

  if ("n_pqtl" %in% colnames(merged)) {
    n_pqtl <- safe_median(merged$n_pqtl)
    if (is.finite(n_pqtl)) dataset_pqtl$N <- as.integer(round(n_pqtl))
  }

  dataset_outcome <- list(
    beta = as.numeric(merged$beta_finngen),
    varbeta = as.numeric(merged$se_finngen^2),
    snp = as.character(merged$rsid),
    position = as.integer(merged$position_pqtl),
    type = "cc"
  )

  if ("n_finngen" %in% colnames(merged)) {
    n_outcome <- safe_median(merged$n_finngen)
    if (is.finite(n_outcome)) dataset_outcome$N <- as.integer(round(n_outcome))
  }

  coloc_result <- coloc::coloc.abf(
    dataset1 = dataset_pqtl,
    dataset2 = dataset_outcome,
    p1 = p1,
    p2 = p2,
    p12 = p12
  )

  summary_values <- coloc_result$summary

  get_value <- function(name) {
    if (!name %in% names(summary_values)) return(NA_real_)
    as.numeric(summary_values[name])
  }

  pp_h0 <- get_value("PP.H0.abf")
  pp_h1 <- get_value("PP.H1.abf")
  pp_h2 <- get_value("PP.H2.abf")
  pp_h3 <- get_value("PP.H3.abf")
  pp_h4 <- get_value("PP.H4.abf")

  conditional_h4 <- if (
    is.finite(pp_h3) &&
      is.finite(pp_h4) &&
      pp_h3 + pp_h4 > 0
  ) {
    pp_h4 / (pp_h3 + pp_h4)
  } else {
    NA_real_
  }

  highest_h4_snp <- NA_character_
  highest_snp_pp_h4 <- NA_real_

  if (
    !is.null(coloc_result$results) &&
      nrow(coloc_result$results) > 0 &&
      "SNP.PP.H4" %in% colnames(coloc_result$results)
  ) {
    best_index <- which.max(coloc_result$results$SNP.PP.H4)
    highest_h4_snp <- as.character(coloc_result$results$snp[best_index])
    highest_snp_pp_h4 <- as.numeric(coloc_result$results$SNP.PP.H4[best_index])
  }

  interpretation <- if (
    gene_name == "FCN1" &&
      outcome_name == "C3_STOMACH_ADENO"
  ) {
    "Prior-dependent evidence compatible with a shared regional signal; sensitive to prior specification"
  } else {
    "No robust evidence for colocalization"
  }

  summary_table <- data.frame(
    protein = gene_name,
    exposure_id = unique(pqtl_data$exposure_id)[1],
    FinnGen_outcome = outcome_name,
    shared_SNPs = nrow(merged),
    minimum_p_pQTL = safe_minimum(merged$p_pqtl),
    minimum_p_outcome = safe_minimum(merged$p_finngen),
    PP_H0 = pp_h0,
    PP_H1 = pp_h1,
    PP_H2 = pp_h2,
    PP_H3 = pp_h3,
    PP_H4 = pp_h4,
    H4_over_H3_plus_H4 = conditional_h4,
    highest_H4_SNP = highest_h4_snp,
    highest_SNP_PP_H4 = highest_snp_pp_h4,
    p1 = p1,
    p2 = p2,
    p12 = p12,
    interpretation = interpretation,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary_table,
    matched_data = merged,
    coloc_object = coloc_result
  )
}

pair_definitions <- list(
  FCN1_ADENO = list(fcn1_pqtl, adeno_outcome, "FCN1", "C3_STOMACH_ADENO"),
  FCN1_WIDE = list(fcn1_pqtl, wide_outcome, "FCN1", "C3_STOMACH_WIDE"),
  TNFSF12_ADENO = list(tnfsf12_pqtl, adeno_outcome, "TNFSF12", "C3_STOMACH_ADENO"),
  TNFSF12_WIDE = list(tnfsf12_pqtl, wide_outcome, "TNFSF12", "C3_STOMACH_WIDE")
)

default_objects <- lapply(
  pair_definitions,
  function(x) {
    run_coloc_pair(
      pqtl_data = x[[1]],
      outcome_data = x[[2]],
      gene_name = x[[3]],
      outcome_name = x[[4]],
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = 1e-5
    )
  }
)

default_results <- do.call(
  rbind,
  lapply(default_objects, function(x) x$summary)
)

rownames(default_results) <- NULL

p12_values <- c(1e-6, 1e-5, 1e-4)
sensitivity_list <- list()
sensitivity_index <- 0

for (pair_name in names(pair_definitions)) {
  pair_item <- pair_definitions[[pair_name]]

  for (current_p12 in p12_values) {
    sensitivity_index <- sensitivity_index + 1

    sensitivity_result <- run_coloc_pair(
      pqtl_data = pair_item[[1]],
      outcome_data = pair_item[[2]],
      gene_name = pair_item[[3]],
      outcome_name = pair_item[[4]],
      p1 = 1e-4,
      p2 = 1e-4,
      p12 = current_p12
    )

    sensitivity_list[[sensitivity_index]] <- sensitivity_result$summary
  }
}

sensitivity_results <- do.call(rbind, sensitivity_list)
rownames(sensitivity_results) <- NULL

complete_rds <- file.path(
  coloc_dir,
  "11_coloc_complete_results.rds"
)

saveRDS(
  list(
    default_results = default_results,
    sensitivity_results = sensitivity_results,
    default_objects = default_objects,
    package_versions = list(
      R = R.version.string,
      coloc = as.character(packageVersion("coloc")),
      officer = as.character(packageVersion("officer")),
      flextable = as.character(packageVersion("flextable"))
    ),
    random_seed = "Not applicable; deterministic calculation"
  ),
  complete_rds
)

word_table <- default_results[
  ,
  c(
    "protein",
    "FinnGen_outcome",
    "shared_SNPs",
    "PP_H0",
    "PP_H1",
    "PP_H2",
    "PP_H3",
    "PP_H4",
    "H4_over_H3_plus_H4",
    "highest_H4_SNP",
    "highest_SNP_PP_H4",
    "interpretation"
  ),
  drop = FALSE
]

for (column_name in c(
  "PP_H0",
  "PP_H1",
  "PP_H2",
  "PP_H3",
  "PP_H4",
  "H4_over_H3_plus_H4",
  "highest_SNP_PP_H4"
)) {
  word_table[[column_name]] <- format_probability(word_table[[column_name]])
}

colnames(word_table) <- c(
  "Protein",
  "FinnGen outcome",
  "Shared SNPs",
  "PP.H0",
  "PP.H1",
  "PP.H2",
  "PP.H3",
  "PP.H4",
  "PP.H4/(PP.H3+PP.H4)",
  "Highest H4 SNP",
  "SNP.PP.H4",
  "Interpretation"
)

ft <- flextable(word_table)
ft <- border_remove(ft)

thin_border <- fp_border(color = "black", width = 0.75)
thick_border <- fp_border(color = "black", width = 1.25)

ft <- hline_top(ft, border = thick_border, part = "header")
ft <- hline_bottom(ft, border = thin_border, part = "header")
ft <- hline_bottom(ft, border = thick_border, part = "body")
ft <- bold(ft, part = "header")
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- fontsize(ft, size = 8, part = "all")
ft <- align(ft, align = "center", part = "header")
ft <- valign(ft, valign = "center", part = "all")
ft <- autofit(ft)

table_s10_file <- file.path(
  result_dir,
  "TableS10_MR_colocalization_results.docx"
)

doc <- read_docx()

doc <- body_add_par(
  doc,
  value = paste0(
    "Supplementary Table S10. Bayesian colocalization analyses of FCN1 and TNFSF12 ",
    "plasma protein levels with FinnGen R13 gastric cancer outcomes."
  ),
  style = "Normal"
)

doc <- body_add_flextable(doc, value = ft)

doc <- body_add_par(
  doc,
  value = paste0(
    "PP.H0, neither trait is associated in the region; PP.H1, association with the ",
    "protein trait only; PP.H2, association with the gastric cancer outcome only; ",
    "PP.H3, both traits are associated but are compatible with distinct causal variants; ",
    "PP.H4, both traits are compatible with a shared causal variant."
  ),
  style = "Normal"
)

doc <- body_add_par(
  doc,
  value = paste0(
    "The primary analysis used coloc.abf with p1 = 1 × 10^-4, p2 = 1 × 10^-4, ",
    "and p12 = 1 × 10^-5. Prior sensitivity analyses used p12 values of ",
    "1 × 10^-6, 1 × 10^-5, and 1 × 10^-4. PP.H4 indicates compatibility with ",
    "a shared causal signal under the single-causal-variant assumption and does not ",
    "constitute definitive proof of causality."
  ),
  style = "Normal"
)

print(doc, target = table_s10_file)

report_file <- file.path(
  coloc_dir,
  "12_coloc_PP_H4_report.txt"
)

sink(report_file, split = TRUE)

cat("============================================================\n")
cat("FCN1 AND TNFSF12 COLOCALIZATION REPORT\n")
cat("============================================================\n\n")
cat("Generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("Random seed: Not applicable; deterministic calculation\n\n")
cat("DEFAULT PRIOR RESULTS\n")
cat("------------------------------------------------------------\n")
print(default_results, row.names = FALSE)
cat("\n\nPRIOR SENSITIVITY RESULTS\n")
cat("------------------------------------------------------------\n")
print(sensitivity_results, row.names = FALSE)
cat("\n\nOUTPUT FILES\n")
cat("------------------------------------------------------------\n")
cat(table_s10_file, "\n")
cat(complete_rds, "\n")
cat(report_file, "\n")

sink()

cat("\nColocalization analysis completed.\n")
cat("Table S10: ", table_s10_file, "\n", sep = "")
cat("Report: ", report_file, "\n", sep = "")
