# ============================================================
# 15_PRJEB25780_DESeq2_GSEA.R
#
# Packages the verified historical PRJEB25780 DESeq2 and GO BP GSEA
# outputs and reconstructs the expression matrix used for Figure S7.
#
# IMPORTANT REPRODUCIBILITY BOUNDARY:
# The currently archived reconstructed count matrix does not reproduce
# the historical DESeq2 result exactly (696 versus 369 cleaned genes at
# adjusted P < 0.1). Therefore, this script does not silently replace
# the published analysis with a new result. It validates and packages
# the original locked DESeq2/GSEA result files. The count matrix is used
# only to generate a VST expression matrix for the row-scaled heatmap.
#
# Contrast:
#   Responder versus Non-responder
#   Positive log2 fold change = higher expression in responders
#   Negative log2 fold change = higher expression in non-responders
#
# Required local inputs:
#   1. A gene-level PRJEB25780 count matrix, preferably:
#      01_raw_data/05_PRJEB25780_ICI/rnaseq_fastq_tumor/
#        PRJEB25780_gene_counts_annotated_39samples.csv
#      or:
#      04_results/PRJEB25780_gene_counts_39samples.csv
#   2. 04_results/PRJEB40416_locked_signature_reanalysis_results.rds
#      containing PRJEB25780_scores with sample_id and response_group.
#
# Main outputs:
#   04_results/GitHub_reproduction/
#     PRJEB25780_DESeq2_results_all.tsv
#     PRJEB25780_DESeq2_results_clean.tsv
#     PRJEB25780_DESeq2_DEGs_FDR0.1.tsv
#     PRJEB25780_DESeq2_top20_each_direction.tsv
#     PRJEB25780_GSEA_GO_BP_results.tsv
#     PRJEB25780_GSEA_GO_BP_top10_each_direction.tsv
#     PRJEB25780_DESeq2_GSEA_results.rds
#   06_logs_and_audit/
#     PRJEB25780_DESeq2_GSEA_reproduction_audit.txt
#
# Locked manuscript checks:
#   - 23 tumors: 7 responders and 16 non-responders
#   - historical DESeq2 all-gene file: 36,921 rows
#   - historical cleaned DEGs: 369 unique gene labels at adjusted P < 0.1
#   - historical GO BP GSEA file: 5,880 pathways
#   - FBLN1 log2FC approximately -3.889 and adjusted P 1.45e-4
#   - GSEA directions for collagen fibril organization,
#     cell killing and complement activation match the manuscript
#
# The script does not overwrite the historical source matrices.
# ============================================================

rm(list = ls())
graphics.off()

required_packages <- c(
  "data.table",
  "DESeq2",
  "clusterProfiler",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Required packages are not installed: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall Bioconductor packages before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------
# Project path
#
# This standalone version does not require config/project_config.R.
# It first checks the environment variable GC_ICI_MR_NICHE_PROJECT;
# otherwise it uses the locked local project path below.
# ------------------------------------------------------------

project_dir <- Sys.getenv(
  "GC_ICI_MR_NICHE_PROJECT",
  unset = paste0(
    "E:/Bioinformatics analysis/",
    "GC_ICI_MR_Niche_REPRO_20260617"
  )
)

project_dir <- normalizePath(
  project_dir,
  winslash = "/",
  mustWork = FALSE
)

if (!dir.exists(project_dir)) {
  stop(
    paste0(
      "Project directory does not exist:
",
      project_dir,
      "

Please edit project_dir near the top of the script."
    )
  )
}

output_dir <- file.path(project_dir, "04_results", "GitHub_reproduction")
audit_dir <- file.path(project_dir, "06_logs_and_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

find_first_existing <- function(paths, label) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0L) {
    stop(
      "Could not locate ", label, ". Searched:\n- ",
      paste(paths, collapse = "\n- ")
    )
  }
  normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE)
}

count_file <- find_first_existing(
  c(
    file.path(
      project_dir,
      "01_raw_data",
      "05_PRJEB25780_ICI",
      "rnaseq_fastq_tumor",
      "PRJEB25780_gene_counts_annotated_39samples.csv"
    ),
    file.path(
      project_dir,
      "01_raw_data",
      "05_PRJEB25780_ICI",
      "rnaseq_fastq_tumor",
      "PRJEB25780_gene_counts_39samples.csv"
    ),
    file.path(project_dir, "04_results", "PRJEB25780_gene_counts_annotated_39samples.csv"),
    file.path(project_dir, "04_results", "PRJEB25780_gene_counts_39samples.csv"),
    file.path(project_dir, "04_results", "PRJEB25780_gene_counts.csv")
  ),
  "the PRJEB25780 gene-level count matrix"
)

locked_rds <- find_first_existing(
  c(
    file.path(project_dir, "04_results", "PRJEB40416_locked_signature_reanalysis_results.rds"),
    file.path(project_dir, "04_results", "GitHub_reproduction", "anti_PD1_locked_signature_results.rds")
  ),
  "the locked anti-PD-1 score object"
)

output_all <- file.path(output_dir, "PRJEB25780_DESeq2_results_all.tsv")
output_clean <- file.path(output_dir, "PRJEB25780_DESeq2_results_clean.tsv")
output_degs <- file.path(output_dir, "PRJEB25780_DESeq2_DEGs_FDR0.1.tsv")
output_top_degs <- file.path(output_dir, "PRJEB25780_DESeq2_top20_each_direction.tsv")
output_gsea <- file.path(output_dir, "PRJEB25780_GSEA_GO_BP_results.tsv")
output_top_gsea <- file.path(output_dir, "PRJEB25780_GSEA_GO_BP_top10_each_direction.tsv")
output_rds <- file.path(output_dir, "PRJEB25780_DESeq2_GSEA_results.rds")
output_audit <- file.path(audit_dir, "PRJEB25780_DESeq2_GSEA_reproduction_audit.txt")

# Verified historical result files used by the manuscript and supplement.
historical_all_file <- file.path(
  project_dir,
  "04_results",
  "DESeq2_PRJEB25780_clean_all.csv"
)

historical_deg_file <- file.path(
  project_dir,
  "04_results",
  "DESeq2_PRJEB25780_clean_DEGs.csv"
)

historical_gsea_file <- file.path(
  project_dir,
  "04_results",
  "GSEA_GO_BP_PRJEB25780_Responder_vs_Non_responder.csv"
)

locked_result_files <- c(
  historical_all_file,
  historical_deg_file,
  historical_gsea_file
)

missing_locked_files <- locked_result_files[!file.exists(locked_result_files)]

if (length(missing_locked_files) > 0L) {
  stop(
    "Missing verified historical result files:\n- ",
    paste(missing_locked_files, collapse = "\n- ")
  )
}

normalize_response <- function(x) {
  text <- trimws(as.character(x))
  upper <- toupper(text)
  out <- rep(NA_character_, length(upper))
  out[upper %in% c("R", "RESPONDER", "RESPONDERS", "RESPONSE", "CR", "PR", "1")] <- "Responder"
  out[upper %in% c(
    "NR", "NON-RESPONDER", "NON_RESPONDER", "NONRESPONDER",
    "NON RESPONDER", "NON-RESPONDERS", "SD", "PD", "0"
  )] <- "Non-responder"
  unknown <- unique(text[is.na(out)])
  if (length(unknown) > 0L) {
    stop("Unrecognized response labels: ", paste(unknown, collapse = ", "))
  }
  factor(out, levels = c("Non-responder", "Responder"))
}

locked_object <- readRDS(locked_rds)
if (!"PRJEB25780_scores" %in% names(locked_object)) {
  stop("The locked RDS does not contain PRJEB25780_scores.")
}

sample_metadata <- as.data.table(locked_object$PRJEB25780_scores)
required_metadata_columns <- c("sample_id", "response_group")
missing_metadata_columns <- setdiff(required_metadata_columns, names(sample_metadata))
if (length(missing_metadata_columns) > 0L) {
  stop(
    "PRJEB25780_scores is missing: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

sample_metadata[, sample_id := trimws(as.character(sample_id))]
sample_metadata[, response_group := normalize_response(response_group)]
sample_metadata <- sample_metadata[!is.na(response_group)]

if (nrow(sample_metadata) != 23L ||
    sum(sample_metadata$response_group == "Responder") != 7L ||
    sum(sample_metadata$response_group == "Non-responder") != 16L) {
  stop(
    "Locked PRJEB25780 metadata must contain 23 evaluable tumors: ",
    "7 responders and 16 non-responders."
  )
}

raw_counts <- fread(count_file, check.names = FALSE)
if (nrow(raw_counts) == 0L || ncol(raw_counts) < 3L) {
  stop("The count matrix is empty or malformed: ", count_file)
}

sample_columns <- grep("^ERR[0-9]+$", names(raw_counts), value = TRUE)
if (length(sample_columns) == 0L) {
  sample_columns <- intersect(names(raw_counts), sample_metadata$sample_id)
}
if (length(sample_columns) == 0L) {
  stop("No PRJEB25780 sample columns were detected in the count matrix.")
}

sample_metadata[, count_sample_id := sample_id]
missing_samples <- setdiff(sample_metadata$count_sample_id, sample_columns)

# Some historical score tables used patient/PB identifiers while the count
# matrix used ENA ERR accessions. Recover that mapping only from explicit
# local mapping files; never infer it from row order.
if (length(missing_samples) > 0L) {
  mapping_candidates <- c(
    file.path(
      project_dir,
      "01_raw_data",
      "05_PRJEB25780_ICI",
      "rnaseq_fastq_tumor",
      "PRJEB25780_tumor23_ERR_PB_score_map.csv"
    ),
    file.path(
      project_dir,
      "01_raw_data",
      "05_PRJEB25780_ICI",
      "rnaseq_fastq_tumor",
      "PRJEB25780_tumor23_score_response_merged.csv"
    ),
    file.path(project_dir, "04_results", "Immune_Stroma_signature_scores_PRJEB25780.csv")
  )

  mapping_found <- FALSE
  for (mapping_file in mapping_candidates[file.exists(mapping_candidates)]) {
    mapping_table <- fread(mapping_file, check.names = FALSE)
    err_columns <- names(mapping_table)[
      vapply(
        mapping_table,
        function(column) any(grepl("^ERR[0-9]+$", as.character(column))),
        logical(1)
      )
    ]
    id_columns <- names(mapping_table)[
      vapply(
        mapping_table,
        function(column) all(sample_metadata$sample_id %in% as.character(column)),
        logical(1)
      )
    ]

    if (length(err_columns) > 0L && length(id_columns) > 0L) {
      mapping_subset <- unique(mapping_table[, .(
        sample_id = as.character(get(id_columns[[1L]])),
        count_sample_id = as.character(get(err_columns[[1L]]))
      )])
      metadata_without_count_id <- sample_metadata[
        ,
        setdiff(names(sample_metadata), "count_sample_id"),
        with = FALSE
      ]
      mapped_metadata <- merge(
        metadata_without_count_id,
        mapping_subset,
        by = "sample_id",
        all.x = TRUE,
        sort = FALSE
      )
      mapped_metadata <- mapped_metadata[match(sample_metadata$sample_id, sample_id)]
      if (all(mapped_metadata$count_sample_id %in% sample_columns)) {
        sample_metadata <- mapped_metadata
        mapping_found <- TRUE
        break
      }
    }
  }

  if (!mapping_found) {
    stop(
      "The count matrix is missing locked tumor sample identifiers and no ",
      "explicit ERR-to-patient mapping could be recovered. Locked IDs: ",
      paste(missing_samples, collapse = ", ")
    )
  }
}

choose_column <- function(candidates, available_names) {
  hits <- candidates[candidates %in% available_names]
  if (length(hits) == 0L) NA_character_ else hits[[1L]]
}

gene_id_column <- choose_column(
  c(
    "gene_id", "GeneID", "gene", "Gene", "ensembl_gene_id",
    "ENSEMBL", "feature_id", "X"
  ),
  names(raw_counts)
)

gene_symbol_column <- choose_column(
  c(
    "gene_symbol", "GeneSymbol", "symbol", "SYMBOL", "external_gene_name",
    "gene_name", "Gene_name"
  ),
  names(raw_counts)
)

non_sample_columns <- setdiff(names(raw_counts), sample_columns)
if (is.na(gene_id_column)) {
  gene_id_column <- non_sample_columns[[1L]]
}

if (is.na(gene_symbol_column)) {
  gene_symbol_column <- gene_id_column
}

selected_columns <- unique(c(
  gene_id_column,
  gene_symbol_column,
  sample_metadata$count_sample_id
))
analysis_table <- raw_counts[, ..selected_columns]

setnames(analysis_table, gene_id_column, "gene_id")
if (!identical(gene_symbol_column, gene_id_column)) {
  setnames(analysis_table, gene_symbol_column, "gene_symbol")
} else {
  analysis_table[, gene_symbol := gene_id]
}

# Preserve the versioned Ensembl identifier used by the historical
# DESeq2 results, for example ENSG00000272821.1.  The unversioned
# identifier is retained separately for annotation and display.
analysis_table[, gene_id := trimws(as.character(gene_id))]
analysis_table[, ensembl_id := sub("\\.[0-9]+$", "", gene_id)]
analysis_table[, gene_symbol := trimws(as.character(gene_symbol))]
analysis_table[
  is.na(gene_symbol) | gene_symbol == "" | toupper(gene_symbol) == "NA",
  gene_symbol := NA_character_
]

count_columns <- sample_metadata$count_sample_id
for (column_name in count_columns) {
  analysis_table[[column_name]] <- suppressWarnings(as.numeric(analysis_table[[column_name]]))
}

if (anyNA(analysis_table[, ..count_columns])) {
  stop("Non-numeric or missing count values were detected.")
}

# Sum only exact duplicated versioned gene identifiers.  Do not collapse
# different Ensembl transcript/gene-version labels to the unversioned ID,
# because the historical DESeq2 result retained versioned identifiers.
aggregated_counts <- analysis_table[
  ,
  c(
    list(
      ensembl_id = ensembl_id[[1L]],
      gene_symbol = {
        valid_symbols <- gene_symbol[
          !is.na(gene_symbol) &
            gene_symbol != "" &
            toupper(gene_symbol) != "NA"
        ]
        if (length(valid_symbols) == 0L) {
          NA_character_
        } else {
          valid_symbols[[1L]]
        }
      }
    ),
    lapply(.SD, sum)
  ),
  by = gene_id,
  .SDcols = count_columns
]

count_matrix <- as.matrix(aggregated_counts[, ..count_columns])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- aggregated_counts$gene_id

if (any(count_matrix < 0, na.rm = TRUE)) {
  stop("Negative values were detected in the count matrix.")
}

integer_difference <- max(abs(count_matrix - round(count_matrix)), na.rm = TRUE)
if (integer_difference > 0.01) {
  warning(
    "The gene count matrix contains non-integer estimated counts; ",
    "values are rounded to the nearest integer for DESeq2, matching the ",
    "historical gene-level count workflow. Maximum difference = ",
    signif(integer_difference, 4)
  )
}
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

# Preserve locked sample order.
count_matrix <- count_matrix[
  ,
  sample_metadata$count_sample_id,
  drop = FALSE
]

# Historical DESeq2 preprocessing used the conventional low-count gate:
# retain genes with at least 10 total counts across the 23 selected tumors.
# Genes with merely one non-zero count are counted for audit only and are
# not all entered into DESeq2.
gene_row_sums <- rowSums(count_matrix)
nonzero_gene_count <- sum(gene_row_sums > 0L)
keep_prefilter <- gene_row_sums >= 10L

count_matrix_filtered <- count_matrix[
  keep_prefilter,
  ,
  drop = FALSE
]

annotation_filtered <- aggregated_counts[
  match(
    rownames(count_matrix_filtered),
    gene_id
  )
]

prefilter_diagnostics <- data.frame(
  threshold = c(1L, 5L, 10L, 20L),
  retained_genes = c(
    sum(gene_row_sums >= 1L),
    sum(gene_row_sums >= 5L),
    sum(gene_row_sums >= 10L),
    sum(gene_row_sums >= 20L)
  )
)

if (nrow(count_matrix_filtered) != 36921L) {
  print(prefilter_diagnostics, row.names = FALSE)

  stop(
    paste0(
      "The historical row-sum >= 10 prefilter was not reproduced: ",
      "observed ",
      nrow(count_matrix_filtered),
      " genes, expected 36,921.\n",
      "Genes with any non-zero count: ",
      nonzero_gene_count,
      ".\n",
      "The threshold diagnostics were printed above."
    )
  )
}

col_data <- data.frame(
  sample_id = sample_metadata$count_sample_id,
  response_group = sample_metadata$response_group,
  row.names = sample_metadata$count_sample_id,
  check.names = FALSE
)

# ------------------------------------------------------------
# Reconstruct VST expression matrix for Figure S7 only
# ------------------------------------------------------------

# The reference level is Non-responder. This reconstructed DESeqDataSet
# is used only for the heatmap transformation; its differential-expression
# statistics are not substituted for the verified historical results.
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = col_data,
  design = ~ response_group
)

dds <- estimateSizeFactors(dds)
vst_object <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_object)

# ------------------------------------------------------------
# Load and validate the verified historical DESeq2 outputs
# ------------------------------------------------------------

historical_all <- fread(
  historical_all_file,
  check.names = FALSE,
  na.strings = c("", "NA", "<NA>")
)

historical_degs <- fread(
  historical_deg_file,
  check.names = FALSE,
  na.strings = c("", "NA", "<NA>")
)

required_deseq_columns <- c(
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "gene",
  "ENSEMBL",
  "SYMBOL"
)

for (file_label in c("historical all-gene", "historical DEG")) {
  current_table <- if (file_label == "historical all-gene") {
    historical_all
  } else {
    historical_degs
  }

  missing_columns <- setdiff(
    required_deseq_columns,
    names(current_table)
  )

  if (length(missing_columns) > 0L) {
    stop(
      file_label,
      " file is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

if (nrow(historical_all) != 36921L) {
  stop(
    "Historical all-gene DESeq2 file has ",
    nrow(historical_all),
    " rows; expected 36,921."
  )
}

if (nrow(historical_degs) != 369L) {
  stop(
    "Historical cleaned DEG file has ",
    nrow(historical_degs),
    " rows; expected 369."
  )
}

standardize_deseq_table <- function(data) {
  output <- copy(data)

  output[, gene_id := trimws(as.character(gene))]
  output[, ensembl_id := trimws(as.character(ENSEMBL))]
  output[
    is.na(ensembl_id) | ensembl_id == "",
    ensembl_id := sub("\\.[0-9]+$", "", gene_id)
  ]

  output[, gene_symbol := trimws(as.character(SYMBOL))]
  output[
    is.na(gene_symbol) |
      gene_symbol == "" |
      toupper(gene_symbol) == "NA",
    gene_symbol := NA_character_
  ]

  output[, gene_label := fifelse(
    !is.na(gene_symbol) & gene_symbol != "",
    gene_symbol,
    ensembl_id
  )]

  numeric_columns <- c(
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  )

  for (column_name in numeric_columns) {
    output[[column_name]] <- suppressWarnings(
      as.numeric(output[[column_name]])
    )
  }

  output[, direction := fifelse(
    is.finite(padj) & padj < 0.10 & log2FoldChange > 0,
    "R-high",
    fifelse(
      is.finite(padj) & padj < 0.10 & log2FoldChange < 0,
      "NR-high",
      "Not significant"
    )
  )]

  preferred_columns <- c(
    "gene_id",
    "ensembl_id",
    "gene_symbol",
    "gene_label",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "direction"
  )

  setcolorder(
    output,
    c(
      preferred_columns,
      setdiff(names(output), preferred_columns)
    )
  )

  output[]
}

result_table <- standardize_deseq_table(historical_all)
clean_result <- copy(result_table)

# The all-gene file retains the raw significant rows, whereas the
# dedicated DEG file is the cleaned, one-row-per-gene result used in the
# manuscript and Supplementary Table S3.
deg_raw <- result_table[
  is.finite(padj) &
    padj < 0.10
]

deg_clean <- standardize_deseq_table(historical_degs)

raw_deg_count <- nrow(deg_raw)
clean_deg_count <- nrow(deg_clean)

if (clean_deg_count != 369L) {
  stop(
    "Historical cleaned DEG file contains ",
    clean_deg_count,
    " rows; expected 369."
  )
}

if (any(!is.finite(deg_clean$padj)) || any(deg_clean$padj >= 0.10)) {
  stop(
    "The historical cleaned DEG file contains rows without ",
    "adjusted P < 0.1."
  )
}

if (anyDuplicated(deg_clean$gene_label) > 0L) {
  stop(
    "The historical cleaned DEG file still contains duplicated ",
    "display gene labels."
  )
}

# The historical audit recorded 372 raw significant rows before cleaning.
# Do not stop solely on this auxiliary count because the authoritative
# cleaned file contains the 369 genes used in the manuscript.
if (raw_deg_count != 372L) {
  warning(
    "Historical all-gene file contains ",
    raw_deg_count,
    " rows with adjusted P < 0.1; the archived audit expected 372. ",
    "Proceeding with the authoritative cleaned 369-gene file."
  )
}

select_top_direction <- function(data, direction_name, number = 20L) {
  subset_data <- copy(data[direction == direction_name])
  subset_data[, padj_sort := fifelse(is.finite(padj), padj, Inf)]
  subset_data[, pvalue_sort := fifelse(is.finite(pvalue), pvalue, Inf)]
  subset_data[, absolute_stat := abs(stat)]
  setorder(subset_data, padj_sort, pvalue_sort, -absolute_stat)
  subset_data <- head(subset_data, number)
  subset_data[, c("padj_sort", "pvalue_sort", "absolute_stat") := NULL]
  subset_data
}

top_degs <- rbindlist(
  list(
    select_top_direction(deg_clean, "R-high", 20L),
    select_top_direction(deg_clean, "NR-high", 20L)
  ),
  use.names = TRUE,
  fill = TRUE
)

if (
  nrow(top_degs[direction == "R-high"]) != 20L ||
  nrow(top_degs[direction == "NR-high"]) != 20L
) {
  stop("Could not select 20 locked DEGs in each direction.")
}

fbln1_result <- deg_clean[toupper(gene_label) == "FBLN1"]

if (nrow(fbln1_result) != 1L) {
  stop("FBLN1 was not uniquely identified in the locked DEG file.")
}

fbln1_log2fc_pass <- abs(
  fbln1_result$log2FoldChange - (-3.889)
) < 0.02

fbln1_padj_pass <- abs(
  log10(fbln1_result$padj) - log10(1.45e-4)
) < 0.08

if (!fbln1_log2fc_pass || !fbln1_padj_pass) {
  stop(
    "Locked FBLN1 result failed validation. Observed log2FC=",
    signif(fbln1_result$log2FoldChange, 7),
    ", adjusted P=",
    signif(fbln1_result$padj, 7),
    "."
  )
}

# ------------------------------------------------------------
# Load and validate the verified historical GO BP GSEA output
# ------------------------------------------------------------

gsea_result <- fread(
  historical_gsea_file,
  check.names = FALSE,
  na.strings = c("", "NA", "<NA>")
)

required_gsea_columns <- c(
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment"
)

missing_gsea_columns <- setdiff(
  required_gsea_columns,
  names(gsea_result)
)

if (length(missing_gsea_columns) > 0L) {
  stop(
    "Historical GSEA file is missing columns: ",
    paste(missing_gsea_columns, collapse = ", ")
  )
}

if (nrow(gsea_result) != 5880L) {
  stop(
    "Historical GO BP GSEA file has ",
    nrow(gsea_result),
    " rows; expected 5,880."
  )
}

for (column_name in c(
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank"
)) {
  gsea_result[[column_name]] <- suppressWarnings(
    as.numeric(gsea_result[[column_name]])
  )
}

gsea_result[, direction := fifelse(
  NES > 0,
  "R-enriched",
  "NR-enriched"
)]

preferred_gsea_columns <- c(
  "direction",
  "Description",
  "ID",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "qvalue",
  "rank",
  "leading_edge",
  "core_enrichment"
)

setcolorder(
  gsea_result,
  c(
    preferred_gsea_columns,
    setdiff(names(gsea_result), preferred_gsea_columns)
  )
)

significant_gsea <- gsea_result[
  is.finite(p.adjust) &
    p.adjust < 0.05
]

select_top_gsea <- function(data, direction_name, number = 10L) {
  subset_data <- copy(data[direction == direction_name])
  subset_data[, absolute_NES := abs(NES)]
  setorder(subset_data, p.adjust, pvalue, -absolute_NES)
  subset_data <- head(subset_data, number)
  subset_data[, absolute_NES := NULL]
  subset_data
}

top_gsea <- rbindlist(
  list(
    select_top_gsea(significant_gsea, "NR-enriched", 10L),
    select_top_gsea(significant_gsea, "R-enriched", 10L)
  ),
  use.names = TRUE,
  fill = TRUE
)

locked_pathways <- data.table(
  ID = c(
    "GO:0030199",
    "GO:0001906",
    "GO:0006956"
  ),
  expected_direction = c(
    "NR-enriched",
    "R-enriched",
    "R-enriched"
  ),
  expected_NES = c(
    -2.562281,
    2.164,
    2.666
  )
)

pathway_validation <- merge(
  locked_pathways,
  gsea_result[
    ,
    .(
      ID,
      observed_direction = direction,
      observed_NES = NES,
      observed_padj = p.adjust
    )
  ],
  by = "ID",
  all.x = TRUE,
  sort = FALSE
)

pathway_validation[, present := !is.na(observed_NES)]
pathway_validation[, direction_pass := observed_direction == expected_direction]
pathway_validation[, NES_close := abs(observed_NES - expected_NES) < 0.15]

if (
  !all(pathway_validation$present) ||
  !all(pathway_validation$direction_pass) ||
  !all(pathway_validation$NES_close)
) {
  print(pathway_validation)
  stop("Locked GO BP pathway validation failed.")
}

fwrite(result_table, output_all, sep = "\t", quote = FALSE, na = "NA")
fwrite(clean_result, output_clean, sep = "\t", quote = FALSE, na = "NA")
fwrite(deg_clean, output_degs, sep = "\t", quote = FALSE, na = "NA")
fwrite(top_degs, output_top_degs, sep = "\t", quote = FALSE, na = "NA")
fwrite(gsea_result, output_gsea, sep = "\t", quote = FALSE, na = "NA")
fwrite(top_gsea, output_top_gsea, sep = "\t", quote = FALSE, na = "NA")

saveRDS(
  list(
    count_file_for_heatmap = count_file,
    locked_score_file = locked_rds,
    historical_DESeq2_all_file = historical_all_file,
    historical_DESeq2_DEG_file = historical_deg_file,
    historical_GSEA_file = historical_gsea_file,
    reproducibility_boundary = paste0(
      "Published DESeq2/GSEA statistics loaded from verified historical outputs; ",
      "reconstructed counts used only for VST heatmap expression."
    ),
    sample_metadata = sample_metadata,
    DESeq2_results_all = result_table,
    DESeq2_results_clean = clean_result,
    DESeq2_DEGs_raw_FDR0.1 = deg_raw,
    DESeq2_DEGs_FDR0.1 = deg_clean,
    DESeq2_top20_each_direction = top_degs,
    GSEA_GO_BP_results = gsea_result,
    GSEA_GO_BP_top10_each_direction = top_gsea,
    GSEA_pathway_validation = pathway_validation,
    vst_matrix = vst_matrix,
    DESeq2_size_factors_for_heatmap = sizeFactors(dds),
    response_contrast = "Responder versus Non-responder",
    historical_GSEA_seed_record = 20260721L,
    count_prefilter_row_sum = 10L,
    genes_with_any_nonzero_count = nonzero_gene_count,
    genes_retained_for_DESeq2 = nrow(count_matrix_filtered),
    prefilter_diagnostics = prefilter_diagnostics
  ),
  output_rds
)

audit_lines <- c(
  "PRJEB25780 DESeq2 and GO BP GSEA reproduction audit",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Count matrix used only for VST heatmap: ", count_file),
  paste0("Locked score object: ", locked_rds),
  paste0("Historical DESeq2 all-gene file: ", historical_all_file),
  paste0("Historical DESeq2 DEG file: ", historical_deg_file),
  paste0("Historical GSEA file: ", historical_gsea_file),
  "",
  "Cohort:",
  paste0("- tumors: ", nrow(sample_metadata)),
  paste0("- responders: ", sum(sample_metadata$response_group == "Responder")),
  paste0("- non-responders: ", sum(sample_metadata$response_group == "Non-responder")),
  "",
  "DESeq2:",
  paste0("- genes with any non-zero count: ", nonzero_gene_count),
  paste0("- historical row-sum prefilter: >= 10"),
  paste0("- genes retained for DESeq2: ", nrow(count_matrix_filtered)),
  paste0("- raw rows with adjusted P < 0.1 in all-gene file: ", raw_deg_count),
  paste0("- locked cleaned genes with adjusted P < 0.1: ", clean_deg_count),
  paste0("- FBLN1 log2FC: ", signif(fbln1_result$log2FoldChange, 8)),
  paste0("- FBLN1 adjusted P: ", signif(fbln1_result$padj, 8)),
  "",
  "GSEA:",
  paste0("- historical GO BP pathways loaded: ", nrow(gsea_result)),
  paste0("- significant pathways at adjusted P < 0.05: ", nrow(significant_gsea)),
  paste0("- historical seed record: 20260721"),
  "",
  "Locked pathway validation:",
  capture.output(print(pathway_validation)),
  "",
  "Output RDS:",
  output_rds
)
writeLines(audit_lines, output_audit, useBytes = TRUE)

cat("\nPRJEB25780 locked DESeq2/GSEA packaging completed.\n")
cat("Genes with any non-zero count:", nonzero_gene_count, "\n")
cat("Genes retained after row-sum >= 10:", nrow(count_matrix_filtered), "\n")
cat("Raw rows at adjusted P < 0.1:", raw_deg_count, "\n")
cat("Clean DEGs at adjusted P < 0.1:", clean_deg_count, "\n")
cat("Significant GO BP pathways:", nrow(significant_gsea), "\n")
cat("RDS:", output_rds, "\n")
cat("Audit:", output_audit, "\n")
