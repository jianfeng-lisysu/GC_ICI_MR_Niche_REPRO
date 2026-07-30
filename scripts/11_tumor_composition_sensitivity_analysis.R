# ============================================================
# 11_tumor_composition_sensitivity_analysis.R
#
# Reproduces the expression-derived composition sensitivity
# analyses for the association between the FCN1/TNFSF12 score
# and the Stroma-ECM score in PRJEB25780 and PRJEB40416.
#
# Required input:
#   04_results/PRJEB40416_locked_signature_reanalysis_results.rds
#
# The two score tables in that RDS must contain:
#   sample_id, FCN1_TNFSF12_score, Stroma_ECM_score,
#   ESTIMATEScore, ImmuneScore and StromalScore.
#
# Optional input for leave-one-Stroma-gene-out analysis:
#   04_results/PRJEB25780_gene_TPM.csv
#   04_results/PRJEB40416_gene_TPM.csv
#
# Outputs:
#   04_results/GitHub_reproduction/
#     tumor_composition_sensitivity_results.rds
#     tumor_composition_main_results.tsv
#     tumor_composition_leave_one_sample_out.tsv
#     tumor_composition_leave_one_gene_out.tsv
#   06_logs_and_audit/
#     tumor_composition_sensitivity_reproduction_audit.txt
# ============================================================

rm(list = ls())

required_packages <- c("data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg)
  }
}
library(data.table)

config_file <- file.path("config", "project_config.R")
if (!file.exists(config_file)) {
  stop("Missing config/project_config.R. Run from the repository root.")
}
source(config_file)

input_rds <- file.path(
  project_dir,
  "04_results",
  "PRJEB40416_locked_signature_reanalysis_results.rds"
)
output_dir <- file.path(project_dir, "04_results", "GitHub_reproduction")
audit_dir <- file.path(project_dir, "06_logs_and_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) {
  stop("Missing locked score object: ", input_rds)
}

locked <- readRDS(input_rds)
required_objects <- c("PRJEB25780_scores", "PRJEB40416_scores")
missing_objects <- setdiff(required_objects, names(locked))
if (length(missing_objects) > 0L) {
  stop("Locked RDS is missing: ", paste(missing_objects, collapse = ", "))
}

resolve_column <- function(data, aliases, label) {
  exact <- aliases[aliases %in% names(data)]
  if (length(exact) > 0L) return(exact[1L])
  normalized_names <- tolower(gsub("[^a-z0-9]", "", names(data)))
  normalized_aliases <- tolower(gsub("[^a-z0-9]", "", aliases))
  hit <- which(normalized_names %in% normalized_aliases)
  if (length(hit) > 0L) return(names(data)[hit[1L]])
  stop("Could not locate ", label, ". Accepted aliases: ", paste(aliases, collapse = ", "))
}

prepare_scores <- function(x, cohort) {
  x <- as.data.frame(x, check.names = FALSE)
  sample_col <- resolve_column(x, c("sample_id", "sample", "Sample_ID"), "sample identifier")
  combined_col <- resolve_column(x, c("FCN1_TNFSF12_score", "FCN1.TNFSF12.score"), "FCN1/TNFSF12 score")
  stroma_col <- resolve_column(x, c("Stroma_ECM_score", "Stroma.ECM.score"), "Stroma-ECM score")
  estimate_col <- resolve_column(x, c("ESTIMATEScore", "ESTIMATE_score", "EstimateScore"), "ESTIMATEScore")
  immune_col <- resolve_column(x, c("ImmuneScore", "Immune_score"), "ImmuneScore")
  stromal_col <- resolve_column(x, c("StromalScore", "Stromal_score"), "StromalScore")

  out <- data.frame(
    sample_id = as.character(x[[sample_col]]),
    FCN1_TNFSF12_score = as.numeric(x[[combined_col]]),
    Stroma_ECM_score = as.numeric(x[[stroma_col]]),
    ESTIMATEScore = as.numeric(x[[estimate_col]]),
    ImmuneScore = as.numeric(x[[immune_col]]),
    StromalScore = as.numeric(x[[stromal_col]]),
    cohort = cohort,
    stringsAsFactors = FALSE
  )
  out
}

scores_25780 <- prepare_scores(locked$PRJEB25780_scores, "PRJEB25780")
scores_40416 <- prepare_scores(locked$PRJEB40416_scores, "PRJEB40416")

if (nrow(scores_25780) != 23L || nrow(scores_40416) != 16L) {
  stop("Unexpected cohort sizes. Expected PRJEB25780 n=23 and PRJEB40416 n=16.")
}

raw_spearman <- function(data, y_column = "Stroma_ECM_score") {
  complete <- complete.cases(data[, c("FCN1_TNFSF12_score", y_column)])
  test_data <- data[complete, , drop = FALSE]
  test <- suppressWarnings(cor.test(
    test_data$FCN1_TNFSF12_score,
    test_data[[y_column]],
    method = "spearman",
    exact = FALSE,
    alternative = "two.sided"
  ))
  c(n = nrow(test_data), rho = unname(test$estimate), p = test$p.value)
}

partial_spearman <- function(data, covariate, y_column = "Stroma_ECM_score") {
  columns <- c("FCN1_TNFSF12_score", y_column, covariate)
  complete <- complete.cases(data[, columns])
  d <- data[complete, columns, drop = FALSE]
  x_rank <- rank(d$FCN1_TNFSF12_score, ties.method = "average")
  y_rank <- rank(d[[y_column]], ties.method = "average")
  z_rank <- rank(d[[covariate]], ties.method = "average")
  x_residual <- residuals(lm(x_rank ~ z_rank))
  y_residual <- residuals(lm(y_rank ~ z_rank))
  test <- cor.test(x_residual, y_residual, method = "pearson", alternative = "two.sided")
  c(n = nrow(d), rho = unname(test$estimate), p = test$p.value)
}

stratified_spearman <- function(data, lower = TRUE) {
  median_value <- median(data$ESTIMATEScore, na.rm = TRUE)
  if (lower) {
    d <- data[data$ESTIMATEScore < median_value, , drop = FALSE]
  } else {
    d <- data[data$ESTIMATEScore >= median_value, , drop = FALSE]
  }
  raw_spearman(d)
}

calculate_main_results <- function(data) {
  cohort <- unique(data$cohort)
  analyses <- list(
    list(name = "Raw Spearman", value = raw_spearman(data)),
    list(name = "Partial Spearman adjusted for ESTIMATEScore", value = partial_spearman(data, "ESTIMATEScore")),
    list(name = "Partial Spearman adjusted for ImmuneScore", value = partial_spearman(data, "ImmuneScore")),
    list(name = "Partial Spearman adjusted for StromalScore", value = partial_spearman(data, "StromalScore")),
    list(name = "Lower ESTIMATEScore stratum (higher inferred purity)", value = stratified_spearman(data, TRUE)),
    list(name = "Higher ESTIMATEScore stratum (lower inferred purity)", value = stratified_spearman(data, FALSE))
  )
  rbindlist(lapply(analyses, function(item) {
    data.table(
      cohort = cohort,
      analysis = item$name,
      n = as.integer(item$value["n"]),
      spearman_rho = as.numeric(item$value["rho"]),
      p_value = as.numeric(item$value["p"])
    )
  }))
}

main_results <- rbindlist(list(
  calculate_main_results(scores_25780),
  calculate_main_results(scores_40416)
))

calculate_loo_sample <- function(data) {
  methods <- c("Raw", "Partial_ESTIMATEScore", "Partial_ImmuneScore", "Partial_StromalScore")
  out <- vector("list", nrow(data) * length(methods))
  index <- 1L
  for (i in seq_len(nrow(data))) {
    d <- data[-i, , drop = FALSE]
    values <- list(
      Raw = raw_spearman(d),
      Partial_ESTIMATEScore = partial_spearman(d, "ESTIMATEScore"),
      Partial_ImmuneScore = partial_spearman(d, "ImmuneScore"),
      Partial_StromalScore = partial_spearman(d, "StromalScore")
    )
    for (method in methods) {
      out[[index]] <- data.table(
        cohort = unique(data$cohort),
        excluded_sample = data$sample_id[i],
        method = method,
        rho = as.numeric(values[[method]]["rho"])
      )
      index <- index + 1L
    }
  }
  rbindlist(out)
}

loo_sample <- rbindlist(list(
  calculate_loo_sample(scores_25780),
  calculate_loo_sample(scores_40416)
))

loo_sample_summary <- loo_sample[, .(
  iterations = .N,
  rho_min = min(rho),
  rho_median = median(rho),
  rho_max = max(rho),
  positive_direction_n = sum(rho > 0),
  negative_direction_n = sum(rho < 0)
), by = .(cohort, method)]

read_tpm_matrix <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- fread(path, check.names = FALSE)
  gene_candidates <- c("gene", "Gene", "gene_symbol", "GeneSymbol", "symbol", "Symbol")
  gene_col <- gene_candidates[gene_candidates %in% names(x)]
  if (length(gene_col) == 0L) gene_col <- names(x)[1L] else gene_col <- gene_col[1L]
  genes <- as.character(x[[gene_col]])
  expression <- as.matrix(x[, setdiff(names(x), gene_col), with = FALSE])
  storage.mode(expression) <- "numeric"
  rownames(expression) <- genes
  expression
}

normalize_sample_name <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

score_from_expression <- function(expression, genes, sample_ids) {
  available <- intersect(genes, rownames(expression))
  if (length(available) < 2L) {
    stop("Too few signature genes in TPM matrix: ", paste(available, collapse = ", "))
  }
  matrix_sample_names <- normalize_sample_name(colnames(expression))
  target_names <- normalize_sample_name(sample_ids)
  match_index <- match(target_names, matrix_sample_names)
  if (anyNA(match_index)) {
    stop("TPM matrix could not be matched to all locked samples: ", paste(sample_ids[is.na(match_index)], collapse = ", "))
  }
  values <- log2(expression[available, match_index, drop = FALSE] + 1)
  z <- t(scale(t(values)))
  colMeans(z, na.rm = TRUE)
}

calculate_loo_gene <- function(data, tpm_path) {
  expression <- read_tpm_matrix(tpm_path)
  if (is.null(expression)) return(NULL)
  genes <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FBLN1", "FN1")
  out <- list()
  index <- 1L
  for (excluded_gene in genes) {
    retained <- setdiff(genes, excluded_gene)
    d <- data
    d$Stroma_ECM_score_LOO <- score_from_expression(expression, retained, d$sample_id)
    values <- list(
      Raw = raw_spearman(d, "Stroma_ECM_score_LOO"),
      Partial_ESTIMATEScore = partial_spearman(d, "ESTIMATEScore", "Stroma_ECM_score_LOO"),
      Partial_ImmuneScore = partial_spearman(d, "ImmuneScore", "Stroma_ECM_score_LOO"),
      Partial_StromalScore = partial_spearman(d, "StromalScore", "Stroma_ECM_score_LOO")
    )
    for (method in names(values)) {
      out[[index]] <- data.table(
        cohort = unique(data$cohort),
        excluded_gene = excluded_gene,
        method = method,
        rho = as.numeric(values[[method]]["rho"])
      )
      index <- index + 1L
    }
  }
  rbindlist(out)
}

tpm_25780 <- file.path(project_dir, "04_results", "PRJEB25780_gene_TPM.csv")
tpm_40416 <- file.path(project_dir, "04_results", "PRJEB40416_gene_TPM.csv")
loo_gene_parts <- Filter(Negate(is.null), list(
  calculate_loo_gene(scores_25780, tpm_25780),
  calculate_loo_gene(scores_40416, tpm_40416)
))

if (length(loo_gene_parts) > 0L) {
  loo_gene <- rbindlist(loo_gene_parts)
  loo_gene_summary <- loo_gene[, .(
    gene_exclusions = .N,
    rho_min = min(rho),
    rho_median = median(rho),
    rho_max = max(rho),
    positive_direction_n = sum(rho > 0),
    negative_direction_n = sum(rho < 0)
  ), by = .(cohort, method)]
} else {
  loo_gene <- data.table()
  loo_gene_summary <- data.table()
}

expected_main <- data.table(
  cohort = rep(c("PRJEB25780", "PRJEB40416"), each = 6L),
  analysis = rep(c(
    "Raw Spearman",
    "Partial Spearman adjusted for ESTIMATEScore",
    "Partial Spearman adjusted for ImmuneScore",
    "Partial Spearman adjusted for StromalScore",
    "Lower ESTIMATEScore stratum (higher inferred purity)",
    "Higher ESTIMATEScore stratum (lower inferred purity)"
  ), 2L),
  expected_rho = c(0.569, 0.443, 0.708, 0.084, 0.527, 0.462, 0.391, 0.299, 0.553, 0.058, 0.071, 0.429),
  expected_p = c(0.005, 0.039, 0.000226, 0.710, 0.096, 0.131, 0.134, 0.279, 0.033, 0.836, 0.867, 0.289)
)
validation <- merge(main_results, expected_main, by = c("cohort", "analysis"), sort = FALSE)
validation[, rho_pass := abs(spearman_rho - expected_rho) < 0.002]
validation[, p_pass := abs(p_value - expected_p) < 0.002]
if (!all(validation$rho_pass) || !all(validation$p_pass)) {
  print(validation)
  stop("Composition sensitivity results did not reproduce the locked values.")
}

result_object <- list(
  main_results = main_results,
  leave_one_sample_out = loo_sample,
  leave_one_sample_out_summary = loo_sample_summary,
  leave_one_gene_out = loo_gene,
  leave_one_gene_out_summary = loo_gene_summary,
  validation = validation
)

saveRDS(result_object, file.path(output_dir, "tumor_composition_sensitivity_results.rds"))
fwrite(main_results, file.path(output_dir, "tumor_composition_main_results.tsv"), sep = "\t")
fwrite(loo_sample_summary, file.path(output_dir, "tumor_composition_leave_one_sample_out.tsv"), sep = "\t")
if (nrow(loo_gene_summary) > 0L) {
  fwrite(loo_gene_summary, file.path(output_dir, "tumor_composition_leave_one_gene_out.tsv"), sep = "\t")
}

report <- c(
  "TUMOR-COMPOSITION SENSITIVITY REPRODUCTION AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "MAIN RESULTS",
  capture.output(print(main_results)),
  "",
  "LEAVE-ONE-SAMPLE-OUT SUMMARY",
  capture.output(print(loo_sample_summary)),
  "",
  "LEAVE-ONE-STROMA-GENE-OUT SUMMARY",
  if (nrow(loo_gene_summary) > 0L) capture.output(print(loo_gene_summary)) else "TPM matrices were not available; this optional component was skipped.",
  "",
  "VALIDATION",
  capture.output(print(validation))
)
writeLines(report, file.path(audit_dir, "tumor_composition_sensitivity_reproduction_audit.txt"), useBytes = TRUE)

cat("Tumor-composition sensitivity analysis reproduced successfully.\n")
print(main_results)
