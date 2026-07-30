# ============================================================
# 07_anti_PD1_locked_signature_analysis.R
#
# Locked anti-PD-1 cohort analysis for PRJEB25780 and PRJEB40416.
#
# This script reads the final locked score object and reproduces:
#   1. FCN1/TNFSF12 correlations with Stroma-ECM, CD8 T and
#      Cytotoxic scores in both cohorts.
#   2. Response-group comparisons for four locked signatures.
#   3. Within-cohort Benjamini-Hochberg correction.
#   4. Validation against the manuscript-locked numerical results.
#
# Required input:
#   04_results/PRJEB40416_locked_signature_reanalysis_results.rds
#
# Required objects inside the RDS:
#   PRJEB25780_scores
#   PRJEB40416_scores
#
# Outputs:
#   04_results/GitHub_reproduction/
#     anti_PD1_locked_signature_results.rds
#     anti_PD1_signature_correlations.tsv
#     anti_PD1_response_group_comparisons.tsv
#   06_logs_and_audit/
#     anti_PD1_locked_signature_reproduction_audit.txt
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
  stop("Missing config/project_config.R. Run from the project root.")
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

output_rds <- file.path(output_dir, "anti_PD1_locked_signature_results.rds")
output_correlations <- file.path(output_dir, "anti_PD1_signature_correlations.tsv")
output_comparisons <- file.path(output_dir, "anti_PD1_response_group_comparisons.tsv")
output_audit <- file.path(audit_dir, "anti_PD1_locked_signature_reproduction_audit.txt")

if (!file.exists(input_rds)) {
  stop("Missing locked score object: ", input_rds)
}

locked_object <- readRDS(input_rds)
required_objects <- c("PRJEB25780_scores", "PRJEB40416_scores")
missing_objects <- setdiff(required_objects, names(locked_object))
if (length(missing_objects) > 0L) {
  stop("Locked RDS is missing: ", paste(missing_objects, collapse = ", "))
}

required_columns <- c(
  "sample_id",
  "response_group",
  "FCN1_TNFSF12_score",
  "Stroma_ECM_score",
  "CD8_T_score",
  "Cytotoxic_score"
)

normalize_response <- function(x) {
  text <- trimws(as.character(x))
  upper <- toupper(text)
  out <- rep(NA_character_, length(upper))
  out[upper %in% c("R", "RESPONDER", "RESPONDERS", "RESPONSE", "CR", "PR", "1")] <- "Responder"
  out[upper %in% c("NR", "NON-RESPONDER", "NON_RESPONDER", "NONRESPONDER", "NON RESPONDER", "NON-RESPONDERS", "SD", "PD", "0")] <- "Non-responder"
  out[upper %in% c("NE", "NOT EVALUABLE", "NOT-EVALUABLE", "NOT_EVALUABLE", "NA", "") | is.na(upper)] <- "Not evaluable"
  unknown <- unique(text[is.na(out)])
  if (length(unknown) > 0L) {
    stop("Unrecognized response labels: ", paste(unknown, collapse = ", "))
  }
  factor(out, levels = c("Responder", "Non-responder", "Not evaluable"))
}

prepare_scores <- function(x, cohort) {
  x <- as.data.frame(x, check.names = FALSE)
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns) > 0L) {
    stop(cohort, " score table is missing: ", paste(missing_columns, collapse = ", "))
  }
  score_columns <- setdiff(required_columns, c("sample_id", "response_group"))
  for (column_name in score_columns) {
    x[[column_name]] <- suppressWarnings(as.numeric(x[[column_name]]))
  }
  x$response_group <- normalize_response(x$response_group)
  x$cohort <- cohort
  x
}

scores_25780 <- prepare_scores(locked_object$PRJEB25780_scores, "PRJEB25780")
scores_40416 <- prepare_scores(locked_object$PRJEB40416_scores, "PRJEB40416")

validate_sample_counts <- function(scores_25780, scores_40416) {
  counts_25780 <- table(scores_25780$response_group, useNA = "ifany")
  counts_40416 <- table(scores_40416$response_group, useNA = "ifany")

  checks <- c(
    PRJEB25780_baseline_n = nrow(scores_25780) == 23L,
    PRJEB25780_responders = as.integer(counts_25780["Responder"]) == 7L,
    PRJEB25780_nonresponders = as.integer(counts_25780["Non-responder"]) == 16L,
    PRJEB40416_baseline_n = nrow(scores_40416) == 16L,
    PRJEB40416_responders = as.integer(counts_40416["Responder"]) == 7L,
    PRJEB40416_nonresponders = as.integer(counts_40416["Non-responder"]) == 8L,
    PRJEB40416_not_evaluable = as.integer(counts_40416["Not evaluable"]) == 1L
  )

  if (!all(checks)) {
    print(data.frame(item = names(checks), pass = checks, row.names = NULL))
    stop("Sample composition does not match the locked analysis.")
  }
  checks
}

sample_checks <- validate_sample_counts(scores_25780, scores_40416)

correlation_targets <- c(
  Stroma_ECM_score = "Stroma-ECM",
  CD8_T_score = "CD8 T",
  Cytotoxic_score = "Cytotoxic"
)

calculate_correlations <- function(data) {
  result_list <- lapply(names(correlation_targets), function(target_column) {
    complete <- is.finite(data$FCN1_TNFSF12_score) & is.finite(data[[target_column]])
    test_data <- data[complete, , drop = FALSE]
    test <- suppressWarnings(cor.test(
      test_data$FCN1_TNFSF12_score,
      test_data[[target_column]],
      method = "spearman",
      exact = FALSE,
      alternative = "two.sided"
    ))
    data.table(
      cohort = unique(test_data$cohort),
      comparison = correlation_targets[[target_column]],
      n = nrow(test_data),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value
    )
  })
  result <- rbindlist(result_list)
  result[, BH_adjusted_P_within_cohort := p.adjust(p_value, method = "BH")]
  result
}

correlation_results <- rbindlist(list(
  calculate_correlations(scores_25780),
  calculate_correlations(scores_40416)
))

score_targets <- c(
  FCN1_TNFSF12_score = "FCN1/TNFSF12",
  Stroma_ECM_score = "Stroma-ECM",
  CD8_T_score = "CD8 T",
  Cytotoxic_score = "Cytotoxic"
)

calculate_response_comparisons <- function(data) {
  evaluable <- data[data$response_group %in% c("Responder", "Non-responder"), , drop = FALSE]
  evaluable$response_group <- droplevels(evaluable$response_group)
  result_list <- lapply(names(score_targets), function(score_column) {
    complete <- is.finite(evaluable[[score_column]]) & !is.na(evaluable$response_group)
    test_data <- evaluable[complete, , drop = FALSE]
    test <- suppressWarnings(wilcox.test(
      test_data[[score_column]] ~ test_data$response_group,
      exact = FALSE,
      alternative = "two.sided"
    ))
    data.table(
      cohort = unique(test_data$cohort),
      score = score_targets[[score_column]],
      n = nrow(test_data),
      responders = sum(test_data$response_group == "Responder"),
      non_responders = sum(test_data$response_group == "Non-responder"),
      responder_median = median(test_data[[score_column]][test_data$response_group == "Responder"], na.rm = TRUE),
      non_responder_median = median(test_data[[score_column]][test_data$response_group == "Non-responder"], na.rm = TRUE),
      wilcoxon_W = unname(test$statistic),
      p_value = test$p.value
    )
  })
  result <- rbindlist(result_list)
  result[, BH_adjusted_P_within_cohort := p.adjust(p_value, method = "BH")]
  result
}

response_results <- rbindlist(list(
  calculate_response_comparisons(scores_25780),
  calculate_response_comparisons(scores_40416)
))

expected_correlations <- data.table(
  cohort = rep(c("PRJEB25780", "PRJEB40416"), each = 3L),
  comparison = rep(c("Stroma-ECM", "CD8 T", "Cytotoxic"), times = 2L),
  expected_rho = c(0.569169960474308, 0.3893281, 0.3675889, 0.391176470588235, 0.1294118, 0.2323529),
  expected_p = c(0.00458974107409376, 0.066320621, 0.084420641, 0.13407128619725525, 0.632882332, 0.386511350)
)

expected_response <- data.table(
  cohort = rep(c("PRJEB25780", "PRJEB40416"), each = 4L),
  score = rep(c("FCN1/TNFSF12", "Stroma-ECM", "CD8 T", "Cytotoxic"), times = 2L),
  expected_p = c(
    0.3670533131621528,
    0.0147378001763191,
    0.0147378001763192,
    0.0487179672469735,
    0.0561971119187939,
    0.5244497218060844,
    0.6025243523612969,
    0.9538571530616036
  )
)

correlation_validation <- merge(
  correlation_results,
  expected_correlations,
  by = c("cohort", "comparison"),
  all.x = TRUE,
  sort = FALSE
)
correlation_validation[, rho_pass := abs(spearman_rho - expected_rho) < 1e-6]
correlation_validation[, p_pass := abs(p_value - expected_p) < 1e-6]

response_validation <- merge(
  response_results,
  expected_response,
  by = c("cohort", "score"),
  all.x = TRUE,
  sort = FALSE
)
response_validation[, p_pass := abs(p_value - expected_p) < 1e-6]

if (!all(correlation_validation$rho_pass) ||
    !all(correlation_validation$p_pass) ||
    !all(response_validation$p_pass)) {
  print(correlation_validation)
  print(response_validation)
  stop("The locked score object did not reproduce the manuscript-locked results.")
}

reproduction_object <- list(
  signature_definitions = list(
    FCN1_TNFSF12 = c("FCN1", "TNFSF12"),
    Stroma_ECM = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FBLN1", "FN1"),
    CD8_T = c("CD8A", "CD8B", "GZMK", "CXCL13"),
    Cytotoxic = c("GZMB", "GZMA", "PRF1", "NKG7", "GNLY")
  ),
  PRJEB25780_scores = scores_25780,
  PRJEB40416_scores = scores_40416,
  correlations = correlation_results,
  response_group_comparisons = response_results,
  validation = list(
    sample_checks = sample_checks,
    correlations = correlation_validation,
    response = response_validation
  )
)

saveRDS(reproduction_object, output_rds)
fwrite(correlation_results, output_correlations, sep = "\t")
fwrite(response_results, output_comparisons, sep = "\t")

report <- c(
  "ANTI-PD-1 LOCKED SIGNATURE REPRODUCTION AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Input: ", normalizePath(input_rds, winslash = "/", mustWork = FALSE)),
  "",
  "SAMPLE CHECKS",
  capture.output(print(data.frame(item = names(sample_checks), pass = sample_checks, row.names = NULL))),
  "",
  "CORRELATIONS",
  capture.output(print(correlation_results)),
  "",
  "RESPONSE-GROUP COMPARISONS",
  capture.output(print(response_results)),
  "",
  "VALIDATION",
  capture.output(print(correlation_validation)),
  capture.output(print(response_validation)),
  "",
  paste0("RDS output: ", normalizePath(output_rds, winslash = "/", mustWork = FALSE))
)
writeLines(report, output_audit, useBytes = TRUE)

cat("Locked anti-PD-1 analysis reproduced successfully.\n")
print(correlation_results)
print(response_results)
