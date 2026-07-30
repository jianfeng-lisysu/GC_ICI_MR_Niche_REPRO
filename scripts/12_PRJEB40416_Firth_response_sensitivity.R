# ============================================================
# 12_PRJEB40416_Firth_response_sensitivity.R
#
# Reproduces the exploratory Firth logistic-regression models
# for response in the PRJEB40416 cohort.
#
# Response is coded as the event (Responder = 1).
# The FCN1/TNFSF12 score is standardized within the 15
# response-evaluable patients. Adjusted models contain one
# clinical covariate at a time because only 7 response events
# are available.
#
# Required input:
#   04_results/PRJEB40416_locked_signature_reanalysis_results.rds
#
# Clinical covariates may be present in PRJEB40416_scores or in:
#   04_results/GitHub_reproduction/PRJEB40416_clinical_covariates.tsv
# Required covariates:
#   poor differentiation status and tumor location.
#
# Outputs:
#   04_results/GitHub_reproduction/
#     PRJEB40416_Firth_response_sensitivity.tsv
#     PRJEB40416_Firth_response_sensitivity.rds
#   06_logs_and_audit/
#     PRJEB40416_Firth_response_sensitivity_audit.txt
# ============================================================

rm(list = ls())

required_packages <- c("data.table", "logistf")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg)
  }
}
library(data.table)
library(logistf)

config_file <- file.path("config", "project_config.R")
if (!file.exists(config_file)) stop("Missing config/project_config.R. Run from the repository root.")
source(config_file)

input_rds <- file.path(project_dir, "04_results", "PRJEB40416_locked_signature_reanalysis_results.rds")
clinical_tsv <- file.path(project_dir, "04_results", "GitHub_reproduction", "PRJEB40416_clinical_covariates.tsv")
output_dir <- file.path(project_dir, "04_results", "GitHub_reproduction")
audit_dir <- file.path(project_dir, "06_logs_and_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) stop("Missing locked score object: ", input_rds)
locked <- readRDS(input_rds)
if (!"PRJEB40416_scores" %in% names(locked)) stop("Locked RDS has no PRJEB40416_scores object.")
data <- as.data.frame(locked$PRJEB40416_scores, check.names = FALSE)

find_column <- function(data, aliases, required = TRUE) {
  exact <- aliases[aliases %in% names(data)]
  if (length(exact) > 0L) return(exact[1L])
  normalized_names <- tolower(gsub("[^a-z0-9]", "", names(data)))
  normalized_aliases <- tolower(gsub("[^a-z0-9]", "", aliases))
  hit <- which(normalized_names %in% normalized_aliases)
  if (length(hit) > 0L) return(names(data)[hit[1L]])
  if (required) stop("Could not find a required column. Accepted aliases: ", paste(aliases, collapse = ", "))
  NA_character_
}

sample_col <- find_column(data, c("sample_id", "sample", "Sample_ID"))
response_col <- find_column(data, c("response_group", "response", "Response"))
score_col <- find_column(data, c("FCN1_TNFSF12_score", "FCN1.TNFSF12.score"))

if (file.exists(clinical_tsv)) {
  clinical <- fread(clinical_tsv, check.names = FALSE)
  clinical_sample_col <- find_column(clinical, c("sample_id", "sample", "Sample_ID"))
  data <- merge(
    data,
    clinical,
    by.x = sample_col,
    by.y = clinical_sample_col,
    all.x = TRUE,
    suffixes = c("", ".clinical"),
    sort = FALSE
  )
  sample_col <- find_column(data, c(sample_col, "sample_id", "sample", "Sample_ID"))
}

poor_col <- find_column(
  data,
  c("poor_differentiation", "poor_differentiation_status", "poorly_differentiated", "differentiation_status", "Differentiation"),
  required = FALSE
)
location_col <- find_column(
  data,
  c("tumor_location", "Tumor_location", "location", "Location", "primary_location"),
  required = FALSE
)

if (is.na(poor_col) || is.na(location_col)) {
  template <- data.table(
    sample_id = as.character(data[[sample_col]]),
    poor_differentiation = NA_integer_,
    tumor_location = NA_character_
  )
  template_path <- file.path(output_dir, "PRJEB40416_clinical_covariates_TEMPLATE.tsv")
  fwrite(template, template_path, sep = "\t")
  stop(
    "Clinical covariates were not found. Complete the generated template and save it as: ",
    clinical_tsv
  )
}

normalize_response <- function(x) {
  upper <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(upper))
  out[upper %in% c("R", "RESPONDER", "RESPONDERS", "CR", "PR", "1")] <- 1L
  out[upper %in% c("NR", "NON-RESPONDER", "NON_RESPONDER", "NONRESPONDER", "SD", "PD", "0")] <- 0L
  out
}

normalize_binary <- function(x, positive_patterns, negative_patterns, label) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) {
    values <- unique(x[is.finite(x)])
    if (all(values %in% c(0, 1))) return(as.integer(x))
  }
  text <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(text))
  for (pattern in positive_patterns) out[grepl(pattern, text)] <- 1L
  for (pattern in negative_patterns) out[grepl(pattern, text)] <- 0L
  if (any(is.na(out) & !is.na(x))) {
    stop("Unrecognized ", label, " values: ", paste(unique(as.character(x)[is.na(out) & !is.na(x)]), collapse = ", "))
  }
  out
}

analysis_data <- data.frame(
  sample_id = as.character(data[[sample_col]]),
  response = normalize_response(data[[response_col]]),
  score = as.numeric(data[[score_col]]),
  poor_differentiation = normalize_binary(
    data[[poor_col]],
    c("poor", "undiffer", "yes", "^1$"),
    c("well", "moderate", "other", "no", "^0$"),
    "differentiation"
  ),
  tumor_body = normalize_binary(
    data[[location_col]],
    c("body", "corpus", "yes", "^1$"),
    c("antrum", "distal", "no", "^0$"),
    "tumor location"
  ),
  stringsAsFactors = FALSE
)

analysis_data <- analysis_data[complete.cases(analysis_data[, c("response", "score")]), , drop = FALSE]
analysis_data$score_sd <- as.numeric(scale(analysis_data$score))

if (nrow(analysis_data) != 15L || sum(analysis_data$response == 1L) != 7L || sum(analysis_data$response == 0L) != 8L) {
  stop("Response-evaluable composition does not match the locked analysis (n=15; R=7; NR=8).")
}

fit_model <- function(formula, description) {
  variables <- all.vars(formula)
  d <- analysis_data[complete.cases(analysis_data[, variables, drop = FALSE]), , drop = FALSE]
  fit <- logistf(formula, data = d, pl = TRUE)
  term <- "score_sd"
  coefficient <- unname(fit$coefficients[term])
  lower <- unname(fit$ci.lower[term])
  upper <- unname(fit$ci.upper[term])
  p_value <- unname(fit$prob[term])
  data.table(
    model = description,
    n = nrow(d),
    responders = sum(d$response == 1L),
    nonresponders = sum(d$response == 0L),
    score_beta = coefficient,
    OR_per_1SD = exp(coefficient),
    CI_lower = exp(lower),
    CI_upper = exp(upper),
    p_value = p_value,
    fit = list(fit)
  )
}

models <- rbindlist(list(
  fit_model(response ~ score_sd, "Unadjusted"),
  fit_model(response ~ score_sd + poor_differentiation, "Adjusted for poor differentiation status"),
  fit_model(response ~ score_sd + tumor_body, "Adjusted for tumor location")
), fill = TRUE)

result_table <- copy(models)
result_table[, fit := NULL]

expected <- data.table(
  model = c("Unadjusted", "Adjusted for poor differentiation status", "Adjusted for tumor location"),
  expected_beta = c(-1.288, -1.188, -1.236),
  expected_OR = c(0.276, 0.305, 0.291),
  expected_lower = c(0.040, 0.048, 0.044),
  expected_upper = c(0.951, 1.038, 0.955),
  expected_p = c(0.040, 0.058, 0.041)
)
validation <- merge(result_table, expected, by = "model", sort = FALSE)
validation[, beta_pass := abs(score_beta - expected_beta) < 0.01]
validation[, OR_pass := abs(OR_per_1SD - expected_OR) < 0.01]
validation[, p_pass := abs(p_value - expected_p) < 0.01]
if (!all(validation$beta_pass) || !all(validation$OR_pass) || !all(validation$p_pass)) {
  print(validation)
  stop("Firth models did not reproduce the locked results. Check clinical coding.")
}

fwrite(result_table, file.path(output_dir, "PRJEB40416_Firth_response_sensitivity.tsv"), sep = "\t")
saveRDS(list(data = analysis_data, results = result_table, validation = validation), file.path(output_dir, "PRJEB40416_Firth_response_sensitivity.rds"))

report <- c(
  "PRJEB40416 FIRTH RESPONSE SENSITIVITY AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Outcome: response (Responder=1; Non-responder=0)",
  "Score scale: per 1 SD increase within the 15 evaluable patients",
  "",
  "RESULTS",
  capture.output(print(result_table)),
  "",
  "VALIDATION",
  capture.output(print(validation))
)
writeLines(report, file.path(audit_dir, "PRJEB40416_Firth_response_sensitivity_audit.txt"), useBytes = TRUE)

cat("PRJEB40416 Firth sensitivity models reproduced successfully.\n")
print(result_table)
