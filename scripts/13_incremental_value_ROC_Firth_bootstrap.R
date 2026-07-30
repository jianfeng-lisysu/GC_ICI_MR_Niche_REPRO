# ============================================================
# 13_incremental_value_ROC_Firth_bootstrap.R
#
# Reproduces the exploratory incremental-value analyses:
#   - direct-score ROC curves and bootstrap AUC confidence limits;
#   - paired DeLong comparisons;
#   - nested Firth models;
#   - leave-one-out cross-validated AUC;
#   - 1,000-replicate bootstrap optimism correction.
#
# Non-response is coded as the positive event.
# Marker directions are prespecified:
#   FCN1/TNFSF12: higher predicts non-response
#   Stroma-ECM: higher predicts non-response
#   CD8 T: lower original score predicts non-response
#
# Required input:
#   04_results/PRJEB40416_locked_signature_reanalysis_results.rds
#
# Outputs:
#   04_results/GitHub_reproduction/incremental_value_results.rds
#   04_results/GitHub_reproduction/incremental_value_*.tsv
#   05_figures/FigS12_exploratory_ROC_discrimination.png/.tiff
#   06_logs_and_audit/incremental_value_reproduction_audit.txt
# ============================================================

rm(list = ls())

required_packages <- c("data.table", "logistf", "pROC", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg)
  }
}
library(data.table)
library(logistf)
library(pROC)
library(ggplot2)

config_file <- file.path("config", "project_config.R")
if (!file.exists(config_file)) stop("Missing config/project_config.R. Run from the repository root.")
source(config_file)

input_rds <- file.path(project_dir, "04_results", "PRJEB40416_locked_signature_reanalysis_results.rds")
output_dir <- file.path(project_dir, "04_results", "GitHub_reproduction")
figure_dir <- file.path(project_dir, "05_figures")
audit_dir <- file.path(project_dir, "06_logs_and_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) stop("Missing locked score object: ", input_rds)
locked <- readRDS(input_rds)

normalize_response <- function(x) {
  upper <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(upper))
  out[upper %in% c("NR", "NON-RESPONDER", "NON_RESPONDER", "NONRESPONDER", "SD", "PD", "0")] <- 1L
  out[upper %in% c("R", "RESPONDER", "RESPONDERS", "CR", "PR", "1")] <- 0L
  out
}

prepare_cohort <- function(x, cohort) {
  x <- as.data.frame(x, check.names = FALSE)
  required <- c("response_group", "FCN1_TNFSF12_score", "Stroma_ECM_score", "CD8_T_score")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) stop(cohort, " score table is missing: ", paste(missing, collapse = ", "))
  d <- data.frame(
    cohort = cohort,
    outcome = normalize_response(x$response_group),
    FCN1_TNFSF12 = as.numeric(x$FCN1_TNFSF12_score),
    Stroma_ECM = as.numeric(x$Stroma_ECM_score),
    CD8_inverted = -as.numeric(x$CD8_T_score),
    stringsAsFactors = FALSE
  )
  d <- d[complete.cases(d), , drop = FALSE]
  d$FCN1_TNFSF12_sd <- as.numeric(scale(d$FCN1_TNFSF12))
  d$Stroma_ECM_sd <- as.numeric(scale(d$Stroma_ECM))
  d$CD8_inverted_sd <- as.numeric(scale(d$CD8_inverted))
  d
}

data_25780 <- prepare_cohort(locked$PRJEB25780_scores, "PRJEB25780")
data_40416 <- prepare_cohort(locked$PRJEB40416_scores, "PRJEB40416")
if (nrow(data_25780) != 23L || sum(data_25780$outcome == 0L) != 7L || sum(data_25780$outcome == 1L) != 16L) stop("PRJEB25780 composition mismatch.")
if (nrow(data_40416) != 15L || sum(data_40416$outcome == 0L) != 7L || sum(data_40416$outcome == 1L) != 8L) stop("PRJEB40416 composition mismatch.")

roc_from_marker <- function(outcome, marker) {
  roc(outcome, marker, levels = c(0, 1), direction = "<", quiet = TRUE)
}

set.seed(20260721)
calculate_direct_auc <- function(data) {
  markers <- c(FCN1_TNFSF12 = "FCN1/TNFSF12 score", Stroma_ECM = "Stroma-ECM score", CD8_inverted = "CD8 T score (inverted)")
  result <- list()
  rocs <- list()
  for (column in names(markers)) {
    roc_object <- roc_from_marker(data$outcome, data[[column]])
    ci <- ci.auc(roc_object, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
    result[[column]] <- data.table(
      cohort = unique(data$cohort),
      n = nrow(data),
      responders = sum(data$outcome == 0L),
      nonresponders = sum(data$outcome == 1L),
      marker = markers[[column]],
      AUC = as.numeric(auc(roc_object)),
      CI_lower = as.numeric(ci[1L]),
      CI_upper = as.numeric(ci[3L])
    )
    rocs[[column]] <- roc_object
  }
  list(table = rbindlist(result), rocs = rocs)
}

direct_25780 <- calculate_direct_auc(data_25780)
direct_40416 <- calculate_direct_auc(data_40416)
direct_results <- rbindlist(list(direct_25780$table, direct_40416$table))
roc_objects <- list(PRJEB25780 = direct_25780$rocs, PRJEB40416 = direct_40416$rocs)

calculate_delong <- function(data, comparator_column, comparator_label) {
  roc_fcn <- roc_from_marker(data$outcome, data$FCN1_TNFSF12)
  roc_comparator <- roc_from_marker(data$outcome, data[[comparator_column]])
  test <- roc.test(roc_fcn, roc_comparator, paired = TRUE, method = "delong", conf.int = TRUE)
  data.table(
    cohort = unique(data$cohort),
    comparison = paste0("FCN1/TNFSF12 vs ", comparator_label),
    AUC_FCN1_TNFSF12 = as.numeric(auc(roc_fcn)),
    AUC_comparator = as.numeric(auc(roc_comparator)),
    Delta_AUC = as.numeric(auc(roc_fcn) - auc(roc_comparator)),
    Delta_CI_lower = as.numeric(test$conf.int[1L]),
    Delta_CI_upper = as.numeric(test$conf.int[2L]),
    paired_DeLong_P = test$p.value
  )
}

delong_results <- rbindlist(list(
  calculate_delong(data_25780, "Stroma_ECM", "Stroma-ECM"),
  calculate_delong(data_25780, "CD8_inverted", "CD8 T"),
  calculate_delong(data_40416, "Stroma_ECM", "Stroma-ECM"),
  calculate_delong(data_40416, "CD8_inverted", "CD8 T")
))

fit_firth <- function(data, predictor_columns) {
  formula <- reformulate(predictor_columns, response = "outcome")
  logistf(formula, data = data, pl = TRUE)
}

predict_firth <- function(fit, newdata) {
  as.numeric(predict(fit, newdata = newdata, type = "response"))
}

auc_numeric <- function(outcome, prediction) {
  if (length(unique(outcome)) < 2L || any(!is.finite(prediction))) return(NA_real_)
  as.numeric(auc(roc_from_marker(outcome, prediction)))
}

calculate_nested <- function(data, base_column, base_label) {
  base_fit <- fit_firth(data, base_column)
  augmented_fit <- fit_firth(data, c(base_column, "FCN1_TNFSF12_sd"))
  term <- "FCN1_TNFSF12_sd"
  coefficient <- unname(augmented_fit$coefficients[term])
  lower <- unname(augmented_fit$ci.lower[term])
  upper <- unname(augmented_fit$ci.upper[term])
  p_value <- unname(augmented_fit$prob[term])
  data.table(
    cohort = unique(data$cohort),
    base_model = base_label,
    n = nrow(data),
    responders = sum(data$outcome == 0L),
    nonresponders = sum(data$outcome == 1L),
    penalized_LR_chisq = qchisq(1 - p_value, df = 1),
    LR_df = 1L,
    nested_LRT_P = p_value,
    added_FCN1_TNFSF12_OR_per_SD = exp(coefficient),
    added_OR_CI_lower = exp(lower),
    added_OR_CI_upper = exp(upper),
    added_coefficient_P = p_value,
    base_fit = list(base_fit),
    augmented_fit = list(augmented_fit)
  )
}

nested_full <- rbindlist(list(
  calculate_nested(data_25780, "Stroma_ECM_sd", "Stroma-ECM"),
  calculate_nested(data_25780, "CD8_inverted_sd", "CD8 T"),
  calculate_nested(data_40416, "Stroma_ECM_sd", "Stroma-ECM"),
  calculate_nested(data_40416, "CD8_inverted_sd", "CD8 T")
), fill = TRUE)
nested_results <- copy(nested_full)
nested_results[, c("base_fit", "augmented_fit") := NULL]

loocv_predictions <- function(data, predictor_columns) {
  predictions <- rep(NA_real_, nrow(data))
  for (i in seq_len(nrow(data))) {
    training <- data[-i, , drop = FALSE]
    testing <- data[i, , drop = FALSE]
    fit <- tryCatch(fit_firth(training, predictor_columns), error = function(e) NULL)
    if (!is.null(fit)) predictions[i] <- tryCatch(predict_firth(fit, testing), error = function(e) NA_real_)
  }
  predictions
}

bootstrap_performance <- function(data, base_column, replicates = 1000L) {
  base_fit_original <- fit_firth(data, base_column)
  augmented_fit_original <- fit_firth(data, c(base_column, "FCN1_TNFSF12_sd"))
  original_base_prediction <- predict_firth(base_fit_original, data)
  original_augmented_prediction <- predict_firth(augmented_fit_original, data)
  apparent_base <- auc_numeric(data$outcome, original_base_prediction)
  apparent_augmented <- auc_numeric(data$outcome, original_augmented_prediction)

  optimism_base <- rep(NA_real_, replicates)
  optimism_augmented <- rep(NA_real_, replicates)
  test_delta <- rep(NA_real_, replicates)
  success <- logical(replicates)

  for (b in seq_len(replicates)) {
    index <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
    bootstrap_data <- data[index, , drop = FALSE]
    if (length(unique(bootstrap_data$outcome)) < 2L) next
    base_fit <- tryCatch(fit_firth(bootstrap_data, base_column), error = function(e) NULL)
    augmented_fit <- tryCatch(fit_firth(bootstrap_data, c(base_column, "FCN1_TNFSF12_sd")), error = function(e) NULL)
    if (is.null(base_fit) || is.null(augmented_fit)) next

    base_train <- tryCatch(predict_firth(base_fit, bootstrap_data), error = function(e) NULL)
    augmented_train <- tryCatch(predict_firth(augmented_fit, bootstrap_data), error = function(e) NULL)
    base_test <- tryCatch(predict_firth(base_fit, data), error = function(e) NULL)
    augmented_test <- tryCatch(predict_firth(augmented_fit, data), error = function(e) NULL)
    if (is.null(base_train) || is.null(augmented_train) || is.null(base_test) || is.null(augmented_test)) next

    train_base_auc <- auc_numeric(bootstrap_data$outcome, base_train)
    train_augmented_auc <- auc_numeric(bootstrap_data$outcome, augmented_train)
    test_base_auc <- auc_numeric(data$outcome, base_test)
    test_augmented_auc <- auc_numeric(data$outcome, augmented_test)
    if (any(!is.finite(c(train_base_auc, train_augmented_auc, test_base_auc, test_augmented_auc)))) next

    optimism_base[b] <- train_base_auc - test_base_auc
    optimism_augmented[b] <- train_augmented_auc - test_augmented_auc
    test_delta[b] <- test_augmented_auc - test_base_auc
    success[b] <- TRUE
  }

  successful_delta <- test_delta[success]
  list(
    apparent_base = apparent_base,
    apparent_augmented = apparent_augmented,
    corrected_base = apparent_base - mean(optimism_base[success]),
    corrected_augmented = apparent_augmented - mean(optimism_augmented[success]),
    delta_median = median(successful_delta),
    delta_lower = unname(quantile(successful_delta, 0.025, type = 1)),
    delta_upper = unname(quantile(successful_delta, 0.975, type = 1)),
    successful = sum(success),
    requested = replicates
  )
}

calculate_performance <- function(data, base_column, base_label, seed_offset) {
  base_fit <- fit_firth(data, base_column)
  augmented_fit <- fit_firth(data, c(base_column, "FCN1_TNFSF12_sd"))
  apparent_base <- auc_numeric(data$outcome, predict_firth(base_fit, data))
  apparent_augmented <- auc_numeric(data$outcome, predict_firth(augmented_fit, data))
  loocv_base <- auc_numeric(data$outcome, loocv_predictions(data, base_column))
  loocv_augmented <- auc_numeric(data$outcome, loocv_predictions(data, c(base_column, "FCN1_TNFSF12_sd")))
  boot <- bootstrap_performance(data, base_column, 1000L)
  data.table(
    cohort = unique(data$cohort),
    base_model = base_label,
    apparent_AUC_base = apparent_base,
    apparent_AUC_augmented = apparent_augmented,
    apparent_Delta_AUC = apparent_augmented - apparent_base,
    LOOCV_AUC_base = loocv_base,
    LOOCV_AUC_augmented = loocv_augmented,
    LOOCV_Delta_AUC = loocv_augmented - loocv_base,
    optimism_corrected_AUC_base = boot$corrected_base,
    optimism_corrected_AUC_augmented = boot$corrected_augmented,
    optimism_corrected_Delta_AUC = boot$corrected_augmented - boot$corrected_base,
    bootstrap_test_Delta_median = boot$delta_median,
    bootstrap_test_Delta_CI_lower = boot$delta_lower,
    bootstrap_test_Delta_CI_upper = boot$delta_upper,
    successful_bootstrap_replicates = boot$successful,
    requested_bootstrap_replicates = boot$requested
  )
}

set.seed(20260721)
performance_results <- rbindlist(list(
  calculate_performance(data_25780, "Stroma_ECM_sd", "Stroma-ECM", 0L),
  calculate_performance(data_25780, "CD8_inverted_sd", "CD8 T", 1L),
  calculate_performance(data_40416, "Stroma_ECM_sd", "Stroma-ECM", 2L),
  calculate_performance(data_40416, "CD8_inverted_sd", "CD8 T", 3L)
))

expected_auc <- data.table(
  cohort = rep(c("PRJEB25780", "PRJEB40416"), each = 3L),
  marker = rep(c("FCN1/TNFSF12 score", "Stroma-ECM score", "CD8 T score (inverted)"), 2L),
  expected_AUC = c(0.6250000, 0.8303571, 0.8303571, 0.8035714, 0.6071429, 0.5892857)
)
auc_validation <- merge(direct_results, expected_auc, by = c("cohort", "marker"), sort = FALSE)
auc_validation[, pass := abs(AUC - expected_AUC) < 1e-6]
if (!all(auc_validation$pass)) {
  print(auc_validation)
  stop("Direct-score AUC values did not reproduce the locked results.")
}

result_object <- list(
  direct_AUC = direct_results,
  paired_DeLong = delong_results,
  nested_Firth = nested_results,
  internal_performance = performance_results,
  validation = auc_validation,
  settings = list(AUC_CI_replicates = 2000L, optimism_replicates = 1000L, seed = 20260721L)
)
saveRDS(result_object, file.path(output_dir, "incremental_value_results.rds"))
fwrite(direct_results, file.path(output_dir, "incremental_value_direct_AUC.tsv"), sep = "\t")
fwrite(delong_results, file.path(output_dir, "incremental_value_paired_DeLong.tsv"), sep = "\t")
fwrite(nested_results, file.path(output_dir, "incremental_value_nested_Firth.tsv"), sep = "\t")
fwrite(performance_results, file.path(output_dir, "incremental_value_internal_performance.tsv"), sep = "\t")

roc_to_data <- function(roc_object, cohort, marker) {
  coordinates <- coords(roc_object, x = "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  data.table(
    cohort = cohort,
    marker = marker,
    false_positive_rate = 1 - coordinates$specificity,
    true_positive_rate = coordinates$sensitivity
  )
}

roc_plot_data <- rbindlist(list(
  roc_to_data(roc_objects$PRJEB25780$FCN1_TNFSF12, "PRJEB25780", "FCN1/TNFSF12"),
  roc_to_data(roc_objects$PRJEB25780$Stroma_ECM, "PRJEB25780", "Stroma-ECM"),
  roc_to_data(roc_objects$PRJEB25780$CD8_inverted, "PRJEB25780", "CD8 T (inverted)"),
  roc_to_data(roc_objects$PRJEB40416$FCN1_TNFSF12, "PRJEB40416", "FCN1/TNFSF12"),
  roc_to_data(roc_objects$PRJEB40416$Stroma_ECM, "PRJEB40416", "Stroma-ECM"),
  roc_to_data(roc_objects$PRJEB40416$CD8_inverted, "PRJEB40416", "CD8 T (inverted)")
))

legend_labels <- direct_results[, .(label = paste0(sub(" score", "", marker), " (AUC = ", sprintf("%.3f", AUC), ")")), by = .(cohort, marker)]
roc_plot_data <- merge(roc_plot_data, legend_labels, by = c("cohort", "marker"), all.x = TRUE)

plot_object <- ggplot(roc_plot_data, aes(false_positive_rate, true_positive_rate, linetype = label)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.45, color = "grey65") +
  geom_step(linewidth = 0.85) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = "1 - Specificity", y = "Sensitivity", linetype = NULL) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    panel.spacing.x = grid::unit(2.0, "lines"),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(face = "bold")
  )

ggsave(file.path(figure_dir, "FigS12_exploratory_ROC_discrimination.png"), plot_object, width = 8.5, height = 4.8, units = "in", dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "FigS12_exploratory_ROC_discrimination.tiff"), plot_object, width = 8.5, height = 4.8, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")

report <- c(
  "INCREMENTAL-VALUE REPRODUCTION AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Positive event: non-response",
  "AUC CI bootstrap replicates: 2000",
  "Optimism-correction replicates: 1000",
  "Base random seed: 20260721",
  "",
  "DIRECT AUC",
  capture.output(print(direct_results)),
  "",
  "PAIRED DELONG",
  capture.output(print(delong_results)),
  "",
  "NESTED FIRTH",
  capture.output(print(nested_results)),
  "",
  "INTERNAL PERFORMANCE",
  capture.output(print(performance_results)),
  "",
  "AUC VALIDATION",
  capture.output(print(auc_validation))
)
writeLines(report, file.path(audit_dir, "incremental_value_reproduction_audit.txt"), useBytes = TRUE)

cat("Incremental-value analyses completed.\n")
print(direct_results)
print(delong_results)
print(nested_results)
print(performance_results)
