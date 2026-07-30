# ============================================================
# 08_anti_PD1_supplementary_figures_S3_S5.R
#
# Rebuilds final supplementary figures:
#   S3: FCN1/TNFSF12 versus Stroma-ECM in both anti-PD-1 cohorts.
#   S5: Cohort-specific correlations with Stroma-ECM, CD8 T and
#       Cytotoxic programs.
#
# Required input:
#   04_results/PRJEB40416_locked_signature_reanalysis_results.rds
#
# Outputs (600 dpi PNG and LZW-compressed TIFF):
#   05_figures/FigS3_ICI_FCN1_TNFSF12_vs_Stroma_ECM_scatter.*
#   05_figures/FigS5_ICI_cohort_signature_correlation_summary.*
# ============================================================

rm(list = ls())
graphics.off()

required_packages <- c("ggplot2", "cowplot", "data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg)
  }
}
library(ggplot2)
library(cowplot)
library(data.table)

config_file <- file.path("config", "project_config.R")
if (!file.exists(config_file)) {
  stop("Missing config/project_config.R. Run from the project root.")
}
source(config_file)

input_rds <- file.path(project_dir, "04_results", "PRJEB40416_locked_signature_reanalysis_results.rds")
figure_dir <- file.path(project_dir, "05_figures")
audit_dir <- file.path(project_dir, "06_logs_and_audit")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) stop("Missing locked score object: ", input_rds)
locked_object <- readRDS(input_rds)
if (!all(c("PRJEB25780_scores", "PRJEB40416_scores") %in% names(locked_object))) {
  stop("Locked RDS does not contain both cohort score tables.")
}

normalize_response <- function(x) {
  text <- trimws(as.character(x))
  upper <- toupper(text)
  out <- rep(NA_character_, length(upper))
  out[upper %in% c("R", "RESPONDER", "RESPONDERS", "RESPONSE", "CR", "PR", "1")] <- "Responder"
  out[upper %in% c("NR", "NON-RESPONDER", "NON_RESPONDER", "NONRESPONDER", "NON RESPONDER", "NON-RESPONDERS", "SD", "PD", "0")] <- "Non-responder"
  out[upper %in% c("NE", "NOT EVALUABLE", "NOT-EVALUABLE", "NOT_EVALUABLE", "NA", "") | is.na(upper)] <- "Not evaluable"
  factor(out, levels = c("Responder", "Non-responder", "Not evaluable"))
}

prepare_scores <- function(x, cohort) {
  x <- as.data.frame(x, check.names = FALSE)
  needed <- c("sample_id", "response_group", "FCN1_TNFSF12_score", "Stroma_ECM_score", "CD8_T_score", "Cytotoxic_score")
  missing <- setdiff(needed, names(x))
  if (length(missing) > 0L) stop(cohort, " score table is missing: ", paste(missing, collapse = ", "))
  for (column_name in setdiff(needed, c("sample_id", "response_group"))) {
    x[[column_name]] <- suppressWarnings(as.numeric(x[[column_name]]))
  }
  x$response_group <- normalize_response(x$response_group)
  x$cohort <- cohort
  x
}

scores_25780 <- prepare_scores(locked_object$PRJEB25780_scores, "PRJEB25780")
scores_40416 <- prepare_scores(locked_object$PRJEB40416_scores, "PRJEB40416")

calculate_spearman <- function(data, outcome_column) {
  complete <- is.finite(data$FCN1_TNFSF12_score) & is.finite(data[[outcome_column]])
  data <- data[complete, , drop = FALSE]
  test <- suppressWarnings(cor.test(
    data$FCN1_TNFSF12_score,
    data[[outcome_column]],
    method = "spearman",
    exact = FALSE,
    alternative = "two.sided"
  ))
  list(n = nrow(data), rho = unname(test$estimate), p = test$p.value)
}

format_p <- function(x) {
  if (!is.finite(x)) return("NA")
  if (x < 0.001) return(formatC(x, format = "e", digits = 2))
  sprintf("%.3f", x)
}

response_colors <- c(
  "Responder" = "#4E79A7",
  "Non-responder" = "#E15759",
  "Not evaluable" = "#9B9B9B"
)

make_scatter_panel <- function(data, cohort_title) {
  statistics <- calculate_spearman(data, "Stroma_ECM_score")
  x_range <- diff(range(data$FCN1_TNFSF12_score, na.rm = TRUE))
  y_range <- diff(range(data$Stroma_ECM_score, na.rm = TRUE))
  if (!is.finite(x_range) || x_range == 0) x_range <- 1
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  x_min <- min(data$FCN1_TNFSF12_score, na.rm = TRUE)
  y_max <- max(data$Stroma_ECM_score, na.rm = TRUE)
  label <- paste0(
    "Spearman rho = ", sprintf("%.3f", statistics$rho),
    "\nP = ", format_p(statistics$p),
    "\nn = ", statistics$n
  )

  ggplot(data, aes(x = FCN1_TNFSF12_score, y = Stroma_ECM_score)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.85, color = "black") +
    geom_point(
      aes(fill = response_group),
      shape = 21,
      size = 3.0,
      stroke = 0.40,
      color = "black",
      alpha = 0.90,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = x_min + 0.05 * x_range,
      y = y_max - 0.05 * y_range,
      label = label,
      hjust = 0,
      vjust = 1,
      family = "Times New Roman",
      size = 4.0,
      lineheight = 1.08
    ) +
    scale_fill_manual(values = response_colors, limits = names(response_colors), drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.08))) +
    labs(title = cohort_title, x = "FCN1/TNFSF12 score", y = "Stroma-ECM score") +
    theme_classic(base_size = 12, base_family = "Times New Roman") +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 7)),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11, color = "black"),
      axis.line = element_line(linewidth = 0.75, color = "black"),
      axis.ticks = element_line(linewidth = 0.60, color = "black"),
      legend.position = "none",
      plot.margin = margin(t = 8, r = 10, b = 5, l = 8)
    )
}

panel_25780 <- make_scatter_panel(scores_25780, "PRJEB25780")
panel_40416 <- make_scatter_panel(scores_40416, "PRJEB40416")

legend_data <- data.frame(
  response_group = factor(names(response_colors), levels = names(response_colors)),
  x = seq_along(response_colors),
  y = 1
)
legend_source <- ggplot(legend_data, aes(x = x, y = y, fill = response_group)) +
  geom_point(shape = 21, size = 3.2, stroke = 0.45, color = "black") +
  scale_fill_manual(values = response_colors, limits = names(response_colors), drop = FALSE, name = "Response") +
  guides(fill = guide_legend(
    title.position = "left",
    nrow = 1,
    byrow = TRUE,
    override.aes = list(shape = 21, size = 3.2, alpha = 1, color = "black")
  )) +
  theme_void(base_family = "Times New Roman") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 11.5, face = "bold", margin = margin(r = 7)),
    legend.text = element_text(size = 11),
    legend.key = element_blank()
  )
shared_legend <- cowplot::get_legend(legend_source)

figure_s3 <- cowplot::plot_grid(
  cowplot::plot_grid(panel_25780, panel_40416, ncol = 2, align = "hv", axis = "tblr"),
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.105)
)

s3_png <- file.path(figure_dir, "FigS3_ICI_FCN1_TNFSF12_vs_Stroma_ECM_scatter.png")
s3_tiff <- file.path(figure_dir, "FigS3_ICI_FCN1_TNFSF12_vs_Stroma_ECM_scatter.tiff")
ggsave(s3_png, figure_s3, width = 10.4, height = 5.2, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(s3_tiff, figure_s3, width = 10.4, height = 5.2, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white", limitsize = FALSE)

correlation_targets <- c(
  Stroma_ECM_score = "Stroma-ECM",
  CD8_T_score = "CD8 T",
  Cytotoxic_score = "Cytotoxic"
)

build_heatmap_data <- function(data, cohort_name) {
  rbindlist(lapply(names(correlation_targets), function(column_name) {
    result <- calculate_spearman(data, column_name)
    data.table(
      Cohort = cohort_name,
      Program = correlation_targets[[column_name]],
      n = result$n,
      rho = result$rho,
      p_value = result$p
    )
  }))
}

heatmap_data <- rbindlist(list(
  build_heatmap_data(scores_25780, "PRJEB25780"),
  build_heatmap_data(scores_40416, "PRJEB40416")
))
heatmap_data[, label := paste0("rho = ", sprintf("%.2f", rho), "\nP = ", format_p(p_value))]
heatmap_data[, Cohort := factor(Cohort, levels = c("PRJEB25780", "PRJEB40416"))]
heatmap_data[, Program := factor(Program, levels = c("Stroma-ECM", "CD8 T", "Cytotoxic"))]

figure_s5 <- ggplot(heatmap_data, aes(x = Program, y = Cohort, fill = rho)) +
  geom_tile(color = "white", linewidth = 1.0) +
  geom_text(aes(label = label), family = "Times New Roman", size = 4.3, lineheight = 1.02) +
  scale_fill_gradient2(low = "#4E79A7", mid = "white", high = "#E15759", midpoint = 0, name = "Spearman\nrho") +
  labs(x = NULL, y = NULL) +
  coord_equal() +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(face = "bold"),
    plot.margin = margin(8, 8, 8, 8)
  )

s5_png <- file.path(figure_dir, "FigS5_ICI_cohort_signature_correlation_summary.png")
s5_tiff <- file.path(figure_dir, "FigS5_ICI_cohort_signature_correlation_summary.tiff")
ggsave(s5_png, figure_s5, width = 6.8, height = 3.4, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(s5_tiff, figure_s5, width = 6.8, height = 3.4, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white", limitsize = FALSE)

expected <- data.table(
  Cohort = rep(c("PRJEB25780", "PRJEB40416"), each = 3L),
  Program = rep(c("Stroma-ECM", "CD8 T", "Cytotoxic"), times = 2L),
  expected_rho = c(0.569169960474308, 0.3893281, 0.3675889, 0.391176470588235, 0.1294118, 0.2323529),
  expected_p = c(0.00458974107409376, 0.066320621, 0.084420641, 0.13407128619725525, 0.632882332, 0.386511350)
)
validation <- merge(heatmap_data, expected, by = c("Cohort", "Program"), sort = FALSE)
validation[, rho_pass := abs(rho - expected_rho) < 1e-6]
validation[, p_pass := abs(p_value - expected_p) < 1e-6]
if (!all(validation$rho_pass) || !all(validation$p_pass)) stop("S3/S5 locked statistics were not reproduced.")

report <- c(
  "SUPPLEMENTARY FIGURES S3 AND S5 REPRODUCTION AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  capture.output(print(validation)),
  "",
  paste0("S3 PNG: ", normalizePath(s3_png, winslash = "/", mustWork = FALSE)),
  paste0("S3 TIFF: ", normalizePath(s3_tiff, winslash = "/", mustWork = FALSE)),
  paste0("S5 PNG: ", normalizePath(s5_png, winslash = "/", mustWork = FALSE)),
  paste0("S5 TIFF: ", normalizePath(s5_tiff, winslash = "/", mustWork = FALSE))
)
writeLines(report, file.path(audit_dir, "FiguresS3_S5_reproduction_audit.txt"), useBytes = TRUE)

cat("Figures S3 and S5 generated successfully.\n")
print(validation)
