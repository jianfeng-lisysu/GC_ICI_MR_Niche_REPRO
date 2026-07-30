# ============================================================
# 09_anti_PD1_supplementary_figure_S4_response_groups.R
#
# Rebuilds final Supplementary Figure S4:
# response-group comparisons of FCN1/TNFSF12, Stroma-ECM,
# CD8 T and Cytotoxic scores in PRJEB25780 and PRJEB40416.
#
# Required input:
#   04_results/PRJEB40416_locked_signature_reanalysis_results.rds
#
# Outputs:
#   05_figures/FigS4_ICI_response_group_comparisons.png
#   05_figures/FigS4_ICI_response_group_comparisons.tiff
# ============================================================

rm(list = ls())
graphics.off()

required_packages <- c("ggplot2", "data.table")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package is not installed: ", pkg)
}
library(ggplot2)
library(data.table)

config_file <- file.path("config", "project_config.R")
if (!file.exists(config_file)) stop("Missing config/project_config.R. Run from the project root.")
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
  factor(out, levels = c("Responder", "Non-responder"))
}

required_columns <- c(
  "sample_id",
  "response_group",
  "FCN1_TNFSF12_score",
  "Stroma_ECM_score",
  "CD8_T_score",
  "Cytotoxic_score"
)

prepare_scores <- function(x, cohort) {
  x <- as.data.frame(x, check.names = FALSE)
  missing <- setdiff(required_columns, names(x))
  if (length(missing) > 0L) stop(cohort, " score table is missing: ", paste(missing, collapse = ", "))
  for (column_name in setdiff(required_columns, c("sample_id", "response_group"))) {
    x[[column_name]] <- suppressWarnings(as.numeric(x[[column_name]]))
  }
  x$response_group <- normalize_response(x$response_group)
  x <- x[!is.na(x$response_group), , drop = FALSE]
  x$Cohort <- cohort
  x
}

scores_25780 <- prepare_scores(locked_object$PRJEB25780_scores, "PRJEB25780")
scores_40416 <- prepare_scores(locked_object$PRJEB40416_scores, "PRJEB40416")

counts_25780 <- table(scores_25780$response_group)
counts_40416 <- table(scores_40416$response_group)
if (nrow(scores_25780) != 23L || as.integer(counts_25780["Responder"]) != 7L || as.integer(counts_25780["Non-responder"]) != 16L) {
  stop("PRJEB25780 response-evaluable composition does not match the locked analysis.")
}
if (nrow(scores_40416) != 15L || as.integer(counts_40416["Responder"]) != 7L || as.integer(counts_40416["Non-responder"]) != 8L) {
  stop("PRJEB40416 response-evaluable composition does not match the locked analysis.")
}

signature_lookup <- c(
  FCN1_TNFSF12_score = "FCN1/TNFSF12",
  Stroma_ECM_score = "Stroma-ECM",
  CD8_T_score = "CD8 T",
  Cytotoxic_score = "Cytotoxic"
)

build_long <- function(data) {
  rbindlist(lapply(names(signature_lookup), function(column_name) {
    data.table(
      sample_id = data$sample_id,
      Cohort = data$Cohort,
      Response = data$response_group,
      Signature = signature_lookup[[column_name]],
      Score = data[[column_name]]
    )
  }))
}

plot_data <- rbindlist(list(build_long(scores_25780), build_long(scores_40416)))
plot_data <- plot_data[is.finite(Score) & !is.na(Response)]
plot_data[, Cohort := factor(Cohort, levels = c("PRJEB25780", "PRJEB40416"))]
plot_data[, Signature := factor(Signature, levels = c("FCN1/TNFSF12", "Stroma-ECM", "CD8 T", "Cytotoxic"))]
plot_data[, Response := factor(Response, levels = c("Responder", "Non-responder"))]

statistics <- plot_data[, {
  test <- suppressWarnings(wilcox.test(Score ~ Response, exact = FALSE, alternative = "two.sided"))
  value_range <- range(Score, na.rm = TRUE)
  span <- diff(value_range)
  if (!is.finite(span) || span == 0) span <- 1
  .(
    n = .N,
    Responders = sum(Response == "Responder"),
    Nonresponders = sum(Response == "Non-responder"),
    W = unname(test$statistic),
    P_value = test$p.value,
    annotation_y = value_range[2] + 0.13 * span
  )
}, by = .(Cohort, Signature)]
statistics[, BH_q := p.adjust(P_value, method = "BH"), by = Cohort]
statistics[, Annotation := paste0("Wilcoxon P = ", sprintf("%.3f", P_value), "\nn = ", n)]
statistics[, annotation_x := 1.5]

expected <- data.table(
  Cohort = rep(c("PRJEB25780", "PRJEB40416"), each = 4L),
  Signature = rep(c("FCN1/TNFSF12", "Stroma-ECM", "CD8 T", "Cytotoxic"), times = 2L),
  Expected_P = c(
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
validation <- merge(statistics, expected, by = c("Cohort", "Signature"), sort = FALSE)
validation[, pass := abs(P_value - Expected_P) < 1e-6]
if (!all(validation$pass)) {
  print(validation)
  stop("Locked response-group P values were not reproduced.")
}

response_colors <- c("Responder" = "#4E79A7", "Non-responder" = "#E15759")

figure_s4 <- ggplot(plot_data, aes(x = Response, y = Score, fill = Response)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.35, color = "grey70") +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.38, linewidth = 0.70, color = "black") +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 20260721),
    shape = 21,
    size = 2.20,
    stroke = 0.38,
    color = "black",
    alpha = 0.92
  ) +
  geom_text(
    data = statistics,
    aes(x = annotation_x, y = annotation_y, label = Annotation),
    inherit.aes = FALSE,
    family = "Times New Roman",
    size = 3.15,
    lineheight = 1.04,
    hjust = 0.5,
    vjust = 0
  ) +
  facet_grid(rows = vars(Cohort), cols = vars(Signature), scales = "free_y", switch = "y") +
  scale_fill_manual(values = response_colors, breaks = names(response_colors), drop = FALSE) +
  scale_x_discrete(labels = c("Responder" = "Responder", "Non-responder" = "Non-\nresponder")) +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.23))) +
  labs(x = NULL, y = "Standardized signature score", fill = "Response") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12.5, base_family = "Times New Roman") +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(0.85, "lines"),
    panel.spacing.y = grid::unit(0.75, "lines"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.70),
    strip.text.x = element_text(face = "bold", size = 12.5, color = "black", margin = margin(t = 5, r = 4, b = 5, l = 4)),
    strip.text.y.left = element_text(face = "bold", size = 12.5, color = "black", angle = 90, margin = margin(t = 4, r = 5, b = 4, l = 5)),
    strip.placement = "outside",
    axis.title.y = element_text(face = "bold", size = 13.5, margin = margin(r = 9)),
    axis.text.x = element_text(size = 9.5, color = "black", lineheight = 0.92, margin = margin(t = 4)),
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.line = element_line(linewidth = 0.65, color = "black"),
    axis.ticks = element_line(linewidth = 0.50, color = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 11.5),
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    plot.margin = margin(t = 10, r = 12, b = 6, l = 10)
  ) +
  guides(fill = guide_legend(
    title.position = "left",
    nrow = 1,
    byrow = TRUE,
    override.aes = list(shape = 21, size = 3.0, alpha = 1, color = "black")
  ))

output_png <- file.path(figure_dir, "FigS4_ICI_response_group_comparisons.png")
output_tiff <- file.path(figure_dir, "FigS4_ICI_response_group_comparisons.tiff")
ggsave(output_png, figure_s4, width = 11.2, height = 6.6, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(output_tiff, figure_s4, width = 11.2, height = 6.6, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white", limitsize = FALSE)

report <- c(
  "SUPPLEMENTARY FIGURE S4 REPRODUCTION AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  capture.output(print(validation)),
  "",
  paste0("PNG: ", normalizePath(output_png, winslash = "/", mustWork = FALSE)),
  paste0("TIFF: ", normalizePath(output_tiff, winslash = "/", mustWork = FALSE))
)
writeLines(report, file.path(audit_dir, "FigureS4_reproduction_audit.txt"), useBytes = TRUE)

cat("Figure S4 generated successfully.\n")
print(validation)
