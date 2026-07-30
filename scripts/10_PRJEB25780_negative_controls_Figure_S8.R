# ============================================================
# 10_PRJEB25780_negative_controls_Figure_S8.R
#
# Rebuilds final Supplementary Figure S8:
# negative-control and comparator correlations of the locked
# FCN1/TNFSF12 score in PRJEB25780.
#
# The final published graphic was redrawn from the locked summary
# statistics in Supplementary Table S15. The patient-level score
# vectors were not preserved in the repository package, so this
# script intentionally reconstructs the figure from those audited
# summary statistics and does not claim to rerun the correlations.
#
# Outputs:
#   05_figures/FigS8_PRJEB25780_negative_control_correlations.png
#   05_figures/FigS8_PRJEB25780_negative_control_correlations.tiff
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

figure_dir <- file.path(project_dir, "05_figures")
audit_dir <- file.path(project_dir, "06_logs_and_audit")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

locked_results <- data.table(
  Program = c("IFN-gamma", "Cell cycle", "Housekeeping"),
  Classification = c("Inflammatory comparator", "Negative control", "Negative control"),
  Spearman_rho = c(0.13, -0.33, -0.40),
  P_value = c(0.544, 0.128, 0.060),
  BH_adjusted_P = c(0.544, 0.191, 0.180),
  N = c(23L, 23L, 23L)
)
locked_results[, Program := factor(Program, levels = c("IFN-gamma", "Cell cycle", "Housekeeping"))]
locked_results[, Direction := ifelse(Spearman_rho >= 0, "Positive", "Negative")]
locked_results[, label_y := ifelse(Spearman_rho >= 0, Spearman_rho + 0.035, Spearman_rho - 0.035)]
locked_results[, label_vjust := ifelse(Spearman_rho >= 0, 0, 1)]
locked_results[, rho_label := paste0("rho = ", sprintf("%.2f", Spearman_rho))]

checks <- c(
  exactly_three_programs = nrow(locked_results) == 3L,
  all_n_23 = all(locked_results$N == 23L),
  all_adjusted_p_above_0_05 = all(locked_results$BH_adjusted_P > 0.05),
  IFN_gamma_rho = abs(locked_results[Program == "IFN-gamma", Spearman_rho] - 0.13) < 1e-12,
  Cell_cycle_rho = abs(locked_results[Program == "Cell cycle", Spearman_rho] + 0.33) < 1e-12,
  Housekeeping_rho = abs(locked_results[Program == "Housekeeping", Spearman_rho] + 0.40) < 1e-12
)
if (!all(checks)) stop("Locked negative-control summary statistics failed validation.")

figure_s8 <- ggplot(locked_results, aes(x = Program, y = Spearman_rho, fill = Direction)) +
  geom_hline(yintercept = 0, linewidth = 0.55, color = "black") +
  geom_col(width = 0.58, color = "grey35", linewidth = 0.55) +
  geom_text(
    aes(y = label_y, label = rho_label, vjust = label_vjust),
    family = "Times New Roman",
    size = 4.4,
    color = "black"
  ) +
  scale_fill_manual(values = c("Positive" = "#E15759", "Negative" = "#2CA6AE"), guide = "none") +
  scale_y_continuous(
    limits = c(-0.52, 0.28),
    breaks = seq(-0.4, 0.2, by = 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Spearman rho",
    caption = "Spearman correlation; all Benjamini-Hochberg-adjusted P > 0.05; n = 23 samples."
  ) +
  theme_classic(base_size = 13, base_family = "Times New Roman") +
  theme(
    axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 8)),
    axis.text.x = element_text(face = "bold", size = 12, color = "black", margin = margin(t = 7)),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.75, color = "black"),
    axis.ticks.x = element_blank(),
    plot.caption = element_text(size = 10.5, hjust = 0.5, margin = margin(t = 8)),
    plot.margin = margin(12, 14, 8, 12)
  )

output_png <- file.path(figure_dir, "FigS8_PRJEB25780_negative_control_correlations.png")
output_tiff <- file.path(figure_dir, "FigS8_PRJEB25780_negative_control_correlations.tiff")
ggsave(output_png, figure_s8, width = 5.8, height = 4.4, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(output_tiff, figure_s8, width = 5.8, height = 4.4, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white", limitsize = FALSE)

report <- c(
  "SUPPLEMENTARY FIGURE S8 REPRODUCTION AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Source: locked summary statistics from Supplementary Table S15.",
  "Patient-level vectors were not available in the repository package; statistics were not recomputed.",
  "",
  capture.output(print(locked_results)),
  "",
  capture.output(print(data.frame(item = names(checks), pass = checks, row.names = NULL))),
  "",
  paste0("PNG: ", normalizePath(output_png, winslash = "/", mustWork = FALSE)),
  paste0("TIFF: ", normalizePath(output_tiff, winslash = "/", mustWork = FALSE))
)
writeLines(report, file.path(audit_dir, "FigureS8_reproduction_audit.txt"), useBytes = TRUE)

cat("Figure S8 reconstructed from locked summary statistics.\n")
print(locked_results)
