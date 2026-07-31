# ============================================================================
# Supplementary Figure S14
# Independent tissue-proteomic assessment of FCN1 protein in CPTAC-STAD
#
# Final submission version:
# 1. Combines the three locked PDC000614 panels without recalculating statistics.
# 2. Uses lowercase panel labels a, b and c.
# 3. Centers panel c at approximately 78% of the total width.
# 4. Writes the formal Figure S14 PNG and TIFF files.
# ============================================================================

options(stringsAsFactors = FALSE, digits = 15)
graphics.off()

required_packages <- c("ggplot2", "patchwork", "png", "grid")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    type = "binary"
  )
}
still_missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing) > 0L) {
  stop("以下R包安装或加载失败：", paste(still_missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
})

project_dir <- Sys.getenv(
  "GC_ICI_MR_NICHE_PROJECT",
  unset = paste0(
    "E:/Bioinformatics analysis/",
    "GC_ICI_MR_Niche_REPRO_20260617"
  )
)
figure_dir <- file.path(project_dir, "05_figures")
audit_dir <- file.path(project_dir, "06_logs_and_audit")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

resolve_panel_file <- function(candidate_paths, panel_name) {
  existing_paths <- candidate_paths[file.exists(candidate_paths)]
  if (length(existing_paths) == 0L) {
    stop(
      "找不到", panel_name, "的任何候选文件：\n",
      paste(candidate_paths, collapse = "\n")
    )
  }
  existing_paths[1]
}

panel_a_file <- resolve_panel_file(
  c(
    file.path(figure_dir, "FigS14A_PDC000614_FCN1_neutrophil_scatter.png"),
    file.path(figure_dir, "FigS13A_PDC000614_FCN1_neutrophil_scatter.png")
  ),
  "Panel a"
)
panel_b_file <- resolve_panel_file(
  c(
    file.path(figure_dir, "FigS14B_PDC000614_FCN1_monocyte_scatter.png"),
    file.path(figure_dir, "FigS13B_PDC000614_FCN1_monocyte_scatter.png")
  ),
  "Panel b"
)
panel_c_file <- resolve_panel_file(
  c(
    file.path(figure_dir, "FigS14C_PDC000614_FCN1_key_correlations_summary.png"),
    file.path(figure_dir, "FigS13C_PDC000614_FCN1_key_correlations_summary.png")
  ),
  "Panel c"
)

output_png <- file.path(
  figure_dir,
  "FigureS14_PDC000614_FCN1_proteomic_validation.png"
)
output_tiff <- file.path(
  figure_dir,
  "FigureS14_PDC000614_FCN1_proteomic_validation.tiff"
)
audit_file <- file.path(
  audit_dir,
  "FigureS14_PDC000614_FCN1_proteomic_assessment_audit.txt"
)

create_raster_panel <- function(image_path, panel_tag) {
  image_array <- png::readPNG(image_path)
  image_height <- dim(image_array)[1]
  image_width <- dim(image_array)[2]
  if (!is.finite(image_height) || !is.finite(image_width) ||
      image_height <= 0 || image_width <= 0) {
    stop("无法读取子图尺寸：", image_path)
  }
  image_grob <- grid::rasterGrob(
    image = image_array,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    interpolate = TRUE
  )
  ggplot() +
    annotation_custom(
      grob = image_grob,
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1
    ) +
    annotate(
      geom = "text",
      x = 0.010,
      y = 0.990,
      label = panel_tag,
      hjust = 0,
      vjust = 1,
      family = "Times New Roman",
      fontface = "bold",
      size = 7.0,
      color = "black"
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void(base_family = "Times New Roman") +
    theme(
      aspect.ratio = image_height / image_width,
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    )
}

panel_a <- create_raster_panel(panel_a_file, "a")
panel_b <- create_raster_panel(panel_b_file, "b")
panel_c <- create_raster_panel(panel_c_file, "c")

top_row <- panel_a + panel_b +
  plot_layout(ncol = 2, widths = c(1, 1))

bottom_row <- (plot_spacer() + panel_c + plot_spacer()) +
  plot_layout(ncol = 3, widths = c(0.11, 0.78, 0.11))

figure_s14 <- top_row / bottom_row +
  plot_layout(ncol = 1, heights = c(1.00, 1.25)) &
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
  )

ggsave(
  filename = output_png,
  plot = figure_s14,
  width = 13.2,
  height = 12.4,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  filename = output_tiff,
  plot = figure_s14,
  width = 13.2,
  height = 12.4,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)

if (!file.exists(output_png) || file.info(output_png)$size <= 10000) {
  stop("Supplementary Figure S14 PNG生成失败。")
}
if (!file.exists(output_tiff) || file.info(output_tiff)$size <= 10000) {
  stop("Supplementary Figure S14 TIFF生成失败。")
}

audit_lines <- c(
  "SUPPLEMENTARY FIGURE S14 COMPOSITION AUDIT",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "INPUT PANELS",
  paste0("Panel a: ", panel_a_file),
  paste0("Panel b: ", panel_b_file),
  paste0("Panel c: ", panel_c_file),
  "",
  "MODIFICATIONS",
  "Panel labels are lowercase a, b and c.",
  "Panel c is centered at approximately 78% of total width.",
  "No statistical result was recalculated.",
  "",
  "OUTPUT FILES",
  output_png,
  output_tiff
)
writeLines(audit_lines, audit_file, useBytes = TRUE)

cat("\nSupplementary Figure S14已重新合成。\n")
cat("PNG：", output_png, "\n", sep = "")
cat("TIFF：", output_tiff, "\n", sep = "")
cat("审计报告：", audit_file, "\n", sep = "")
print(figure_s14)
