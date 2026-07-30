# ============================================================
# 16_PRJEB25780_supplementary_figure_S7_heatmap.R
#
# Generates the locked Supplementary Figure S7:
# row-scaled VST heatmap of the top 20 responder-high and
# top 20 non-responder-high genes in PRJEB25780.
#
# Figure S6 is intentionally not regenerated here. The archived all-gene
# DESeq2 table does not reproduce the exact historical volcano-plot layer,
# so the already locked formal Figure S6 must be retained unchanged.
#
# Required input, searched in this order:
#   1. <project>/04_results/GitHub_reproduction/
#        PRJEB25780_DESeq2_GSEA_results.rds
#   2. <repository>/results_summary/
#        PRJEB25780_DESeq2_GSEA_results.rds
#
# Formal outputs, overwritten in place:
#   <project>/05_figures/FigS7_PRJEB25780_top_DEG_heatmap.png
#   <project>/05_figures/FigS7_PRJEB25780_top_DEG_heatmap.tiff
#
# PNG and TIFF are saved at 600 dpi.
# ============================================================

rm(list = ls())
graphics.off()

required_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Required packages are not installed: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# Locate the current script and repository root when possible
# ------------------------------------------------------------

get_current_script_path <- function() {
  command_arguments <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_arguments, value = TRUE)

  if (length(file_argument) > 0L) {
    candidate <- sub("^--file=", "", file_argument[[1L]])

    if (file.exists(candidate)) {
      return(
        normalizePath(
          candidate,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  if (
    requireNamespace("rstudioapi", quietly = TRUE) &&
    isTRUE(
      tryCatch(
        rstudioapi::isAvailable(),
        error = function(e) FALSE
      )
    )
  ) {
    candidate <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )

    if (nzchar(candidate) && file.exists(candidate)) {
      return(
        normalizePath(
          candidate,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  NA_character_
}

current_script_path <- get_current_script_path()

repository_root <- if (!is.na(current_script_path)) {
  normalizePath(
    file.path(dirname(current_script_path), ".."),
    winslash = "/",
    mustWork = FALSE
  )
} else {
  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = FALSE
  )
}

# ------------------------------------------------------------
# Project path
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

project_input_rds <- file.path(
  project_dir,
  "04_results",
  "GitHub_reproduction",
  "PRJEB25780_DESeq2_GSEA_results.rds"
)

repository_input_rds <- file.path(
  repository_root,
  "results_summary",
  "PRJEB25780_DESeq2_GSEA_results.rds"
)

input_candidates <- c(
  project_input_rds,
  repository_input_rds
)

existing_inputs <- input_candidates[file.exists(input_candidates)]

if (length(existing_inputs) == 0L) {
  stop(
    "Could not locate PRJEB25780_DESeq2_GSEA_results.rds. Searched:\n- ",
    paste(input_candidates, collapse = "\n- ")
  )
}

input_rds <- normalizePath(
  existing_inputs[[1L]],
  winslash = "/",
  mustWork = TRUE
)

if (dir.exists(project_dir)) {
  figure_dir <- file.path(project_dir, "05_figures")
  audit_dir <- file.path(project_dir, "06_logs_and_audit")
} else {
  figure_dir <- file.path(repository_root, "reproduced_figures")
  audit_dir <- file.path(repository_root, "reproduction_audit")
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load and validate the locked object
# ------------------------------------------------------------

analysis_object <- readRDS(input_rds)

required_objects <- c(
  "sample_metadata",
  "DESeq2_top20_each_direction",
  "vst_matrix"
)

missing_objects <- setdiff(
  required_objects,
  names(analysis_object)
)

if (length(missing_objects) > 0L) {
  stop(
    "The locked RDS is missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

top_degs <- as.data.table(
  analysis_object$DESeq2_top20_each_direction
)

sample_metadata <- as.data.table(
  analysis_object$sample_metadata
)

vst_matrix <- analysis_object$vst_matrix

if (
  nrow(top_degs[direction == "R-high"]) != 20L ||
  nrow(top_degs[direction == "NR-high"]) != 20L
) {
  stop(
    "Figure S7 requires exactly 20 R-high and 20 NR-high genes."
  )
}

if (nrow(sample_metadata) != 23L) {
  stop(
    "Figure S7 requires the locked 23-patient PRJEB25780 cohort."
  )
}

# ------------------------------------------------------------
# Match the 40 locked genes to the VST matrix
# ------------------------------------------------------------

top_gene_ids <- as.character(top_degs$gene_id)
vst_row_ids <- rownames(vst_matrix)

vst_match <- match(
  top_gene_ids,
  vst_row_ids
)

missing_positions <- which(is.na(vst_match))

if (length(missing_positions) > 0L) {
  top_unversioned <- sub(
    "\\.[0-9]+$",
    "",
    top_gene_ids[missing_positions]
  )

  vst_unversioned <- sub(
    "\\.[0-9]+$",
    "",
    vst_row_ids
  )

  for (position_index in seq_along(missing_positions)) {
    target_position <- missing_positions[[position_index]]

    candidate_rows <- which(
      vst_unversioned == top_unversioned[[position_index]]
    )

    if (length(candidate_rows) == 1L) {
      vst_match[[target_position]] <- candidate_rows[[1L]]
    }
  }
}

if (anyNA(vst_match)) {
  missing_vst_genes <- top_gene_ids[is.na(vst_match)]

  stop(
    "Top DEG rows are missing or ambiguous in the VST matrix: ",
    paste(missing_vst_genes, collapse = ", ")
  )
}

if (anyDuplicated(vst_match) > 0L) {
  stop(
    "Two locked genes mapped to the same VST row."
  )
}

matched_vst_gene_ids <- vst_row_ids[vst_match]

# ------------------------------------------------------------
# Lock sample order
# ------------------------------------------------------------

if (!"count_sample_id" %in% names(sample_metadata)) {
  sample_metadata[, count_sample_id := sample_id]
}

sample_metadata[, response_group := factor(
  as.character(response_group),
  levels = c(
    "Responder",
    "Non-responder"
  )
)]

setorder(
  sample_metadata,
  response_group,
  sample_id
)

sample_order <- as.character(
  sample_metadata$count_sample_id
)

sample_labels <- setNames(
  as.character(sample_metadata$sample_id),
  sample_order
)

missing_vst_samples <- setdiff(
  sample_order,
  colnames(vst_matrix)
)

if (length(missing_vst_samples) > 0L) {
  stop(
    "Locked samples are missing from the VST matrix: ",
    paste(missing_vst_samples, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Row-scaled heatmap matrix
# ------------------------------------------------------------

heatmap_matrix <- vst_matrix[
  matched_vst_gene_ids,
  sample_order,
  drop = FALSE
]

rownames(heatmap_matrix) <- top_gene_ids

row_means <- rowMeans(heatmap_matrix)
row_sds <- apply(heatmap_matrix, 1L, sd)
row_sds[!is.finite(row_sds) | row_sds == 0] <- 1

heatmap_z <- sweep(
  heatmap_matrix,
  1L,
  row_means,
  "-"
)

heatmap_z <- sweep(
  heatmap_z,
  1L,
  row_sds,
  "/"
)

heatmap_z[heatmap_z > 2.5] <- 2.5
heatmap_z[heatmap_z < -2.5] <- -2.5

heatmap_long <- as.data.table(
  as.table(heatmap_z)
)

setnames(
  heatmap_long,
  c(
    "gene_id",
    "sample_id",
    "row_Z"
  )
)

heatmap_long <- merge(
  heatmap_long,
  top_degs[
    ,
    .(
      gene_id,
      gene_label,
      direction
    )
  ],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(heatmap_long$gene_label)) {
  stop(
    "One or more heatmap rows could not be assigned a display label."
  )
}

# Preserve the locked top-gene order: R-high block first, then NR-high.
gene_display_order <- rev(
  as.character(top_degs$gene_label)
)

heatmap_long[, gene_label := factor(
  gene_label,
  levels = gene_display_order
)]

heatmap_long[, sample_id := factor(
  sample_id,
  levels = sample_order
)]

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

response_colors <- c(
  "Responder" = "#D95F59",
  "Non-responder" = "#4C78A8"
)

annotation_data <- sample_metadata[
  ,
  .(
    sample_id = factor(
      count_sample_id,
      levels = sample_order
    ),
    response_group
  )
]

annotation_panel <- ggplot(
  annotation_data,
  aes(
    x = sample_id,
    y = "Response",
    fill = response_group
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.20
  ) +
  scale_fill_manual(
    values = response_colors,
    breaks = c(
      "Responder",
      "Non-responder"
    ),
    name = "Response",
    drop = FALSE
  ) +
  scale_x_discrete(
    labels = sample_labels,
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  theme_void(
    base_family = "Times New Roman"
  ) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(
      size = 9.5
    ),
    plot.margin = margin(
      0,
      8,
      0,
      82
    )
  )

heatmap_panel <- ggplot(
  heatmap_long,
  aes(
    x = sample_id,
    y = gene_label,
    fill = row_Z
  )
) +
  geom_tile(
    color = NA
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -2.5,
      2.5
    ),
    breaks = c(
      -2,
      0,
      2
    ),
    name = "Row Z-score"
  ) +
  scale_x_discrete(
    labels = sample_labels,
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 10,
    base_family = "Times New Roman"
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 7.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 7.7,
      color = "black"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 9.5
    ),
    legend.text = element_text(
      size = 8.5
    ),
    plot.margin = margin(
      0,
      8,
      6,
      6
    )
  )

figure_s7 <- annotation_panel / heatmap_panel +
  plot_layout(
    heights = c(
      0.65,
      10.0
    ),
    guides = "keep"
  )

output_s7_png <- file.path(
  figure_dir,
  "FigS7_PRJEB25780_top_DEG_heatmap.png"
)

output_s7_tiff <- file.path(
  figure_dir,
  "FigS7_PRJEB25780_top_DEG_heatmap.tiff"
)

ggsave(
  filename = output_s7_png,
  plot = figure_s7,
  width = 8.8,
  height = 9.6,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = output_s7_tiff,
  plot = figure_s7,
  width = 8.8,
  height = 9.6,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)

if (!file.exists(output_s7_png)) {
  stop("Figure S7 PNG was not generated.")
}

if (!file.exists(output_s7_tiff)) {
  stop("Figure S7 TIFF was not generated.")
}

# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------

output_audit <- file.path(
  audit_dir,
  "PRJEB25780_FigureS7_reproduction_audit.txt"
)

audit_lines <- c(
  "PRJEB25780 Supplementary Figure S7 reproduction audit",
  "============================================================",
  paste0(
    "Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S %Z"
    )
  ),
  paste0("Input RDS: ", input_rds),
  "",
  "Figure S7:",
  paste0("- genes: ", length(top_gene_ids)),
  paste0("- samples: ", length(sample_order)),
  paste0(
    "- responders: ",
    sum(sample_metadata$response_group == "Responder")
  ),
  paste0(
    "- non-responders: ",
    sum(sample_metadata$response_group == "Non-responder")
  ),
  paste0("- PNG: ", output_s7_png),
  paste0("- TIFF: ", output_s7_tiff),
  "- resolution: 600 dpi",
  "",
  paste0(
    "Figure S6 was intentionally not regenerated; ",
    "the locked formal S6 remains unchanged."
  )
)

writeLines(
  audit_lines,
  output_audit,
  useBytes = TRUE
)

cat("\nSupplementary Figure S7 generated.\n")
cat("S7 PNG:", output_s7_png, "\n")
cat("S7 TIFF:", output_s7_tiff, "\n")
cat("Audit:", output_audit, "\n")
cat("Figure S6 was not modified.\n")
