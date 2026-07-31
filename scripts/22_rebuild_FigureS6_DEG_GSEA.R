# ============================================================================
# Supplementary Figure S6
# PRJEB25780 differential expression and pathway enrichment overview
#
# Final submission version:
# 1. Rebuilds the figure from locked DESeq2 and GSEA result files.
# 2. Labels five genes per response direction (10 total).
# 3. Uses italic HGNC gene symbols and upright Ensembl identifiers.
# 4. Prevents adjusted-P numerical underflow from expanding the y axis.
# 5. Determines a compact y-axis limit from the observed distribution.
# 6. Uses lowercase panel labels a and b.
# 7. Overwrites the formal PNG and TIFF files; no version suffixes are created.
# ============================================================================

options(stringsAsFactors = FALSE, digits = 15)
graphics.off()

required_packages <- c("data.table", "ggplot2", "ggrepel", "patchwork")
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
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

project_root <- Sys.getenv(
  "GC_ICI_MR_NICHE_PROJECT",
  unset = paste0(
    "E:/Bioinformatics analysis/",
    "GC_ICI_MR_Niche_REPRO_20260617"
  )
)
result_dir <- file.path(project_root, "04_results")
figure_dir <- file.path(project_root, "05_figures")
audit_dir <- file.path(project_root, "06_logs_and_audit")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

deseq_all_file <- file.path(result_dir, "DESeq2_PRJEB25780_clean_all.csv")
deseq_deg_file <- file.path(result_dir, "DESeq2_PRJEB25780_clean_DEGs.csv")
gsea_file <- file.path(
  result_dir,
  "GSEA_GO_BP_PRJEB25780_Responder_vs_Non_responder.csv"
)
output_png <- file.path(
  figure_dir,
  "FigS6_PRJEB25780_DEG_GSEA_overview.png"
)
output_tiff <- file.path(
  figure_dir,
  "FigS6_PRJEB25780_DEG_GSEA_overview.tiff"
)
audit_file <- file.path(
  audit_dir,
  "FigS6_PRJEB25780_DEG_GSEA_rebuild_audit.txt"
)

input_files <- c(deseq_all_file, deseq_deg_file, gsea_file)
missing_input_files <- input_files[!file.exists(input_files)]
if (length(missing_input_files) > 0L) {
  stop("缺少以下正式输入文件：\n", paste(missing_input_files, collapse = "\n"))
}

deseq_all <- fread(deseq_all_file, data.table = TRUE, showProgress = FALSE)
deseq_degs <- fread(deseq_deg_file, data.table = TRUE, showProgress = FALSE)
gsea_results <- fread(gsea_file, data.table = TRUE, showProgress = FALSE)

required_deseq_columns <- c(
  "baseMean", "log2FoldChange", "pvalue", "padj", "gene", "ENSEMBL", "SYMBOL"
)
required_gsea_columns <- c("ID", "Description", "NES", "p.adjust")

missing_all_columns <- setdiff(required_deseq_columns, names(deseq_all))
missing_deg_columns <- setdiff(required_deseq_columns, names(deseq_degs))
missing_gsea_columns <- setdiff(required_gsea_columns, names(gsea_results))
if (length(missing_all_columns) > 0L) {
  stop("DESeq2全基因结果缺少字段：", paste(missing_all_columns, collapse = ", "))
}
if (length(missing_deg_columns) > 0L) {
  stop("DESeq2 DEG结果缺少字段：", paste(missing_deg_columns, collapse = ", "))
}
if (length(missing_gsea_columns) > 0L) {
  stop("GSEA结果缺少字段：", paste(missing_gsea_columns, collapse = ", "))
}
if (nrow(deseq_degs) != 369L) {
  stop("清理后DEG结果应为369行，当前为：", nrow(deseq_degs))
}

standardize_deseq <- function(data) {
  output <- copy(data)
  for (current_column in c("baseMean", "log2FoldChange", "pvalue", "padj")) {
    output[, (current_column) := suppressWarnings(as.numeric(get(current_column)))]
  }
  output[, ensembl_id := trimws(as.character(ENSEMBL))]
  output[
    is.na(ensembl_id) | ensembl_id == "" | toupper(ensembl_id) == "NA",
    ensembl_id := sub("\\.[0-9]+$", "", trimws(as.character(gene)))
  ]
  output[, gene_symbol := trimws(as.character(SYMBOL))]
  output[
    is.na(gene_symbol) | gene_symbol == "" | toupper(gene_symbol) == "NA",
    gene_symbol := NA_character_
  ]
  output[, has_gene_symbol := !is.na(gene_symbol) & gene_symbol != ""]
  output[, gene_label := fifelse(has_gene_symbol == TRUE, gene_symbol, ensembl_id)]
  output[]
}

deseq_all <- standardize_deseq(deseq_all)
deseq_degs <- standardize_deseq(deseq_degs)

fdr_threshold <- 0.10
minimum_display_padj <- 1e-300
x_display_cap <- 16

deseq_all[, direction := fifelse(
  is.finite(padj) & padj < fdr_threshold & is.finite(log2FoldChange) & log2FoldChange > 0,
  "R-high",
  fifelse(
    is.finite(padj) & padj < fdr_threshold & is.finite(log2FoldChange) & log2FoldChange < 0,
    "NR-high",
    "Other genes"
  )
)]
deseq_all[, direction := factor(
  direction,
  levels = c("R-high", "NR-high", "Other genes")
)]

deseq_degs[, direction := fifelse(
  is.finite(log2FoldChange) & log2FoldChange > 0,
  "R-high",
  fifelse(
    is.finite(log2FoldChange) & log2FoldChange < 0,
    "NR-high",
    NA_character_
  )
)]
deseq_degs[, direction := factor(direction, levels = c("R-high", "NR-high"))]

add_plot_coordinates <- function(data, y_limit = NULL) {
  output <- copy(data)
  output[, adjusted_p_for_plot := fifelse(
    is.finite(padj) & padj > 0,
    pmax(padj, minimum_display_padj),
    fifelse(padj == 0, minimum_display_padj, 1)
  )]
  output[, raw_neg_log10_padj := -log10(adjusted_p_for_plot)]
  if (is.null(y_limit)) {
    finite_y <- output[is.finite(raw_neg_log10_padj), raw_neg_log10_padj]
    if (length(finite_y) == 0L) {
      stop("没有可用于确定火山图y轴范围的数值。")
    }
    robust_reference <- as.numeric(
      quantile(finite_y, probs = 0.999, na.rm = TRUE, names = FALSE, type = 7)
    )
    y_limit <- max(12, min(20, ceiling(robust_reference + 2)))
  }
  output[, x_was_truncated := is.finite(log2FoldChange) & abs(log2FoldChange) > x_display_cap]
  output[, y_was_truncated := is.finite(raw_neg_log10_padj) & raw_neg_log10_padj > y_limit]
  output[, plot_x := pmax(pmin(log2FoldChange, x_display_cap), -x_display_cap)]
  output[, plot_y := pmin(raw_neg_log10_padj, y_limit)]
  list(data = output, y_limit = y_limit)
}

all_plot_result <- add_plot_coordinates(deseq_all)
deseq_all <- all_plot_result$data
y_display_cap <- all_plot_result$y_limit
deseq_degs <- add_plot_coordinates(deseq_degs, y_limit = y_display_cap)$data

volcano_data <- deseq_all[is.finite(plot_x) & is.finite(plot_y)]
if (nrow(volcano_data) == 0L) {
  stop("没有可用于绘制火山图的数据。")
}

select_direction_labels <- function(
    data,
    direction_name,
    number_to_label = 5L,
    forced_symbols = character(0)
) {
  candidate_data <- copy(data[
    direction == direction_name &
      !is.na(gene_label) & gene_label != "" &
      is.finite(plot_x) & is.finite(plot_y)
  ])
  if (nrow(candidate_data) == 0L) {
    stop("没有找到可标注的", direction_name, "基因。")
  }
  candidate_data[, symbol_priority := fifelse(has_gene_symbol == TRUE, 0L, 1L)]
  candidate_data[, absolute_log2fc := abs(log2FoldChange)]
  candidate_data[, padj_sort := fifelse(is.finite(padj), padj, Inf)]
  setorder(candidate_data, symbol_priority, -absolute_log2fc, padj_sort)
  candidate_data <- candidate_data[!duplicated(gene_label)]
  forced_rows <- candidate_data[
    has_gene_symbol == TRUE & gene_symbol %in% forced_symbols
  ]
  other_rows <- candidate_data[!(gene_label %in% forced_rows$gene_label)]
  number_remaining <- max(0L, number_to_label - nrow(forced_rows))
  selected_rows <- rbindlist(
    list(forced_rows, head(other_rows, number_remaining)),
    use.names = TRUE,
    fill = TRUE
  )
  selected_rows <- head(selected_rows, number_to_label)
  selected_rows[, c("symbol_priority", "absolute_log2fc", "padj_sort") := NULL]
  selected_rows[]
}

labels_r <- select_direction_labels(
  deseq_degs,
  direction_name = "R-high",
  number_to_label = 5L
)
labels_nr <- select_direction_labels(
  deseq_degs,
  direction_name = "NR-high",
  number_to_label = 5L,
  forced_symbols = "FBLN1"
)
volcano_labels <- rbindlist(list(labels_r, labels_nr), use.names = TRUE, fill = TRUE)
if (nrow(volcano_labels) != 10L || anyDuplicated(volcano_labels$gene_label) > 0L) {
  stop("火山图标签数量或唯一性核验失败。")
}

symbol_labels <- volcano_labels[has_gene_symbol == TRUE]
ensembl_labels <- volcano_labels[has_gene_symbol == FALSE | is.na(has_gene_symbol)]

volcano_colors <- c(
  "R-high" = "#E56B63",
  "NR-high" = "#4C78A8",
  "Other genes" = "#C9C9C9"
)
caption_text <- NULL
if (any(volcano_data$y_was_truncated, na.rm = TRUE)) {
  caption_text <- paste0(
    "Values exceeding the displayed y-axis limit (",
    y_display_cap,
    ") are shown at the upper boundary."
  )
}

panel_a <- ggplot(volcano_data, aes(x = plot_x, y = plot_y)) +
  geom_hline(
    yintercept = -log10(fdr_threshold),
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey45"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.45,
    color = "grey45"
  ) +
  geom_point(
    data = volcano_data[direction == "Other genes" & y_was_truncated == FALSE],
    aes(color = direction),
    shape = 16,
    size = 1.10,
    alpha = 0.55,
    show.legend = TRUE
  ) +
  geom_point(
    data = volcano_data[direction != "Other genes" & y_was_truncated == FALSE],
    aes(color = direction),
    shape = 16,
    size = 1.50,
    alpha = 0.84,
    show.legend = TRUE
  ) +
  geom_point(
    data = volcano_data[y_was_truncated == TRUE],
    aes(color = direction),
    shape = 17,
    size = 2.10,
    alpha = 0.90,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = symbol_labels,
    aes(label = gene_label, color = direction),
    family = "Times New Roman",
    fontface = "italic",
    size = 3.20,
    box.padding = 0.80,
    point.padding = 0.40,
    min.segment.length = 0,
    segment.size = 0.35,
    segment.alpha = 0.75,
    force = 2.8,
    force_pull = 0.65,
    max.overlaps = Inf,
    seed = 20260731,
    direction = "both",
    xlim = c(-15.5, 15.5),
    ylim = c(0.8, y_display_cap - 0.25),
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = ensembl_labels,
    aes(label = gene_label, color = direction),
    family = "Times New Roman",
    fontface = "plain",
    size = 2.90,
    box.padding = 0.80,
    point.padding = 0.40,
    min.segment.length = 0,
    segment.size = 0.35,
    segment.alpha = 0.75,
    force = 2.8,
    force_pull = 0.65,
    max.overlaps = Inf,
    seed = 20260731,
    direction = "both",
    xlim = c(-15.5, 15.5),
    ylim = c(0.8, y_display_cap - 0.25),
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = volcano_colors,
    breaks = c("R-high", "NR-high", "Other genes"),
    drop = FALSE,
    name = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 16, size = 3, alpha = 1)
    )
  ) +
  scale_x_continuous(
    limits = c(-16.5, 16.5),
    breaks = seq(-15, 15, by = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, y_display_cap + 0.25),
    breaks = pretty(c(0, y_display_cap), n = 6),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = expression(log[2] * " fold change (Responder vs Non-responder)"),
    y = expression(-log[10] * " BH-adjusted P"),
    caption = caption_text
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Times New Roman", base_size = 11) +
  theme(
    axis.title = element_text(size = 11.5, color = "black"),
    axis.text = element_text(size = 9.5, color = "black"),
    axis.line = element_line(linewidth = 0.55, color = "black"),
    axis.ticks = element_line(linewidth = 0.50, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 9.5, color = "black"),
    legend.key.width = grid::unit(0.45, "cm"),
    plot.caption = element_text(
      size = 8,
      color = "grey30",
      hjust = 0,
      margin = margin(t = 7)
    ),
    plot.margin = margin(t = 10, r = 16, b = 8, l = 8)
  )

gsea_pathway_table <- data.table(
  description_key = c(
    "complement activation",
    "positive regulation of leukocyte mediated cytotoxicity",
    "positive regulation of cell killing",
    "natural killer cell mediated immunity",
    "cell killing",
    "dna replication",
    "dna-templated dna replication",
    "defense response to virus",
    "lymphocyte mediated immunity",
    "response to virus",
    "skeletal system development",
    "axon guidance",
    "action potential",
    "heart development",
    "heart morphogenesis",
    "homophilic cell-cell adhesion",
    "extracellular matrix organization",
    "extracellular structure organization",
    "external encapsulating structure organization",
    "collagen fibril organization"
  ),
  display_label = c(
    "complement activation",
    "positive regulation of leukocyte mediated cytotoxicity",
    "positive regulation of cell killing",
    "natural killer cell mediated immunity",
    "cell killing",
    "DNA replication",
    "DNA-templated DNA replication",
    "defense response to virus",
    "lymphocyte mediated immunity",
    "response to virus",
    "skeletal system development",
    "axon guidance",
    "action potential",
    "heart development",
    "heart morphogenesis",
    "homophilic cell-cell adhesion",
    "extracellular matrix organization",
    "extracellular structure organization",
    "external encapsulating structure organization",
    "collagen fibril organization"
  )
)

gsea_results[, Description := trimws(as.character(Description))]
gsea_results[, description_key := tolower(Description)]
gsea_results[, NES := suppressWarnings(as.numeric(NES))]
gsea_results[, p.adjust := suppressWarnings(as.numeric(p.adjust))]

locked_gsea_hits <- gsea_results[
  description_key %in% gsea_pathway_table$description_key
]
setorder(locked_gsea_hits, description_key, p.adjust)
locked_gsea_hits <- locked_gsea_hits[!duplicated(description_key)]

gsea_plot_data <- merge(
  gsea_pathway_table,
  locked_gsea_hits,
  by = "description_key",
  all.x = TRUE,
  sort = FALSE
)
gsea_plot_data[, display_order := match(
  description_key,
  gsea_pathway_table$description_key
)]
setorder(gsea_plot_data, display_order)

missing_locked_pathways <- gsea_plot_data[is.na(NES), display_label]
if (length(missing_locked_pathways) > 0L) {
  stop(
    "GSEA文件缺少以下锁定通路：\n",
    paste(missing_locked_pathways, collapse = "\n")
  )
}
if (nrow(gsea_plot_data) != 20L || any(!is.finite(gsea_plot_data$NES))) {
  stop("Figure S6b通路数量或NES核验失败。")
}

gsea_plot_data[, enrichment_direction := fifelse(
  NES > 0,
  "R-enriched",
  "NR-enriched"
)]
gsea_plot_data[, enrichment_direction := factor(
  enrichment_direction,
  levels = c("R-enriched", "NR-enriched")
)]
gsea_plot_data[, display_label := factor(
  display_label,
  levels = rev(gsea_pathway_table$display_label)
)]

panel_b <- ggplot(
  gsea_plot_data,
  aes(x = NES, y = display_label, fill = enrichment_direction)
) +
  geom_vline(xintercept = 0, linewidth = 0.55, color = "black") +
  geom_col(width = 0.72) +
  scale_fill_manual(
    values = c("R-enriched" = "#D95F59", "NR-enriched" = "#4C78A8"),
    breaks = c("R-enriched", "NR-enriched"),
    drop = FALSE,
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(-2.85, 2.85),
    breaks = c(-2, -1, 0, 1, 2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Normalized enrichment score", y = NULL) +
  theme_classic(base_family = "Times New Roman", base_size = 11) +
  theme(
    axis.title.x = element_text(size = 11.5, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9.5, color = "black"),
    axis.text.y = element_text(size = 9.3, color = "black"),
    axis.line = element_line(linewidth = 0.55, color = "black"),
    axis.ticks = element_line(linewidth = 0.50, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 9.5, color = "black"),
    plot.margin = margin(t = 10, r = 12, b = 8, l = 8)
  )

figure_s6 <- (panel_a | panel_b) +
  plot_layout(widths = c(0.95, 1.35)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(
        family = "Times New Roman",
        face = "bold",
        color = "black",
        size = 14
      ),
      plot.tag.position = c(0.01, 0.99)
    )
  )

ggsave(
  filename = output_png,
  plot = figure_s6,
  width = 13.4,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
ggsave(
  filename = output_tiff,
  plot = figure_s6,
  width = 13.4,
  height = 7.2,
  units = "in",
  dpi = 600,
  device = "tiff",
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)

if (!file.exists(output_png) || file.info(output_png)$size <= 10000) {
  stop("Figure S6 PNG生成失败。")
}
if (!file.exists(output_tiff) || file.info(output_tiff)$size <= 10000) {
  stop("Figure S6 TIFF生成失败。")
}

audit_lines <- c(
  "SUPPLEMENTARY FIGURE S6 REBUILD AUDIT",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0("Automatic y-axis display limit: ", y_display_cap),
  paste0("Y-truncated points: ", sum(volcano_data$y_was_truncated, na.rm = TRUE)),
  paste0("X-truncated points: ", sum(volcano_data$x_was_truncated, na.rm = TRUE)),
  "",
  "VOLCANO LABELS",
  capture.output(print(
    volcano_labels[, .(
      direction,
      gene_label,
      has_gene_symbol,
      log2FoldChange,
      padj
    )],
    row.names = FALSE
  )),
  "",
  "OUTPUT FILES",
  output_png,
  output_tiff
)
writeLines(audit_lines, audit_file, useBytes = TRUE)

cat("\nSupplementary Figure S6已重新生成。\n")
cat("自动确定的y轴上限：", y_display_cap, "\n", sep = "")
cat("PNG：", output_png, "\n", sep = "")
cat("TIFF：", output_tiff, "\n", sep = "")
cat("审计报告：", audit_file, "\n", sep = "")
print(figure_s6)
