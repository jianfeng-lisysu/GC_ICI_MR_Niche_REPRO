# 04_GSE167297_singlecell_expression.R
# Coarse single-cell expression validation for MR-prioritized candidates in GSE167297.

library(data.table)
library(ggplot2)

config_file <- file.path("config", "project_config.R")
if (!file.exists(config_file)) {
  stop("Run this script from the repository/project root, or copy config/project_config.R into the project root.")
}
source(config_file)
raw_dir <- file.path(project_dir, "01_raw_data/03_GSE167297_scRNA/raw_count_matrices")
sc_res_dir <- file.path(project_dir, "04_results/02_singlecell_results")
fig_main_dir <- file.path(project_dir, "05_figures/01_main_figures")
fig_supp_dir <- file.path(project_dir, "05_figures/02_supplementary_figures")
dir.create(sc_res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp_dir, recursive = TRUE, showWarnings = FALSE)

raw_files <- list.files(raw_dir, pattern = "CountMatrix.txt.gz$", full.names = TRUE)
if (length(raw_files) == 0) stop("No GSE167297 CountMatrix.txt.gz files found.")

candidate_genes <- c("FCN1", "SCARF2", "TNFSF12", "CD1D", "DLL1", "CD274", "SIGLEC7")

marker_sets <- list(
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  t_cell = c("CD3D", "CD3E", "TRAC"),
  cd8_t = c("CD8A", "CD8B", "GZMK"),
  exhausted_t = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4"),
  myeloid = c("LYZ", "LST1", "TYROBP", "AIF1", "C1QA", "C1QB", "FCGR3A", "CD68", "MS4A7", "S100A8", "S100A9"),
  dc = c("ITGAX", "LAMP3", "CCR7", "CLEC9A", "CD1C", "FCER1A", "HLA-DRA"),
  nk = c("NKG7", "GNLY", "PRF1", "GZMB", "KLRD1"),
  caf = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "TAGLN", "FAP"),
  endothelial = c("PECAM1", "VWF", "KDR", "ENG", "PLVAP"),
  b_plasma = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN")
)

all_genes <- unique(c(candidate_genes, unlist(marker_sets)))

annotate_one_sample <- function(f) {
  message("Reading: ", basename(f))
  dt <- fread(f)
  setnames(dt, 1, "gene")
  sample_name <- sub("_CountMatrix.txt.gz$", "", basename(f))
  region <- ifelse(grepl("Normal", sample_name), "Normal",
                   ifelse(grepl("Superficial", sample_name), "Superficial",
                          ifelse(grepl("Deep", sample_name), "Deep", "Unknown")))
  x <- dt[gene %in% all_genes]
  missing_genes <- setdiff(all_genes, x$gene)
  if (length(missing_genes) > 0) {
    zero_dt <- data.table(gene = missing_genes)
    for (cc in setdiff(names(dt), "gene")) zero_dt[[cc]] <- 0L
    x <- rbind(x, zero_dt, fill = TRUE)
  }
  mat <- as.data.frame(x)
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  mat <- as.matrix(mat)
  logmat <- log1p(mat)
  score_dt <- data.table(cell = colnames(logmat))
  for (ct in names(marker_sets)) {
    markers <- intersect(marker_sets[[ct]], rownames(logmat))
    score_dt[[ct]] <- colMeans(logmat[markers, , drop = FALSE])
  }
  score_mat <- as.matrix(score_dt[, names(marker_sets), with = FALSE])
  max_score <- apply(score_mat, 1, max)
  assigned <- names(marker_sets)[max.col(score_mat, ties.method = "first")]
  assigned[max_score == 0] <- "Unknown"
  anno <- data.table(cell = colnames(logmat), sample = sample_name, region = region, assigned_celltype = assigned)
  cand <- logmat[candidate_genes, , drop = FALSE]
  cand_long <- as.data.table(as.table(cand))
  setnames(cand_long, c("gene", "cell", "log_expr"))
  cand_long[, count_positive := log_expr > 0]
  merge(cand_long, anno, by = "cell")
}

cell_expr <- rbindlist(lapply(raw_files, annotate_one_sample), fill = TRUE)

candidate_celltype_summary <- cell_expr[, .(
  n_cells = .N,
  n_expr = sum(count_positive),
  pct_expr = 100 * mean(count_positive),
  mean_log_expr_all = mean(log_expr),
  mean_log_expr_positive = ifelse(sum(count_positive) > 0, mean(log_expr[count_positive]), NA_real_)
), by = .(gene, region, assigned_celltype)][order(gene, region, -pct_expr)]

fwrite(candidate_celltype_summary, file.path(sc_res_dir, "GSE167297_candidate_expr_by_coarse_celltype.tsv"), sep = "\t")

plot_dt <- candidate_celltype_summary[
  assigned_celltype %in% c("epithelial", "myeloid", "dc", "t_cell", "cd8_t", "exhausted_t", "nk", "b_plasma", "caf", "endothelial")
]

celltype_labels <- c(
  epithelial = "Epithelial",
  myeloid = "Myeloid",
  dc = "DC",
  t_cell = "T cell",
  cd8_t = "CD8 T",
  exhausted_t = "Exhausted T",
  nk = "NK",
  b_plasma = "B/plasma",
  caf = "CAF",
  endothelial = "Endothelial"
)

plot_dt[, assigned_celltype := factor(assigned_celltype, levels = names(celltype_labels))]
plot_dt[, celltype_label := factor(celltype_labels[as.character(assigned_celltype)], levels = celltype_labels)]
plot_dt[, region := factor(region, levels = c("Normal", "Superficial", "Deep"))]

make_expr_plot <- function(dt, genes, title_text) {
  pdat <- copy(dt[gene %in% genes])
  pdat[, gene := factor(gene, levels = genes)]
  ggplot(pdat, aes(x = celltype_label, y = pct_expr, fill = region)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.68) +
    facet_wrap(~ gene, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(Normal = "#7A7A7A", Superficial = "#2F80ED", Deep = "#D64541")) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "#F0F0F0", color = NA),
      strip.text = element_text(face = "bold", size = 13),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.spacing = unit(1.1, "lines"),
      plot.title = element_text(face = "bold", size = 15, hjust = 0)
    ) +
    labs(x = "Coarse cell type", y = "% expressing cells", fill = "Region", title = title_text)
}

p_main <- make_expr_plot(
  plot_dt,
  c("FCN1", "TNFSF12", "SCARF2", "DLL1"),
  "Core MR-prioritized targets across coarse cell types in GSE167297"
)
p_supp <- make_expr_plot(
  plot_dt,
  c("CD274", "SIGLEC7", "CD1D"),
  "Additional immune-checkpoint and antigen-presentation candidates in GSE167297"
)

ggsave(file.path(fig_main_dir, "Figure3_core_targets_expression_by_celltype_region.pdf"), p_main, width = 10, height = 7)
ggsave(file.path(fig_main_dir, "Figure3_core_targets_expression_by_celltype_region.png"), p_main, width = 10, height = 7, dpi = 600)
ggsave(file.path(fig_supp_dir, "FigureS_candidate_secondary_targets_expression_by_celltype_region.pdf"), p_supp, width = 10, height = 5.5)
ggsave(file.path(fig_supp_dir, "FigureS_candidate_secondary_targets_expression_by_celltype_region.png"), p_supp, width = 10, height = 5.5, dpi = 600)

print(candidate_celltype_summary[gene %in% c("FCN1", "TNFSF12", "SCARF2", "DLL1")][order(gene, -pct_expr)][1:40])
