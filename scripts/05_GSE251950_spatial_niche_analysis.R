# 05_GSE251950_spatial_niche_analysis.R
# Spatial validation of MR-prioritized targets in GSE251950.

library(data.table)
library(Matrix)
library(ggplot2)

PROJECT <- path.expand("~/Desktop/Scientific Research/Bioinformatics analysis/GC_ICI_MR_Niche_REPRO_20260617")
sp_dir <- file.path(PROJECT, "01_raw_data/04_GSE251950_spatial/extracted")
res_dir <- file.path(PROJECT, "04_results/03_spatial_results")
fig_dir <- file.path(PROJECT, "05_figures/03_spatial_figures")

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# NOTE:
# This script assumes GSE251950_RAW.tar has been downloaded and all GSM sample tar.gz
# files have been extracted into 10 Visium sample folders under sp_dir.

sample_dirs <- list.dirs(sp_dir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[file.exists(file.path(sample_dirs, paste0(basename(sample_dirs), "_matrix.mtx.gz")))]

target_genes <- c("FCN1", "TNFSF12", "SCARF2", "DLL1", "CD274", "SIGLEC7", "CD1D")

marker_sets <- list(
  myeloid = c("LYZ", "LST1", "TYROBP", "AIF1", "C1QA", "C1QB", "FCGR3A", "CD68", "MS4A7", "S100A8", "S100A9"),
  caf = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "TAGLN", "FAP"),
  endothelial = c("PECAM1", "VWF", "KDR", "ENG", "PLVAP"),
  cd8_t = c("CD8A", "CD8B", "GZMK", "GZMB"),
  exhausted_t = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4"),
  dc_lamp3 = c("LAMP3", "CCR7", "ITGAX", "HLA-DRA", "CD1C", "FCER1A", "CLEC9A"),
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1")
)

all_needed_genes <- unique(c(target_genes, unlist(marker_sets)))

read_one_visium_light <- function(sdir) {
  sample_id <- basename(sdir)
  message("Reading sample: ", sample_id)

  mat <- readMM(file.path(sdir, paste0(sample_id, "_matrix.mtx.gz")))
  barcodes <- fread(file.path(sdir, paste0(sample_id, "_barcodes.tsv.gz")), header = FALSE)$V1
  features <- fread(file.path(sdir, paste0(sample_id, "_features.tsv.gz")), header = FALSE)
  gene_symbol <- if (ncol(features) >= 2) features[[2]] else features[[1]]

  colnames(mat) <- barcodes
  lib_size <- Matrix::colSums(mat)
  lib_size[lib_size == 0] <- NA

  keep_idx <- which(gene_symbol %in% all_needed_genes)
  sub_mat <- mat[keep_idx, , drop = FALSE]
  sub_gene <- gene_symbol[keep_idx]

  expr_mat <- do.call(rbind, lapply(all_needed_genes, function(g) {
    idx <- which(sub_gene == g)
    if (length(idx) == 0) rep(0, ncol(mat)) else Matrix::colSums(sub_mat[idx, , drop = FALSE])
  }))

  rownames(expr_mat) <- all_needed_genes
  colnames(expr_mat) <- barcodes
  norm_mat <- log1p(t(t(expr_mat) / lib_size * 10000))

  spot_dt <- data.table(sample = sample_id, barcode = barcodes)

  for (g in target_genes) spot_dt[[g]] <- as.numeric(norm_mat[g, ])

  for (nm in names(marker_sets)) {
    genes <- intersect(marker_sets[[nm]], rownames(norm_mat))
    spot_dt[[paste0(nm, "_score")]] <- colMeans(norm_mat[genes, , drop = FALSE])
  }

  pos_file <- file.path(sdir, paste0(sample_id, "_tissue_positions_list.csv"))
  if (file.exists(pos_file)) {
    pos <- fread(pos_file, header = FALSE)
    setnames(pos, c("barcode", "in_tissue", "array_row", "array_col", "pxl_col", "pxl_row"))
    spot_dt <- merge(spot_dt, pos, by = "barcode", all.x = TRUE)
    spot_dt <- spot_dt[in_tissue == 1 | is.na(in_tissue)]
  }

  spot_dt[]
}

spatial_spots <- rbindlist(lapply(sample_dirs, read_one_visium_light), fill = TRUE)
fwrite(spatial_spots, file.path(res_dir, "GSE251950_spatial_spot_scores_targets.tsv"), sep = "\t")

target_for_cor <- target_genes
score_cols <- c("myeloid_score", "caf_score", "endothelial_score", "cd8_t_score", "exhausted_t_score", "dc_lamp3_score", "epithelial_score")

cor_one <- function(dt, target, score) {
  x <- dt[[target]]
  y <- dt[[score]]
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  data.table(target = target, score = score, rho = unname(ct$estimate), pval = ct$p.value, n_spots = sum(ok))
}

cor_res <- rbindlist(lapply(target_for_cor, function(tg) {
  rbindlist(lapply(score_cols, function(sc) cor_one(spatial_spots, tg, sc)))
}))

cor_res[, FDR := p.adjust(pval, method = "BH")]
fwrite(cor_res, file.path(res_dir, "GSE251950_target_niche_score_spearman.tsv"), sep = "\t")

plot_cor <- copy(cor_res)
plot_cor[, target := factor(target, levels = rev(c("FCN1", "TNFSF12", "SCARF2", "DLL1", "CD274", "SIGLEC7", "CD1D")))]
plot_cor[, score := factor(score,
  levels = c("myeloid_score", "caf_score", "endothelial_score", "dc_lamp3_score", "cd8_t_score", "exhausted_t_score", "epithelial_score"),
  labels = c("Myeloid", "CAF", "Endothelial", "DC/LAMP3", "CD8 T", "Exhausted T", "Epithelial")
)]

p <- ggplot(plot_cor, aes(x = score, y = target, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 4, color = "black") +
  scale_fill_gradient2(low = "#2F80ED", mid = "white", high = "#D64541", midpoint = 0, limits = c(-0.15, 0.35)) +
  coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black", face = "bold"),
        axis.title = element_blank(),
        legend.title = element_text(face = "bold"),
        plot.title = element_blank())

ggsave(file.path(fig_dir, "Figure4_spatial_niche_correlation_600dpi.tiff"), p, width = 7.5, height = 5.2, dpi = 600, compression = "lzw")
