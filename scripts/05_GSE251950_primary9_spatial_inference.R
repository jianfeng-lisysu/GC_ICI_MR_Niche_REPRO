# ============================================================
# GSE251950 primary-nine spatial inference
#
# Reproduces the patient-level analysis described in the final
# manuscript. The paired metastatic section GC6-PM is excluded
# from primary inference. Spot-pooled results are descriptive.
#
# Primary tests:
#   FCN1 or TNFSF12 versus CAF, myeloid and endothelial scores
#   within each of nine independent primary-tumour sections.
#   The nine section-level Spearman coefficients are tested
#   against zero with a two-sided exact Wilcoxon signed-rank test.
#   BH correction is applied across the six candidate-program pairs.
#
# Run from the project root.
# ============================================================

required_packages <- c("data.table", "Matrix", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg)
  }
}

library(data.table)
library(Matrix)
library(ggplot2)

config_file <- file.path("config", "project_config.R")
if (!file.exists(config_file)) {
  stop("Missing config/project_config.R. Run from the project root.")
}
source(config_file)

spatial_dir <- file.path(
  project_dir,
  "01_raw_data",
  "04_GSE251950_spatial",
  "extracted"
)
result_dir <- file.path(project_dir, "04_results", "GSE251950_spatial")
figure_dir <- file.path(project_dir, "05_figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

primary_sections <- paste0("GC", 1:9)
excluded_section <- "GC6-PM"
target_genes <- c("FCN1", "TNFSF12")

marker_sets <- list(
  CAF = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "TAGLN", "FAP"),
  Myeloid = c("LYZ", "LST1", "TYROBP", "AIF1", "C1QA", "C1QB", "FCGR3A", "CD68", "MS4A7", "S100A8", "S100A9"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ENG", "PLVAP")
)

all_needed_genes <- unique(c(target_genes, unlist(marker_sets, use.names = FALSE)))

find_sample_label <- function(path) {
  name <- basename(path)
  hits <- c(primary_sections, excluded_section)
  matched <- hits[vapply(hits, function(x) grepl(x, name, fixed = TRUE), logical(1))]
  if (length(matched) != 1L) {
    stop("Could not assign a unique GC section label to directory: ", path)
  }
  matched
}

find_one_file <- function(directory, patterns, label) {
  files <- list.files(directory, full.names = TRUE, recursive = FALSE)
  keep <- rep(FALSE, length(files))
  for (pattern in patterns) {
    keep <- keep | grepl(pattern, basename(files), ignore.case = TRUE)
  }
  matches <- files[keep]
  if (length(matches) != 1L) {
    stop("Expected one ", label, " file in ", directory, "; found ", length(matches))
  }
  matches
}

read_one_section <- function(section_dir) {
  section <- find_sample_label(section_dir)
  message("Reading ", section, " from ", section_dir)

  matrix_file <- find_one_file(section_dir, c("matrix\\.mtx", "matrix\\.mtx\\.gz"), "matrix")
  barcode_file <- find_one_file(section_dir, c("barcodes\\.tsv", "barcodes\\.tsv\\.gz"), "barcode")
  feature_file <- find_one_file(section_dir, c("features\\.tsv", "features\\.tsv\\.gz"), "feature")

  counts <- Matrix::readMM(matrix_file)
  barcodes <- data.table::fread(barcode_file, header = FALSE)[[1]]
  features <- data.table::fread(feature_file, header = FALSE)
  gene_symbols <- if (ncol(features) >= 2L) features[[2]] else features[[1]]

  if (ncol(counts) != length(barcodes)) {
    stop("Barcode count does not match matrix columns for ", section)
  }

  colnames(counts) <- barcodes
  library_size <- Matrix::colSums(counts)
  library_size[library_size <= 0] <- NA_real_

  keep <- which(gene_symbols %in% all_needed_genes)
  sub_counts <- counts[keep, , drop = FALSE]
  sub_symbols <- gene_symbols[keep]

  collapsed <- do.call(
    rbind,
    lapply(all_needed_genes, function(gene) {
      indices <- which(sub_symbols == gene)
      if (length(indices) == 0L) {
        rep(0, ncol(counts))
      } else {
        Matrix::colSums(sub_counts[indices, , drop = FALSE])
      }
    })
  )
  rownames(collapsed) <- all_needed_genes
  colnames(collapsed) <- barcodes

  normalized <- log1p(t(t(collapsed) / library_size * 10000))
  section_data <- data.table(section = section, barcode = barcodes)

  for (gene in target_genes) {
    section_data[[gene]] <- as.numeric(normalized[gene, ])
  }

  for (program in names(marker_sets)) {
    present <- intersect(marker_sets[[program]], rownames(normalized))
    if (length(present) == 0L) {
      stop("No genes from the ", program, " marker set were present in ", section)
    }
    section_data[[paste0(program, "_score")]] <- colMeans(
      normalized[present, , drop = FALSE]
    )
    section_data[[paste0(program, "_genes_used")]] <- length(present)
  }

  position_files <- list.files(
    section_dir,
    pattern = "tissue_positions.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(position_files) == 1L) {
    positions <- fread(position_files)
    if (ncol(positions) >= 6L) {
      if (!"barcode" %in% names(positions)) {
        setnames(
          positions,
          old = names(positions)[1:6],
          new = c("barcode", "in_tissue", "array_row", "array_col", "pxl_col", "pxl_row")
        )
      }
      section_data <- merge(section_data, positions, by = "barcode", all.x = TRUE)
      if ("in_tissue" %in% names(section_data)) {
        section_data <- section_data[is.na(in_tissue) | in_tissue == 1]
      }
    }
  }

  section_data[]
}

sample_dirs <- list.dirs(spatial_dir, recursive = FALSE, full.names = TRUE)
if (length(sample_dirs) == 0L) {
  stop("No extracted GSE251950 section directories found in: ", spatial_dir)
}

labels <- vapply(sample_dirs, find_sample_label, character(1))
if (anyDuplicated(labels)) {
  stop("Duplicate section labels detected: ", paste(labels[duplicated(labels)], collapse = ", "))
}

all_spots <- rbindlist(lapply(sample_dirs, read_one_section), fill = TRUE)
fwrite(all_spots, file.path(result_dir, "GSE251950_spot_scores.tsv.gz"), sep = "\t")

primary_spots <- all_spots[section %in% primary_sections]
if (!setequal(unique(primary_spots$section), primary_sections)) {
  missing <- setdiff(primary_sections, unique(primary_spots$section))
  stop("Missing primary sections: ", paste(missing, collapse = ", "))
}

correlate_one <- function(data, target, program) {
  score_name <- paste0(program, "_score")
  valid <- is.finite(data[[target]]) & is.finite(data[[score_name]])
  if (sum(valid) < 10L) {
    return(data.table(rho = NA_real_, p_value = NA_real_, n_spots = sum(valid)))
  }
  test <- suppressWarnings(cor.test(
    data[[target]][valid],
    data[[score_name]][valid],
    method = "spearman",
    exact = FALSE
  ))
  data.table(
    rho = unname(test$estimate),
    p_value = test$p.value,
    n_spots = sum(valid)
  )
}

section_results <- rbindlist(lapply(primary_sections, function(section_name) {
  section_data <- primary_spots[section == section_name]
  rbindlist(lapply(target_genes, function(target) {
    rbindlist(lapply(names(marker_sets), function(program) {
      result <- correlate_one(section_data, target, program)
      result[, `:=`(section = section_name, target = target, program = program)]
      result
    }))
  }))
}))

setcolorder(section_results, c("section", "target", "program", "rho", "p_value", "n_spots"))
section_results[, section := factor(section, levels = primary_sections)]
section_results[, target := factor(target, levels = target_genes)]
section_results[, program := factor(program, levels = names(marker_sets))]
setorder(section_results, target, program, section)

summary_results <- section_results[, {
  values <- rho[is.finite(rho)]
  signed_rank <- suppressWarnings(wilcox.test(
    values,
    mu = 0,
    alternative = "two.sided",
    exact = TRUE,
    paired = FALSE
  ))
  .(
    n_sections = length(values),
    positive_sections = sum(values > 0),
    median_rho = median(values),
    q1_rho = unname(quantile(values, probs = 0.25, names = FALSE)),
    q3_rho = unname(quantile(values, probs = 0.75, names = FALSE)),
    exact_signed_rank_p = signed_rank$p.value
  )
}, by = .(target, program)]
summary_results[, BH_adjusted_p := p.adjust(exact_signed_rank_p, method = "BH")]

fwrite(section_results, file.path(result_dir, "GSE251950_primary9_section_correlations.tsv"), sep = "\t")
fwrite(summary_results, file.path(result_dir, "GSE251950_primary9_summary.tsv"), sep = "\t")

# Descriptive pooled-spot correlations; not used for primary inference.
pooled_results <- rbindlist(lapply(target_genes, function(target) {
  rbindlist(lapply(names(marker_sets), function(program) {
    result <- correlate_one(primary_spots, target, program)
    result[, `:=`(target = target, program = program)]
    result
  }))
}))
pooled_results[, BH_adjusted_p := p.adjust(p_value, method = "BH")]
fwrite(pooled_results, file.path(result_dir, "GSE251950_primary9_pooled_spot_correlations_descriptive.tsv"), sep = "\t")

# Patient-level heatmap (Supplementary Figure S9 concept).
heatmap_data <- copy(section_results)
heatmap_data[, section := factor(section, levels = rev(primary_sections))]
p_heatmap <- ggplot(heatmap_data, aes(x = program, y = section, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 3.2) +
  facet_grid(. ~ target) +
  scale_fill_gradient2(low = "white", mid = "white", high = "firebrick", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Spearman rho") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold")
  )

ggsave(
  file.path(figure_dir, "FigS9_GSE251950_primary9_patient_level_heatmap.tiff"),
  p_heatmap,
  width = 7.2,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)

cat("Primary-nine spatial analysis completed.\n")
print(summary_results)
