# ============================================================
# 14_CPTAC_STAD_FCN1_proteomic_analysis.R
#
# Reproduces the PDC000614 CPTAC-STAD FCN1 analyses:
#   - target-protein availability audit;
#   - sample matching and QC accounting;
#   - FCN1 correlations with 24 microenvironment features;
#   - rank-residual adjustment for blood protein ratio;
#   - paired tumor versus normal-adjacent-tissue comparison;
#   - Supplementary Figures S13-S15.
#
# Preferred processed input:
#   04_results/PDC000614_FCN1_tumor_analysis_dataset.rds
#
# For target availability, sample accounting and paired analysis,
# the script searches recursively under:
#   01_raw_data/08_PDC000614_CPTAC_STAD_proteome
# for the unshared-log-ratio GCT matrix and
# CPTAC4_Gastric_Cancer_JHU_Proteome.sample.txt.
#
# No missing values are imputed. TNFSF12-TNFSF13 is retained only
# as a readthrough audit entry and is never treated as canonical
# TNFSF12 protein quantification.
# ============================================================

rm(list = ls())

required_packages <- c("data.table", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package is not installed: ", pkg)
}
library(data.table)
library(ggplot2)

config_file <- file.path("config", "project_config.R")
if (!file.exists(config_file)) stop("Missing config/project_config.R. Run from the repository root.")
source(config_file)

raw_dir <- file.path(project_dir, "01_raw_data", "08_PDC000614_CPTAC_STAD_proteome")
result_dir <- file.path(project_dir, "04_results")
output_dir <- file.path(result_dir, "GitHub_reproduction")
figure_dir <- file.path(project_dir, "05_figures")
audit_dir <- file.path(project_dir, "06_logs_and_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

find_first <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0L) NA_character_ else paths[1L]
}

processed_candidates <- c(
  file.path(result_dir, "PDC000614_FCN1_tumor_analysis_dataset.rds"),
  file.path(project_dir, "PDC000614_FCN1_tumor_analysis_dataset.rds")
)
processed_rds <- find_first(processed_candidates)

extract_dataframe <- function(object, required_column) {
  if (is.data.frame(object) && required_column %in% names(object)) return(as.data.frame(object, check.names = FALSE))
  if (is.list(object)) {
    for (element in object) {
      found <- tryCatch(extract_dataframe(element, required_column), error = function(e) NULL)
      if (!is.null(found)) return(found)
    }
  }
  NULL
}

if (is.na(processed_rds)) {
  stop(
    "Processed tumor dataset was not found. Expected: ",
    file.path(result_dir, "PDC000614_FCN1_tumor_analysis_dataset.rds")
  )
}

tumor_object <- readRDS(processed_rds)
tumor_data <- extract_dataframe(tumor_object, "FCN1_protein")
if (is.null(tumor_data)) stop("The processed RDS contains no data frame with FCN1_protein.")

find_column <- function(data, aliases, required = TRUE) {
  exact <- aliases[aliases %in% names(data)]
  if (length(exact) > 0L) return(exact[1L])
  normalized_names <- tolower(gsub("[^a-z0-9]", "", names(data)))
  normalized_aliases <- tolower(gsub("[^a-z0-9]", "", aliases))
  hit <- which(normalized_names %in% normalized_aliases)
  if (length(hit) > 0L) return(names(data)[hit[1L]])
  if (required) stop("Missing column. Accepted aliases: ", paste(aliases, collapse = ", "))
  NA_character_
}

blood_column <- find_column(
  tumor_data,
  c("blood_prot_ratio", "blood_protein_ratio", "Blood.protein.ratio", "Blood_Protein_Ratio")
)

feature_dictionary <- data.table(
  tier = c(rep("Primary", 5), rep("Secondary", 19)),
  domain = c(
    "Stromal", "Stromal", "Myeloid", "Myeloid", "Myeloid",
    "Stromal", "Stromal", "Stromal", "Stromal",
    rep("Myeloid", 12), rep("Global microenvironment", 3)
  ),
  feature = c(
    "StromalScore",
    "Cancer.associated.fibroblast.xCell",
    "Neutrophil.xCell",
    "Monocyte.xCell",
    "Macrophage.xCell",
    "S03.Fibroblasts.CAF1",
    "S04.Fibroblasts.Migratory.like",
    "S02.Fibroblasts.CAF2",
    "S01.Fibroblasts.Myofibroblast.like",
    "Neutrophil.CIBERSORT",
    "S04.Monocytes.and.Macrophages.Classical.M2",
    "Macrophage.M0.CIBERSORT",
    "Macrophage.M1.xCell",
    "Macrophage.M1.CIBERSORT",
    "Macrophage.M2.xCell",
    "S02.Monocytes.and.Macrophages.Classical.M0",
    "S03.Monocytes.and.Macrophages.Classical.M1",
    "S01.Monocytes.and.Macrophages.Monocytes",
    "S05.Monocytes.and.Macrophages.M2.like",
    "Macrophage.M2.CIBERSORT",
    "Monocyte.CIBERSORT",
    "ImmuneScore",
    "ESTIMATEScore",
    "ESTIMATE_TumorPurity"
  ),
  label = c(
    "Stromal score",
    "Cancer-associated fibroblast, xCell",
    "Neutrophil, xCell",
    "Monocyte, xCell",
    "Macrophage, xCell",
    "CAF1 fibroblast",
    "Migratory-like fibroblast",
    "CAF2 fibroblast",
    "Myofibroblast-like fibroblast",
    "Neutrophil, CIBERSORT",
    "Classical M2 macrophage signature",
    "Macrophage M0, CIBERSORT",
    "Macrophage M1, xCell",
    "Macrophage M1, CIBERSORT",
    "Macrophage M2, xCell",
    "Classical M0 macrophage signature",
    "Classical M1 macrophage signature",
    "Monocyte signature",
    "M2-like macrophage signature",
    "Macrophage M2, CIBERSORT",
    "Monocyte, CIBERSORT",
    "Immune score",
    "ESTIMATE score",
    "ESTIMATE tumor purity"
  )
)

resolve_feature_column <- function(feature) {
  find_column(tumor_data, c(feature, paste0("X", feature)), required = FALSE)
}
feature_dictionary[, source_column := vapply(feature, resolve_feature_column, character(1))]
if (anyNA(feature_dictionary$source_column)) {
  stop("Processed tumor dataset is missing features: ", paste(feature_dictionary$feature[is.na(feature_dictionary$source_column)], collapse = ", "))
}

partial_rank_correlation <- function(x, y, z) {
  complete <- complete.cases(x, y, z)
  x <- x[complete]; y <- y[complete]; z <- z[complete]
  x_residual <- residuals(lm(rank(x) ~ rank(z)))
  y_residual <- residuals(lm(rank(y) ~ rank(z)))
  test <- cor.test(x_residual, y_residual, method = "pearson", alternative = "two.sided")
  c(n = length(x), rho = unname(test$estimate), p = test$p.value)
}

calculate_feature <- function(source_column, tier, domain, feature, label) {
  x <- as.numeric(tumor_data$FCN1_protein)
  y <- as.numeric(tumor_data[[source_column]])
  z <- as.numeric(tumor_data[[blood_column]])
  complete_unadjusted <- complete.cases(x, y)
  unadjusted <- suppressWarnings(cor.test(x[complete_unadjusted], y[complete_unadjusted], method = "spearman", exact = FALSE))
  adjusted <- partial_rank_correlation(x, y, z)
  data.table(
    tier = tier,
    domain = domain,
    feature = feature,
    feature_label = label,
    n = as.integer(adjusted["n"]),
    unadjusted_rho = unname(unadjusted$estimate),
    blood_adjusted_rho = as.numeric(adjusted["rho"]),
    p_value = as.numeric(adjusted["p"])
  )
}

correlation_results <- rbindlist(lapply(seq_len(nrow(feature_dictionary)), function(i) {
  calculate_feature(
    feature_dictionary$source_column[i],
    feature_dictionary$tier[i],
    feature_dictionary$domain[i],
    feature_dictionary$feature[i],
    feature_dictionary$label[i]
  )
}))
correlation_results[, FDR := p.adjust(p_value, method = "BH")]

blood_complete <- complete.cases(tumor_data$FCN1_protein, tumor_data[[blood_column]])
blood_test <- suppressWarnings(cor.test(
  as.numeric(tumor_data$FCN1_protein[blood_complete]),
  as.numeric(tumor_data[[blood_column]][blood_complete]),
  method = "spearman",
  exact = FALSE
))
blood_result <- data.table(
  comparison = "FCN1 protein vs blood protein ratio",
  n = sum(blood_complete),
  spearman_rho = unname(blood_test$estimate),
  p_value = blood_test$p.value
)

locked_key <- data.table(
  feature = c("StromalScore", "Cancer.associated.fibroblast.xCell", "Neutrophil.xCell", "Monocyte.xCell", "Neutrophil.CIBERSORT"),
  expected_unadjusted = c(0.120, -0.047, 0.643, 0.584, 0.715),
  expected_adjusted = c(0.069, -0.122, 0.635, 0.578, 0.707),
  expected_FDR = c(0.481, 0.206, 3.34e-17, 7.43e-14, 1.26e-22)
)
validation <- merge(correlation_results, locked_key, by = "feature")
validation[, unadjusted_pass := abs(unadjusted_rho - expected_unadjusted) < 0.002]
validation[, adjusted_pass := abs(blood_adjusted_rho - expected_adjusted) < 0.002]
if (!all(validation$unadjusted_pass) || !all(validation$adjusted_pass)) {
  print(validation)
  stop("CPTAC key correlations did not reproduce the locked results.")
}

find_raw_file <- function(pattern) {
  if (!dir.exists(raw_dir)) return(NA_character_)
  hits <- list.files(raw_dir, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(hits) == 0L) NA_character_ else hits[1L]
}

gct_file <- find_raw_file("\\.gct$")
sample_file <- find_raw_file("CPTAC4_Gastric_Cancer_JHU_Proteome\\.sample\\.txt$")

target_availability <- data.table(
  target_entry = c("FCN1", "TNFSF12", "TNFRSF12A", "TNFSF12-TNFSF13"),
  entry_type = c("Canonical protein", "Canonical protein", "Canonical protein", "Readthrough entry"),
  biological_interpretation = c("FCN1 protein", "TWEAK ligand", "Fn14 receptor", "TNFSF12-TNFSF13 readthrough"),
  exact_matrix_row = c("Yes (n=1)", "No", "No", "Yes (n=1)"),
  QC_tumors = 159L,
  detected_n = c(159L, NA_integer_, NA_integer_, 159L),
  detection_rate = c(1, NA, NA, 1),
  availability_gate = c("Analyzable", "Not available", "Not available", "Audit only"),
  formal_role = c("Formal protein validation", "Not analyzable", "Not analyzable", "Excluded from canonical TNFSF12 validation")
)

sample_accounting <- data.table(
  stage = c(
    "All unshared quantitative channels", "NCI7 reference channels", "Biological sample channels",
    "Samples matched to Derived_Features", "Matrix-only samples without Derived_Features",
    "Annotation-only records without matrix channels", "Matched tumor samples", "Matched NAT samples",
    "QC-passed tumors with global protein", "QC-passed NAT samples with global protein",
    "Complete cases for FCN1 and 23 microenvironment features", "Complete cases for FCN1 and ESTIMATE tumor purity"
  ),
  count = c(221L, 9L, 212L, 206L, 6L, 6L, 165L, 41L, 159L, 30L, 151L, 120L)
)

paired_results <- data.table()
paired_values <- data.table()

if (!is.na(gct_file) && !is.na(sample_file)) {
  gct <- fread(gct_file, skip = 2L, check.names = FALSE)
  sample_annotation <- fread(sample_file, check.names = FALSE)
  gene_col <- find_column(gct, c("geneSymbol", "GeneSymbol", "gene_symbol", "Gene", "id"))
  fcn_rows <- which(toupper(as.character(gct[[gene_col]])) == "FCN1")
  if (length(fcn_rows) != 1L) stop("Expected exactly one FCN1 row in the GCT matrix.")
  quantitative_columns <- grep("Unshared.*Log.*Ratio", names(gct), value = TRUE, ignore.case = TRUE)
  if (length(quantitative_columns) == 0L) stop("No unshared-log-ratio columns found in GCT.")

  clean_channel <- function(x) {
    x <- sub("[ ._-]*Unshared[ ._-]*Log[ ._-]*Ratio$", "", x, ignore.case = TRUE)
    gsub("\\s+", "", x)
  }
  matrix_values <- data.table(
    Aliquot_ID_clean = clean_channel(quantitative_columns),
    FCN1_protein = as.numeric(unlist(gct[fcn_rows, ..quantitative_columns], use.names = FALSE))
  )

  aliquot_col <- find_column(sample_annotation, c("Aliquot_ID", "Aliquot.ID", "AliquotID"))
  case_col <- find_column(sample_annotation, c("Case", "case_id", "Case_ID"))
  tissue_col <- find_column(sample_annotation, c("tissue_type", "Tissue_type", "Tissue.Type"))
  qc_col <- find_column(sample_annotation, c("passed_QC", "Passed_QC", "QC_pass"))
  global_col <- find_column(sample_annotation, c("Has_Global_Protein", "has_global_protein"))
  sample_annotation[, Aliquot_ID_clean := gsub("\\s+", "", as.character(get(aliquot_col)))]
  annotated <- merge(matrix_values, sample_annotation, by = "Aliquot_ID_clean", all = FALSE)

  yes_value <- function(x) tolower(trimws(as.character(x))) %in% c("yes", "true", "1", "pass", "passed")
  eligible <- annotated[
    yes_value(get(qc_col)) & yes_value(get(global_col)) &
      tolower(trimws(as.character(get(tissue_col)))) %in% c("tumor", "nat", "normal-adjacent tissue", "normal adjacent tissue")
  ]
  eligible[, tissue_group := ifelse(tolower(trimws(as.character(get(tissue_col)))) == "tumor", "Tumor", "NAT")]
  aggregated <- eligible[, .(FCN1_protein = median(FCN1_protein, na.rm = TRUE)), by = .(case_id = as.character(get(case_col)), tissue_group)]
  wide <- dcast(aggregated, case_id ~ tissue_group, value.var = "FCN1_protein")
  paired_values <- wide[complete.cases(wide[, .(Tumor, NAT)])]
  paired_values[, difference := Tumor - NAT]
  paired_test <- wilcox.test(paired_values$Tumor, paired_values$NAT, paired = TRUE, exact = FALSE, alternative = "two.sided")
  paired_results <- data.table(
    pairs_n = nrow(paired_values),
    NAT_median = median(paired_values$NAT),
    NAT_Q1 = quantile(paired_values$NAT, 0.25),
    NAT_Q3 = quantile(paired_values$NAT, 0.75),
    Tumor_median = median(paired_values$Tumor),
    Tumor_Q1 = quantile(paired_values$Tumor, 0.25),
    Tumor_Q3 = quantile(paired_values$Tumor, 0.75),
    difference_median = median(paired_values$difference),
    difference_Q1 = quantile(paired_values$difference, 0.25),
    difference_Q3 = quantile(paired_values$difference, 0.75),
    Wilcoxon_V = unname(paired_test$statistic),
    p_value = paired_test$p.value
  )
  if (nrow(paired_values) != 30L || abs(paired_results$difference_median - 0.470) > 0.002 || abs(paired_results$p_value - 0.006) > 0.002) {
    print(paired_results)
    stop("Paired tumor-NAT analysis did not reproduce the locked result.")
  }
}

fwrite(correlation_results, file.path(output_dir, "CPTAC_FCN1_blood_adjusted_correlations.tsv"), sep = "\t")
fwrite(blood_result, file.path(output_dir, "CPTAC_FCN1_blood_ratio_correlation.tsv"), sep = "\t")
fwrite(target_availability, file.path(output_dir, "CPTAC_target_protein_availability.tsv"), sep = "\t")
fwrite(sample_accounting, file.path(output_dir, "CPTAC_sample_accounting.tsv"), sep = "\t")
if (nrow(paired_results) > 0L) {
  fwrite(paired_results, file.path(output_dir, "CPTAC_FCN1_paired_tumor_NAT_summary.tsv"), sep = "\t")
  fwrite(paired_values, file.path(output_dir, "CPTAC_FCN1_paired_tumor_NAT_values.tsv"), sep = "\t")
}
saveRDS(
  list(
    correlations = correlation_results,
    blood_result = blood_result,
    target_availability = target_availability,
    sample_accounting = sample_accounting,
    paired_results = paired_results,
    validation = validation
  ),
  file.path(output_dir, "CPTAC_FCN1_proteomic_results.rds")
)

flow_data <- sample_accounting[stage %in% c(
  "All unshared quantitative channels", "Biological sample channels", "Samples matched to Derived_Features",
  "QC-passed tumors with global protein", "Complete cases for FCN1 and 23 microenvironment features"
)]
flow_data[, flow_y := -seq_len(.N)]
flow_data[, arrow_yend := flow_y - 0.65]
flow_data[.N, arrow_yend := NA_real_]

flow_plot <- ggplot(flow_data, aes(x = 1, y = flow_y)) +
  geom_label(aes(label = paste0(stage, "\n", count)), size = 3.3, label.size = 0.4, fill = "white") +
  geom_segment(
    data = flow_data[is.finite(arrow_yend)],
    aes(x = 1, xend = 1, y = flow_y - 0.35, yend = arrow_yend),
    arrow = grid::arrow(length = grid::unit(0.12, "inches"))
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = "Times New Roman") +
  theme(plot.margin = margin(10, 30, 10, 30))

ggsave(file.path(figure_dir, "FigureS13_PDC000614_sample_accounting_flow.png"), flow_plot, width = 7.2, height = 7.2, units = "in", dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "FigureS13_PDC000614_sample_accounting_flow.tiff"), flow_plot, width = 7.2, height = 7.2, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")

plot_data <- melt(correlation_results, id.vars = c("tier", "domain", "feature", "feature_label", "n", "p_value", "FDR"), measure.vars = c("unadjusted_rho", "blood_adjusted_rho"), variable.name = "estimate_type", value.name = "rho")
plot_data[, estimate_type := factor(estimate_type, levels = c("unadjusted_rho", "blood_adjusted_rho"), labels = c("Unadjusted", "Blood-adjusted"))]
plot_data[, feature_label := factor(feature_label, levels = rev(unique(correlation_results$feature_label)))]
correlation_plot <- ggplot(plot_data, aes(rho, feature_label, shape = estimate_type)) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 2.2, position = position_dodge(width = 0.45)) +
  facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Spearman correlation coefficient", y = NULL, shape = NULL) +
  theme_classic(base_size = 10.5, base_family = "Times New Roman") +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(figure_dir, "FigureS14_PDC000614_FCN1_proteomic_validation.png"), correlation_plot, width = 8.2, height = 9.0, units = "in", dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "FigureS14_PDC000614_FCN1_proteomic_validation.tiff"), correlation_plot, width = 8.2, height = 9.0, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")

if (nrow(paired_values) > 0L) {
  paired_long <- melt(paired_values, id.vars = "case_id", measure.vars = c("NAT", "Tumor"), variable.name = "tissue", value.name = "FCN1_protein")
  paired_long[, tissue := factor(tissue, levels = c("NAT", "Tumor"))]
  pair_plot <- ggplot(paired_long, aes(tissue, FCN1_protein, group = case_id)) +
    geom_line(linewidth = 0.45, alpha = 0.55) +
    geom_point(size = 1.8) +
    labs(x = NULL, y = "FCN1 protein abundance\n(unshared log2 ratio)") +
    theme_classic(base_size = 12, base_family = "Times New Roman") +
    theme(axis.title.y = element_text(face = "bold"))
  ggsave(file.path(figure_dir, "FigureS15_PDC000614_FCN1_paired_tumor_NAT.png"), pair_plot, width = 5.8, height = 5.2, units = "in", dpi = 600, bg = "white")
  ggsave(file.path(figure_dir, "FigureS15_PDC000614_FCN1_paired_tumor_NAT.tiff"), pair_plot, width = 5.8, height = 5.2, units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")
}

report <- c(
  "CPTAC-STAD FCN1 PROTEOMIC REPRODUCTION AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Processed tumor input: ", processed_rds),
  paste0("Raw GCT located: ", !is.na(gct_file)),
  paste0("Sample annotation located: ", !is.na(sample_file)),
  "",
  "BLOOD-RATIO ASSOCIATION",
  capture.output(print(blood_result)),
  "",
  "MICROENVIRONMENT CORRELATIONS",
  capture.output(print(correlation_results)),
  "",
  "TARGET AVAILABILITY",
  capture.output(print(target_availability)),
  "",
  "SAMPLE ACCOUNTING",
  capture.output(print(sample_accounting)),
  "",
  "PAIRED TUMOR-NAT",
  if (nrow(paired_results) > 0L) capture.output(print(paired_results)) else "Raw GCT/sample annotation not available; paired analysis was not rerun.",
  "",
  "KEY VALIDATION",
  capture.output(print(validation))
)
writeLines(report, file.path(audit_dir, "CPTAC_FCN1_proteomic_reproduction_audit.txt"), useBytes = TRUE)

cat("CPTAC FCN1 proteomic analysis completed.\n")
print(blood_result)
print(correlation_results)
if (nrow(paired_results) > 0L) print(paired_results)
