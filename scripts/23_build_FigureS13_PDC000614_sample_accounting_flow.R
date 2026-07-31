# ============================================================
# Supplementary Figure S13
# PDC000614 CPTAC-STAD蛋白质组中样本匹配、
# 质量控制及FCN1蛋白分析队列构建流程
#
# 本次修复：
# 1. 删除Matched到Tumor的重复箭头
# 2. 将蛋白矩阵与Derived_Features设为两条独立来源
# 3. 正确表达：
#    212个矩阵样本 - 6个matrix-only
#    212个注释记录 - 6个annotation-only
#    共同形成206个匹配样本
# 4. 优化节点间距和流程对称性
# 5. 直接覆盖正式Supplementary Figure S13
#
# 输入：
# 03_intermediate/
# PDC000614_sample_accounting_audit.rds
#
# 输出：
# 05_figures/
# FigureS13_PDC000614_sample_accounting_flow.png
#
# 05_figures/
# FigureS13_PDC000614_sample_accounting_flow.tiff
# ============================================================

# ------------------------------------------------------------
# 1. 清理旧图形设备
# ------------------------------------------------------------

graphics.off()

# ------------------------------------------------------------
# 2. 项目和文件路径
# ------------------------------------------------------------

project_dir <- paste0(
    "E:/Bioinformatics analysis/",
    "GC_ICI_MR_Niche_REPRO_20260617"
)

input_rds <- file.path(
    project_dir,
    "03_intermediate",
    "PDC000614_sample_accounting_audit.rds"
)

output_png <- file.path(
    project_dir,
    "05_figures",
    "FigureS13_PDC000614_sample_accounting_flow.png"
)

output_tiff <- file.path(
    project_dir,
    "05_figures",
    "FigureS13_PDC000614_sample_accounting_flow.tiff"
)

if (!file.exists(input_rds)) {

    stop(
        "找不到样本审计RDS：\n",
        input_rds
    )
}

if (!dir.exists(dirname(output_png))) {

    dir.create(
        dirname(output_png),
        recursive = TRUE
    )
}

# ------------------------------------------------------------
# 3. 安装并加载程序包
# ------------------------------------------------------------

required_packages <- c(
    "ggplot2",
    "data.table",
    "grid"
)

for (pkg in required_packages) {

    if (!requireNamespace(pkg, quietly = TRUE)) {

        install.packages(
            pkg
        )
    }
}

library(ggplot2)
library(data.table)
library(grid)

# ------------------------------------------------------------
# 4. 读取Table S17样本审计结果
# ------------------------------------------------------------

audit_object <- readRDS(
    input_rds
)

if (
    !is.list(audit_object) ||
        !"audit_table" %in% names(audit_object)
) {

    stop(
        "样本审计RDS结构不正确，缺少audit_table。"
    )
}

audit_table <- as.data.table(
    audit_object$audit_table
)

required_columns <- c(
    "Stage",
    "Count"
)

missing_columns <- setdiff(
    required_columns,
    colnames(audit_table)
)

if (length(missing_columns) > 0) {

    stop(
        "audit_table缺少字段：",
        paste(
            missing_columns,
            collapse = ", "
        )
    )
}

audit_table[
    ,
    Count := suppressWarnings(
        as.integer(
            Count
        )
    )
]

# ------------------------------------------------------------
# 5. 提取样本数量
# ------------------------------------------------------------

get_count_exact <- function(stage_name) {

    target_value <- audit_table[
        Stage == stage_name,
        Count
    ]

    if (length(target_value) != 1) {

        stop(
            "无法唯一匹配审计阶段：",
            stage_name
        )
    }

    as.integer(
        target_value
    )
}

get_count_pattern <- function(stage_pattern) {

    target_value <- audit_table[
        grepl(
            pattern = stage_pattern,
            x = Stage,
            ignore.case = TRUE
        ),
        Count
    ]

    if (length(target_value) != 1) {

        stop(
            "无法唯一匹配审计阶段模式：",
            stage_pattern
        )
    }

    as.integer(
        target_value
    )
}

n_all_unshared <- get_count_exact(
    "All unshared quantitative channels"
)

n_reference <- get_count_exact(
    "NCI7 reference channels"
)

n_biological <- get_count_exact(
    "Biological sample channels"
)

n_matched <- get_count_exact(
    "Samples matched to Derived_Features"
)

n_matrix_only <- get_count_exact(
    "Matrix-only samples without Derived_Features"
)

n_annotation_only <- get_count_exact(
    "Annotation-only records without matrix channels"
)

n_tumor <- get_count_exact(
    "Matched tumor samples"
)

n_nat <- get_count_exact(
    "Matched NAT samples"
)

n_qc_tumor <- get_count_exact(
    "QC-passed tumors with global protein"
)

n_qc_nat <- get_count_exact(
    "QC-passed NAT samples with global protein"
)

n_main_complete <- get_count_pattern(
    "^Complete cases for FCN1 and [0-9]+ microenvironment features$"
)

n_purity_complete <- get_count_exact(
    "Complete cases for FCN1 and ESTIMATE tumor purity"
)

n_annotation_total <- n_matched +
    n_annotation_only

n_tumor_qc_failed <- n_tumor -
    n_qc_tumor

n_nat_qc_failed <- n_nat -
    n_qc_nat

# ------------------------------------------------------------
# 6. 检查关键数量
# ------------------------------------------------------------

expected_counts <- c(
    n_all_unshared = 221,
    n_reference = 9,
    n_biological = 212,
    n_annotation_total = 212,
    n_matched = 206,
    n_matrix_only = 6,
    n_annotation_only = 6,
    n_tumor = 165,
    n_nat = 41,
    n_qc_tumor = 159,
    n_qc_nat = 30,
    n_main_complete = 151,
    n_purity_complete = 120
)

observed_counts <- c(
    n_all_unshared = n_all_unshared,
    n_reference = n_reference,
    n_biological = n_biological,
    n_annotation_total = n_annotation_total,
    n_matched = n_matched,
    n_matrix_only = n_matrix_only,
    n_annotation_only = n_annotation_only,
    n_tumor = n_tumor,
    n_nat = n_nat,
    n_qc_tumor = n_qc_tumor,
    n_qc_nat = n_qc_nat,
    n_main_complete = n_main_complete,
    n_purity_complete = n_purity_complete
)

if (any(
    expected_counts != observed_counts
)) {

    warning(
        "部分样本数量与预期不一致，请检查输入RDS。"
    )
}

# ------------------------------------------------------------
# 7. 定义流程图节点
# ------------------------------------------------------------

node_data <- data.frame(

    node_id = c(
        "all_channels",
        "reference_channels",
        "biological_channels",
        "annotation_records",
        "matrix_only",
        "annotation_only",
        "matched",
        "tumor",
        "nat",
        "tumor_failed",
        "nat_failed",
        "qc_tumor",
        "qc_nat",
        "main_complete",
        "purity_complete"
    ),

    x = c(
        -2.20,
        -4.55,
        -2.20,
        2.20,
        -4.55,
        4.55,
        0.00,
        -1.90,
        1.90,
        -4.20,
        4.20,
        -1.90,
        1.90,
        -1.90,
        1.90
    ),

    y = c(
        10.20,
        9.10,
        8.55,
        8.55,
        7.30,
        7.30,
        6.20,
        4.85,
        4.85,
        3.85,
        3.85,
        2.95,
        2.95,
        1.20,
        1.20
    ),

    label = c(
        paste0(
            "All unshared quantitative channels\n",
            "n = ",
            n_all_unshared
        ),

        paste0(
            "NCI7 reference channels excluded\n",
            "n = ",
            n_reference
        ),

        paste0(
            "Biological matrix channels retained\n",
            "n = ",
            n_biological
        ),

        paste0(
            "Derived_Features annotation records\n",
            "n = ",
            n_annotation_total
        ),

        paste0(
            "Matrix-only samples\n",
            "without Derived_Features\n",
            "n = ",
            n_matrix_only
        ),

        paste0(
            "Annotation-only records\n",
            "without matrix channels\n",
            "n = ",
            n_annotation_only
        ),

        paste0(
            "Matched matrix and annotation records\n",
            "n = ",
            n_matched
        ),

        paste0(
            "Matched tumor samples\n",
            "n = ",
            n_tumor
        ),

        paste0(
            "Matched NAT samples\n",
            "n = ",
            n_nat
        ),

        paste0(
            "Tumor samples not passing QC\n",
            "n = ",
            n_tumor_qc_failed
        ),

        paste0(
            "NAT samples not passing QC\n",
            "n = ",
            n_nat_qc_failed
        ),

        paste0(
            "QC-passed tumors with global protein\n",
            "Primary FCN1 protein-analysis cohort\n",
            "n = ",
            n_qc_tumor
        ),

        paste0(
            "QC-passed NAT samples\n",
            "Auxiliary reference cohort\n",
            "n = ",
            n_qc_nat
        ),

        paste0(
            "Most FCN1–microenvironment correlations\n",
            "Complete cases, n = ",
            n_main_complete
        ),

        paste0(
            "FCN1–ESTIMATE tumor-purity analysis\n",
            "Complete cases, n = ",
            n_purity_complete
        )
    ),

    node_type = c(
        "starting",
        "excluded",
        "retained",
        "annotation",
        "excluded",
        "excluded",
        "matched",
        "tumor",
        "auxiliary",
        "excluded",
        "excluded",
        "formal",
        "auxiliary",
        "analysis",
        "analysis"
    ),

    stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 8. 定义主要流程箭头
# ------------------------------------------------------------

main_edges <- data.frame(

    x = c(
        -2.20,
        -1.90,
        1.90,
        -0.20,
        0.20,
        -1.90,
        1.90,
        -2.05,
        -1.75
    ),

    y = c(
        9.88,
        8.20,
        8.20,
        5.88,
        5.88,
        4.52,
        4.52,
        2.62,
        2.62
    ),

    xend = c(
        -2.20,
        -0.25,
        0.25,
        -1.60,
        1.60,
        -1.90,
        1.90,
        -1.90,
        1.60
    ),

    yend = c(
        8.90,
        6.52,
        6.52,
        5.18,
        5.18,
        3.28,
        3.28,
        1.55,
        1.55
    ),

    stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 9. 定义排除流程箭头
# ------------------------------------------------------------

excluded_edges <- data.frame(

    x = c(
        -2.55,
        -2.55,
        2.55,
        -2.20,
        2.20
    ),

    y = c(
        9.92,
        8.25,
        8.25,
        4.58,
        4.58
    ),

    xend = c(
        -4.20,
        -4.15,
        4.15,
        -3.85,
        3.85
    ),

    yend = c(
        9.30,
        7.60,
        7.60,
        4.05,
        4.05
    ),

    stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 10. 设置节点颜色
# ------------------------------------------------------------

node_fill_values <- c(
    "starting" = "#E2E3E5",
    "retained" = "#D9EAF7",
    "annotation" = "#D9EAF7",
    "matched" = "#D9EAF7",
    "tumor" = "#D9EAF7",
    "formal" = "#DDEBD2",
    "analysis" = "#DDEBD2",
    "auxiliary" = "#F7E6D5",
    "excluded" = "#EEEEEE"
)

# ------------------------------------------------------------
# 11. 绘制流程图
# ------------------------------------------------------------

figure_plot <- ggplot() +

    geom_segment(
        data = main_edges,
        mapping = aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend
        ),
        linewidth = 0.78,
        color = "black",
        lineend = "round",
        arrow = grid::arrow(
            length = grid::unit(
                0.13,
                "inches"
            ),
            type = "closed"
        )
    ) +

    geom_segment(
        data = excluded_edges,
        mapping = aes(
            x = x,
            y = y,
            xend = xend,
            yend = yend
        ),
        linewidth = 0.65,
        color = "grey45",
        lineend = "round",
        arrow = grid::arrow(
            length = grid::unit(
                0.11,
                "inches"
            ),
            type = "closed"
        )
    ) +

    geom_label(
        data = node_data,
        mapping = aes(
            x = x,
            y = y,
            label = label,
            fill = node_type
        ),
        family = "Times New Roman",
        size = 3.25,
        lineheight = 1.02,
        label.size = 0.55,
        label.padding = grid::unit(
            0.15,
            "lines"
        ),
        label.r = grid::unit(
            0.08,
            "lines"
        ),
        color = "black"
    ) +

    scale_fill_manual(
        values = node_fill_values,
        guide = "none"
    ) +

    coord_cartesian(
        xlim = c(
            -5.70,
            5.70
        ),
        ylim = c(
            0.45,
            10.65
        ),
        clip = "off"
    ) +

    theme_void(
        base_family = "Times New Roman"
    ) +

    theme(
        plot.margin = margin(
            t = 18,
            r = 24,
            b = 18,
            l = 24
        )
    )

# ------------------------------------------------------------
# 12. 保存600 dpi PNG
# ------------------------------------------------------------

ggsave(
    filename = output_png,
    plot = figure_plot,
    width = 9.5,
    height = 9.2,
    units = "in",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
)

# ------------------------------------------------------------
# 13. 保存600 dpi TIFF
# ------------------------------------------------------------

ggsave(
    filename = output_tiff,
    plot = figure_plot,
    width = 9.5,
    height = 9.2,
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
)

# ------------------------------------------------------------
# 14. 检查输出
# ------------------------------------------------------------

if (!file.exists(output_png)) {

    stop(
        "Supplementary Figure S13 PNG生成失败。"
    )
}

if (!file.exists(output_tiff)) {

    stop(
        "Supplementary Figure S13 TIFF生成失败。"
    )
}

# ------------------------------------------------------------
# 15. 控制台输出
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("Supplementary Figure S13已重新生成并覆盖原文件\n")
cat("============================================================\n")

cat(
    "蛋白矩阵生物学通道：",
    n_biological,
    "\n"
)

cat(
    "Derived_Features注释记录：",
    n_annotation_total,
    "\n"
)

cat(
    "成功匹配样本：",
    n_matched,
    "\n"
)

cat(
    "质控合格肿瘤：",
    n_qc_tumor,
    "\n"
)

cat(
    "主要相关分析完整病例：",
    n_main_complete,
    "\n"
)

cat(
    "肿瘤纯度分析完整病例：",
    n_purity_complete,
    "\n\n"
)

cat(
    "PNG：\n",
    normalizePath(
        output_png,
        winslash = "/",
        mustWork = FALSE
    ),
    "\n\n",
    sep = ""
)

cat(
    "TIFF：\n",
    normalizePath(
        output_tiff,
        winslash = "/",
        mustWork = FALSE
    ),
    "\n",
    sep = ""
)

