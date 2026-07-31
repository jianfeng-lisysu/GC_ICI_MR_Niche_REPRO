## ============================================================================
## Figure 4：GSE251950九位独立原发胃癌患者的空间关联
##
## 修订：
##   Panel A/B：每个空间程序上方标注BH校正P值
##   Panel C：仅标注方向一致性，不重复P值；基因符号斜体
##
## 输出覆盖：
##   05_figures/Figure4_GSE251950_primary9_spatial_inference.png
##   05_figures/Figure4_GSE251950_primary9_spatial_inference.tiff
## ============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
})

PROJECT <- "E:/Bioinformatics analysis/GC_ICI_MR_Niche_REPRO_20260617"

SECTION_FILE <- file.path(
    PROJECT,
    "04_results/03_spatial_results",
    "GSE251950_section_level_candidate_program_spearman.tsv"
)

INFERENCE_FILE <- file.path(
    PROJECT,
    "04_results/03_spatial_results",
    "GSE251950_primary9_section_level_inference.tsv"
)

FIG_DIR <- file.path(
    PROJECT,
    "05_figures"
)

FIGURE4_PNG <- file.path(
    FIG_DIR,
    "Figure4_GSE251950_primary9_spatial_inference.png"
)

FIGURE4_TIFF <- file.path(
    FIG_DIR,
    "Figure4_GSE251950_primary9_spatial_inference.tiff"
)

dir.create(
    FIG_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(SECTION_FILE)) {
    stop("逐切片相关结果不存在：", SECTION_FILE)
}

if (!file.exists(INFERENCE_FILE)) {
    stop("九患者正式推断结果不存在：", INFERENCE_FILE)
}

## --------------------------------------------------------------------------
## 1. 读取数据
## --------------------------------------------------------------------------

section_results <- fread(
    SECTION_FILE,
    data.table = TRUE,
    showProgress = FALSE
)

summary_results <- fread(
    INFERENCE_FILE,
    data.table = TRUE,
    showProgress = FALSE
)

required_section_columns <- c(
    "section",
    "patient",
    "tissue",
    "candidate",
    "program",
    "rho"
)

required_summary_columns <- c(
    "candidate",
    "program",
    "n_patients",
    "median_rho",
    "q1_rho",
    "q3_rho",
    "n_positive",
    "exact_p",
    "bh_adjusted_p"
)

missing_section_columns <- setdiff(
    required_section_columns,
    colnames(section_results)
)

missing_summary_columns <- setdiff(
    required_summary_columns,
    colnames(summary_results)
)

if (length(missing_section_columns) > 0L) {
    stop(
        "逐切片结果缺少列：",
        paste(missing_section_columns, collapse = ", ")
    )
}

if (length(missing_summary_columns) > 0L) {
    stop(
        "正式推断结果缺少列：",
        paste(missing_summary_columns, collapse = ", ")
    )
}

## --------------------------------------------------------------------------
## 2. 保留9张独立原发切片
## --------------------------------------------------------------------------

primary_data <- section_results[
    tolower(as.character(tissue)) == "primary"
]

if (uniqueN(primary_data$patient) != 9L) {
    stop(
        "独立患者数不是9，实际为：",
        uniqueN(primary_data$patient)
    )
}

if (uniqueN(primary_data$section) != 9L) {
    stop(
        "原发切片数不是9，实际为：",
        uniqueN(primary_data$section)
    )
}

if (nrow(primary_data) != 54L) {
    stop(
        "患者级结果应为54行，实际为：",
        nrow(primary_data)
    )
}

if (nrow(summary_results) != 6L) {
    stop(
        "正式汇总结果应为6行，实际为：",
        nrow(summary_results)
    )
}

## --------------------------------------------------------------------------
## 3. 固定显示顺序
## --------------------------------------------------------------------------

candidate_levels <- c(
    "FCN1",
    "TNFSF12"
)

program_levels <- c(
    "CAF",
    "Myeloid",
    "Endothelial"
)

patient_levels <- paste0(
    "P",
    1:9
)

primary_data[
    ,
    candidate := factor(
        as.character(candidate),
        levels = candidate_levels
    )
]

primary_data[
    ,
    program := factor(
        as.character(program),
        levels = program_levels
    )
]

primary_data[
    ,
    patient := factor(
        as.character(patient),
        levels = patient_levels
    )
]

summary_results[
    ,
    candidate := factor(
        as.character(candidate),
        levels = candidate_levels
    )
]

summary_results[
    ,
    program := factor(
        as.character(program),
        levels = program_levels
    )
]

setorder(
    primary_data,
    candidate,
    program,
    patient
)

setorder(
    summary_results,
    candidate,
    program
)

## --------------------------------------------------------------------------
## 4. 生成标签
## --------------------------------------------------------------------------

summary_results[
    ,
    bh_label := sprintf(
        "BH P = %.4f",
        bh_adjusted_p
    )
]

summary_results[
    ,
    direction_label := sprintf(
        "%d/%d positive",
        n_positive,
        n_patients
    )
]

summary_results[
    ,
    pair_label := paste(
        as.character(candidate),
        as.character(program),
        sep = " – "
    )
]

pair_levels <- c(
    "FCN1 – CAF",
    "FCN1 – Myeloid",
    "FCN1 – Endothelial",
    "TNFSF12 – CAF",
    "TNFSF12 – Myeloid",
    "TNFSF12 – Endothelial"
)

summary_results[
    ,
    pair_label := factor(
        pair_label,
        levels = rev(pair_levels)
    )
]

summary_results[
    ,
    y_position := as.numeric(pair_label)
]

## Panel c：仅将人类基因符号设为斜体；空间程序名称保持正体
panel_c_axis_labels <- expression(
    italic(TNFSF12) ~ "\u2013" ~ "Endothelial",
    italic(TNFSF12) ~ "\u2013" ~ "Myeloid",
    italic(TNFSF12) ~ "\u2013" ~ "CAF",
    italic(FCN1) ~ "\u2013" ~ "Endothelial",
    italic(FCN1) ~ "\u2013" ~ "Myeloid",
    italic(FCN1) ~ "\u2013" ~ "CAF"
)

## --------------------------------------------------------------------------
## 5. 统一主题
## --------------------------------------------------------------------------

base_theme <- theme_classic(
    base_size = 11,
    base_family = "Arial"
) +
    theme(
        axis.text = element_text(
            color = "black",
            size = 9
        ),
        axis.title = element_text(
            color = "black",
            size = 10
        ),
        axis.line = element_line(
            linewidth = 0.45,
            color = "black"
        ),
        axis.ticks = element_line(
            linewidth = 0.40,
            color = "black"
        ),
        plot.margin = margin(
            8,
            8,
            8,
            8
        )
    )

## --------------------------------------------------------------------------
## 6. 创建Panel A/B
## --------------------------------------------------------------------------

make_candidate_panel <- function(
    candidate_name,
    point_fill,
    lower_limit,
    upper_limit,
    label_y
) {
    current_data <- primary_data[
        candidate == candidate_name
    ]

    current_summary <- summary_results[
        candidate == candidate_name
    ]

    if (nrow(current_data) != 27L) {
        stop(
            candidate_name,
            "患者级结果不是27行。"
        )
    }

    if (nrow(current_summary) != 3L) {
        stop(
            candidate_name,
            "汇总结果不是3行。"
        )
    }

    ggplot(
        current_data,
        aes(
            x = program,
            y = rho
        )
    ) +
        geom_hline(
            yintercept = 0,
            linetype = "dashed",
            linewidth = 0.45,
            color = "grey45"
        ) +
        geom_point(
            position = position_jitter(
                width = 0.085,
                height = 0,
                seed = 20260727
            ),
            shape = 21,
            size = 2.7,
            stroke = 0.45,
            fill = point_fill,
            color = "black",
            alpha = 0.88
        ) +
        geom_errorbar(
            data = current_summary,
            mapping = aes(
                x = program,
                ymin = q1_rho,
                ymax = q3_rho
            ),
            inherit.aes = FALSE,
            width = 0.18,
            linewidth = 0.75,
            color = "black"
        ) +
        geom_point(
            data = current_summary,
            mapping = aes(
                x = program,
                y = median_rho
            ),
            inherit.aes = FALSE,
            shape = 23,
            size = 3.6,
            stroke = 0.75,
            fill = "white",
            color = "black"
        ) +
        geom_text(
            data = current_summary,
            mapping = aes(
                x = program,
                y = label_y,
                label = bh_label
            ),
            inherit.aes = FALSE,
            size = 2.7,
            color = "black"
        ) +
        coord_cartesian(
            ylim = c(
                lower_limit,
                upper_limit
            ),
            clip = "off"
        ) +
        labs(
            x = NULL,
            y = expression(
                "Section-level Spearman " * rho
            ),
            title = candidate_name
        ) +
        base_theme +
        theme(
            plot.title = element_text(
                face = "italic",
                size = 12,
                hjust = 0.5,
                color = "black"
            )
        )
}

## --------------------------------------------------------------------------
## 7. Panel A：FCN1
## --------------------------------------------------------------------------

panel_a <- make_candidate_panel(
    candidate_name = "FCN1",
    point_fill = "#4C78A8",
    lower_limit = -0.03,
    upper_limit = 0.235,
    label_y = 0.222
)

## --------------------------------------------------------------------------
## 8. Panel B：TNFSF12
## --------------------------------------------------------------------------

panel_b <- make_candidate_panel(
    candidate_name = "TNFSF12",
    point_fill = "#E45756",
    lower_limit = -0.03,
    upper_limit = 0.365,
    label_y = 0.348
)

## --------------------------------------------------------------------------
## 9. Panel C：中位rho和IQR，仅标方向一致性
## --------------------------------------------------------------------------

panel_c <- ggplot(
    summary_results,
    aes(
        x = median_rho,
        y = y_position
    )
) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.45,
        color = "grey45"
    ) +
    geom_segment(
        aes(
            x = q1_rho,
            xend = q3_rho,
            y = y_position,
            yend = y_position
        ),
        linewidth = 0.75,
        color = "black"
    ) +
    geom_segment(
        aes(
            x = q1_rho,
            xend = q1_rho,
            y = y_position - 0.10,
            yend = y_position + 0.10
        ),
        linewidth = 0.65,
        color = "black"
    ) +
    geom_segment(
        aes(
            x = q3_rho,
            xend = q3_rho,
            y = y_position - 0.10,
            yend = y_position + 0.10
        ),
        linewidth = 0.65,
        color = "black"
    ) +
    geom_point(
        aes(
            fill = candidate
        ),
        shape = 21,
        size = 3.4,
        stroke = 0.55,
        color = "black"
    ) +
    geom_text(
        aes(
            x = 0.315,
            label = direction_label
        ),
        hjust = 0,
        size = 2.85,
        color = "black"
    ) +
    scale_fill_manual(
        values = c(
            "FCN1" = "#4C78A8",
            "TNFSF12" = "#E45756"
        ),
        drop = FALSE
    ) +
    scale_y_continuous(
        breaks = 1:6,
        labels = panel_c_axis_labels,
        expand = expansion(
            mult = c(
                0.08,
                0.08
            )
        )
    ) +
    scale_x_continuous(
        breaks = seq(
            0,
            0.4,
            by = 0.1
        )
    ) +
    coord_cartesian(
        xlim = c(
            -0.02,
            0.46
        ),
        clip = "off"
    ) +
    labs(
        x = expression(
            "Median section-level Spearman " * rho
        ),
        y = NULL
    ) +
    base_theme +
    theme(
        legend.position = "none",
        axis.text.y = element_text(
            color = "black",
            size = 8.7,
            face = "plain"
        ),
        plot.margin = margin(
            8,
            12,
            8,
            8
        )
    )

## --------------------------------------------------------------------------
## 10. 合并Figure 4
## --------------------------------------------------------------------------

figure4 <- (
    panel_a |
        panel_b |
        panel_c
) +
    plot_layout(
        widths = c(
            1,
            1,
            1.28
        )
    ) +
    plot_annotation(
        tag_levels = "a",
        theme = theme(
            plot.tag = element_text(
                face = "bold",
                size = 13,
                color = "black"
            )
        )
    )

## --------------------------------------------------------------------------
## 11. 覆盖保存
## --------------------------------------------------------------------------

ggsave(
    filename = FIGURE4_PNG,
    plot = figure4,
    width = 12.0,
    height = 5.3,
    units = "in",
    dpi = 600,
    bg = "white"
)

ggsave(
    filename = FIGURE4_TIFF,
    plot = figure4,
    width = 12.0,
    height = 5.3,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
)

## --------------------------------------------------------------------------
## 12. 最终检查
## --------------------------------------------------------------------------

output_files <- c(
    FIGURE4_PNG,
    FIGURE4_TIFF
)

if (any(!file.exists(output_files))) {
    stop("Figure 4的PNG或TIFF未成功生成。")
}

if (any(file.info(output_files)$size <= 0L)) {
    stop("Figure 4输出文件大小异常。")
}

cat("\n============================================================\n")
cat("Figure 4修订完成\n")
cat("============================================================\n")

for (current_file in output_files) {
    cat(
        basename(current_file),
        " | ",
        round(
            file.info(current_file)$size / 1024 / 1024,
            2
        ),
        " MB\n",
        sep = ""
    )
}

cat(
    "\n输出目录：\n",
    FIG_DIR,
    "\n",
    sep = ""
)

sessionInfo()

