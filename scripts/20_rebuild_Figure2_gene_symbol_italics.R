# ============================================================
# Figure 2
# GSE167297单细胞数据中FCN1和TNFSF12的细胞来源
#
# Figure 2A：
# GSE167297胃癌单细胞数据中注释细胞类型的UMAP分布
#
# Figure 2B：
# FCN1和TNFSF12在不同细胞类型中的表达水平与表达比例
#
# Figure 2C：
# FCN1和TNFSF12在单细胞UMAP空间中的表达分布
#
# 变量定义及赋值
# 数据来源：
# 03_intermediate/scRNA_GSE167297_processed_annotated.rds
#
# celltype：
# 原始细胞类型注释
#
# celltype_display：
# B_cell      → B cell
# Other       → Other cells
# T_NK        → T/NK
# Epithelial  → Epithelial
# Fibroblast  → Fibroblast
# Endothelial → Endothelial
# Myeloid     → Myeloid
#
# Figure 2B：
# 点颜色表示平均表达水平
# 点大小表示表达该基因的细胞比例
#
# Figure 2C：
# 颜色表示Seurat标准化表达值
# 每个基因使用独立色阶
# 不用于比较FCN1和TNFSF12之间的绝对表达强度
#
# 统计学方法：
# 本图为单细胞表达分布的描述性可视化，
# 不进行新的假设检验。
# ============================================================


# ------------------------------------------------------------
# 1. 检查并载入R包
# ------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "ggplot2",
  "patchwork",
  "scales",
  "grid"
)

for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org",
      type = "binary"
    )
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("程序包安装失败：", pkg)
  }
}

library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)


# ------------------------------------------------------------
# 2. 项目路径
# ------------------------------------------------------------

proj <- paste0(
  "E:/Bioinformatics analysis/",
  "GC_ICI_MR_Niche_REPRO_20260617"
)

fig_dir <- file.path(
  proj,
  "05_figures"
)

mid_dir <- file.path(
  proj,
  "03_intermediate"
)

input_file <- file.path(
  mid_dir,
  "scRNA_GSE167297_processed_annotated.rds"
)

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop("找不到单细胞RDS文件：", input_file)
}


# ------------------------------------------------------------
# 3. 读取单细胞对象
# ------------------------------------------------------------

sce <- readRDS(input_file)

if (!inherits(sce, "Seurat")) {
  stop("输入文件不是Seurat对象。")
}

if (!"celltype" %in% colnames(sce@meta.data)) {
  stop("Seurat对象的meta.data中不存在celltype变量。")
}

required_genes <- c(
  "FCN1",
  "TNFSF12"
)

missing_genes <- setdiff(
  required_genes,
  rownames(sce)
)

if (length(missing_genes) > 0) {

  stop(
    "Seurat对象中缺少以下基因：",
    paste(missing_genes, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 4. 规范细胞类型显示名称
#
# 不修改原始celltype变量，
# 仅新增celltype_display用于作图
# ------------------------------------------------------------

original_celltype_levels <- c(
  "B_cell",
  "Other",
  "T_NK",
  "Epithelial",
  "Fibroblast",
  "Endothelial",
  "Myeloid"
)

display_celltype_levels <- c(
  "B cell",
  "Other cells",
  "T/NK",
  "Epithelial",
  "Fibroblast",
  "Endothelial",
  "Myeloid"
)

unexpected_celltypes <- setdiff(
  unique(as.character(sce$celltype)),
  original_celltype_levels
)

if (length(unexpected_celltypes) > 0) {

  warning(
    "发现未列入预设顺序的细胞类型：",
    paste(unexpected_celltypes, collapse = ", ")
  )
}

sce$celltype_display <- factor(
  as.character(sce$celltype),
  levels = original_celltype_levels,
  labels = display_celltype_levels
)

if (any(is.na(sce$celltype_display))) {

  stop(
    "部分细胞无法映射到celltype_display，",
    "请检查原始celltype取值。"
  )
}

cat("细胞总数：", ncol(sce), "\n\n")

cat("各细胞类型数量：\n")

print(
  table(
    sce$celltype_display,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 5. Figure 2A
# 注释细胞类型UMAP
# ------------------------------------------------------------

p2a <- DimPlot(
  object = sce,
  reduction = "umap",
  group.by = "celltype_display",
  label = FALSE,
  pt.size = 0.25
) +
  labs(
    x = "UMAP_1",
    y = "UMAP_2",
    color = NULL
  ) +
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  ) +
  theme_classic(
    base_size = 13,
    base_family = "sans"
  ) +
  theme(
    plot.title = element_blank(),

    axis.title = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),

    axis.text = element_text(
      color = "black",
      size = 11
    ),

    axis.line = element_line(
      linewidth = 0.6,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.6,
      color = "black"
    ),

    legend.title = element_blank(),

    legend.text = element_text(
      size = 11,
      color = "black"
    ),

    legend.position = "right",

    plot.margin = margin(
      t = 6,
      r = 8,
      b = 6,
      l = 6
    )
  )


# ------------------------------------------------------------
# 6. 保存Figure 2A
# ------------------------------------------------------------

out_2a_png <- file.path(
  fig_dir,
  "Fig2A_scRNA_UMAP_celltype.png"
)

out_2a_tiff <- file.path(
  fig_dir,
  "Fig2A_scRNA_UMAP_celltype.tiff"
)

ggsave(
  filename = out_2a_png,
  plot = p2a,
  width = 7.8,
  height = 5.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_2a_tiff,
  plot = p2a,
  width = 7.8,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 7. Figure 2B
# FCN1和TNFSF12 DotPlot
# ------------------------------------------------------------

p2b <- DotPlot(
  object = sce,
  features = c(
    "FCN1",
    "TNFSF12"
  ),
  group.by = "celltype_display"
) +
  labs(
    x = NULL,
    y = NULL
  ) +
  guides(
    color = guide_colorbar(
      title = "Average\nexpression",
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(
        3.2,
        "cm"
      ),
      barwidth = unit(
        0.5,
        "cm"
      )
    ),

    size = guide_legend(
      title = "Percent\nexpressed",
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme_classic(
    base_size = 10,
    base_family = "sans"
  ) +
  theme(
    plot.title = element_blank(),

    axis.text.x = element_text(
      color = "black",
      size = 11,
      face = "italic"
    ),

    axis.text.y = element_text(
      color = "black",
      size = 10
    ),

    axis.line = element_line(
      linewidth = 0.5,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.5,
      color = "black"
    ),

    legend.title = element_text(
      size = 7,
      color = "black"
    ),

    legend.text = element_text(
      size = 7,
      color = "black"
    ),

    legend.key.size = unit(
      0.35,
      "cm"
    ),

    legend.spacing.y = unit(
      0.05,
      "cm"
    ),

    plot.margin = margin(
      t = 6,
      r = 6,
      b = 6,
      l = 6
    )
  )


# ------------------------------------------------------------
# 8. 保存Figure 2B
# ------------------------------------------------------------

out_2b_png <- file.path(
  fig_dir,
  "Fig2B_scRNA_FCN1_TNFSF12_DotPlot.png"
)

out_2b_tiff <- file.path(
  fig_dir,
  "Fig2B_scRNA_FCN1_TNFSF12_DotPlot.tiff"
)

ggsave(
  filename = out_2b_png,
  plot = p2b,
  width = 4.8,
  height = 3.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_2b_tiff,
  plot = p2b,
  width = 4.8,
  height = 3.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 9. 删除FeaturePlot原有颜色刻度的辅助函数
#
# 目的是重新添加统一标题和统一小数格式，
# 避免出现“Scale for colour is already present”警告
# ------------------------------------------------------------

remove_existing_colour_scale <- function(plot_object) {

  if (length(plot_object$scales$scales) == 0) {
    return(plot_object)
  }

  keep_scale <- vapply(
    plot_object$scales$scales,
    function(scale_object) {

      scale_aesthetics <- scale_object$aesthetics

      !any(
        scale_aesthetics %in% c(
          "colour",
          "color"
        )
      )
    },
    logical(1)
  )

  plot_object$scales$scales <-
    plot_object$scales$scales[keep_scale]

  plot_object
}


# ------------------------------------------------------------
# 10. Figure 2C
# 分别生成两个FeaturePlot
# ------------------------------------------------------------

feature_plot_list <- FeaturePlot(
  object = sce,
  features = c(
    "FCN1",
    "TNFSF12"
  ),
  reduction = "umap",
  combine = FALSE,
  pt.size = 0.18,
  order = TRUE,
  min.cutoff = 0,
  max.cutoff = "q95",
  keep.scale = "feature"
)


# ------------------------------------------------------------
# 11. 统一两个FeaturePlot的图例标题和刻度格式
# ------------------------------------------------------------

feature_plot_list <- lapply(
  seq_along(feature_plot_list),
  function(i) {

    current_plot <- feature_plot_list[[i]]

    current_plot <- remove_existing_colour_scale(
      current_plot
    )

    current_plot +
      scale_color_gradient(
        low = "grey85",
        high = "blue",
        name = "Normalized\nexpression",
        breaks = scales::breaks_pretty(
          n = 5
        ),
        labels = scales::label_number(
          accuracy = 0.1,
          trim = TRUE
        ),
        na.value = "grey85",
        guide = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barheight = unit(
            2.2,
            "cm"
          ),
          barwidth = unit(
            0.35,
            "cm"
          )
        )
      ) +
      labs(
        x = "UMAP_1",
        y = "UMAP_2"
      ) +
      theme_classic(
        base_size = 10,
        base_family = "sans"
      ) +
      theme(
        plot.title = element_text(
          face = "bold.italic",
          hjust = 0.5,
          size = 12,
          color = "black"
        ),

        axis.title = element_text(
          size = 11,
          color = "black"
        ),

        axis.text = element_text(
          color = "black",
          size = 9
        ),

        axis.line = element_line(
          linewidth = 0.5,
          color = "black"
        ),

        axis.ticks = element_line(
          linewidth = 0.5,
          color = "black"
        ),

        legend.title = element_text(
          size = 7.5,
          color = "black"
        ),

        legend.text = element_text(
          size = 8,
          color = "black"
        ),

        plot.margin = margin(
          t = 5,
          r = 6,
          b = 5,
          l = 6
        )
      )
  }
)


# ------------------------------------------------------------
# 12. 合并两个FeaturePlot
#
# keep.scale = "feature"意味着两个基因各自使用独立范围。
# 每幅图保留自己的Normalized expression图例。
# ------------------------------------------------------------

p2c <- wrap_plots(
  feature_plot_list,
  ncol = 2,
  guides = "keep"
)


# ------------------------------------------------------------
# 13. 保存Figure 2C
# ------------------------------------------------------------

out_2c_png <- file.path(
  fig_dir,
  "Fig2C_scRNA_FCN1_TNFSF12_FeaturePlot.png"
)

out_2c_tiff <- file.path(
  fig_dir,
  "Fig2C_scRNA_FCN1_TNFSF12_FeaturePlot.tiff"
)

ggsave(
  filename = out_2c_png,
  plot = p2c,
  width = 7.2,
  height = 3.6,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_2c_tiff,
  plot = p2c,
  width = 7.2,
  height = 3.6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 14. 调整Figure 2B在合成图中的图例尺寸
# ------------------------------------------------------------

p2b_combined <- p2b +
  theme(
    legend.title = element_text(
      size = 7
    ),

    legend.text = element_text(
      size = 7
    ),

    legend.key.size = unit(
      0.35,
      "cm"
    ),

    legend.spacing.y = unit(
      0.05,
      "cm"
    )
  )


# ------------------------------------------------------------
# 15. 合成完整Figure 2
#
# a和b为上排；
# c的两个FeaturePlot作为一个整体面板处理。
# ------------------------------------------------------------

top_row <- p2a +
  p2b_combined +
  plot_layout(
    widths = c(
      1.35,
      0.85
    )
  )

panel_c_wrapped <- wrap_elements(
  full = p2c
)

p_fig2 <- (
  top_row /
    panel_c_wrapped
) +
  plot_layout(
    heights = c(
      1.05,
      1.00
    )
  ) +
  plot_annotation(
    tag_levels = "a"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 14,
      color = "black"
    ),

    plot.tag.position = c(
      0.01,
      0.99
    )
  )


# ------------------------------------------------------------
# 16. 保存完整Figure 2
# 覆盖原有同名PNG和TIFF
# ------------------------------------------------------------

out_fig2_png <- file.path(
  fig_dir,
  "Fig2_scRNA_FCN1_TNFSF12_cellular_source.png"
)

out_fig2_tiff <- file.path(
  fig_dir,
  "Fig2_scRNA_FCN1_TNFSF12_cellular_source.tiff"
)

ggsave(
  filename = out_fig2_png,
  plot = p_fig2,
  width = 11.2,
  height = 8.0,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_fig2_tiff,
  plot = p_fig2,
  width = 11.2,
  height = 8.0,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ------------------------------------------------------------
# 17. 输出检查
# ------------------------------------------------------------

output_files <- c(
  out_2a_png,
  out_2a_tiff,
  out_2b_png,
  out_2b_tiff,
  out_2c_png,
  out_2c_tiff,
  out_fig2_png,
  out_fig2_tiff
)

output_check <- data.frame(
  file_name = basename(output_files),
  exists = file.exists(output_files),
  size_MB = round(
    file.info(output_files)$size /
      1024 /
      1024,
    2
  ),
  stringsAsFactors = FALSE
)

cat("Figure 2全部生成完成。\n\n")

print(
  output_check,
  row.names = FALSE
)

cat(
  "\n完整Figure 2 PNG：\n",
  out_fig2_png,
  "\n",
  sep = ""
)

cat(
  "\n完整Figure 2 TIFF：\n",
  out_fig2_tiff,
  "\n",
  sep = ""
)

p_fig2

