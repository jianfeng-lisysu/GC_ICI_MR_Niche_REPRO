# ============================================================
# Figure 1：修正图内措辞并将面板字母改为小写
#
# 修正内容：
# 1. A/B 改为 a/b
# 2. 空间队列描述改为九个独立原发肿瘤的主推断口径
# 3. 总结横幅降级为探索性方向关系
#
# 输入并覆盖：
# 05_figures/Figure1_study_design_and_MR_candidates.tiff
#
# 同步输出：
# 05_figures/Figure1_study_design_and_MR_candidates.png
#
# 说明：
# 该脚本只修改图像内文字，不改变统计结果、森林图或其他图形元素。
# ============================================================

# ------------------------------------------------------------
# 1. 安装并载入 magick
# ------------------------------------------------------------

if (!requireNamespace("magick", quietly = TRUE)) {
  install.packages(
    "magick",
    repos = "https://cloud.r-project.org",
    type = "binary"
  )
}

library(magick)

# ------------------------------------------------------------
# 2. 项目路径
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

if (!dir.exists(project_dir)) {
  stop(
    "项目目录不存在：",
    project_dir
  )
}

figure_dir <- file.path(
  project_dir,
  "05_figures"
)

input_tiff <- file.path(
  figure_dir,
  "Figure1_study_design_and_MR_candidates.tiff"
)

output_png <- file.path(
  figure_dir,
  "Figure1_study_design_and_MR_candidates.png"
)

if (!file.exists(input_tiff)) {
  stop(
    "找不到 Figure 1 TIFF：",
    input_tiff
  )
}

# ------------------------------------------------------------
# 3. 辅助函数
# ------------------------------------------------------------

make_text_patch <- function(
    width,
    height,
    background,
    lines,
    font_size,
    font_color,
    font_weight = 700,
    font_family = "Arial",
    line_height = NULL
) {

  patch <- image_blank(
    width = width,
    height = height,
    color = background
  )

  number_of_lines <- length(lines)

  if (is.null(line_height)) {
    line_height <- floor(
      height /
        number_of_lines
    )
  }

  total_text_height <- line_height *
    number_of_lines

  top_offset <- floor(
    (
      height -
        total_text_height
    ) /
      2
  )

  for (
    line_index in seq_along(lines)
  ) {

    line_canvas <- image_blank(
      width = width,
      height = line_height,
      color = "none"
    )

    line_canvas <- image_annotate(
      image = line_canvas,
      text = lines[[line_index]],
      gravity = "center",
      size = font_size,
      font = font_family,
      weight = font_weight,
      color = font_color
    )

    line_y <- top_offset +
      (
        line_index -
          1L
      ) *
      line_height

    patch <- image_composite(
      image = patch,
      composite_image = line_canvas,
      operator = "over",
      offset = paste0(
        "+0+",
        line_y
      )
    )
  }

  patch
}

place_patch <- function(
    image,
    patch,
    x,
    y
) {

  image_composite(
    image = image,
    composite_image = patch,
    operator = "over",
    offset = paste0(
      "+",
      as.integer(x),
      "+",
      as.integer(y)
    )
  )
}

# ------------------------------------------------------------
# 4. 读取并检查正式 Figure 1
# ------------------------------------------------------------

figure1 <- image_read(
  input_tiff
) |>
  image_background(
    color = "white",
    flatten = TRUE
  ) |>
  image_convert(
    colorspace = "sRGB"
  )

figure1_information <- image_info(
  figure1
)

expected_width <- 9600L
expected_height <- 7341L

if (
  figure1_information$width[[1L]] != expected_width ||
    figure1_information$height[[1L]] != expected_height
) {
  stop(
    paste0(
      "Figure 1 尺寸与锁定版不一致。\n",
      "观察到：",
      figure1_information$width[[1L]],
      " × ",
      figure1_information$height[[1L]],
      "\n预期：",
      expected_width,
      " × ",
      expected_height
    )
  )
}

# ------------------------------------------------------------
# 5. 面板字母 A/B 改为 a/b
# ------------------------------------------------------------

label_a_patch <- make_text_patch(
  width = 240L,
  height = 210L,
  background = "white",
  lines = "a",
  font_size = 170L,
  font_color = "black",
  font_weight = 700L,
  line_height = 210L
)

figure1 <- place_patch(
  image = figure1,
  patch = label_a_patch,
  x = 40L,
  y = 20L
)

# b 标签覆盖区域从横幅下边框以下开始，避免擦除边框。
label_b_patch <- make_text_patch(
  width = 260L,
  height = 250L,
  background = "white",
  lines = "b",
  font_size = 170L,
  font_color = "black",
  font_weight = 700L,
  line_height = 250L
)

figure1 <- place_patch(
  image = figure1,
  patch = label_b_patch,
  x = 950L,
  y = 3750L
)

# ------------------------------------------------------------
# 6. 修正空间队列主推断口径
# ------------------------------------------------------------

spatial_patch <- make_text_patch(
  width = 1725L,
  height = 315L,
  background = "white",
  lines = c(
    "GSE251950 spatial transcriptomics",
    paste0(
      "9 primary tumors from 9 patients ",
      "(29,808 spots);"
    ),
    paste0(
      "1 paired metastatic section excluded ",
      "from primary inference"
    )
  ),
  font_size = 60L,
  font_color = "#2D2D2D",
  font_weight = 400L,
  line_height = 95L
)

figure1 <- place_patch(
  image = figure1,
  patch = spatial_patch,
  x = 5770L,
  y = 1900L
)

# ------------------------------------------------------------
# 7. 修正总结横幅中的过度宣称
# ------------------------------------------------------------

banner_patch <- make_text_patch(
  width = 5900L,
  height = 365L,
  background = "#F7F8FA",
  lines = c(
    paste0(
      "FCN1 and TNFSF12 mark a stroma–myeloid ",
      "microenvironmental program"
    ),
    paste0(
      "showing an exploratory directional ",
      "relationship"
    ),
    paste0(
      "with anti–PD-1 non-response in ",
      "gastric cancer."
    )
  ),
  font_size = 104L,
  font_color = "#27313D",
  font_weight = 700L,
  line_height = 115L
)

figure1 <- place_patch(
  image = figure1,
  patch = banner_patch,
  x = 2600L,
  y = 3240L
)

# ------------------------------------------------------------
# 8. 覆盖保存 600 dpi TIFF
# ------------------------------------------------------------

image_write(
  image = figure1,
  path = input_tiff,
  format = "tiff",
  compression = "lzw",
  density = "600x600"
)

# ------------------------------------------------------------
# 9. 覆盖保存 300 dpi PNG
# ------------------------------------------------------------

figure1_png <- image_resize(
  figure1,
  "4800x"
)

image_write(
  image = figure1_png,
  path = output_png,
  format = "png",
  density = "300x300"
)

# ------------------------------------------------------------
# 10. 输出检查
# ------------------------------------------------------------

if (!file.exists(input_tiff)) {
  stop(
    "Figure 1 TIFF 保存失败。"
  )
}

if (!file.exists(output_png)) {
  stop(
    "Figure 1 PNG 保存失败。"
  )
}

final_tiff_information <- image_info(
  image_read(
    input_tiff
  )
)

final_png_information <- image_info(
  image_read(
    output_png
  )
)

cat(
  "\nFigure 1 修正完成。\n"
)

cat(
  "TIFF：",
  input_tiff,
  "\n"
)

cat(
  "TIFF 尺寸：",
  final_tiff_information$width[[1L]],
  " × ",
  final_tiff_information$height[[1L]],
  " pixels\n"
)

cat(
  "PNG：",
  output_png,
  "\n"
)

cat(
  "PNG 尺寸：",
  final_png_information$width[[1L]],
  " × ",
  final_png_information$height[[1L]],
  " pixels\n"
)

figure1_png
