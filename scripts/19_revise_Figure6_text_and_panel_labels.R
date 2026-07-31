# ============================================================
# Figure 6：修正图内措辞并将面板字母改为小写
#
# 修正内容：
# 1. A/B/C 改为 a/b/c
# 2. Panel c 红色文字降级为探索性方向关联
#
# 输入并覆盖：
# 05_figures/Figure6.tiff
#
# 同步输出：
# 05_figures/Figure6.png
#
# 说明：
# 该脚本只修改合成图中的文字，不改变三个面板的统计结果或图形元素。
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
  "Figure6.tiff"
)

output_png <- file.path(
  figure_dir,
  "Figure6.png"
)

if (!file.exists(input_tiff)) {
  stop(
    "找不到 Figure 6 TIFF：",
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
# 4. 读取并检查正式 Figure 6
# ------------------------------------------------------------

figure6 <- image_read(
  input_tiff
) |>
  image_background(
    color = "white",
    flatten = TRUE
  ) |>
  image_convert(
    colorspace = "sRGB"
  )

figure6_information <- image_info(
  figure6
)

expected_width <- 3600L
expected_height <- 4630L

if (
  figure6_information$width[[1L]] != expected_width ||
    figure6_information$height[[1L]] != expected_height
) {
  stop(
    paste0(
      "Figure 6 尺寸与锁定版不一致。\n",
      "观察到：",
      figure6_information$width[[1L]],
      " × ",
      figure6_information$height[[1L]],
      "\n预期：",
      expected_width,
      " × ",
      expected_height
    )
  )
}

# ------------------------------------------------------------
# 5. 面板字母 A/B/C 改为 a/b/c
# ------------------------------------------------------------

label_a_patch <- make_text_patch(
  width = 155L,
  height = 140L,
  background = "white",
  lines = "a",
  font_size = 95L,
  font_color = "black",
  font_weight = 700L,
  line_height = 140L
)

label_b_patch <- make_text_patch(
  width = 155L,
  height = 160L,
  background = "white",
  lines = "b",
  font_size = 95L,
  font_color = "black",
  font_weight = 700L,
  line_height = 160L
)

label_c_patch <- make_text_patch(
  width = 155L,
  height = 165L,
  background = "white",
  lines = "c",
  font_size = 95L,
  font_color = "black",
  font_weight = 700L,
  line_height = 165L
)

figure6 <- place_patch(
  image = figure6,
  patch = label_a_patch,
  x = 45L,
  y = 50L
)

figure6 <- place_patch(
  image = figure6,
  patch = label_b_patch,
  x = 45L,
  y = 2070L
)

figure6 <- place_patch(
  image = figure6,
  patch = label_c_patch,
  x = 45L,
  y = 3740L
)

# ------------------------------------------------------------
# 6. 修正 Panel c 红色过度宣称文字
# ------------------------------------------------------------

panel_c_text_patch <- make_text_patch(
  width = 715L,
  height = 163L,
  background = "#FCFCFC",
  lines = c(
    "Exploratory directional",
    "association with anti–PD-1",
    "non-response"
  ),
  font_size = 52L,
  font_color = "#D21C39",
  font_weight = 700L,
  line_height = 51L
)

figure6 <- place_patch(
  image = figure6,
  patch = panel_c_text_patch,
  x = 2580L,
  y = 4225L
)

# ------------------------------------------------------------
# 7. 覆盖保存 600 dpi TIFF
# ------------------------------------------------------------

image_write(
  image = figure6,
  path = input_tiff,
  format = "tiff",
  compression = "lzw",
  density = "600x600"
)

# ------------------------------------------------------------
# 8. 覆盖保存 300 dpi PNG
# ------------------------------------------------------------

figure6_png <- image_resize(
  figure6,
  "1800x"
)

image_write(
  image = figure6_png,
  path = output_png,
  format = "png",
  density = "300x300"
)

# ------------------------------------------------------------
# 9. 输出检查
# ------------------------------------------------------------

if (!file.exists(input_tiff)) {
  stop(
    "Figure 6 TIFF 保存失败。"
  )
}

if (!file.exists(output_png)) {
  stop(
    "Figure 6 PNG 保存失败。"
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
  "\nFigure 6 修正完成。\n"
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

figure6_png
