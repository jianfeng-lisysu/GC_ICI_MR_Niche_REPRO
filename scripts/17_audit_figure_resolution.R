# ============================================================
# 17_audit_figure_resolution.R
#
# Audits raster source figures in 05_figures without modifying them.
# The default requirement is a conservative effective resolution of
# at least 600 dpi at a 180-mm final print width.
#
# Optional width map:
#   06_logs_and_audit/figure_width_targets.tsv
# with columns:
#   filename    target_width_mm
#
# Outputs:
#   06_logs_and_audit/FIGURE_RESOLUTION_AUDIT.tsv
#   06_logs_and_audit/FIGURE_RESOLUTION_AUDIT.txt
#
# This script does not inspect Word-compressed embedded copies. It audits
# the source PNG/TIFF/JPEG files that should be submitted separately.
# ============================================================

rm(list = ls())
graphics.off()

required_packages <- c("data.table")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Required package is not installed: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
})

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
    "Project directory does not exist:\n",
    project_dir
  )
}

figure_dir <- file.path(project_dir, "05_figures")
audit_dir <- file.path(project_dir, "06_logs_and_audit")

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(figure_dir)) {
  stop("Figure directory does not exist: ", figure_dir)
}

output_tsv <- file.path(
  audit_dir,
  "FIGURE_RESOLUTION_AUDIT.tsv"
)

output_txt <- file.path(
  audit_dir,
  "FIGURE_RESOLUTION_AUDIT.txt"
)

width_map_file <- file.path(
  audit_dir,
  "figure_width_targets.tsv"
)

default_target_width_mm <- 180
required_effective_dpi <- 600

raster_files <- list.files(
  figure_dir,
  pattern = "\\.(png|tif|tiff|jpg|jpeg)$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

if (length(raster_files) == 0L) {
  stop("No PNG, TIFF or JPEG files were found in: ", figure_dir)
}

read_png_dimensions <- function(path) {
  if (!requireNamespace("png", quietly = TRUE)) {
    return(c(NA_integer_, NA_integer_))
  }

  image <- png::readPNG(path, native = TRUE, info = TRUE)
  dimensions <- dim(image)

  if (length(dimensions) < 2L) {
    return(c(NA_integer_, NA_integer_))
  }

  c(
    width = as.integer(dimensions[[2L]]),
    height = as.integer(dimensions[[1L]])
  )
}

read_tiff_dimensions <- function(path) {
  if (!requireNamespace("tiff", quietly = TRUE)) {
    return(c(NA_integer_, NA_integer_))
  }

  image <- tiff::readTIFF(
    path,
    native = TRUE,
    info = TRUE,
    all = FALSE
  )

  dimensions <- dim(image)

  if (length(dimensions) < 2L) {
    return(c(NA_integer_, NA_integer_))
  }

  c(
    width = as.integer(dimensions[[2L]]),
    height = as.integer(dimensions[[1L]])
  )
}

read_jpeg_dimensions <- function(path) {
  if (!requireNamespace("jpeg", quietly = TRUE)) {
    return(c(NA_integer_, NA_integer_))
  }

  image <- jpeg::readJPEG(path, native = TRUE)
  dimensions <- dim(image)

  if (length(dimensions) < 2L) {
    return(c(NA_integer_, NA_integer_))
  }

  c(
    width = as.integer(dimensions[[2L]]),
    height = as.integer(dimensions[[1L]])
  )
}

read_dimensions <- function(path) {
  extension <- tolower(tools::file_ext(path))

  if (extension == "png") {
    return(read_png_dimensions(path))
  }

  if (extension %in% c("tif", "tiff")) {
    return(read_tiff_dimensions(path))
  }

  if (extension %in% c("jpg", "jpeg")) {
    return(read_jpeg_dimensions(path))
  }

  c(width = NA_integer_, height = NA_integer_)
}

if (file.exists(width_map_file)) {
  width_map <- fread(width_map_file, check.names = FALSE)

  required_width_columns <- c(
    "filename",
    "target_width_mm"
  )

  missing_width_columns <- setdiff(
    required_width_columns,
    names(width_map)
  )

  if (length(missing_width_columns) > 0L) {
    stop(
      "Width-map file is missing columns: ",
      paste(missing_width_columns, collapse = ", ")
    )
  }

  width_map[, filename := basename(as.character(filename))]
  width_map[, target_width_mm := as.numeric(target_width_mm)]
} else {
  width_map <- data.table(
    filename = character(),
    target_width_mm = numeric()
  )
}

audit_rows <- lapply(
  raster_files,
  function(path) {
    dimensions <- read_dimensions(path)
    file_name <- basename(path)

    mapped_width <- width_map[
      filename == file_name,
      target_width_mm
    ]

    target_width_mm <- if (
      length(mapped_width) == 1L &&
      is.finite(mapped_width[[1L]]) &&
      mapped_width[[1L]] > 0
    ) {
      mapped_width[[1L]]
    } else {
      default_target_width_mm
    }

    width_px <- unname(dimensions[["width"]])
    height_px <- unname(dimensions[["height"]])

    effective_dpi <- if (
      is.finite(width_px) &&
      width_px > 0
    ) {
      width_px / (target_width_mm / 25.4)
    } else {
      NA_real_
    }

    maximum_width_mm_at_600dpi <- if (
      is.finite(width_px) &&
      width_px > 0
    ) {
      width_px / required_effective_dpi * 25.4
    } else {
      NA_real_
    }

    status <- if (!is.finite(effective_dpi)) {
      "UNREADABLE"
    } else if (effective_dpi >= required_effective_dpi) {
      "PASS"
    } else {
      "FAIL"
    }

    data.table(
      filename = file_name,
      format = toupper(tools::file_ext(path)),
      width_px = width_px,
      height_px = height_px,
      target_width_mm = target_width_mm,
      effective_dpi_at_target_width = effective_dpi,
      maximum_width_mm_at_600dpi = maximum_width_mm_at_600dpi,
      required_dpi = required_effective_dpi,
      status = status,
      file_size_bytes = file.info(path)$size,
      full_path = normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }
)

audit_table <- rbindlist(
  audit_rows,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  audit_table,
  status,
  filename
)

fwrite(
  audit_table,
  output_tsv,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

status_counts <- audit_table[, .N, by = status]

potential_duplicate_files <- audit_table[
  grepl(
    "\\([0-9]+\\)|FINAL|REVISED|_v[0-9]+|v[0-9]+$",
    filename,
    ignore.case = TRUE
  ),
  filename
]

report_lines <- c(
  "FIGURE RESOLUTION AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Figure directory: ", figure_dir),
  paste0("Default final width: ", default_target_width_mm, " mm"),
  paste0("Required effective resolution: ", required_effective_dpi, " dpi"),
  "",
  "STATUS COUNTS",
  capture.output(print(status_counts, row.names = FALSE)),
  "",
  "FAILED OR UNREADABLE FILES",
  if (nrow(audit_table[status != "PASS"]) == 0L) {
    "None"
  } else {
    capture.output(
      print(
        audit_table[
          status != "PASS",
          .(
            filename,
            width_px,
            height_px,
            target_width_mm,
            effective_dpi_at_target_width,
            maximum_width_mm_at_600dpi,
            status
          )
        ],
        row.names = FALSE
      )
    )
  },
  "",
  "POTENTIAL DUPLICATE OR VERSIONED FILENAMES",
  if (length(potential_duplicate_files) == 0L) {
    "None"
  } else {
    paste0("- ", potential_duplicate_files)
  },
  "",
  "INTERPRETATION",
  paste0(
    "PASS means the raster width supports at least ",
    required_effective_dpi,
    " dpi when printed at the configured target width."
  ),
  paste0(
    "The audit is intentionally conservative. A vector PDF is not ",
    "assessed by this raster script."
  ),
  "",
  paste0("Detailed TSV: ", output_tsv)
)

writeLines(
  report_lines,
  output_txt,
  useBytes = TRUE
)

cat("\nFigure-resolution audit completed.\n")
cat("Files audited:", nrow(audit_table), "\n")
cat("PASS:", sum(audit_table$status == "PASS"), "\n")
cat("FAIL:", sum(audit_table$status == "FAIL"), "\n")
cat("UNREADABLE:", sum(audit_table$status == "UNREADABLE"), "\n")
cat("TSV:", output_tsv, "\n")
cat("Report:", output_txt, "\n")
