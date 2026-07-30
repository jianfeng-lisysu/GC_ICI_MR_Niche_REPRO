# ============================================================
# Export software versions for manuscript and GitHub audit
#
# Project:
# GC_ICI_MR_Niche_REPRO_20260617
#
# Purpose:
# 1. Record the exact R, Bioconductor, operating-system, and
#    package versions in the environment used for analysis.
# 2. Create machine-readable and human-readable audit files.
# 3. Do not install, update, or modify any package.
#
# Output:
# 06_logs_and_audit/SOFTWARE_VERSIONS.tsv
# 06_logs_and_audit/SOFTWARE_VERSIONS.txt
# 06_logs_and_audit/SESSION_INFO.txt
# ============================================================

rm(list = ls())

project_dir <- paste0(
    "E:/Bioinformatics analysis/",
    "GC_ICI_MR_Niche_REPRO_20260617"
)

audit_dir <- file.path(
    project_dir,
    "06_logs_and_audit"
)

if (!dir.exists(audit_dir)) {
    dir.create(
        audit_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )
}

output_tsv <- file.path(
    audit_dir,
    "SOFTWARE_VERSIONS.tsv"
)

output_txt <- file.path(
    audit_dir,
    "SOFTWARE_VERSIONS.txt"
)

session_txt <- file.path(
    audit_dir,
    "SESSION_INFO.txt"
)

# Packages directly named in the manuscript plus packages used
# in scoring, modelling, plotting, Word output, and auditing.
packages_to_audit <- unique(c(
    "TwoSampleMR",
    "ieugwasr",
    "coloc",
    "data.table",
    "dplyr",
    "tidyr",
    "readr",
    "readxl",
    "openxlsx",
    "Seurat",
    "SeuratObject",
    "Matrix",
    "DESeq2",
    "tximport",
    "biomaRt",
    "clusterProfiler",
    "AnnotationDbi",
    "enrichplot",
    "fgsea",
    "org.Hs.eg.db",
    "GSVA",
    "limma",
    "edgeR",
    "estimate",
    "preprocessCore",
    "logistf",
    "pROC",
    "boot",
    "survival",
    "ggplot2",
    "ggrepel",
    "patchwork",
    "cowplot",
    "scales",
    "officer",
    "flextable",
    "DiagrammeR",
    "DiagrammeRsvg",
    "rsvg",
    "png",
    "tiff",
    "digest"
))

get_package_version <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        return(NA_character_)
    }

    as.character(
        utils::packageVersion(pkg)
    )
}

package_versions <- data.frame(
    package = packages_to_audit,
    version = vapply(
        packages_to_audit,
        get_package_version,
        character(1)
    ),
    installed = vapply(
        packages_to_audit,
        requireNamespace,
        logical(1),
        quietly = TRUE
    ),
    stringsAsFactors = FALSE
)

package_versions <- package_versions[
    order(
        package_versions$package
    ),
    ,
    drop = FALSE
]

r_version <- paste0(
    R.version$version.string
)

platform <- R.version$platform

operating_system <- paste(
    Sys.info()[c(
        "sysname",
        "release",
        "version",
        "machine"
    )],
    collapse = " | "
)

bioconductor_version <- NA_character_

if (requireNamespace("BiocManager", quietly = TRUE)) {
    bioconductor_version <- as.character(
        BiocManager::version()
    )
}

metadata <- data.frame(
    item = c(
        "Generated",
        "R version",
        "Platform",
        "Operating system",
        "Bioconductor version",
        "Working library paths"
    ),
    value = c(
        format(
            Sys.time(),
            "%Y-%m-%d %H:%M:%S %Z"
        ),
        r_version,
        platform,
        operating_system,
        bioconductor_version,
        paste(
            .libPaths(),
            collapse = " | "
        )
    ),
    stringsAsFactors = FALSE
)

utils::write.table(
    package_versions,
    file = output_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NOT_INSTALLED",
    fileEncoding = "UTF-8"
)

human_readable <- c(
    "SOFTWARE VERSION AUDIT",
    "============================================================",
    "",
    paste0(
        metadata$item,
        ": ",
        metadata$value
    ),
    "",
    "PACKAGE VERSIONS",
    "------------------------------------------------------------",
    capture.output(
        print(
            package_versions,
            row.names = FALSE
        )
    ),
    "",
    "MANUSCRIPT-CORE PACKAGE SUMMARY",
    "------------------------------------------------------------"
)

core_packages <- c(
    "TwoSampleMR",
    "coloc",
    "DESeq2",
    "clusterProfiler",
    "Seurat",
    "logistf",
    "pROC",
    "data.table",
    "ggplot2",
    "ggrepel",
    "patchwork",
    "officer",
    "flextable"
)

core_table <- package_versions[
    match(
        core_packages,
        package_versions$package
    ),
    ,
    drop = FALSE
]

human_readable <- c(
    human_readable,
    capture.output(
        print(
            core_table,
            row.names = FALSE
        )
    )
)

writeLines(
    human_readable,
    con = output_txt,
    useBytes = TRUE
)

session_output <- c(
    paste0(
        "Generated: ",
        format(
            Sys.time(),
            "%Y-%m-%d %H:%M:%S %Z"
        )
    ),
    "",
    capture.output(
        sessionInfo()
    )
)

writeLines(
    session_output,
    con = session_txt,
    useBytes = TRUE
)

if (
    !file.exists(output_tsv) ||
    !file.exists(output_txt) ||
    !file.exists(session_txt)
) {
    stop(
        "软件版本审计文件生成失败。"
    )
}

cat("\n")
cat("============================================================\n")
cat("软件版本审计完成\n")
cat("============================================================\n")
cat("TSV：\n", normalizePath(output_tsv, winslash = "/"), "\n\n", sep = "")
cat("TXT：\n", normalizePath(output_txt, winslash = "/"), "\n\n", sep = "")
cat("Session info：\n", normalizePath(session_txt, winslash = "/"), "\n\n", sep = "")

cat("核心程序包版本：\n")
print(
    core_table,
    row.names = FALSE
)

missing_core <- core_table$package[
    is.na(
        core_table$version
    )
]

if (length(missing_core) > 0) {
    cat(
        "\n注意：以下核心程序包在当前R库中未检测到：\n",
        paste(
            missing_core,
            collapse = ", "
        ),
        "\n",
        sep = ""
    )
}
