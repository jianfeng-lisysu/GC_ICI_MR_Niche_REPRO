# Shared project-path configuration.
#
# Recommended use:
#   1. Run scripts from the project/repository root; or
#   2. Define the environment variable GC_ICI_MR_NICHE_PROJECT.

project_dir <- Sys.getenv(
  "GC_ICI_MR_NICHE_PROJECT",
  unset = normalizePath(".", winslash = "/", mustWork = FALSE)
)

project_dir <- normalizePath(
  project_dir,
  winslash = "/",
  mustWork = FALSE
)

if (!dir.exists(project_dir)) {
  stop("Project directory does not exist: ", project_dir)
}
