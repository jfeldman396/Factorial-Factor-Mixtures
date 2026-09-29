#!/usr/bin/env Rscript

# Export the compact, analysis-ready release for the primary threshold-1
# row-wise few-shot IFEval study. Large RDS fits and optimization diagnostics
# remain under results/full and are reproducible from REPRODUCE.md.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1L]])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

copy_required <- function(source_root, filenames, destination) {
  sources <- file.path(source_root, filenames)
  missing <- sources[!file.exists(sources)]
  if (length(missing)) {
    stop("Missing release inputs: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  targets <- file.path(destination, filenames)
  ok <- file.copy(sources, targets, overwrite = TRUE)
  if (!all(ok)) stop("Could not copy: ", paste(sources[!ok], collapse = ", "))
  targets
}

rowwise_root <- get_env(
  "IFEVAL_ROWWISE_ROOT",
  file.path(repo_root, "results", "full", "ifeval_rowwise_fewshot_cv")
)
fixed_root <- get_env(
  "IFEVAL_FIXED_ROOT",
  file.path(rowwise_root, "alternative_H4_G2-2-3-1_fixed_orientation_lambda4")
)
gaussian_root <- get_env(
  "IFEVAL_GAUSSIAN_ROOT",
  file.path(rowwise_root, "gaussian_H4_lambda4_rowwise_fewshot_cv")
)
saved_root <- get_env(
  "IFEVAL_ROWWISE_SAVED_ROOT",
  file.path(repo_root, "results", "saved", "ifeval", "threshold_1", "rowwise_fewshot")
)
plot_root <- get_env(
  "IFEVAL_ROWWISE_PLOT_ROOT",
  file.path(repo_root, "results", "selected_plots", "ifeval", "threshold_1", "rowwise_fewshot")
)

copied <- character(0)
copied <- c(copied, copy_required(
  rowwise_root,
  c(
    "rowwise_cv_results.csv",
    "rowwise_cv_selected_model.csv",
    "rowwise_cv_fold_scores.csv",
    "rowwise_cv_fold_assignments.csv",
    "rowwise_cv_reveal_plan_index.csv"
  ),
  file.path(saved_root, "model_selection")
))
copied <- c(copied, copy_required(
  fixed_root,
  c(
    "fixed_orientation_factor_summary.csv",
    "fixed_orientation_mixture_parameters.csv",
    "fixed_orientation_item_loadings_metadata.csv",
    "fixed_orientation_family_loading_summary.csv",
    "fixed_orientation_instruction_loading_summary.csv"
  ),
  file.path(saved_root, "interpretation_model")
))
copied <- c(copied, copy_required(
  gaussian_root,
  c(
    "gaussian_rowwise_settings.csv",
    "gaussian_rowwise_summary.csv",
    "gaussian_rowwise_fold_scores.csv",
    "gaussian_rowwise_unit_scores.csv",
    "pifa_gaussian_rowwise_comparison.csv",
    "pifa_gaussian_paired_fold_scores.csv"
  ),
  file.path(saved_root, "gaussian_comparison")
))
copied <- c(copied, copy_required(
  fixed_root,
  "fixed_orientation_factor_marginals.png",
  plot_root
))

absolute_paths <- normalizePath(copied, mustWork = TRUE)
repo_prefix <- paste0(repo_root, .Platform$file.sep)
if (!all(startsWith(absolute_paths, repo_prefix))) {
  stop("Every released artifact must be inside the repository.")
}
manifest <- data.frame(
  path = substring(absolute_paths, nchar(repo_prefix) + 1L),
  bytes = file.info(copied)$size,
  md5 = unname(tools::md5sum(copied)),
  stringsAsFactors = FALSE
)
dir.create(file.path(repo_root, "results", "saved", "ifeval"), recursive = TRUE, showWarnings = FALSE)
write.csv(
  manifest,
  file.path(repo_root, "results", "saved", "ifeval", "artifact_manifest.csv"),
  row.names = FALSE
)
cat("Exported ", nrow(manifest), " row-wise IFEval release artifacts.\n", sep = "")
