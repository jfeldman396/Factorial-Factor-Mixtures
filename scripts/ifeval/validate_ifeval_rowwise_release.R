#!/usr/bin/env Rscript

# Validate the compact threshold-1 IFEval release without loading large fits.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
release_root <- Sys.getenv(
  "IFEVAL_ROWWISE_SAVED_ROOT",
  file.path(repo_root, "results", "saved", "ifeval", "threshold_1", "rowwise_fewshot")
)
manifest_path <- file.path(repo_root, "results", "saved", "ifeval", "artifact_manifest.csv")

read_release <- function(...) {
  read.csv(file.path(release_root, ...), check.names = FALSE)
}

results <- read_release("model_selection", "rowwise_cv_results.csv")
fold_scores <- read_release("model_selection", "rowwise_cv_fold_scores.csv")
fold_assignments <- read_release("model_selection", "rowwise_cv_fold_assignments.csv")
reveal_plan <- read_release("model_selection", "rowwise_cv_reveal_plan_index.csv")
item_loadings <- read_release(
  "interpretation_model", "fixed_orientation_item_loadings_metadata.csv"
)
gaussian_units <- read_release(
  "gaussian_comparison", "gaussian_rowwise_unit_scores.csv"
)

stopifnot(
  nrow(results) == 462L,
  nrow(fold_scores) == 2310L,
  nrow(fold_assignments) == 122L,
  nrow(reveal_plan) == 366L,
  nrow(item_loadings) == 534L,
  nrow(gaussian_units) == 366L,
  all(results$n_completed_folds == 5L),
  all(nzchar(reveal_plan$reveal_index)),
  all(nzchar(reveal_plan$prediction_index))
)

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
manifest_files <- file.path(repo_root, manifest$path)
stopifnot(
  nrow(manifest) == 17L,
  all(file.exists(manifest_files)),
  identical(unname(tools::md5sum(manifest_files)), manifest$md5)
)

cat("IFEval row-wise release validation passed.\n")
