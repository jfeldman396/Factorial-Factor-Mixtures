#!/usr/bin/env Rscript

# Build the IFEval analysis matrix used by the factor-mixture scripts.
#
# Input:
#   A full OpenEval model-by-item binary matrix, usually
#   openeval_binary_matrix_raw.csv from
#   scripts/data/format_openeval_binary_matrix.py.  The first column is the
#   model name and the remaining columns are item ids from one or more
#   OpenEval benchmarks.
#
# Output:
#   data/ifeval/openeval_ifeval_only_binary_matrix.csv
#   data/ifeval/openeval_item_metadata.csv      if metadata is supplied
#   data/ifeval/openeval_item_instruction_metadata_long.csv if supplied
#   data/ifeval/openeval_model_metadata.csv
#   data/ifeval/ifeval_analysis_matrix_build_summary.csv
#
# The filter is deliberately simple and mirrors the fitted analysis:
# keep IFEval columns; iteratively remove low-coverage models/items; require
# complete item columns for the retained models; then retain nonconstant items,
# i.e. items with at least one 0 and at least one 1.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || identical(x, "")) y else x
}

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else ""
repo_root <- normalizePath(file.path(dirname(script_path %||% "."), "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "scripts"))) {
  repo_root <- normalizePath(getwd(), mustWork = FALSE)
}

full_matrix_path <- Sys.getenv(
  "OPENEVAL_FULL_MATRIX",
  file.path(repo_root, "data", "openeval_ifeval_formatted_uncapped", "openeval_binary_matrix_raw.csv")
)
item_metadata_path <- Sys.getenv(
  "OPENEVAL_ITEM_METADATA",
  file.path(repo_root, "data", "openeval", "openeval_item_metadata.csv")
)
instruction_metadata_path <- Sys.getenv(
  "OPENEVAL_ITEM_INSTRUCTION_METADATA",
  file.path(dirname(item_metadata_path), "openeval_item_instruction_metadata_long.csv")
)
out_dir <- Sys.getenv("OUT_DIR", file.path(repo_root, "data", "ifeval"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
min_item_response_prop <- as.numeric(Sys.getenv("MIN_ITEM_RESPONSE_PROP", "0.25"))
min_model_response_prop <- as.numeric(Sys.getenv("MIN_MODEL_RESPONSE_PROP", "0.25"))
max_filter_iter <- as.integer(Sys.getenv("MAX_COVERAGE_FILTER_ITER", "10"))

if (!file.exists(full_matrix_path)) {
  stop(
    "Full OpenEval matrix not found: ", full_matrix_path, "\n",
    "Create it first with scripts/data/format_openeval_binary_matrix.py, or set OPENEVAL_FULL_MATRIX."
  )
}

raw <- read.csv(full_matrix_path, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(raw) < 2L) stop("Expected a first model-name column plus item columns.")

model_name <- raw[[1L]]
X_all <- as.matrix(raw[, -1L, drop = FALSE])
storage.mode(X_all) <- "numeric"
rownames(X_all) <- model_name

ifeval_cols <- grepl("^ifeval_", colnames(X_all))
X_ifeval_raw <- X_all[, ifeval_cols, drop = FALSE]
if (ncol(X_ifeval_raw) == 0L) stop("No item columns with prefix 'ifeval_' were found.")

X_coverage <- X_ifeval_raw
for (iter in seq_len(max_filter_iter)) {
  old_dim <- dim(X_coverage)
  item_keep <- colMeans(!is.na(X_coverage)) >= min_item_response_prop
  X_coverage <- X_coverage[, item_keep, drop = FALSE]
  model_keep <- rowMeans(!is.na(X_coverage)) >= min_model_response_prop
  X_coverage <- X_coverage[model_keep, , drop = FALSE]
  if (identical(dim(X_coverage), old_dim)) break
}

models_dropped_by_coverage <- setdiff(rownames(X_ifeval_raw), rownames(X_coverage))
items_dropped_by_coverage <- setdiff(colnames(X_ifeval_raw), colnames(X_coverage))

complete_items <- colSums(is.na(X_coverage)) == 0L
X_complete <- X_coverage[, complete_items, drop = FALSE]
items_dropped_for_missing_after_coverage <- setdiff(colnames(X_coverage), colnames(X_complete))

nonconstant_items <- colSums(X_complete == 1) > 0L & colSums(X_complete == 0) > 0L
X_final <- X_complete[, nonconstant_items, drop = FALSE]
items_dropped_constant <- setdiff(colnames(X_complete), colnames(X_final))

analysis_matrix <- data.frame(model_name = rownames(X_final), X_final, check.names = FALSE)
matrix_out <- file.path(out_dir, "openeval_ifeval_only_binary_matrix.csv")
write.csv(analysis_matrix, matrix_out, row.names = FALSE)

model_meta <- data.frame(
  model_name = rownames(X_final),
  n_ifeval_observed_items_before_filter = rowSums(!is.na(X_ifeval_raw[rownames(X_final), , drop = FALSE])),
  n_ifeval_items_after_filter = ncol(X_final),
  mean_correct_ifeval_after_filter = rowMeans(X_final),
  check.names = FALSE
)
write.csv(model_meta, file.path(out_dir, "openeval_model_metadata.csv"), row.names = FALSE)

metadata_rows_written <- NA_integer_
instruction_metadata_rows_written <- NA_integer_
if (file.exists(item_metadata_path)) {
  item_meta <- read.csv(item_metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  if ("item_id" %in% names(item_meta)) {
    item_meta <- item_meta[item_meta$item_id %in% colnames(X_final), , drop = FALSE]
    item_meta <- item_meta[match(colnames(X_final), item_meta$item_id), , drop = FALSE]
    write.csv(item_meta, file.path(out_dir, "openeval_item_metadata.csv"), row.names = FALSE)
    metadata_rows_written <- nrow(item_meta)
  }
}
if (file.exists(instruction_metadata_path)) {
  instruction_meta <- read.csv(instruction_metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  if ("item_id" %in% names(instruction_meta)) {
    instruction_meta <- instruction_meta[instruction_meta$item_id %in% colnames(X_final), , drop = FALSE]
    item_order <- match(instruction_meta$item_id, colnames(X_final))
    instruction_meta <- instruction_meta[order(item_order, instruction_meta$instruction_position), , drop = FALSE]
    write.csv(
      instruction_meta,
      file.path(out_dir, "openeval_item_instruction_metadata_long.csv"),
      row.names = FALSE
    )
    instruction_metadata_rows_written <- nrow(instruction_meta)
  }
}

summary <- data.frame(
  input_matrix = full_matrix_path,
  input_rows = nrow(X_all),
  input_item_columns = ncol(X_all),
  ifeval_item_columns_before_filter = ncol(X_ifeval_raw),
  min_item_response_prop = min_item_response_prop,
  min_model_response_prop = min_model_response_prop,
  dropped_models_by_coverage = length(models_dropped_by_coverage),
  dropped_model_names_by_coverage = paste(models_dropped_by_coverage, collapse = ";"),
  dropped_ifeval_items_by_coverage = length(items_dropped_by_coverage),
  ifeval_item_columns_after_coverage = ncol(X_coverage),
  dropped_ifeval_items_with_missing = length(items_dropped_for_missing_after_coverage),
  ifeval_item_columns_complete = ncol(X_complete),
  dropped_complete_ifeval_items_constant = length(items_dropped_constant),
  dropped_constant_item_ids = paste(items_dropped_constant, collapse = ";"),
  final_ifeval_items = ncol(X_final),
  final_models = nrow(X_final),
  metadata_rows_written = metadata_rows_written,
  instruction_metadata_rows_written = instruction_metadata_rows_written,
  check.names = FALSE
)
write.csv(summary, file.path(out_dir, "ifeval_analysis_matrix_build_summary.csv"), row.names = FALSE)

print(summary)
cat("Wrote IFEval analysis matrix to:", matrix_out, "\n")
