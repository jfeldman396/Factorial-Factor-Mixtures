#!/usr/bin/env Rscript

# Regenerate focused loading-interpretation tables from the fitted IFEval
# independent-mixture probit lambda matrix.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

loadings_path <- Sys.getenv(
  "LOADINGS_PATH",
  unset = file.path(
    bundle_root,
    "results",
    "reproduced_openeval_ifeval_H3_G3_interpretation",
    "openeval_item_intercepts_loadings_metadata.csv"
  )
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(bundle_root, "results", "loadings_crossloadings")
)
threshold <- as.numeric(Sys.getenv("LOADING_THRESHOLD", unset = "0.50"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(loadings_path, check.names = FALSE, stringsAsFactors = FALSE)
lambda_cols <- paste0("loading_factor_", 1:3)
lambda <- as.matrix(d[, lambda_cols, drop = FALSE])
colnames(lambda) <- paste0("F", 1:3)

abs_lambda <- abs(lambda)
active <- abs_lambda >= threshold
lambda_sparse <- lambda
lambda_sparse[!active] <- 0
strongest <- max.col(abs_lambda, ties.method = "first")
active_count <- rowSums(active)

all_items <- cbind(
  d[, intersect(
    c("item_id", "empirical_accuracy", "benchmark", "primary_semantic", "question_snippet"),
    names(d)
  ), drop = FALSE],
  data.frame(lambda, check.names = FALSE),
  strongest_factor = paste0("F", strongest),
  max_abs_loading = apply(abs_lambda, 1L, max),
  n_active_factors = active_count,
  active_pattern = apply(active, 1L, function(x) paste0(ifelse(x, "1", "0"), collapse = "")),
  stringsAsFactors = FALSE
)
all_items <- all_items[order(all_items$strongest_factor, -all_items$max_abs_loading), ]
write.csv(all_items, file.path(out_dir, "ifeval_sparse_learned_lambda_all_items_ordered.csv"), row.names = FALSE)

sparse_items <- cbind(
  all_items[, intersect(
    c("item_id", "empirical_accuracy", "benchmark", "primary_semantic", "question_snippet"),
    names(all_items)
  ), drop = FALSE],
  data.frame(lambda_sparse[match(all_items$item_id, d$item_id), , drop = FALSE], check.names = FALSE),
  strongest_factor = all_items$strongest_factor,
  n_active_factors = all_items$n_active_factors,
  active_pattern = all_items$active_pattern,
  stringsAsFactors = FALSE
)
write.csv(
  sparse_items,
  file.path(out_dir, sprintf("ifeval_sparse_lambda_threshold_%s_ordered.csv", gsub("[.]", "p", format(threshold, trim = TRUE)))),
  row.names = FALSE
)

sparse_matrix <- data.frame(
  item_id = d$item_id,
  lambda_sparse,
  check.names = FALSE
)
write.csv(
  sparse_matrix,
  file.path(out_dir, sprintf("ifeval_sparse_lambda_matrix_threshold_%s.csv", gsub("[.]", "p", format(threshold, trim = TRUE)))),
  row.names = FALSE
)

cross <- all_items[all_items$n_active_factors >= 2L, ]
write.csv(cross, file.path(out_dir, "ifeval_G3_cross_loading_items.csv"), row.names = FALSE)

factor_only <- all_items[all_items$n_active_factors == 1L, ]
write.csv(factor_only, file.path(out_dir, "ifeval_G3_factor_only_exact_items.csv"), row.names = FALSE)

plot_loading_heatmap <- function(items, filename, title) {
  if (nrow(items) == 0L) return(invisible(NULL))
  mat <- as.matrix(items[, paste0("F", 1:3), drop = FALSE])
  rownames(mat) <- items$item_id
  max_abs <- max(abs(mat), na.rm = TRUE)
  pal <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)

  png(file.path(out_dir, filename), width = 1100, height = 1800, res = 170)
  op <- par(mar = c(4, 8, 4, 2))
  on.exit(par(op), add = TRUE)
  image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1L, , drop = FALSE]),
    col = pal,
    breaks = seq(-max_abs, max_abs, length.out = 102),
    axes = FALSE,
    xlab = "factor",
    ylab = "item",
    main = title
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
  axis(2, labels = FALSE)
  box()
  dev.off()
}

plot_loading_heatmap(
  all_items,
  "ifeval_sparse_learned_lambda_all_items_ordered_heatmap.png",
  "IFEval learned lambda, ordered by strongest factor"
)
plot_loading_heatmap(
  cross,
  "ifeval_G3_cross_loading_items_heatmap.png",
  paste0("IFEval cross-loading items; |lambda| >= ", threshold)
)

sparse_plot_items <- sparse_items
names(sparse_plot_items)[names(sparse_plot_items) %in% paste0("F", 1:3)] <- paste0("F", 1:3)
plot_loading_heatmap(
  sparse_plot_items,
  sprintf("ifeval_sparse_lambda_threshold_%s_heatmap.png", gsub("[.]", "p", format(threshold, trim = TRUE))),
  paste0("IFEval sparse reported Lambda; |lambda| >= ", threshold)
)

cat("Saved loading summaries in: ", normalizePath(out_dir), "\n", sep = "")
cat("Threshold: ", threshold, "\n", sep = "")
cat("Cross-loading items: ", nrow(cross), "\n", sep = "")
cat("Factor-only items: ", nrow(factor_only), "\n", sep = "")
cat("Sparse Lambda nonzero entries: ", sum(lambda_sparse != 0), " of ", length(lambda_sparse), "\n", sep = "")
