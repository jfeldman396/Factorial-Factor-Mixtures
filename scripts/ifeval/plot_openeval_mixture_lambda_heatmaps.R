#!/usr/bin/env Rscript

# Visualize the fitted OpenEval independent-mixture probit loading matrix.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

fit_dir <- Sys.getenv(
  "FIT_DIR",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ifeval_H3_G3_interpretation")
)
load_path <- Sys.getenv(
  "LOAD_PATH",
  unset = file.path(fit_dir, "openeval_item_intercepts_loadings_metadata.csv")
)
out_dir <- Sys.getenv("OUT_DIR", unset = fit_dir)

d <- read.csv(load_path, check.names = FALSE, stringsAsFactors = FALSE)
lambda_cols <- grep("^loading_factor_[0-9]+$", names(d), value = TRUE)
if (!length(lambda_cols)) stop("No loading_factor_* columns found in: ", load_path)
L <- as.matrix(d[, lambda_cols])
rownames(L) <- d$item_id
colnames(L) <- paste0("F", seq_len(ncol(L)))

plot_lambda_heatmap <- function(L, main, filename, row_labels = NULL) {
  png(file.path(out_dir, filename), width = 1450, height = 1800, res = 170)
  op <- par(mar = c(4, 8, 4, 7))
  on.exit(par(op), add = TRUE)
  max_abs <- quantile(abs(L), 0.995, na.rm = TRUE)
  pal <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)
  image(
    x = seq_len(ncol(L)),
    y = seq_len(nrow(L)),
    z = t(L[nrow(L):1L, , drop = FALSE]),
    col = pal,
    breaks = seq(-max_abs, max_abs, length.out = 102),
    axes = FALSE,
    xlab = "factor",
    ylab = "items",
    main = main
  )
  axis(1, at = seq_len(ncol(L)), labels = colnames(L))
  axis(2, labels = FALSE)
  box()
  if (!is.null(row_labels)) {
    y_pos <- nrow(L) - row_labels$mid + 1L
    axis(4, at = y_pos, labels = row_labels$label, las = 1, cex.axis = 0.65, tick = FALSE)
    abline(h = nrow(L) - row_labels$end[-length(row_labels$end)] + 0.5, col = "gray45", lty = 3)
  }
  dev.off()
}

bench_order <- order(d$benchmark, d$item_id)
bench_tab <- split(seq_along(bench_order), d$benchmark[bench_order])
bench_labels <- do.call(rbind, lapply(names(bench_tab), function(b) {
  idx <- bench_tab[[b]]
  data.frame(
    label = b,
    start = min(idx),
    end = max(idx),
    mid = round(mean(range(idx)))
  )
}))

plot_lambda_heatmap(
  L[bench_order, , drop = FALSE],
  "OpenEval fitted Lambda, ordered by benchmark",
  "openeval_lambda_full_ordered_by_benchmark.png",
  bench_labels
)

strongest <- apply(abs(L), 1L, which.max)
max_abs <- apply(abs(L), 1L, max)
ord_strength <- order(strongest, -max_abs)
strength_tab <- split(seq_along(ord_strength), strongest[ord_strength])
strength_labels <- do.call(rbind, lapply(names(strength_tab), function(h) {
  idx <- strength_tab[[h]]
  data.frame(
    label = paste0("F", h, " strongest"),
    start = min(idx),
    end = max(idx),
    mid = round(mean(range(idx)))
  )
}))

plot_lambda_heatmap(
  L[ord_strength, , drop = FALSE],
  "OpenEval fitted Lambda, ordered by strongest factor",
  "openeval_lambda_full_ordered_by_strongest_factor.png",
  strength_labels
)

summary <- data.frame(
  factor = colnames(L),
  min = apply(L, 2L, min),
  q05 = apply(L, 2L, quantile, 0.05),
  median = apply(L, 2L, median),
  q95 = apply(L, 2L, quantile, 0.95),
  max = apply(L, 2L, max),
  mean_abs = colMeans(abs(L)),
  n_abs_gt_0_1 = colSums(abs(L) > 0.1),
  n_abs_gt_0_5 = colSums(abs(L) > 0.5),
  n_abs_gt_1_0 = colSums(abs(L) > 1.0)
)
write.csv(summary, file.path(out_dir, "openeval_lambda_numeric_summary.csv"), row.names = FALSE)

active <- data.frame(
  benchmark = d$benchmark,
  strongest_factor = paste0("F", strongest),
  max_abs = max_abs
)
active <- active[active$max_abs > 0.1, ]
tab <- as.data.frame.matrix(table(active$benchmark, active$strongest_factor))
tab$benchmark <- rownames(tab)
tab <- tab[, c("benchmark", setdiff(names(tab), "benchmark"))]
write.csv(tab, file.path(out_dir, "openeval_lambda_strongest_factor_by_benchmark.csv"), row.names = FALSE)

cat("\nLambda numeric summary:\n")
print(summary)
cat("\nActive item counts by benchmark and strongest factor:\n")
print(tab)
