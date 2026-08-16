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

lambda_palette <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(201)

lambda_raster <- function(L, limit = NULL) {
  if (is.null(limit)) {
    limit <- as.numeric(quantile(abs(L), 0.995, na.rm = TRUE))
  }
  if (!is.finite(limit) || limit <= 0) limit <- 1
  idx <- floor((as.vector(L) + limit) / (2 * limit) * (length(lambda_palette) - 1)) + 1L
  idx <- pmax(1L, pmin(length(lambda_palette), idx))
  as.raster(matrix(lambda_palette[idx], nrow = nrow(L), ncol = ncol(L)))
}

plot_lambda_heatmap <- function(L, main, filename, row_labels = NULL,
                                width = 1450, height = 1800, cex_axis = 0.65) {
  png(file.path(out_dir, filename), width = width, height = height, res = 170)
  op <- par(mar = c(4, 8, 4, 7))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  plot(
    NA,
    xlim = c(0.5, ncol(L) + 0.5),
    ylim = c(nrow(L) + 0.5, 0.5),
    axes = FALSE,
    xlab = "factor",
    ylab = "items",
    main = main
  )
  rasterImage(lambda_raster(L), 0.5, nrow(L) + 0.5, ncol(L) + 0.5, 0.5,
              interpolate = FALSE)
  axis(1, at = seq_len(ncol(L)), labels = colnames(L))
  axis(2, labels = FALSE)
  box()
  if (!is.null(row_labels)) {
    axis(4, at = row_labels$mid, labels = row_labels$label, las = 1,
         cex.axis = cex_axis, tick = FALSE)
    abline(h = row_labels$end[-length(row_labels$end)] + 0.5, col = "gray45", lty = 3)
  }
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

if (ncol(L) >= 2L && any(strongest == 2L)) {
  f2_idx <- which(strongest == 2L)
  f2_idx <- f2_idx[order(-max_abs[f2_idx])]
  L_f2 <- L[f2_idx, , drop = FALSE]
  f2_labels <- paste0(
    sub(".*Z_", "", rownames(L_f2)),
    " max=",
    sprintf("%.2f", max_abs[f2_idx])
  )

  png(file.path(out_dir, "openeval_lambda_F2_block_zoom.png"),
      width = 1600, height = 1900, res = 220)
  op <- par(mar = c(5, 6, 4, 8))
  plot(
    NA,
    xlim = c(0.5, ncol(L_f2) + 0.5),
    ylim = c(nrow(L_f2) + 0.5, 0.5),
    axes = FALSE,
    xlab = "factor",
    ylab = "F2-strongest items",
    main = "F2 strongest block, ordered by max |loading|"
  )
  rasterImage(lambda_raster(L_f2), 0.5, nrow(L_f2) + 0.5, ncol(L_f2) + 0.5, 0.5,
              interpolate = FALSE)
  axis(1, at = seq_len(ncol(L_f2)), labels = colnames(L_f2))
  axis(2, labels = FALSE)
  axis(4, at = seq_along(f2_labels), labels = f2_labels, las = 1,
       cex.axis = 0.38, tick = FALSE)
  box()
  par(op)
  dev.off()

  clean_idx <- head(f2_idx, 4L)
  weak_idx <- tail(f2_idx, 4L)
  example_idx <- c(clean_idx, weak_idx)
  L_examples <- L[example_idx, , drop = FALSE]
  rownames(L_examples) <- c(
    paste0("clean ", sub(".*Z_", "", rownames(L)[clean_idx])),
    paste0("weak ", sub(".*Z_", "", rownames(L)[weak_idx]))
  )
  x_limit <- max(abs(L_examples), na.rm = TRUE)
  x_limit <- if (is.finite(x_limit)) max(0.5, x_limit) else 1
  x_lim <- c(min(-0.5, min(L_examples, na.rm = TRUE) - 0.1), x_limit + 0.2)
  dot_cols <- c("#B23A48", "#355C9A", "#2E8B57", "#6A3D9A")
  dot_cols <- rep(dot_cols, length.out = ncol(L_examples))

  png(file.path(out_dir, "openeval_lambda_F2_clean_vs_weak_dotplot.png"),
      width = 1550, height = 950, res = 180)
  op <- par(mar = c(5, 12, 4, 2))
  y_pos <- nrow(L_examples):1L
  plot(
    NA,
    xlim = x_lim,
    ylim = c(0.5, nrow(L_examples) + 0.5),
    yaxt = "n",
    xlab = "loading",
    ylab = "",
    main = "F2 examples: exact loading values",
    bty = "n"
  )
  abline(v = 0, col = "gray75")
  for (j in seq_len(ncol(L_examples))) {
    points(L_examples[, j], y_pos, pch = 19, col = dot_cols[j], cex = 1.1)
    text(L_examples[, j], y_pos + 0.18, labels = sprintf("%.2f", L_examples[, j]),
         col = dot_cols[j], cex = 0.55)
  }
  axis(2, at = y_pos, labels = rownames(L_examples), las = 1, cex.axis = 0.65)
  legend("bottomright", legend = colnames(L_examples), col = dot_cols,
         pch = 19, bty = "n", horiz = TRUE)
  par(op)
  dev.off()
}

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
