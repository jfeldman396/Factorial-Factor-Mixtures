#!/usr/bin/env Rscript

# Side-by-side loading-matrix comparison:
#   independent-mixture probit Lambda vs ordinary probit Lambda.
#
# The ordinary probit factor labels are arbitrary.  We align them to the mixture
# factors by choosing the permutation/signs that maximize factor-score
# correlations, then apply that same permutation/sign convention to Lambda.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  unset = file.path(bundle_root, "data", "openeval_ifeval_only_binary_matrix.csv")
)
mix_path <- Sys.getenv(
  "MIXTURE_LOADINGS_PATH",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ifeval_H3_G3_interpretation", "openeval_item_intercepts_loadings_metadata.csv")
)
ordinary_path <- Sys.getenv(
  "ORDINARY_LAMBDA_PATH",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ifeval_ordinary_probit_H3_visualization", "ordinary_probit_lambda.csv")
)
cor_path <- Sys.getenv(
  "ORDINARY_MIXTURE_COR_PATH",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ordinary_vs_mixture_H3", "ordinary_mixture_factor_correlation_matrix.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ordinary_vs_mixture_H3")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_matrix_item_ids <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, nrows = 1)
  names(raw)[-1L]
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  out <- do.call(rbind, lapply(seq_along(x), function(i) {
    rest <- x[-i]
    cbind(x[i], all_permutations(rest))
  }))
  out
}

align_ordinary_to_mixture <- function(L_ord, factor_cor) {
  H <- ncol(factor_cor)
  perms <- all_permutations(seq_len(H))
  scores <- apply(perms, 1L, function(perm) {
    sum(abs(factor_cor[cbind(perm, seq_len(H))]))
  })
  perm <- perms[which.max(scores), ]
  signs <- sign(factor_cor[cbind(perm, seq_len(H))])
  signs[signs == 0] <- 1
  L_aligned <- sweep(L_ord[, perm, drop = FALSE], 2L, signs, "*")
  colnames(L_aligned) <- paste0("F", seq_len(H))
  list(Lambda = L_aligned, permutation = perm, signs = signs)
}

plot_three_panel <- function(L_mix, L_ord, meta, ord, filename, title, side_labels = NULL) {
  L_mix <- L_mix[ord, , drop = FALSE]
  L_ord <- L_ord[ord, , drop = FALSE]
  L_diff <- L_mix - L_ord

  png(file.path(out_dir, filename), width = 2300, height = 1500, res = 170)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 1))
  op <- par(mar = c(4, 5, 5, 1), oma = c(0, 0, 3, 0))
  on.exit(par(op), add = TRUE)

  max_abs <- quantile(abs(c(L_mix, L_ord)), 0.995, na.rm = TRUE)
  max_diff <- quantile(abs(L_diff), 0.995, na.rm = TRUE)
  pal <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)

  draw <- function(mat, main, lim, show_y = FALSE, labels = NULL) {
    image(
      x = seq_len(ncol(mat)),
      y = seq_len(nrow(mat)),
      z = t(mat[nrow(mat):1L, , drop = FALSE]),
      col = pal,
      breaks = seq(-lim, lim, length.out = 102),
      axes = FALSE,
      xlab = "factor",
      ylab = if (show_y) "items" else "",
      main = main
    )
    axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
    if (show_y) axis(2, labels = FALSE)
    if (!is.null(labels)) {
      y_pos <- nrow(mat) - labels$mid + 1L
      axis(4, at = y_pos, labels = labels$label, las = 1, cex.axis = 0.60, tick = FALSE)
      abline(h = nrow(mat) - labels$end[-length(labels$end)] + 0.5, col = "gray45", lty = 3)
    }
    box()
  }

  draw(L_mix, "mixture probit Lambda", max_abs, show_y = TRUE)
  draw(L_ord, "ordinary probit Lambda, aligned", max_abs)
  draw(L_diff, "mixture - ordinary", max_diff, labels = side_labels)
  mtext(title, outer = TRUE, cex = 1.1, font = 2)
  dev.off()
}

item_ids <- read_matrix_item_ids(matrix_path)

mix <- read.csv(mix_path, check.names = FALSE, stringsAsFactors = FALSE)
L_mix <- as.matrix(mix[, paste0("loading_factor_", 1:3)])
rownames(L_mix) <- mix$item_id
colnames(L_mix) <- paste0("F", 1:3)

ord_raw <- read.csv(ordinary_path, check.names = FALSE, stringsAsFactors = FALSE)
L_ord <- as.matrix(ord_raw[, paste0("factor_", 1:3)])
rownames(L_ord) <- item_ids
colnames(L_ord) <- paste0("ordinary_F", 1:3)

common_items <- intersect(rownames(L_mix), rownames(L_ord))
L_mix <- L_mix[common_items, , drop = FALSE]
L_ord <- L_ord[common_items, , drop = FALSE]
meta <- mix[match(common_items, mix$item_id), , drop = FALSE]

factor_cor <- as.matrix(read.csv(cor_path, row.names = 1L, check.names = FALSE))
aligned <- align_ordinary_to_mixture(L_ord, factor_cor)
L_ord_aligned <- aligned$Lambda

alignment <- data.frame(
  mixture_factor = paste0("F", 1:3),
  ordinary_factor_used = paste0("ordinary_F", aligned$permutation),
  sign_applied = aligned$signs,
  score_correlation = factor_cor[cbind(aligned$permutation, seq_len(3))]
)
write.csv(alignment, file.path(out_dir, "openeval_lambda_alignment.csv"), row.names = FALSE)

comparison_summary <- data.frame(
  factor = paste0("F", 1:3),
  cor_lambda_mixture_ordinary_aligned = diag(cor(L_mix, L_ord_aligned)),
  rmse_lambda = sqrt(colMeans((L_mix - L_ord_aligned)^2)),
  mean_abs_mixture = colMeans(abs(L_mix)),
  mean_abs_ordinary_aligned = colMeans(abs(L_ord_aligned))
)
write.csv(comparison_summary, file.path(out_dir, "openeval_lambda_side_by_side_summary.csv"), row.names = FALSE)

bench_ord <- order(meta$benchmark, meta$item_id)
bench_split <- split(seq_along(bench_ord), meta$benchmark[bench_ord])
bench_labels <- do.call(rbind, lapply(names(bench_split), function(b) {
  idx <- bench_split[[b]]
  data.frame(label = b, start = min(idx), end = max(idx), mid = round(mean(range(idx))))
}))

plot_three_panel(
  L_mix,
  L_ord_aligned,
  meta,
  bench_ord,
  "openeval_lambda_side_by_side_ordered_by_benchmark.png",
  "OpenEval fitted loading matrices, ordered by benchmark",
  bench_labels
)

strongest <- apply(abs(L_mix), 1L, which.max)
max_abs <- apply(abs(L_mix), 1L, max)
strength_ord <- order(strongest, -max_abs)
strength_split <- split(seq_along(strength_ord), strongest[strength_ord])
strength_labels <- do.call(rbind, lapply(names(strength_split), function(h) {
  idx <- strength_split[[h]]
  data.frame(
    label = paste0("strongest F", h),
    start = min(idx),
    end = max(idx),
    mid = round(mean(range(idx)))
  )
}))

plot_three_panel(
  L_mix,
  L_ord_aligned,
  meta,
  strength_ord,
  "openeval_lambda_side_by_side_ordered_by_mixture_strength.png",
  "OpenEval fitted loading matrices, ordered by strongest mixture loading",
  strength_labels
)

cat("\nLambda alignment:\n")
print(alignment)
cat("\nLambda comparison summary:\n")
print(comparison_summary)
cat("\nOutputs saved in: ", normalizePath(out_dir), "\n", sep = "")
