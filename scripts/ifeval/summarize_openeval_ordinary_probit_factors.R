#!/usr/bin/env Rscript

# OpenEval ordinary binary probit factor summary and comparison to the selected
# independent-mixture probit fit.  The script detects the fitted rank from the
# saved factor-score columns, so it works for H=3, H=4, or another selected H.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- Sys.getenv(
  "IFEVAL_BUNDLE_ROOT",
  unset = normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
)
ordinary_dir <- Sys.getenv(
  "ORDINARY_DIR",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ifeval_ordinary_probit_H3_visualization")
)
mixture_dir <- Sys.getenv(
  "MIXTURE_DIR",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ifeval_H3_G3_interpretation")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ordinary_vs_mixture_H3")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ordinary <- read.csv(file.path(ordinary_dir, "ordinary_probit_factor_scores.csv"), check.names = FALSE)
names(ordinary)[1L] <- "model_id"
ordinary_factor_cols <- grep("^factor_[0-9]+$", names(ordinary), value = TRUE)
if (!length(ordinary_factor_cols)) {
  ordinary_factor_cols <- grep("^[0-9]+$", names(ordinary), value = TRUE)
}
if (!length(ordinary_factor_cols)) {
  stop("No ordinary factor-score columns found in: ", ordinary_dir)
}
H_ordinary <- length(ordinary_factor_cols)
names(ordinary)[match(ordinary_factor_cols, names(ordinary))] <- paste0("ordinary_F", seq_len(H_ordinary))

mixture <- read.csv(file.path(mixture_dir, "openeval_model_factor_scores_profiles.csv"), check.names = FALSE)
mixture_factor_cols <- grep("^factor_[0-9]+$", names(mixture), value = TRUE)
mixture_group_cols <- grep("^group_factor_[0-9]+$", names(mixture), value = TRUE)
H_mixture <- length(mixture_factor_cols)
H_compare <- min(H_ordinary, H_mixture)
if (H_compare < 1L) stop("No common factor dimensions found.")
mixture <- mixture[, c(
  "model_id",
  "accuracy",
  mixture_factor_cols,
  mixture_group_cols,
  "profile_id"
)]
names(mixture)[match(mixture_factor_cols, names(mixture))] <- paste0("mixture_F", seq_len(H_mixture))

dat <- merge(ordinary, mixture, by = c("model_id", "accuracy"), all = FALSE)
dat <- dat[order(dat$accuracy, decreasing = TRUE), ]
write.csv(dat, file.path(out_dir, "openeval_ordinary_and_mixture_factor_scores.csv"), row.names = FALSE)

ordinary_mat <- as.matrix(dat[, paste0("ordinary_F", seq_len(H_compare)), drop = FALSE])
mixture_mat <- as.matrix(dat[, paste0("mixture_F", seq_len(H_compare)), drop = FALSE])

acc_cor <- data.frame(
  factor = paste0("ordinary_F", seq_len(H_compare)),
  cor_with_accuracy = as.numeric(cor(ordinary_mat, dat$accuracy))
)
write.csv(acc_cor, file.path(out_dir, "ordinary_factor_accuracy_correlations.csv"), row.names = FALSE)

factor_cor <- cor(ordinary_mat, mixture_mat)
write.csv(factor_cor, file.path(out_dir, "ordinary_mixture_factor_correlation_matrix.csv"))

plot_pairs <- combn(seq_len(H_compare), 2L, simplify = FALSE)
png(file.path(out_dir, "ordinary_probit_pairwise_factor_scatter.png"), width = 1900, height = 1250, res = 160)
op <- par(mfrow = c(ceiling(length(plot_pairs) / 3), 3), mar = c(5, 5, 3, 1))
pal <- colorRampPalette(c("#355C9A", "#F2C14E", "#B23A48"))(100)
idx <- pmax(1, pmin(100, as.integer(cut(dat$accuracy, breaks = 100, labels = FALSE))))
for (pair in plot_pairs) {
  x <- ordinary_mat[, pair[1]]
  y <- ordinary_mat[, pair[2]]
  plot(
    x, y,
    pch = 19,
    col = pal[idx],
    xlab = paste0("ordinary F", pair[1]),
    ylab = paste0("ordinary F", pair[2]),
    main = paste0("ordinary F", pair[1], " vs F", pair[2])
  )
  abline(h = 0, v = 0, col = "gray70", lty = 3)
}
par(op)
dev.off()

png(file.path(out_dir, "ordinary_factor_scores_by_llm_heatmap.png"), width = 3200, height = 850, res = 180)
op <- par(mar = c(12, 5, 4, 2))
mat <- t(ordinary_mat)
colnames(mat) <- dat$model_id
rownames(mat) <- paste0("F", seq_len(H_compare))
max_abs <- max(abs(mat), na.rm = TRUE)
pal2 <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)
image(
  seq_len(ncol(mat)),
  seq_len(nrow(mat)),
  t(mat),
  col = pal2,
  breaks = seq(-max_abs, max_abs, length.out = 102),
  axes = FALSE,
  xlab = "",
  ylab = "ordinary probit factor",
  main = "Ordinary probit factor scores by LLM, ordered by accuracy"
)
axis(2, at = seq_len(H_compare), labels = rownames(mat), las = 1)
axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.28, tick = FALSE)
box()
par(op)
dev.off()

png(file.path(out_dir, "ordinary_vs_mixture_factor_correlation_heatmap.png"), width = 900, height = 800, res = 160)
op <- par(mar = c(5, 5, 4, 2))
pal3 <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)
image(
  seq_len(H_compare),
  seq_len(H_compare),
  t(factor_cor),
  col = pal3,
  breaks = seq(-1, 1, length.out = 102),
  axes = FALSE,
  xlab = "ordinary probit factor",
  ylab = "mixture factor",
  main = "Correlation: ordinary vs mixture factors"
)
axis(1, at = seq_len(H_compare), labels = paste0("Ord F", seq_len(H_compare)))
axis(2, at = seq_len(H_compare), labels = paste0("Mix F", seq_len(H_compare)), las = 1)
for (i in seq_len(H_compare)) {
  for (j in seq_len(H_compare)) {
    text(i, j, labels = sprintf("%.2f", factor_cor[i, j]), cex = 1.2)
  }
}
box()
par(op)
dev.off()

cat("\nOrdinary factor correlations with accuracy:\n")
print(acc_cor)
cat("\nOrdinary vs mixture factor correlation matrix:\n")
print(round(factor_cor, 3))
cat("\nTop and bottom models by accuracy with ordinary factors:\n")
print(head(dat[, c("model_id", "accuracy", paste0("ordinary_F", seq_len(H_compare)))], 12), row.names = FALSE)
cat("\n")
print(tail(dat[, c("model_id", "accuracy", paste0("ordinary_F", seq_len(H_compare)))], 12), row.names = FALSE)
cat("\nOutputs saved in: ", normalizePath(out_dir), "\n", sep = "")
