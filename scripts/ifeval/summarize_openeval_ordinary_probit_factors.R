#!/usr/bin/env Rscript

# OpenEval ordinary binary probit H=3 factor summary and comparison to the
# selected independent-mixture H=3, G=3 fit.

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
names(ordinary)[3:5] <- paste0("ordinary_F", 1:3)

mixture <- read.csv(file.path(mixture_dir, "openeval_model_factor_scores_profiles.csv"), check.names = FALSE)
mixture <- mixture[, c("model_id", "accuracy", paste0("factor_", 1:3), paste0("group_factor_", 1:3), "profile_id")]
names(mixture)[3:5] <- paste0("mixture_F", 1:3)

dat <- merge(ordinary, mixture, by = c("model_id", "accuracy"), all = FALSE)
dat <- dat[order(dat$accuracy, decreasing = TRUE), ]
write.csv(dat, file.path(out_dir, "openeval_ordinary_and_mixture_factor_scores.csv"), row.names = FALSE)

ordinary_mat <- as.matrix(dat[, paste0("ordinary_F", 1:3)])
mixture_mat <- as.matrix(dat[, paste0("mixture_F", 1:3)])

acc_cor <- data.frame(
  factor = paste0("ordinary_F", 1:3),
  cor_with_accuracy = as.numeric(cor(ordinary_mat, dat$accuracy))
)
write.csv(acc_cor, file.path(out_dir, "ordinary_factor_accuracy_correlations.csv"), row.names = FALSE)

factor_cor <- cor(ordinary_mat, mixture_mat)
write.csv(factor_cor, file.path(out_dir, "ordinary_mixture_factor_correlation_matrix.csv"))

png(file.path(out_dir, "ordinary_probit_pairwise_factor_scatter.png"), width = 1900, height = 650, res = 160)
op <- par(mfrow = c(1, 3), mar = c(5, 5, 3, 1))
pal <- colorRampPalette(c("#355C9A", "#F2C14E", "#B23A48"))(100)
idx <- pmax(1, pmin(100, as.integer(cut(dat$accuracy, breaks = 100, labels = FALSE))))
pairs_to_plot <- list(c(1, 2), c(1, 3), c(2, 3))
for (pair in pairs_to_plot) {
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
rownames(mat) <- paste0("F", 1:3)
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
axis(2, at = seq_len(3), labels = rownames(mat), las = 1)
axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.28, tick = FALSE)
box()
par(op)
dev.off()

png(file.path(out_dir, "ordinary_vs_mixture_factor_correlation_heatmap.png"), width = 900, height = 800, res = 160)
op <- par(mar = c(5, 5, 4, 2))
pal3 <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)
image(
  seq_len(3),
  seq_len(3),
  t(factor_cor),
  col = pal3,
  breaks = seq(-1, 1, length.out = 102),
  axes = FALSE,
  xlab = "ordinary probit factor",
  ylab = "mixture factor",
  main = "Correlation: ordinary vs mixture factors"
)
axis(1, at = seq_len(3), labels = paste0("Ord F", 1:3))
axis(2, at = seq_len(3), labels = paste0("Mix F", 1:3), las = 1)
for (i in 1:3) {
  for (j in 1:3) {
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
print(head(dat[, c("model_id", "accuracy", paste0("ordinary_F", 1:3))], 12), row.names = FALSE)
cat("\n")
print(tail(dat[, c("model_id", "accuracy", paste0("ordinary_F", 1:3))], 12), row.names = FALSE)
cat("\nOutputs saved in: ", normalizePath(out_dir), "\n", sep = "")
