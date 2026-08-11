#!/usr/bin/env Rscript

# Base-R visualizations for the selected IFEval mixture-factor fit.  These are
# dependency-light fallbacks for the Python/Plotly visualizations.

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
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(bundle_root, "results", "ifeval_3d_factor_visualizations")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(file.path(fit_dir, "openeval_model_factor_scores_profiles.csv"), check.names = FALSE)
scores$profile_id <- as.character(scores$profile_id)

profile_levels <- sort(unique(scores$profile_id))
profile_cols <- grDevices::hcl.colors(length(profile_levels), palette = "Dark 3")
point_col <- profile_cols[match(scores$profile_id, profile_levels)]
point_cex <- 0.65 + 1.7 * (scores$accuracy - min(scores$accuracy)) /
  pmax(diff(range(scores$accuracy)), 1e-12)

project3 <- function(mat) {
  data.frame(
    x = mat[, 1L] + 0.42 * mat[, 3L],
    y = mat[, 2L] + 0.24 * mat[, 3L]
  )
}

F_mat <- as.matrix(scores[, paste0("factor_", 1:3)])
xy <- project3(F_mat)

png(file.path(out_dir, "ifeval_3d_factor_scatter_profiles.png"), width = 1400, height = 1000, res = 160)
op <- par(mar = c(5, 5, 4, 9))
plot(
  xy$x, xy$y,
  pch = 19,
  col = point_col,
  cex = point_cex,
  xlab = "projected F1/F3",
  ylab = "projected F2/F3",
  main = "IFEval learned mixture factors"
)
abline(h = 0, v = 0, col = "gray80", lty = 3)
legend("right", inset = c(-0.28, 0), xpd = TRUE, legend = profile_levels,
       col = profile_cols, pch = 19, title = "profile", cex = 0.65, bty = "n")
par(op)
dev.off()

png(file.path(out_dir, "ifeval_pairwise_factor_scatter_profiles.png"), width = 1700, height = 600, res = 160)
op <- par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))
for (pair in list(c(1, 2), c(1, 3), c(2, 3))) {
  plot(
    scores[[paste0("factor_", pair[1L])]],
    scores[[paste0("factor_", pair[2L])]],
    pch = 19,
    col = point_col,
    cex = point_cex,
    xlab = paste0("F", pair[1L]),
    ylab = paste0("F", pair[2L]),
    main = paste0("F", pair[1L], " vs F", pair[2L])
  )
  abline(h = 0, v = 0, col = "gray80", lty = 3)
}
par(op)
dev.off()

ordered <- scores[order(scores$accuracy, decreasing = TRUE), ]
factor_mat <- t(as.matrix(ordered[, paste0("factor_", 1:3)]))
colnames(factor_mat) <- ordered$model_id
rownames(factor_mat) <- paste0("F", 1:3)
max_abs <- max(abs(factor_mat), na.rm = TRUE)
pal <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)

png(file.path(out_dir, "ifeval_factor_score_heatmap_by_llm.png"), width = 3200, height = 850, res = 180)
op <- par(mar = c(12, 5, 4, 2))
image(seq_len(ncol(factor_mat)), seq_len(nrow(factor_mat)), t(factor_mat),
      col = pal, breaks = seq(-max_abs, max_abs, length.out = 102),
      axes = FALSE, xlab = "", ylab = "factor",
      main = "IFEval factor scores by LLM, ordered by accuracy")
axis(2, at = seq_len(3), labels = rownames(factor_mat), las = 1)
axis(1, at = seq_len(ncol(factor_mat)), labels = colnames(factor_mat), las = 2, cex.axis = 0.28, tick = FALSE)
box()
par(op)
dev.off()

group_mat <- t(as.matrix(ordered[, paste0("group_factor_", 1:3)]))
colnames(group_mat) <- ordered$model_id
rownames(group_mat) <- paste0("F", 1:3)
group_cols <- colorRampPalette(c("#3B6EA8", "#F2C14E", "#9E2F44"))(max(group_mat, na.rm = TRUE))

png(file.path(out_dir, "ifeval_map_group_heatmap_by_llm.png"), width = 3200, height = 850, res = 180)
op <- par(mar = c(12, 5, 4, 2))
image(seq_len(ncol(group_mat)), seq_len(nrow(group_mat)), t(group_mat),
      col = group_cols, breaks = seq(0.5, max(group_mat, na.rm = TRUE) + 0.5, by = 1),
      axes = FALSE, xlab = "", ylab = "factor",
      main = "IFEval MAP mixture groups by LLM, ordered by accuracy")
axis(2, at = seq_len(3), labels = rownames(group_mat), las = 1)
axis(1, at = seq_len(ncol(group_mat)), labels = colnames(group_mat), las = 2, cex.axis = 0.28, tick = FALSE)
legend("topright", fill = group_cols, legend = paste0("group ", seq_along(group_cols)), bty = "n")
box()
par(op)
dev.off()

profile_summary <- aggregate(accuracy ~ profile_id, scores, function(z) {
  c(n = length(z), mean = mean(z), min = min(z), max = max(z))
})
profile_summary <- data.frame(
  profile_id = profile_summary$profile_id,
  n_models = profile_summary$accuracy[, "n"],
  mean_accuracy = profile_summary$accuracy[, "mean"],
  min_accuracy = profile_summary$accuracy[, "min"],
  max_accuracy = profile_summary$accuracy[, "max"]
)
profile_summary <- profile_summary[order(profile_summary$mean_accuracy, decreasing = TRUE), ]
write.csv(profile_summary, file.path(out_dir, "ifeval_visual_profile_summary.csv"), row.names = FALSE)

cat("Saved base-R learned factor visualizations in: ", normalizePath(out_dir), "\n", sep = "")
