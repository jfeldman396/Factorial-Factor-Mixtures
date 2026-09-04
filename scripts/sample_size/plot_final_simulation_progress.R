#!/usr/bin/env Rscript

# Progress dashboard for the final product/Viroli simulation.
#
# This script is intentionally read-only with respect to the simulation
# checkpoint.  It can be run while the long simulation is active.

options(stringsAsFactors = FALSE)

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

checkpoint <- get_env(
  "CHECKPOINT",
  file.path(repo_root, "results", "full", "final_cross_ifeval_product_viroli", "comparison_results_checkpoint.csv")
)
out_dir <- get_env("OUT_DIR", dirname(checkpoint))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(checkpoint)) {
  stop("Checkpoint not found: ", checkpoint, call. = FALSE)
}

d <- read.csv(checkpoint, stringsAsFactors = FALSE)
if (!nrow(d)) stop("Checkpoint has zero rows: ", checkpoint, call. = FALSE)

method_labels <- c(
  independent_marginal_mixture = "product MAP",
  viroli_laplace_gibbs = "Viroli Laplace",
  viroli_gaussian_gibbs = "Viroli Gaussian"
)
method_colors <- c(
  independent_marginal_mixture = "#2f6db5",
  viroli_laplace_gibbs = "#c92d32",
  viroli_gaussian_gibbs = "#37a567"
)

if (!"method" %in% names(d)) stop("Checkpoint is missing the method column.", call. = FALSE)
d$method_label <- ifelse(
  d$method %in% names(method_labels),
  method_labels[d$method],
  d$method
)

short_G <- function(x) {
  vapply(strsplit(as.character(x), "-", fixed = TRUE), function(z) {
    z <- z[nzchar(z)]
    if (length(unique(z)) == 1L) unique(z) else paste(z, collapse = "-")
  }, character(1))
}

d$G_short <- short_G(d$G_true)
d$cell <- paste0(
  "n=", d$n,
  ", p=", d$p,
  ", H=", d$H_true,
  ", G=", d$G_short,
  ", sep=", d$separation
)

metric_cols <- c(
  factor_score_rmse = "factor scores",
  lambda_rmse = "loadings",
  alpha_rmse = "intercepts",
  marginal_mu_rmse = "mixture means",
  marginal_var_rmse = "mixture variances",
  marginal_weight_rmse = "mixture weights",
  seconds = "runtime"
)
metric_cols <- metric_cols[names(metric_cols) %in% names(d)]

group_cols <- c("method", "method_label", "n", "p", "H_true", "G_true", "G_short", "separation")
key <- interaction(d[, group_cols, drop = FALSE], drop = TRUE, sep = " | ")

summary_list <- lapply(split(d, key), function(dd) {
  base <- dd[1L, group_cols, drop = FALSE]
  out <- data.frame(n_reps = nrow(dd))
  for (metric in names(metric_cols)) {
    x <- dd[[metric]]
    x <- x[is.finite(x)]
    out[[paste0(metric, "_mean")]] <- if (length(x)) mean(x) else NA_real_
    out[[paste0(metric, "_sd")]] <- if (length(x) > 1L) sd(x) else NA_real_
  }
  cbind(base, out)
})
summary_df <- do.call(rbind, summary_list)
summary_df <- summary_df[order(summary_df$n, summary_df$p, summary_df$separation, summary_df$G_short, summary_df$H_true, summary_df$method), ]

summary_file <- file.path(out_dir, "progress_summary_current.csv")
write.csv(summary_df, summary_file, row.names = FALSE)

latest_table <- summary_df[, c(
  "method_label", "n", "p", "H_true", "G_short", "separation", "n_reps",
  intersect(paste0(names(metric_cols), "_mean"), names(summary_df))
), drop = FALSE]
latest_table_file <- file.path(out_dir, "progress_latest_table.csv")
write.csv(latest_table, latest_table_file, row.names = FALSE)

# Plot the earliest n,p cell with the most completed rows.  This makes the
# dashboard useful early in the run and naturally becomes a comparison plot once
# Gibbs baselines populate the checkpoint.
np_key <- interaction(summary_df$n, summary_df$p, summary_df$G_short, summary_df$separation, drop = TRUE)
np_counts <- sort(tapply(summary_df$n_reps, np_key, sum), decreasing = TRUE)
chosen_key <- names(np_counts)[1L]
plot_df <- summary_df[np_key == chosen_key, , drop = FALSE]

plot_file <- file.path(out_dir, "progress_parameter_rmse_runtime_panel.png")
png(plot_file, width = 1800, height = 1100, res = 150)
op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)
par(mfrow = c(2, 4), mar = c(4.2, 4.5, 3.2, 1), oma = c(0, 0, 4.2, 0), xaxs = "i")

plot_metric <- function(metric, title, ylab) {
  mean_col <- paste0(metric, "_mean")
  sd_col <- paste0(metric, "_sd")
  if (!mean_col %in% names(plot_df)) {
    plot.new()
    title(title)
    return(invisible(NULL))
  }
  vals <- plot_df[[mean_col]]
  sds <- if (sd_col %in% names(plot_df)) plot_df[[sd_col]] else rep(NA_real_, nrow(plot_df))
  ymax <- max(vals + ifelse(is.finite(sds), 2 * sds, 0), na.rm = TRUE)
  ymin <- min(vals - ifelse(is.finite(sds), 2 * sds, 0), na.rm = TRUE)
  if (!is.finite(ymax)) ymax <- 1
  if (!is.finite(ymin)) ymin <- 0
  if (metric != "seconds") ymin <- min(0, ymin)
  ylim <- range(c(ymin, ymax))
  if (diff(ylim) <= 0) ylim <- ylim + c(-0.1, 0.1)

  plot(
    NA,
    xlim = range(plot_df$H_true),
    ylim = ylim,
    xlab = "H",
    ylab = ylab,
    main = title,
    xaxt = "n"
  )
  axis(1, at = sort(unique(plot_df$H_true)))
  grid(col = "#eeeeee")

  for (method in unique(plot_df$method)) {
    dd <- plot_df[plot_df$method == method, , drop = FALSE]
    dd <- dd[order(dd$H_true), , drop = FALSE]
    col <- method_colors[[method]]
    if (is.null(col)) col <- "#555555"
    lines(dd$H_true, dd[[mean_col]], type = "b", pch = 19, lwd = 2, col = col)
    if (sd_col %in% names(dd)) {
      err <- 2 * dd[[sd_col]]
      ok <- is.finite(err)
      arrows(
        dd$H_true[ok], dd[[mean_col]][ok] - err[ok],
        dd$H_true[ok], dd[[mean_col]][ok] + err[ok],
        angle = 90, code = 3, length = 0.04, col = col
      )
    }
  }
}

for (metric in names(metric_cols)) {
  ylab <- if (metric == "seconds") "seconds" else "RMSE"
  plot_metric(metric, metric_cols[[metric]], ylab)
}
if (length(metric_cols) < 8L) {
  for (i in seq_len(8L - length(metric_cols))) plot.new()
}

legend_methods <- unique(plot_df$method)
legend_labels <- ifelse(
  legend_methods %in% names(method_labels),
  method_labels[legend_methods],
  legend_methods
)
legend_cols <- method_colors[legend_methods]
legend_cols[is.na(legend_cols)] <- "#555555"
mtext(
  paste0(
    "Progress dashboard: ",
    unique(plot_df$n), " systems, ",
    unique(plot_df$p), " items, G=", unique(plot_df$G_short),
    ", sep=", unique(plot_df$separation),
    " | checkpoint rows=", nrow(d)
  ),
  outer = TRUE,
  cex = 1.25,
  font = 2
)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend(
  "bottom",
  legend = legend_labels,
  col = legend_cols,
  lwd = 3,
  pch = 19,
  horiz = TRUE,
  bty = "n",
  cex = 1.05
)
dev.off()

cat("Rows:", nrow(d), "\n")
cat("Summary:", summary_file, "\n")
cat("Table:", latest_table_file, "\n")
cat("Plot:", plot_file, "\n")
print(latest_table, row.names = FALSE)
