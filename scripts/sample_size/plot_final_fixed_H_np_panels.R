#!/usr/bin/env Rscript

# Fixed-H progress panels for the final Product MAP / Viroli simulation.
#
# The long simulation appends one row per method/replication/scenario to a
# checkpoint CSV.  This script is read-only with respect to that checkpoint and
# can be run while the simulation is active.  It makes one six-panel figure for
# each completed or partially completed fixed-H, fixed-G, fixed-separation
# setting, with n on the x-axis, method encoded by color, and p encoded by line
# type.

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
  file.path(
    repo_root,
    "results",
    "full",
    "final_cross_ifeval_product_viroli",
    "comparison_results_checkpoint.csv"
  )
)
out_dir <- get_env("OUT_DIR", dirname(checkpoint))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(checkpoint)) {
  stop("Checkpoint not found: ", checkpoint, call. = FALSE)
}

d <- read.csv(checkpoint, stringsAsFactors = FALSE)
if (!nrow(d)) stop("Checkpoint has zero rows: ", checkpoint, call. = FALSE)

required_cols <- c("method", "n", "p", "H_true", "G_true", "separation")
missing_cols <- setdiff(required_cols, names(d))
if (length(missing_cols)) {
  stop("Checkpoint is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

if (!"G_config" %in% names(d)) d$G_config <- d$G_true
if (!"loading_design" %in% names(d)) d$loading_design <- "unknown_loading"
if (!"block_size_mode" %in% names(d)) d$block_size_mode <- "balanced"

method_order <- c(
  "independent_marginal_mixture",
  "viroli_laplace_gibbs",
  "viroli_gaussian_gibbs",
  "joint_mixture_factor_gibbs",
  "viroli_independent_factor_gibbs"
)
method_labels <- c(
  independent_marginal_mixture = "product MAP",
  viroli_laplace_gibbs = "Viroli Laplace",
  viroli_gaussian_gibbs = "Viroli Gaussian",
  joint_mixture_factor_gibbs = "joint Gibbs",
  viroli_independent_factor_gibbs = "Viroli Gibbs"
)
method_colors <- c(
  independent_marginal_mixture = "#2f6db5",
  viroli_laplace_gibbs = "#c92d32",
  viroli_gaussian_gibbs = "#37a567",
  joint_mixture_factor_gibbs = "#111111",
  viroli_independent_factor_gibbs = "#37a567"
)

p_line_types <- c(
  `250` = 3,
  `500` = 1,
  `1000` = 2,
  `2000` = 4,
  `4000` = 5
)

loading_label <- function(x) {
  out <- ifelse(grepl("dense|cross|multi", x), "Cross", "Sparse")
  out[is.na(out) | !nzchar(out)] <- "Unknown"
  out
}

block_label <- function(x) {
  labels <- c(
    balanced = "balanced blocks",
    ifeval_like = "strongly unbalanced blocks",
    moderate_ifeval_like = "unbalanced blocks"
  )
  out <- labels[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

format_G_config <- function(x) {
  x <- as.character(x)
  vapply(strsplit(x, "-", fixed = TRUE), function(parts) {
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(x)
    if (length(unique(parts)) == 1L && length(parts) > 4L) {
      paste0(unique(parts), " x ", length(parts))
    } else {
      paste(parts, collapse = "-")
    }
  }, character(1L))
}

sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  gsub("(^_+|_+$)", "", x)
}

metric_labels <- c(
  alpha_rmse = "intercepts",
  lambda_rmse = "loadings",
  marginal_mu_rmse = "mixture means",
  marginal_var_rmse = "mixture variances",
  marginal_weight_rmse = "mixture weights",
  factor_score_rmse = "factor scores"
)
metric_labels <- metric_labels[names(metric_labels) %in% names(d)]
if (!length(metric_labels)) {
  stop("No RMSE metric columns found in checkpoint.", call. = FALSE)
}

summarize_checkpoint <- function(d, metric_names) {
  group_cols <- c(
    "method",
    "n",
    "p",
    "H_true",
    "G_true",
    "G_config",
    "separation",
    "loading_design",
    "block_size_mode"
  )
  key <- interaction(d[, group_cols, drop = FALSE], drop = TRUE, sep = " | ")
  pieces <- lapply(split(d, key), function(dd) {
    row <- dd[1L, group_cols, drop = FALSE]
    row$n_reps <- length(unique(dd$rep))
    for (metric in metric_names) {
      x <- dd[[metric]]
      x <- x[is.finite(x)]
      row[[paste0(metric, "_mean")]] <- if (length(x)) mean(x) else NA_real_
      row[[paste0(metric, "_sd")]] <- if (length(x) > 1L) sd(x) else NA_real_
      row[[paste0(metric, "_lower")]] <- if (length(x) > 1L) mean(x) - 2 * sd(x) else NA_real_
      row[[paste0(metric, "_upper")]] <- if (length(x) > 1L) mean(x) + 2 * sd(x) else NA_real_
    }
    row
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out[order(
    out$loading_design,
    out$block_size_mode,
    out$H_true,
    out$G_config,
    out$separation,
    out$p,
    out$n,
    out$method
  ), , drop = FALSE]
}

draw_intervals <- function(x, lower, upper, col, ylim) {
  ok <- is.finite(x) & is.finite(lower) & is.finite(upper) & abs(upper - lower) > 1e-10
  if (!any(ok)) return(invisible(NULL))
  lower <- pmax(ylim[1L], lower)
  upper <- pmin(ylim[2L], upper)
  cap <- pmax(0.9, x * 0.025)
  segments(x[ok], lower[ok], x[ok], upper[ok], col = col, lwd = 1.1)
  segments(x[ok] - cap[ok], lower[ok], x[ok] + cap[ok], lower[ok], col = col, lwd = 1.1)
  segments(x[ok] - cap[ok], upper[ok], x[ok] + cap[ok], upper[ok], col = col, lwd = 1.1)
  invisible(NULL)
}

plot_fixed_panel <- function(summary_df, panel_row, out_file) {
  z <- summary_df[
    summary_df$loading_design == panel_row$loading_design &
      summary_df$block_size_mode == panel_row$block_size_mode &
      summary_df$H_true == panel_row$H_true &
      summary_df$G_config == panel_row$G_config &
      summary_df$separation == panel_row$separation,
    ,
    drop = FALSE
  ]
  if (!nrow(z)) return(FALSE)

  present_methods <- intersect(method_order, unique(z$method))
  present_methods <- c(present_methods, setdiff(unique(z$method), present_methods))
  present_p <- sort(unique(z$p))
  p_lty <- p_line_types[as.character(present_p)]
  missing_lty <- is.na(p_lty)
  if (any(missing_lty)) {
    p_lty[missing_lty] <- rep_len(c(1, 2, 3, 4, 5, 6), sum(missing_lty))
  }
  names(p_lty) <- as.character(present_p)

  n_panels <- length(metric_labels)
  png(out_file, width = 1800, height = if (n_panels <= 3L) 760 else 1100, res = 170)
  op <- par(
    mfrow = c(ceiling(n_panels / 3), 3),
    mar = c(4.5, 4.6, 3.0, 1.0),
    oma = c(4.8, 0, 3.5, 0)
  )
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  for (metric in names(metric_labels)) {
    mean_col <- paste0(metric, "_mean")
    lower_col <- paste0(metric, "_lower")
    upper_col <- paste0(metric, "_upper")
    zm <- z[is.finite(z[[mean_col]]), , drop = FALSE]
    if (!nrow(zm)) {
      plot.new()
      title(main = metric_labels[[metric]])
      next
    }
    all_y <- c(0, zm[[mean_col]], zm[[lower_col]], zm[[upper_col]])
    ylim <- range(all_y[is.finite(all_y)], na.rm = TRUE)
    if (!all(is.finite(ylim)) || diff(ylim) <= 0) ylim <- c(0, 1)
    ylim[1L] <- 0
    ylim[2L] <- max(ylim[2L] * 1.08, 0.01)
    ns <- sort(unique(zm$n))

    plot(
      NA,
      xlim = range(ns),
      ylim = ylim,
      xaxt = "n",
      xlab = "n",
      ylab = "RMSE",
      main = metric_labels[[metric]]
    )
    axis(1, at = ns, labels = ns)
    grid(col = "#e9e9e9")
    abline(h = 0, col = "#cccccc")
    box()

    for (method in present_methods) {
      for (p_value in present_p) {
        zz <- zm[zm$method == method & zm$p == p_value, , drop = FALSE]
        if (!nrow(zz)) next
        zz <- zz[order(zz$n), , drop = FALSE]
        col <- method_colors[[method]]
        if (is.null(col) || is.na(col)) col <- "#555555"
        lty <- p_lty[[as.character(p_value)]]
        lines(zz$n, zz[[mean_col]], col = col, lwd = 2.4, lty = lty)
        points(zz$n, zz[[mean_col]], col = col, pch = 16, cex = 0.9)
        draw_intervals(zz$n, zz[[lower_col]], zz[[upper_col]], col, ylim)
      }
    }
  }

  method_names <- method_labels[present_methods]
  names(method_names) <- NULL
  method_cols <- method_colors[present_methods]
  method_cols[is.na(method_cols)] <- "#555555"
  p_labels <- paste0("p=", present_p)

  legend(
    "bottom",
    inset = c(0, -0.18),
    xpd = NA,
    horiz = TRUE,
    legend = c(method_names, p_labels),
    col = c(method_cols, rep("#222222", length(present_p))),
    lwd = c(rep(3.0, length(present_methods)), rep(2.6, length(present_p))),
    lty = c(rep(1, length(present_methods)), unname(p_lty)),
    pch = c(rep(16, length(present_methods)), rep(NA_integer_, length(present_p))),
    pt.cex = c(rep(1.15, length(present_methods)), rep(1, length(present_p))),
    bty = "o",
    bg = "white",
    box.col = "#555555",
    cex = 0.95
  )
  mtext(
    sprintf(
      "Loadings = \"%s\" | %s | H=%s, G=%s, sep=%s: mean recovery metrics +/- 2 sd",
      unique(loading_label(panel_row$loading_design)),
      unique(block_label(panel_row$block_size_mode)),
      panel_row$H_true,
      format_G_config(panel_row$G_config),
      panel_row$separation
    ),
    outer = TRUE,
    cex = 1.04
  )
  TRUE
}

summary_df <- summarize_checkpoint(d, names(metric_labels))
summary_file <- file.path(out_dir, "progress_fixed_H_np_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)

panel_cols <- c("loading_design", "block_size_mode", "H_true", "G_config", "separation")
panels <- unique(summary_df[, panel_cols, drop = FALSE])
panels <- panels[order(panels$loading_design, panels$block_size_mode, panels$H_true, panels$G_config, panels$separation), , drop = FALSE]

plot_dir <- file.path(out_dir, "fixed_H_np_panels")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
out_files <- character(0)
for (i in seq_len(nrow(panels))) {
  panel <- panels[i, , drop = FALSE]
  file_tag <- paste(
    sanitize(loading_label(panel$loading_design)),
    sanitize(panel$block_size_mode),
    paste0("H", panel$H_true),
    paste0("G", sanitize(panel$G_config)),
    paste0("sep", sanitize(panel$separation)),
    sep = "_"
  )
  out_file <- file.path(plot_dir, paste0("progress_fixed_H_np_rmse_panel_", file_tag, ".png"))
  made <- plot_fixed_panel(summary_df, panel, out_file)
  if (isTRUE(made)) out_files <- c(out_files, out_file)
}

cat("Rows:", nrow(d), "\n")
cat("Summary:", normalizePath(summary_file, mustWork = FALSE), "\n")
cat("Fixed-H plots:\n")
cat(paste(normalizePath(out_files, mustWork = FALSE), collapse = "\n"), "\n")
