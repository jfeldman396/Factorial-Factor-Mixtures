#!/usr/bin/env Rscript

# Plot timing summaries from the sample-size comparison. Each output figure is
# one simulation setting, with mean log(seconds) +/- 2 sd across repetitions.

options(stringsAsFactors = FALSE)

get_env <- function(name, default, FUN = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  FUN(value)
}

out_dir <- get_env(
  "OUT_DIR",
  file.path(
    "..",
    "results",
    "moderate_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP_intercepts_centered"
  ),
  as.character
)
checkpoint_file <- file.path(out_dir, "comparison_results_checkpoint.csv")
if (!file.exists(checkpoint_file)) {
  stop("Checkpoint file not found: ", checkpoint_file)
}

results <- read.csv(checkpoint_file)
if (!nrow(results)) stop("Checkpoint file has no rows.")
if (!("seconds" %in% names(results))) stop("Checkpoint file has no seconds column.")

method_colors <- c(
  independent_marginal_mixture = "#2B6CB0",
  joint_mixture_factor_gibbs = "#C53030"
)
method_labels <- c(
  independent_marginal_mixture = "product mixture",
  joint_mixture_factor_gibbs = "joint Gibbs"
)
loading_labels <- c(
  balanced_moderate_few_positive_cross = 'Loadings = "Sparse"',
  balanced_moderate_dense_signed_cross = 'Loadings = "Cross"'
)
block_size_labels <- c(
  balanced = "balanced blocks",
  ifeval_like = "strongly unbalanced blocks",
  moderate_ifeval_like = "unbalanced blocks"
)

summarize_timing <- function(results) {
  keep <- is.finite(results$seconds) & results$seconds > 0
  results <- results[keep, , drop = FALSE]
  results$log_seconds <- log(results$seconds)

  group_cols <- c("loading_design", "method", "H_true", "G_true", "n", "p")
  if ("block_size_mode" %in% names(results)) {
    group_cols <- c("block_size_mode", group_cols)
  }
  group_key <- interaction(results[, group_cols, drop = FALSE], drop = TRUE, sep = " | ")
  out <- lapply(split(results, group_key), function(d) {
    row <- d[1L, group_cols, drop = FALSE]
    row$n_reps <- length(unique(d$rep))
    row$mean_seconds <- mean(d$seconds)
    row$median_seconds <- median(d$seconds)
    row$mean_log_seconds <- mean(d$log_seconds)
    row$sd_log_seconds <- if (nrow(d) > 1L) sd(d$log_seconds) else NA_real_
    row$lower_log_seconds <- if (nrow(d) > 1L) row$mean_log_seconds - 2 * row$sd_log_seconds else NA_real_
    row$upper_log_seconds <- if (nrow(d) > 1L) row$mean_log_seconds + 2 * row$sd_log_seconds else NA_real_
    row$geometric_mean_seconds <- exp(row$mean_log_seconds)
    row
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$loading_design, out$H_true, out$G_true, out$n, out$method), , drop = FALSE]
}

sanitize_file_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

plot_timing_lines <- function(summary, loading_design, H_value, G_value, block_size_mode, out_file) {
  G_value <- as.character(G_value)
  d0 <- summary[
    summary$loading_design == loading_design &
      summary$H_true == H_value &
      as.character(summary$G_true) == G_value &
      (!("block_size_mode" %in% names(summary)) | summary$block_size_mode == block_size_mode),
    ,
    drop = FALSE
  ]
  if (!nrow(d0)) return(invisible(FALSE))

  loading_label <- loading_labels[[loading_design]]
  if (is.null(loading_label)) loading_label <- loading_design
  block_label <- if (!is.na(block_size_mode) && block_size_mode %in% names(block_size_labels)) {
    block_size_labels[[block_size_mode]]
  } else if (!is.na(block_size_mode) && is.finite(nchar(block_size_mode))) {
    block_size_mode
  } else {
    "balanced blocks"
  }

  y_values <- c(d0$mean_log_seconds, d0$lower_log_seconds, d0$upper_log_seconds)
  y_values <- y_values[is.finite(y_values)]
  ylim <- range(y_values, na.rm = TRUE)
  if (!all(is.finite(ylim)) || diff(ylim) == 0) ylim <- ylim + c(-0.5, 0.5)
  ylim <- ylim + c(-0.05, 0.05) * diff(ylim)

  ns <- sort(unique(d0$n))
  png(out_file, width = 1200, height = 800, res = 170)
  op <- par(mar = c(4.6, 5.1, 3.5, 1.2), oma = c(3.5, 0, 0, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  plot(
    NA,
    xlim = range(ns),
    ylim = ylim,
    log = if (length(ns) > 1L) "x" else "",
    xaxt = "n",
    xlab = "n",
    ylab = "log(seconds)",
    main = sprintf("%s | %s | H=%d, G=%s: runtime", loading_label, block_label, H_value, G_value)
  )
  axis(1, at = ns, labels = ns)
  grid(col = "#E2E2E2")
  box()

  for (method in unique(d0$method)) {
    d <- d0[d0$method == method & is.finite(d0$mean_log_seconds), , drop = FALSE]
    if (!nrow(d)) next
    d <- d[order(d$n), , drop = FALSE]
    col <- method_colors[[method]]
    lines(d$n, d$mean_log_seconds, col = col, lwd = 2.5)
    points(d$n, d$mean_log_seconds, col = col, pch = 16, cex = 1.05)

    has_interval <- is.finite(d$lower_log_seconds) & is.finite(d$upper_log_seconds) &
      abs(d$upper_log_seconds - d$lower_log_seconds) > 1e-10
    if (any(has_interval)) {
      x <- d$n[has_interval]
      cap_width <- x * 0.035
      segments(x, d$lower_log_seconds[has_interval], x, d$upper_log_seconds[has_interval],
               col = col, lwd = 1.2)
      segments(x - cap_width, d$lower_log_seconds[has_interval],
               x + cap_width, d$lower_log_seconds[has_interval], col = col, lwd = 1.2)
      segments(x - cap_width, d$upper_log_seconds[has_interval],
               x + cap_width, d$upper_log_seconds[has_interval], col = col, lwd = 1.2)
    }
  }

  used_methods <- unique(d0$method)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend(
    x = 0.5,
    y = 0.045,
    xjust = 0.5,
    yjust = 0.5,
    horiz = TRUE,
    legend = method_labels[used_methods],
    col = method_colors[used_methods],
    lwd = 3,
    pch = 16,
    pt.cex = 1.25,
    bty = "o",
    bg = "white",
    box.col = "#555555",
    cex = 1
  )

  invisible(TRUE)
}

summary <- summarize_timing(results)
write.csv(summary, file.path(out_dir, "checkpoint_timing_log_seconds_summary.csv"), row.names = FALSE)

panel_cols <- c("loading_design", "H_true", "G_true")
if ("block_size_mode" %in% names(summary)) panel_cols <- c("block_size_mode", panel_cols)
panels <- unique(summary[, panel_cols, drop = FALSE])
panels <- panels[order(panels$loading_design, panels$H_true, panels$G_true), , drop = FALSE]
out_files <- character(0)
for (i in seq_len(nrow(panels))) {
  loading_design <- panels$loading_design[i]
  H_value <- panels$H_true[i]
  G_value <- as.character(panels$G_true[i])
  block_size_mode <- if ("block_size_mode" %in% names(panels)) panels$block_size_mode[i] else NA_character_
  tag_prefix <- if (is.na(block_size_mode)) "" else paste0(sanitize_file_tag(block_size_mode), "_")
  tag <- sprintf("%s%s_H%d_G%s", tag_prefix, sanitize_file_tag(loading_design), H_value, sanitize_file_tag(G_value))
  out_file <- file.path(out_dir, paste0("checkpoint_timing_log_seconds_", tag, ".png"))
  plot_timing_lines(summary, loading_design, H_value, G_value, block_size_mode, out_file)
  out_files <- c(out_files, out_file)
}

cat("Wrote timing summary to:", normalizePath(file.path(out_dir, "checkpoint_timing_log_seconds_summary.csv")), "\n")
cat("Wrote timing plots:\n")
cat(paste(normalizePath(out_files, mustWork = FALSE), collapse = "\n"), "\n")
