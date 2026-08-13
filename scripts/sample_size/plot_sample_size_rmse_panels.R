#!/usr/bin/env Rscript

# Plot parameter recovery metrics from the sample-size comparison. Most panels
# are RMSE/error metrics; the factor-score panels report correlations.
# The input is the resumable comparison_results_checkpoint.csv produced by
# compare_original_simulation_joint_mfa_gibbs.R.

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

compute_weight_rmse_from_recovery_tables <- function(out_dir) {
  files <- list.files(out_dir, pattern = "^joint_parameter_recovery_.*[.]csv$", full.names = TRUE)
  if (!length(files)) return(data.frame())

  rows <- lapply(files, function(path) {
    tab <- tryCatch(read.csv(path), error = function(e) NULL)
    if (is.null(tab) || !all(c("method", "true_weight", "est_weight") %in% names(tab))) {
      return(NULL)
    }
    scenario <- basename(path)
    scenario <- sub("^joint_parameter_recovery_(ours|gibbs)_", "", scenario)
    scenario <- sub("[.]csv$", "", scenario)
    x <- tab$est_weight - tab$true_weight
    x <- x[is.finite(x)]
    data.frame(
      scenario = scenario,
      method = tab$method[1L],
      joint_weight_rmse = if (length(x)) sqrt(mean(x^2)) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) data.frame() else out
}

weight_rmse <- compute_weight_rmse_from_recovery_tables(out_dir)
if (nrow(weight_rmse)) {
  results <- merge(
    results,
    weight_rmse,
    by = c("scenario", "method"),
    all.x = TRUE,
    sort = FALSE
  )
} else if (!("joint_weight_rmse" %in% names(results)) &&
           all(c("joint_weight_l1", "K_joint") %in% names(results))) {
  # Fallback for old output folders that do not retain the per-profile recovery
  # tables.  Prefer the exact vectorized RMSE above whenever possible.
  results$joint_weight_rmse <- results$joint_weight_l1 / sqrt(pmax(results$K_joint, 1))
}

metric_labels <- c(
  alpha_rmse = "intercepts",
  lambda_rmse = "loadings",
  joint_mu_rmse = "mixture means",
  joint_var_rmse = "mixture variances",
  joint_weight_rmse = "mixture weights",
  flat_parameter_corr = "flat factor correlation"
)
metric_labels <- metric_labels[names(metric_labels) %in% names(results)]

correlation_metrics <- c("flat_parameter_corr")

summarize_metrics <- function(results, metric_names) {
  group_cols <- c("loading_design", "method", "H_true", "G_true", "n", "p")
  group_key <- interaction(results[, group_cols, drop = FALSE], drop = TRUE, sep = " | ")
  out <- lapply(split(results, group_key), function(d) {
    row <- d[1L, group_cols, drop = FALSE]
    row$n_reps <- length(unique(d$rep))
    for (metric in metric_names) {
      x <- d[[metric]]
      x <- x[is.finite(x)]
      row[[paste0("mean_", metric)]] <- if (length(x)) mean(x) else NA_real_
      row[[paste0("sd_", metric)]] <- if (length(x) > 1L) sd(x) else NA_real_
      row[[paste0("lower_", metric)]] <- if (length(x) > 1L) mean(x) - 2 * sd(x) else NA_real_
      row[[paste0("upper_", metric)]] <- if (length(x) > 1L) mean(x) + 2 * sd(x) else NA_real_
    }
    row
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$loading_design, out$H_true, out$G_true, out$n, out$method), , drop = FALSE]
}

plot_recovery_panel <- function(summary, loading_design, H_value, G_value, out_file) {
  d0 <- summary[
    summary$loading_design == loading_design &
      summary$H_true == H_value &
      summary$G_true == G_value,
    ,
    drop = FALSE
  ]
  if (!nrow(d0)) return(invisible(FALSE))
  loading_label <- loading_labels[[loading_design]]
  if (is.null(loading_label)) loading_label <- loading_design

  n_col <- 3L
  n_row <- ceiling(length(metric_labels) / n_col)
  png(out_file, width = 1800, height = max(900, 510 * n_row), res = 170)
  op <- par(mfrow = c(n_row, n_col), mar = c(4.4, 4.5, 3.1, 1), oma = c(4.8, 0, 3.2, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  for (metric in names(metric_labels)) {
    mean_col <- paste0("mean_", metric)
    lower_col <- paste0("lower_", metric)
    upper_col <- paste0("upper_", metric)
    finite_y <- is.finite(d0[[mean_col]])
    if (!any(finite_y)) {
      plot.new()
      title(main = metric_labels[[metric]])
      next
    }

    y_values <- c(0, d0[[mean_col]], d0[[lower_col]], d0[[upper_col]])
    if (metric %in% correlation_metrics) {
      ylim <- range(c(0, 1, y_values[is.finite(y_values)]), na.rm = TRUE)
      ylim <- c(max(-0.05, ylim[1L]), min(1.02, ylim[2L]))
    } else {
      ylim <- range(y_values[is.finite(y_values)], na.rm = TRUE)
      ylim[1L] <- 0
      if (!is.finite(ylim[2L]) || ylim[2L] <= 0) ylim[2L] <- 1
      ylim[2L] <- ylim[2L] * 1.08
    }
    ns <- sort(unique(d0$n))

    plot(
      NA,
      xlim = range(ns),
      ylim = ylim,
      log = if (length(ns) > 1L) "x" else "",
      xaxt = "n",
      xlab = "n",
      ylab = if (metric %in% correlation_metrics) "correlation" else "RMSE / error",
      main = metric_labels[[metric]]
    )
    axis(1, at = ns, labels = ns)
    if (metric %in% correlation_metrics) {
      abline(h = c(0, 0.9), col = c("#CFCFCF", "#777777"), lty = c(1, 3))
    } else {
      abline(h = 0, col = "#CFCFCF")
    }
    box()

    for (method in unique(d0$method)) {
      d <- d0[d0$method == method & is.finite(d0[[mean_col]]), , drop = FALSE]
      if (!nrow(d)) next
      d <- d[order(d$n), , drop = FALSE]
      col <- method_colors[[method]]
      lines(d$n, d[[mean_col]], col = col, lwd = 2)
      points(d$n, d[[mean_col]], col = col, pch = 16, cex = 0.95)

      y_lower <- if (metric %in% correlation_metrics) {
        pmax(ylim[1L], d[[lower_col]])
      } else {
        pmax(0, d[[lower_col]])
      }
      y_upper <- pmin(ylim[2L], d[[upper_col]])
      has_interval <- is.finite(y_lower) & is.finite(y_upper) &
        abs(y_upper - y_lower) > 1e-10
      if (any(has_interval)) {
        x <- d$n[has_interval]
        cap_width <- x * 0.035
        segments(x, y_lower[has_interval], x, y_upper[has_interval], col = col, lwd = 1.1)
        segments(x - cap_width, y_lower[has_interval], x + cap_width, y_lower[has_interval],
                 col = col, lwd = 1.1)
        segments(x - cap_width, y_upper[has_interval], x + cap_width, y_upper[has_interval],
                 col = col, lwd = 1.1)
      }
    }
  }

  used_methods <- unique(d0$method)
  mtext(
    sprintf(
      "%s | H=%d, G=%d: mean recovery metrics +/- 2 sd",
      loading_label,
      H_value,
      G_value
    ),
    outer = TRUE,
    cex = 1.02
  )
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend(
    x = 0.5,
    y = 0.035,
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
    cex = 1.08
  )
  invisible(TRUE)
}

sanitize_file_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

summary <- summarize_metrics(results, names(metric_labels))
write.csv(summary, file.path(out_dir, "checkpoint_parameter_recovery_summary.csv"), row.names = FALSE)

panels <- unique(summary[, c("loading_design", "H_true", "G_true"), drop = FALSE])
panels <- panels[order(panels$loading_design, panels$H_true, panels$G_true), , drop = FALSE]
out_files <- character(0)
for (i in seq_len(nrow(panels))) {
  loading_design <- panels$loading_design[i]
  H_value <- panels$H_true[i]
  G_value <- panels$G_true[i]
  tag <- sprintf("%s_H%d_G%d", sanitize_file_tag(loading_design), H_value, G_value)
  out_file <- file.path(out_dir, paste0("checkpoint_parameter_recovery_panel_", tag, ".png"))
  plot_recovery_panel(summary, loading_design, H_value, G_value, out_file)
  out_files <- c(out_files, out_file)
}

cat("Wrote parameter recovery summary to:", normalizePath(file.path(out_dir, "checkpoint_parameter_recovery_summary.csv")), "\n")
cat("Wrote parameter recovery panels:\n")
cat(paste(normalizePath(out_files, mustWork = FALSE), collapse = "\n"), "\n")
