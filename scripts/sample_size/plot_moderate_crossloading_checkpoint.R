#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

get_env <- function(name, default, FUN = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  FUN(value)
}

out_dir <- get_env(
  "OUT_DIR",
  "results/moderate_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP",
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

summarize_metric <- function(results, metric) {
  group_cols <- c("loading_design", "method", "H_true", "G_true", "n", "p")
  group_key <- interaction(results[, group_cols, drop = FALSE], drop = TRUE, sep = " | ")
  out <- lapply(split(results, group_key), function(d) {
    x <- d[[metric]]
    x <- x[is.finite(x)]
    data.frame(
      loading_design = d$loading_design[1L],
      method = d$method[1L],
      H_true = d$H_true[1L],
      G_true = d$G_true[1L],
      n = d$n[1L],
      p = d$p[1L],
      n_reps = length(x),
      mean = if (length(x)) mean(x) else NA_real_,
      sd = if (length(x) > 1L) sd(x) else NA_real_,
      lower = if (length(x) > 1L) mean(x) - 2 * sd(x) else NA_real_,
      upper = if (length(x) > 1L) mean(x) + 2 * sd(x) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$loading_design, out$H_true, out$G_true, out$n, out$method), ]
}

plot_metric <- function(summary, metric_label, out_file, ylim = c(-0.05, 1.02)) {
  designs <- sort(unique(summary$loading_design))
  panels <- unique(summary[, c("loading_design", "H_true", "G_true"), drop = FALSE])
  panels <- panels[order(panels$loading_design, panels$H_true, panels$G_true), , drop = FALSE]
  if (!nrow(panels)) return(invisible(FALSE))

  n_col <- length(unique(panels$G_true))
  n_row <- ceiling(nrow(panels) / n_col)
  png(out_file, width = max(1300, 620 * n_col), height = max(850, 470 * n_row), res = 160)
  op <- par(mfrow = c(n_row, n_col), mar = c(4.3, 4.2, 3.4, 1), oma = c(4, 0, 3, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  for (idx in seq_len(nrow(panels))) {
    panel <- panels[idx, ]
    d <- summary[
      summary$loading_design == panel$loading_design &
        summary$H_true == panel$H_true &
        summary$G_true == panel$G_true,
      ,
      drop = FALSE
    ]
    ns <- sort(unique(d$n))
    plot(
      NA,
      xlim = range(ns),
      ylim = ylim,
      log = if (length(ns) > 1L) "x" else "",
      xaxt = "n",
      xlab = "n",
      ylab = metric_label,
      main = sprintf("%s\nH=%s, G=%s", panel$loading_design, panel$H_true, panel$G_true)
    )
    axis(1, at = ns, labels = ns)
    abline(h = c(0, 0.9), col = c("#D0D0D0", "#777777"), lty = c(1, 3))
    box()

    for (method in unique(d$method)) {
      dm <- d[d$method == method, , drop = FALSE]
      dm <- dm[order(dm$n), , drop = FALSE]
      col <- method_colors[[method]]
      lines(dm$n, dm$mean, col = col, lwd = 2)
      points(dm$n, dm$mean, col = col, pch = 16, cex = 0.9)
      has_interval <- is.finite(dm$lower) & is.finite(dm$upper)
      if (any(has_interval)) {
        arrows(
          x0 = dm$n[has_interval],
          y0 = pmax(ylim[1L], dm$lower[has_interval]),
          x1 = dm$n[has_interval],
          y1 = pmin(ylim[2L], dm$upper[has_interval]),
          angle = 90,
          code = 3,
          length = 0.035,
          col = col,
          lwd = 1.1
        )
      }
      text(
        dm$n,
        pmin(ylim[2L] - 0.04, dm$mean + 0.045),
        labels = dm$n_reps,
        col = col,
        cex = 0.65
      )
    }
    if (idx == 1L) {
      methods_present <- unique(d$method)
      legend(
        "bottomright",
        legend = method_labels[methods_present],
        col = method_colors[methods_present],
        lwd = 2,
        pch = 16,
        bty = "n",
        cex = 0.8
      )
    }
  }

  mtext(
    sprintf("%s: mean +/- 2 sd across completed reps; point labels are completed reps", metric_label),
    outer = TRUE,
    cex = 1.05
  )
  invisible(TRUE)
}

plot_progress <- function(results, out_file) {
  agg <- aggregate(rep ~ loading_design + method + H_true + G_true + n, results, function(x) length(unique(x)))
  names(agg)[names(agg) == "rep"] <- "n_reps"
  agg$cell <- paste0("H", agg$H_true, " G", agg$G_true, " n", agg$n)
  agg$loading_design <- factor(agg$loading_design, levels = sort(unique(agg$loading_design)))
  agg$method <- factor(agg$method, levels = names(method_labels))
  cells <- unique(agg$cell[order(agg$H_true, agg$G_true, agg$n)])
  rows <- interaction(agg$loading_design, agg$method, sep = " | ", drop = TRUE)
  row_levels <- unique(as.character(rows))
  mat <- matrix(0, nrow = length(row_levels), ncol = length(cells), dimnames = list(row_levels, cells))
  for (i in seq_len(nrow(agg))) {
    mat[as.character(rows[i]), agg$cell[i]] <- agg$n_reps[i]
  }

  png(out_file, width = 1700, height = max(700, 170 * nrow(mat)), res = 160)
  par(mar = c(8, 13, 4, 5), xaxs = "i", yaxs = "i")
  pal <- colorRampPalette(c("#F7FBFF", "#9ECAE1", "#08519C"))(26)
  image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = pal,
    zlim = c(0, 25),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = sprintf("Checkpoint progress: %d / 1600 result rows", nrow(results))
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.75)
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 1, cex.axis = 0.75)
  for (i in seq_len(nrow(mat))) {
    for (j in seq_len(ncol(mat))) {
      text(j, nrow(mat) - i + 1L, labels = mat[i, j], cex = 0.7)
    }
  }
  box()
  invisible(dev.off())
}

sanitize_file_tag <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

plot_parameter_panel <- function(results, loading_design, H_value, G_value, out_file, summary_file) {
  d0 <- results[
    results$loading_design == loading_design &
      results$H_true == H_value &
      results$G_true == G_value,
    ,
    drop = FALSE
  ]
  if (!nrow(d0)) return(invisible(FALSE))

  metric_labels <- c(
    flat_parameter_corr = "flat parameters",
    all_parameter_corr = "all parameters",
    alpha_corr = "intercepts",
    mean_factor_abs_cor = "factor scores",
    lambda_corr = "loadings",
    joint_mu_corr = "mixture means",
    joint_var_corr = "mixture variances",
    joint_weight_corr = "mixture weights"
  )
  metric_labels <- metric_labels[names(metric_labels) %in% names(d0)]

  group_cols <- c("method", "n", "p", "H_true", "G_true", "K_joint", "loading_design")
  group_key <- interaction(d0[, group_cols, drop = FALSE], drop = TRUE, sep = " | ")
  summary_rows <- lapply(split(d0, group_key), function(d) {
    out <- d[1L, group_cols, drop = FALSE]
    out$n_reps <- length(unique(d$rep))
    for (metric in names(metric_labels)) {
      x <- d[[metric]]
      x <- x[is.finite(x)]
      m <- if (length(x)) mean(x) else NA_real_
      s <- if (length(x) > 1L) sd(x) else NA_real_
      out[[paste0("mean_", metric)]] <- m
      out[[paste0("sd_", metric)]] <- s
      out[[paste0("lower_", metric)]] <- if (is.finite(s)) m - 2 * s else NA_real_
      out[[paste0("upper_", metric)]] <- if (is.finite(s)) m + 2 * s else NA_real_
    }
    out
  })
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  summary <- summary[order(summary$n, summary$method), , drop = FALSE]
  write.csv(summary, summary_file, row.names = FALSE)

  n_col <- 3L
  n_row <- ceiling(length(metric_labels) / n_col)
  png(out_file, width = 1800, height = max(900, 520 * n_row), res = 170)
  op <- par(mfrow = c(n_row, n_col), mar = c(4.4, 4.5, 3.2, 1), oma = c(0, 0, 3.2, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  for (metric in names(metric_labels)) {
    mean_col <- paste0("mean_", metric)
    lower_col <- paste0("lower_", metric)
    upper_col <- paste0("upper_", metric)
    finite_y <- is.finite(summary[[mean_col]])

    if (!any(finite_y)) {
      plot.new()
      title(main = metric_labels[[metric]])
      next
    }

    ylim <- range(c(summary[[mean_col]], summary[[lower_col]], summary[[upper_col]], 0, 1), na.rm = TRUE)
    ylim <- c(max(-1, ylim[1L]), min(1, ylim[2L]))
    ns <- sort(unique(summary$n))

    plot(
      NA,
      xlim = range(summary$n[finite_y]),
      ylim = ylim,
      log = if (length(ns) > 1L) "x" else "",
      xaxt = "n",
      xlab = "n",
      ylab = "correlation",
      main = metric_labels[[metric]]
    )
    axis(1, at = ns, labels = ns)
    abline(h = c(0, 0.9), col = c("#CFCFCF", "#777777"), lty = c(1, 3), lwd = c(1, 1.2))

    for (method in unique(summary$method)) {
      d <- summary[summary$method == method & is.finite(summary[[mean_col]]), , drop = FALSE]
      if (!nrow(d)) next
      d <- d[order(d$n), , drop = FALSE]
      col <- method_colors[[method]]

      lines(d$n, d[[mean_col]], col = col, lwd = 2)
      points(d$n, d[[mean_col]], col = col, pch = 16, cex = 0.95)

      has_sd <- is.finite(d[[lower_col]]) & is.finite(d[[upper_col]]) &
        abs(d[[upper_col]] - d[[lower_col]]) > 1e-10
      if (any(has_sd)) {
        arrows(
          d$n[has_sd],
          pmax(-1, d[[lower_col]][has_sd]),
          d$n[has_sd],
          pmin(1, d[[upper_col]][has_sd]),
          code = 3,
          angle = 90,
          length = 0.04,
          col = col,
          lwd = 1.1
        )
      }
      text(d$n, d[[mean_col]], labels = d$n_reps, pos = 3, cex = 0.58, col = col)
    }
    box()
  }

  used_methods <- unique(summary$method)
  legend(
    "bottom",
    inset = c(0, -0.04),
    xpd = NA,
    horiz = TRUE,
    legend = method_labels[used_methods],
    col = method_colors[used_methods],
    lwd = 2,
    pch = 16,
    bty = "n",
    cex = 0.9
  )
  mtext(
    sprintf(
      "%s | H=%d, G=%d (K=%d): mean correlations +/- 2 sd; point labels are completed reps",
      loading_design,
      H_value,
      G_value,
      G_value^H_value
    ),
    outer = TRUE,
    cex = 1.02
  )
  invisible(TRUE)
}

plot_all_parameter_panels <- function(results, out_dir) {
  panels <- unique(results[, c("loading_design", "H_true", "G_true"), drop = FALSE])
  panels <- panels[order(panels$loading_design, panels$H_true, panels$G_true), , drop = FALSE]
  out_files <- character(0)
  for (i in seq_len(nrow(panels))) {
    loading_design <- panels$loading_design[i]
    H_value <- panels$H_true[i]
    G_value <- panels$G_true[i]
    tag <- sprintf("%s_H%d_G%d", sanitize_file_tag(loading_design), H_value, G_value)
    plot_file <- file.path(out_dir, paste0("checkpoint_all_parameter_panel_", tag, ".png"))
    summary_file <- file.path(out_dir, paste0("checkpoint_all_parameter_panel_", tag, "_summary.csv"))
    plot_parameter_panel(results, loading_design, H_value, G_value, plot_file, summary_file)
    out_files <- c(out_files, plot_file)
  }
  out_files
}

flat_summary <- summarize_metric(results, "flat_parameter_corr")
alpha_summary <- summarize_metric(results, "alpha_corr")
alpha_rmse_summary <- summarize_metric(results, "alpha_rmse")
alpha_raw_summary <- if ("alpha_raw_corr" %in% names(results)) summarize_metric(results, "alpha_raw_corr") else data.frame()
alpha_raw_rmse_summary <- if ("alpha_raw_rmse" %in% names(results)) summarize_metric(results, "alpha_raw_rmse") else data.frame()
var_summary <- summarize_metric(results, "joint_var_corr")
lambda_summary <- summarize_metric(results, "lambda_corr")
time_summary <- summarize_metric(results, "seconds")

write.csv(flat_summary, file.path(out_dir, "checkpoint_flat_parameter_corr_summary.csv"), row.names = FALSE)
write.csv(alpha_summary, file.path(out_dir, "checkpoint_alpha_corr_summary.csv"), row.names = FALSE)
write.csv(alpha_rmse_summary, file.path(out_dir, "checkpoint_alpha_rmse_summary.csv"), row.names = FALSE)
if (nrow(alpha_raw_summary)) {
  write.csv(alpha_raw_summary, file.path(out_dir, "checkpoint_alpha_raw_corr_summary.csv"), row.names = FALSE)
}
if (nrow(alpha_raw_rmse_summary)) {
  write.csv(alpha_raw_rmse_summary, file.path(out_dir, "checkpoint_alpha_raw_rmse_summary.csv"), row.names = FALSE)
}
write.csv(var_summary, file.path(out_dir, "checkpoint_joint_var_corr_summary.csv"), row.names = FALSE)
write.csv(lambda_summary, file.path(out_dir, "checkpoint_lambda_corr_summary.csv"), row.names = FALSE)
write.csv(time_summary, file.path(out_dir, "checkpoint_seconds_summary.csv"), row.names = FALSE)

plot_metric(flat_summary, "flat parameter correlation", file.path(out_dir, "checkpoint_flat_parameter_corr_lines.png"))
plot_metric(alpha_summary, "alpha correlation", file.path(out_dir, "checkpoint_alpha_corr_lines.png"))
if (nrow(alpha_raw_summary)) {
  plot_metric(alpha_raw_summary, "raw Gibbs alpha correlation", file.path(out_dir, "checkpoint_alpha_raw_corr_lines.png"))
}
plot_metric(
  alpha_rmse_summary,
  "alpha RMSE",
  file.path(out_dir, "checkpoint_alpha_rmse_lines.png"),
  ylim = range(c(0, alpha_rmse_summary$upper, alpha_rmse_summary$mean), na.rm = TRUE)
)
if (nrow(alpha_raw_rmse_summary)) {
  plot_metric(
    alpha_raw_rmse_summary,
    "raw Gibbs alpha RMSE",
    file.path(out_dir, "checkpoint_alpha_raw_rmse_lines.png"),
    ylim = range(c(0, alpha_raw_rmse_summary$upper, alpha_raw_rmse_summary$mean), na.rm = TRUE)
  )
}
plot_metric(var_summary, "joint variance correlation", file.path(out_dir, "checkpoint_joint_var_corr_lines.png"), ylim = c(-0.1, 1.02))
plot_metric(lambda_summary, "lambda correlation", file.path(out_dir, "checkpoint_lambda_corr_lines.png"))
plot_metric(time_summary, "seconds", file.path(out_dir, "checkpoint_seconds_lines.png"), ylim = range(c(0, time_summary$upper, time_summary$mean), na.rm = TRUE))
plot_progress(results, file.path(out_dir, "checkpoint_progress_heatmap.png"))
panel_files <- plot_all_parameter_panels(results, out_dir)

cat("Wrote checkpoint plots to:", normalizePath(out_dir), "\n")
if (length(panel_files)) {
  cat("Wrote all-parameter panels:\n")
  cat(paste(normalizePath(panel_files, mustWork = FALSE), collapse = "\n"), "\n")
}
