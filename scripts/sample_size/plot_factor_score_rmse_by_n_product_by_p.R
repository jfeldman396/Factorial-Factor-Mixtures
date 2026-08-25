#!/usr/bin/env Rscript

# Plot product-mixture factor-score RMSE as n grows, with separate lines for
# each p.  This is the clearest view of factor-score convergence in the full
# p-grid because the proposed product-mixture estimator is fit for every p,
# while the Gibbs comparator is only fit for p in {250, 500}.

options(stringsAsFactors = FALSE)

get_env <- function(name, default, FUN = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  FUN(value)
}

repo_root <- normalizePath(get_env("REPO_ROOT", ".", as.character))
selected_table_dir <- file.path(repo_root, "results", "selected_tables", "sample_size")
selected_plot_dir <- file.path(repo_root, "results", "selected_plots", "sample_size")
dir.create(selected_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(selected_plot_dir, recursive = TRUE, showWarnings = FALSE)

result_files <- c(
  balanced = file.path(
    selected_table_dir,
    "balanced_sampledZ_pgrid_product_allp_gibbs_smallp_comparison_results.csv"
  ),
  unbalanced = file.path(
    selected_table_dir,
    "unbalanced_sampledZ_pgrid_product_allp_gibbs_smallp_comparison_results.csv"
  )
)

missing_files <- result_files[!file.exists(result_files)]
if (length(missing_files)) {
  stop("Missing selected result file(s):\n", paste(missing_files, collapse = "\n"))
}

loading_labels <- c(
  balanced_moderate_few_positive_cross = "Sparse",
  balanced_moderate_dense_signed_cross = "Cross"
)
block_labels <- c(
  balanced = "Balanced blocks",
  unbalanced = "Unbalanced blocks"
)
p_colors <- c(
  `250` = "#2B6CB0",
  `500` = "#2F855A",
  `1000` = "#D69E2E",
  `2000` = "#C53030"
)

summarize_factor_rmse <- function(results, block_label) {
  keep <- results$method == "independent_marginal_mixture" &
    is.finite(results$factor_score_rmse)
  d <- results[keep, , drop = FALSE]
  d$block_label <- block_label

  group_cols <- c("block_label", "loading_design", "H_true", "G_true", "n", "p")
  group_key <- interaction(d[, group_cols, drop = FALSE], drop = TRUE, sep = " | ")
  out <- lapply(split(d, group_key), function(x) {
    row <- x[1L, group_cols, drop = FALSE]
    y <- x$factor_score_rmse[is.finite(x$factor_score_rmse)]
    row$n_reps <- length(unique(x$rep))
    row$mean_factor_score_rmse <- mean(y)
    row$sd_factor_score_rmse <- if (length(y) > 1L) sd(y) else NA_real_
    row$lower_factor_score_rmse <- if (length(y) > 1L) mean(y) - 2 * sd(y) else NA_real_
    row$upper_factor_score_rmse <- if (length(y) > 1L) mean(y) + 2 * sd(y) else NA_real_
    row
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$loading_design, out$H_true, out$G_true, out$p, out$n), , drop = FALSE]
}

panel_order <- function(summary) {
  loading_order <- c(
    "balanced_moderate_few_positive_cross",
    "balanced_moderate_dense_signed_cross"
  )
  panels <- unique(summary[, c("loading_design", "H_true", "G_true"), drop = FALSE])
  panels$loading_rank <- match(panels$loading_design, loading_order)
  panels <- panels[order(panels$loading_rank, panels$H_true, panels$G_true), , drop = FALSE]
  panels$loading_rank <- NULL
  panels
}

plot_factor_rmse_by_n <- function(summary, block_name, out_file) {
  panels <- panel_order(summary)
  n_values <- sort(unique(summary$n))
  p_values <- sort(unique(summary$p))
  plot_cols <- p_colors[as.character(p_values)]
  missing_cols <- is.na(plot_cols)
  if (any(missing_cols)) {
    plot_cols[missing_cols] <- rep_len(c("#2B6CB0", "#2F855A", "#D69E2E", "#C53030"),
                                       sum(missing_cols))
  }
  names(plot_cols) <- as.character(p_values)

  png(out_file, width = 2300, height = 1450, res = 185)
  op <- par(mfrow = c(2, 4), mar = c(4.2, 4.45, 3.1, 1.1), oma = c(5.9, 0, 4.1, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  for (i in seq_len(nrow(panels))) {
    panel <- panels[i, , drop = FALSE]
    d0 <- summary[
      summary$loading_design == panel$loading_design &
        summary$H_true == panel$H_true &
        summary$G_true == panel$G_true,
      ,
      drop = FALSE
    ]
    y_values <- c(
      d0$mean_factor_score_rmse,
      pmax(0, d0$lower_factor_score_rmse),
      d0$upper_factor_score_rmse
    )
    ylim <- range(c(0, y_values[is.finite(y_values)]), na.rm = TRUE)
    if (!all(is.finite(ylim)) || ylim[2L] <= 0) ylim <- c(0, 1)
    ylim[2L] <- ylim[2L] * 1.08

    label <- loading_labels[[panel$loading_design]]
    if (is.null(label)) label <- panel$loading_design
    main <- sprintf("%s | H=%d, G=%s", label, panel$H_true, panel$G_true)

    plot(
      NA,
      xlim = range(n_values),
      ylim = ylim,
      log = "x",
      xaxt = "n",
      xlab = "n",
      ylab = "factor score RMSE",
      main = main
    )
    axis(1, at = n_values, labels = n_values)
    grid(col = "#E9E9E9")
    box()

    for (p_value in p_values) {
      d <- d0[d0$p == p_value, , drop = FALSE]
      if (!nrow(d)) next
      d <- d[order(d$n), , drop = FALSE]
      col <- plot_cols[[as.character(p_value)]]
      lines(d$n, d$mean_factor_score_rmse, col = col, lwd = 2.3)
      points(d$n, d$mean_factor_score_rmse, col = col, pch = 16, cex = 0.95)

      y_lower <- pmax(0, d$lower_factor_score_rmse)
      y_upper <- pmin(ylim[2L], d$upper_factor_score_rmse)
      has_interval <- is.finite(y_lower) & is.finite(y_upper) &
        abs(y_upper - y_lower) > 1e-10
      if (any(has_interval)) {
        x <- d$n[has_interval]
        cap_width <- x * 0.035
        segments(x, y_lower[has_interval], x, y_upper[has_interval], col = col, lwd = 1.05)
        segments(x - cap_width, y_lower[has_interval], x + cap_width, y_lower[has_interval],
                 col = col, lwd = 1.05)
        segments(x - cap_width, y_upper[has_interval], x + cap_width, y_upper[has_interval],
                 col = col, lwd = 1.05)
      }
    }
  }

  mtext(
    sprintf("%s: product-mixture factor score RMSE across n", block_labels[[block_name]]),
    outer = TRUE,
    cex = 1.25,
    font = 2
  )
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend(
    x = 0.5,
    y = 0.04,
    xjust = 0.5,
    yjust = 0.5,
    horiz = TRUE,
    legend = paste0("p=", p_values),
    col = plot_cols[as.character(p_values)],
    lwd = 3.2,
    pch = 16,
    pt.cex = 1.25,
    bty = "o",
    bg = "white",
    box.col = "#555555",
    cex = 1.08,
    x.intersp = 0.85,
    seg.len = 2.1
  )
  invisible(TRUE)
}

all_summaries <- list()
for (block_name in names(result_files)) {
  results <- read.csv(result_files[[block_name]])
  summary <- summarize_factor_rmse(results, block_name)
  all_summaries[[block_name]] <- summary
  write.csv(
    summary,
    file.path(selected_table_dir, paste0(block_name, "_factor_score_rmse_by_n_product_by_p_summary.csv")),
    row.names = FALSE
  )
  plot_factor_rmse_by_n(
    summary,
    block_name,
    file.path(selected_plot_dir, paste0(block_name, "_factor_score_rmse_by_n_product_by_p.png"))
  )
}

combined <- do.call(rbind, all_summaries)
write.csv(
  combined,
  file.path(selected_table_dir, "combined_factor_score_rmse_by_n_product_by_p_summary.csv"),
  row.names = FALSE
)

cat("Wrote factor-score RMSE across-n/by-p plots and summaries to selected results.\n")
