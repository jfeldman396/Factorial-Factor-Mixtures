#!/usr/bin/env Rscript

# Progress plots for the fixed-DGP IFEval-like loading simulation grid.
#
# The script can be run while the simulation is still in progress.  It reads
# the combined results file when available and otherwise collects per-task
# checkpoint files from results/full/<run_label>/chunks.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

split_csv <- function(x) {
  if (!nzchar(x)) return(character(0))
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

safe_token <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "value")
}

rbind_fill <- function(x) {
  x <- x[!vapply(x, is.null, logical(1L))]
  if (!length(x)) return(data.frame())
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) d[[nm]] <- NA
    d[, all_names, drop = FALSE]
  })
  do.call(rbind, x)
}

collect_results <- function(results_dir) {
  chunk_dir <- file.path(results_dir, "chunks")
  task_dirs <- unique(dirname(list.files(
    chunk_dir,
    pattern = "comparison_results(_checkpoint)?\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )))
  if (!length(task_dirs)) {
    final <- file.path(results_dir, "comparison_results.csv")
    if (file.exists(final)) return(read.csv(final, check.names = FALSE))
    return(data.frame())
  }
  files <- vapply(task_dirs, function(d) {
    f <- file.path(d, "comparison_results.csv")
    c <- file.path(d, "comparison_results_checkpoint.csv")
    if (file.exists(f)) f else c
  }, character(1L))
  files <- files[file.exists(files)]
  rbind_fill(lapply(files, function(path) {
    d <- read.csv(path, check.names = FALSE)
    d$source_chunk <- basename(dirname(path))
    d
  }))
}

method_label <- function(x) {
  labels <- c(
    independent_marginal_mixture = "Product MAP",
    viroli_laplace_gibbs = "Viroli Laplace",
    viroli_gaussian_gibbs = "Viroli Gaussian"
  )
  out <- labels[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  unname(out)
}

plot_metric_lines <- function(d, metric, ylab, out_file, title) {
  d <- d[is.finite(d[[metric]]) & is.finite(d$n) & is.finite(d$p), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  methods <- unique(d$method)
  p_values <- sort(unique(d$p))
  colors <- c(
    independent_marginal_mixture = "#2b6db6",
    viroli_laplace_gibbs = "#cf2f34",
    viroli_gaussian_gibbs = "#2f9b57"
  )
  line_types <- setNames(seq_along(p_values), as.character(p_values))
  ylim <- range(d[[metric]], na.rm = TRUE)
  if (!all(is.finite(ylim))) return(invisible(FALSE))
  pad <- diff(ylim)
  if (!is.finite(pad) || pad <= 0) pad <- max(abs(ylim), 1)
  ylim <- ylim + c(-0.08, 0.12) * pad

  png(out_file, width = 1900, height = 1200, res = 160)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(4.5, 5, 4, 1.5), xpd = NA)
  plot(
    NA,
    xlim = range(d$n),
    ylim = ylim,
    xlab = "n",
    ylab = ylab,
    main = title,
    xaxt = "n",
    bty = "l"
  )
  axis(1, at = sort(unique(d$n)))
  grid(col = "#e6e6e6", lty = 3)

  for (method in methods) {
    for (p in p_values) {
      sub <- d[d$method == method & d$p == p, , drop = FALSE]
      if (!nrow(sub)) next
      means <- tapply(sub[[metric]], sub$n, mean, na.rm = TRUE)
      sds <- tapply(sub[[metric]], sub$n, sd, na.rm = TRUE)
      n_rep <- tapply(is.finite(sub[[metric]]), sub$n, sum)
      agg <- data.frame(
        n = as.numeric(names(means)),
        mean = as.numeric(means),
        sd = as.numeric(sds[names(means)]),
        n_rep = as.integer(n_rep[names(means)])
      )
      agg <- agg[order(agg$n), , drop = FALSE]
      col <- if (method %in% names(colors)) colors[[method]] else "black"
      lty <- line_types[[as.character(p)]]
      lines(agg$n, agg$mean, type = "b", pch = 19, lwd = 2.2, col = col, lty = lty)
      if (any(is.finite(agg$sd)) && any(agg$n_rep > 1L)) {
        se <- 2 * agg$sd
        arrows(
          agg$n,
          agg$mean - se,
          agg$n,
          agg$mean + se,
          angle = 90,
          code = 3,
          length = 0.04,
          col = col,
          lwd = 1.2
        )
      }
    }
  }

  legend(
    "topright",
    legend = method_label(methods),
    col = colors[methods],
    lty = 1,
    pch = 19,
    lwd = 2.4,
    bty = "n",
    cex = 0.95
  )
  legend(
    "bottomright",
    legend = paste0("p=", p_values),
    col = "#333333",
    lty = line_types,
    lwd = 2.4,
    bty = "n",
    cex = 0.95
  )
  invisible(TRUE)
}

plot_metric_boxplots <- function(d, metric, ylab, out_file, title) {
  d <- d[is.finite(d[[metric]]) & is.finite(d$n) & is.finite(d$p), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  d$method_label <- method_label(d$method)
  d$cell <- paste0("n=", d$n, "\np=", d$p)
  cells <- unique(d$cell[order(d$n, d$p)])
  methods <- unique(d$method_label)
  colors <- c("Product MAP" = "#2b6db6", "Viroli Laplace" = "#cf2f34", "Viroli Gaussian" = "#2f9b57")

  png(out_file, width = 2200, height = 1200, res = 160)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(6.5, 5, 4, 1.5), xpd = NA)
  vals <- split(d[[metric]], interaction(d$cell, d$method_label, drop = TRUE))
  at <- numeric(0)
  box_vals <- list()
  box_cols <- character(0)
  labels_at <- numeric(0)
  labels <- character(0)
  pos <- 1
  for (cell in cells) {
    cell_positions <- numeric(0)
    for (method in methods) {
      key <- paste(cell, method, sep = ".")
      if (!key %in% names(vals)) next
      box_vals[[length(box_vals) + 1L]] <- vals[[key]]
      at <- c(at, pos)
      cell_positions <- c(cell_positions, pos)
      box_cols <- c(box_cols, if (method %in% names(colors)) colors[[method]] else "#777777")
      pos <- pos + 1
    }
    if (length(cell_positions)) {
      labels_at <- c(labels_at, mean(cell_positions))
      labels <- c(labels, cell)
      pos <- pos + 1
    }
  }
  if (!length(box_vals)) return(invisible(FALSE))
  boxplot(
    box_vals,
    at = at,
    xaxt = "n",
    col = adjustcolor(box_cols, alpha.f = 0.55),
    border = box_cols,
    outline = FALSE,
    ylab = ylab,
    main = title,
    bty = "l"
  )
  axis(1, at = labels_at, labels = labels, las = 2, cex.axis = 0.85)
  grid(nx = NA, ny = NULL, col = "#e6e6e6", lty = 3)
  legend(
    "topright",
    legend = methods,
    fill = adjustcolor(colors[methods], alpha.f = 0.55),
    border = colors[methods],
    bty = "n",
    cex = 0.95
  )
  invisible(TRUE)
}

summarize_by_cell <- function(d, metrics) {
  group_cols <- intersect(
    c("method", "H_true", "G_true", "n", "p", "block_size_mode",
      "loading_strength", "cross_loading_prob"),
    names(d)
  )
  key <- do.call(interaction, c(
    as.data.frame(lapply(d[, group_cols, drop = FALSE], as.character), stringsAsFactors = FALSE),
    list(drop = TRUE, sep = "\r")
  ))
  pieces <- lapply(split(seq_len(nrow(d)), key), function(idx) {
    sub <- d[idx, , drop = FALSE]
    out <- sub[1L, group_cols, drop = FALSE]
    out$n_result_rows <- nrow(sub)
    for (metric in names(metrics)) {
      if (!metric %in% names(sub)) next
      values <- sub[[metric]]
      values <- values[is.finite(values)]
      out[[paste0(metric, "_mean")]] <- if (length(values)) mean(values) else NA_real_
      out[[paste0(metric, "_sd")]] <- if (length(values) > 1L) sd(values) else NA_real_
      out[[paste0(metric, "_n")]] <- length(values)
    }
    out
  })
  out <- rbind_fill(pieces)
  out[order(out$method, out$H_true, out$G_true, out$n, out$p), , drop = FALSE]
}

make_vs_gibbs_table <- function(cell_summary, metrics) {
  product_method <- "independent_marginal_mixture"
  gibbs_methods <- intersect(c("viroli_laplace_gibbs", "viroli_gaussian_gibbs"), unique(cell_summary$method))
  if (!length(gibbs_methods)) return(data.frame())

  match_cols <- intersect(
    c("H_true", "G_true", "n", "p", "block_size_mode", "loading_strength", "cross_loading_prob"),
    names(cell_summary)
  )
  product <- cell_summary[cell_summary$method == product_method, , drop = FALSE]
  out <- list()
  for (gibbs_method in gibbs_methods) {
    gibbs <- cell_summary[cell_summary$method == gibbs_method, , drop = FALSE]
    if (!nrow(product) || !nrow(gibbs)) next
    pkey <- do.call(paste, c(product[, match_cols, drop = FALSE], sep = "\r"))
    gkey <- do.call(paste, c(gibbs[, match_cols, drop = FALSE], sep = "\r"))
    common <- intersect(pkey, gkey)
    if (!length(common)) next
    for (key in common) {
      prow <- product[pkey == key, , drop = FALSE][1L, , drop = FALSE]
      grow <- gibbs[gkey == key, , drop = FALSE][1L, , drop = FALSE]
      row <- prow[, match_cols, drop = FALSE]
      row$gibbs_method <- gibbs_method
      row$product_rows <- prow$n_result_rows
      row$gibbs_rows <- grow$n_result_rows
      for (metric in names(metrics)) {
        mean_col <- paste0(metric, "_mean")
        if (!mean_col %in% names(prow) || !mean_col %in% names(grow)) next
        row[[paste0(metric, "_product")]] <- prow[[mean_col]]
        row[[paste0(metric, "_gibbs")]] <- grow[[mean_col]]
        row[[paste0(metric, "_product_minus_gibbs")]] <- prow[[mean_col]] - grow[[mean_col]]
        if (metric == "seconds" && is.finite(prow[[mean_col]]) && prow[[mean_col]] > 0) {
          row$seconds_gibbs_div_product <- grow[[mean_col]] / prow[[mean_col]]
        }
      }
      out[[length(out) + 1L]] <- row
    }
  }
  rbind_fill(out)
}

matched_method_rows <- function(results, product_method, comparator_method, metrics) {
  if (!"method" %in% names(results)) return(data.frame())
  product <- results[as.character(results$method) == product_method, , drop = FALSE]
  comparator <- results[as.character(results$method) == comparator_method, , drop = FALSE]
  if (!nrow(product) || !nrow(comparator)) return(data.frame())

  candidate_match_cols <- c(
    "scenario", "rep", "n", "p", "H_true", "G_true", "G_config",
    "separation", "mixture_param_mode", "mixture_variance_mode",
    "intercept_mode", "loading_design", "block_size_mode",
    "loading_sign_mode", "alignment_mode", "loading_strength",
    "primary_loading_min", "primary_loading_max",
    "cross_loading_min", "cross_loading_max",
    "cross_loading_prob", "cross_sign_mode",
    "fix_dgp_parameters", "dgp_parameter_seed", "loading_parameter_seed",
    "mixture_parameter_seed", "data_seed", "dgp_p_max"
  )
  match_cols <- intersect(candidate_match_cols, names(results))
  metric_cols <- intersect(names(metrics), names(results))
  if (!length(match_cols) || !length(metric_cols)) return(data.frame())

  make_key <- function(d) {
    parts <- lapply(d[, match_cols, drop = FALSE], function(x) {
      x <- as.character(x)
      x[is.na(x)] <- "<NA>"
      x
    })
    do.call(paste, c(parts, sep = "\r"))
  }

  product$.match_key <- make_key(product)
  comparator$.match_key <- make_key(comparator)
  common <- intersect(unique(product$.match_key), unique(comparator$.match_key))
  if (!length(common)) return(data.frame())

  out <- vector("list", length(common))
  for (i in seq_along(common)) {
    key <- common[[i]]
    prow <- product[product$.match_key == key, , drop = FALSE][1L, , drop = FALSE]
    crow <- comparator[comparator$.match_key == key, , drop = FALSE][1L, , drop = FALSE]
    row <- prow[, match_cols, drop = FALSE]
    row$product_method <- product_method
    row$comparator_method <- comparator_method
    for (metric in metric_cols) {
      row[[paste0(metric, "_product")]] <- prow[[metric]]
      row[[paste0(metric, "_comparator")]] <- crow[[metric]]
      row[[paste0(metric, "_product_minus_comparator")]] <- prow[[metric]] - crow[[metric]]
    }
    out[[i]] <- row
  }
  rbind_fill(out)
}

plot_matched_product_vs_viroli_laplace <- function(results, metrics, out_file, table_file = NULL) {
  matched <- matched_method_rows(
    results = results,
    product_method = "independent_marginal_mixture",
    comparator_method = "viroli_laplace_gibbs",
    metrics = metrics
  )
  if (!nrow(matched)) return(invisible(FALSE))
  if (!is.null(table_file)) write.csv(matched, table_file, row.names = FALSE)

  plot_metrics <- c(
    factor_score_rmse = "factor score RMSE",
    lambda_rmse = "loading RMSE",
    alpha_rmse = "intercept RMSE",
    marginal_mu_rmse = "mixture mean RMSE",
    marginal_var_rmse = "mixture variance RMSE",
    marginal_weight_rmse = "mixture weight RMSE",
    seconds = "runtime seconds"
  )
  plot_metrics <- plot_metrics[
    paste0(names(plot_metrics), "_product") %in% names(matched) &
      paste0(names(plot_metrics), "_comparator") %in% names(matched)
  ]
  if (!length(plot_metrics)) return(invisible(FALSE))

  h_col <- intersect(c("H_true", "H"), names(matched))[1L]
  p_col <- intersect(c("p", "P"), names(matched))[1L]
  n_col <- intersect(c("n", "N"), names(matched))[1L]
  if (is.na(h_col) || is.na(p_col) || is.na(n_col)) return(invisible(FALSE))

  h_values <- sort(unique(matched[[h_col]][is.finite(matched[[h_col]])]))
  p_values <- sort(unique(matched[[p_col]][is.finite(matched[[p_col]])]))
  n_values <- sort(unique(matched[[n_col]][is.finite(matched[[n_col]])]))
  h_colors <- setNames(
    rep_len(c("#2b6db6", "#cf2f34", "#2f9b57", "#7b4cc2", "#d99019"), length(h_values)),
    as.character(h_values)
  )
  p_shapes <- setNames(rep_len(c(16, 17, 15, 18, 8, 3), length(p_values)), as.character(p_values))

  size_for_n <- function(x) {
    x <- as.numeric(x)
    if (length(n_values) <= 1L || diff(range(n_values)) <= 0) return(rep(1.05, length(x)))
    0.65 + 1.15 * (x - min(n_values)) / diff(range(n_values))
  }

  png(out_file, width = 2400, height = 1500, res = 170)
  on.exit(dev.off(), add = TRUE)
  layout(
    matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 2L, byrow = TRUE),
    widths = c(1, 1, 1, 1.02)
  )
  op <- par(mar = c(4.7, 4.9, 3.0, 1.1), oma = c(0, 0, 3.1, 0), xpd = FALSE)
  on.exit(par(op), add = TRUE)

  for (metric in names(plot_metrics)) {
    x <- matched[[paste0(metric, "_product")]]
    y <- matched[[paste0(metric, "_comparator")]]
    ok <- is.finite(x) & is.finite(y)
    if (!any(ok)) {
      plot.new()
      next
    }
    lim <- range(c(x[ok], y[ok]), na.rm = TRUE)
    pad <- diff(lim)
    if (!is.finite(pad) || pad <= 0) pad <- max(abs(lim), 1)
    lim <- lim + c(-0.08, 0.08) * pad
    if (lim[1L] > 0) lim[1L] <- max(0, lim[1L])

    plot(
      NA,
      xlim = lim,
      ylim = lim,
      xlab = "Product MAP",
      ylab = "Viroli Laplace",
      main = unname(plot_metrics[[metric]]),
      bty = "l"
    )
    grid(col = "#e6e6e6", lty = 3)
    abline(0, 1, col = "#777777", lty = 2, lwd = 1.3)
    for (i in which(ok)) {
      h_key <- as.character(matched[[h_col]][i])
      p_key <- as.character(matched[[p_col]][i])
      points(
        x[i],
        y[i],
        pch = p_shapes[[p_key]],
        cex = size_for_n(matched[[n_col]][i]),
        col = adjustcolor(h_colors[[h_key]], alpha.f = 0.48),
        bg = adjustcolor(h_colors[[h_key]], alpha.f = 0.48)
      )
    }
  }

  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend(
    "center",
    legend = c(paste0("H=", h_values), paste0("p=", p_values), paste0("n=", n_values), "diagonal"),
    col = c(
      h_colors[as.character(h_values)],
      rep("#111111", length(p_values)),
      rep("#555555", length(n_values)),
      "#777777"
    ),
    pch = c(
      rep(16, length(h_values)),
      p_shapes[as.character(p_values)],
      rep(16, length(n_values)),
      NA
    ),
    pt.cex = c(
      rep(1.2, length(h_values)),
      rep(1.2, length(p_values)),
      size_for_n(n_values),
      NA
    ),
    lty = c(
      rep(NA, length(h_values) + length(p_values) + length(n_values)),
      2
    ),
    lwd = c(
      rep(NA, length(h_values) + length(p_values) + length(n_values)),
      1.6
    ),
    bty = "n",
    cex = 1.0,
    y.intersp = 1.25
  )
  mtext("Matched-rep Product MAP vs Viroli Laplace Gibbs comparisons", outer = TRUE, cex = 1.25, font = 2)
  invisible(TRUE)
}

plot_across_p_recovery_panel <- function(d, metrics, out_file, title) {
  keep_metrics <- intersect(
    c("factor_score_rmse", "lambda_rmse", "alpha_rmse",
      "marginal_mu_rmse", "marginal_var_rmse", "marginal_weight_rmse"),
    names(metrics)
  )
  keep_metrics <- keep_metrics[keep_metrics %in% names(d)]
  if (!length(keep_metrics)) return(invisible(FALSE))
  d <- d[is.finite(d$p) & is.finite(d$n), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))

  colors <- c(
    independent_marginal_mixture = "#2b6db6",
    viroli_laplace_gibbs = "#cf2f34",
    viroli_gaussian_gibbs = "#2f9b57"
  )
  methods <- intersect(names(colors), unique(as.character(d$method)))
  methods <- c(methods, setdiff(unique(as.character(d$method)), methods))
  n_values <- sort(unique(d$n))
  line_types <- setNames(seq_along(n_values), as.character(n_values))
  p_values <- sort(unique(d$p))

  png(out_file, width = 2600, height = 1600, res = 170)
  on.exit(dev.off(), add = TRUE)
  layout(rbind(matrix(seq_len(6L), 2L, 3L, byrow = TRUE), rep(7L, 3L)), heights = c(1, 1, 0.16))
  op <- par(mar = c(4.4, 4.9, 3, 1.2), oma = c(0, 0, 3.2, 0), xpd = FALSE)
  on.exit(par(op), add = TRUE)
  for (metric in keep_metrics) {
    values <- d[[metric]]
    ylim <- range(values[is.finite(values)], na.rm = TRUE)
    if (!all(is.finite(ylim))) ylim <- c(0, 1)
    if (ylim[1L] > 0) ylim[1L] <- 0
    pad <- diff(ylim)
    if (!is.finite(pad) || pad <= 0) pad <- max(abs(ylim), 1)
    ylim[2L] <- ylim[2L] + 0.16 * pad
    plot(
      NA,
      xlim = range(p_values),
      ylim = ylim,
      xlab = "p",
      ylab = "RMSE",
      main = metrics[[metric]],
      xaxt = "n",
      bty = "l"
    )
    axis(1, at = p_values)
    grid(col = "#e6e6e6", lty = 3)
    for (method in methods) {
      for (n_value in n_values) {
        sub <- d[d$method == method & d$n == n_value, , drop = FALSE]
        sub <- sub[is.finite(sub[[metric]]), , drop = FALSE]
        if (!nrow(sub)) next
        means <- tapply(sub[[metric]], sub$p, mean, na.rm = TRUE)
        sds <- tapply(sub[[metric]], sub$p, sd, na.rm = TRUE)
        n_rep <- tapply(is.finite(sub[[metric]]), sub$p, sum)
        agg <- data.frame(
          p = as.numeric(names(means)),
          mean = as.numeric(means),
          sd = as.numeric(sds[names(means)]),
          n_rep = as.integer(n_rep[names(means)])
        )
        agg <- agg[order(agg$p), , drop = FALSE]
        col <- if (method %in% names(colors)) colors[[method]] else "#555555"
        lty <- line_types[[as.character(n_value)]]
        lines(agg$p, agg$mean, type = "b", pch = 19, lwd = 2.1, col = col, lty = lty)
        if (any(agg$n_rep > 1L) && any(is.finite(agg$sd))) {
          arrows(
            agg$p,
            agg$mean - 2 * agg$sd,
            agg$p,
            agg$mean + 2 * agg$sd,
            angle = 90,
            code = 3,
            length = 0.035,
            col = col,
            lwd = 1.1
          )
        }
      }
    }
  }
  mtext(title, outer = TRUE, cex = 1.15, font = 2)
  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  legend(
    x = 0.5,
    y = 0.68,
    legend = method_label(methods),
    col = colors[methods],
    lty = 1,
    pch = 19,
    lwd = 2.5,
    horiz = TRUE,
    bty = "n",
    cex = 1.0,
    xjust = 0.5,
    yjust = 0.5
  )
  legend(
    x = 0.5,
    y = 0.28,
    legend = paste0("n=", n_values),
    col = "#222222",
    lty = line_types,
    lwd = 2.5,
    horiz = TRUE,
    bty = "n",
    cex = 1.0,
    xjust = 0.5,
    yjust = 0.5
  )
  invisible(TRUE)
}

run_label <- get_env("RUN_LABEL", "fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10")
results_dir <- get_env(
  "RESULTS_DIR",
  file.path(repo_root, "results", "full", run_label)
)
plot_dir <- get_env(
  "PLOT_DIR",
  file.path(repo_root, "results", "selected_plots", "sample_size", run_label)
)
table_dir <- get_env(
  "TABLE_DIR",
  file.path(repo_root, "results", "selected_tables", "sample_size")
)
method_filter <- split_csv(get_env("METHOD_FILTER", ""))
output_tag <- get_env("OUTPUT_TAG", "")
if (!nzchar(output_tag) && identical(method_filter, "independent_marginal_mixture")) {
  output_tag <- "product_map_only"
}
output_prefix <- if (nzchar(output_tag)) paste0(safe_token(output_tag), "_") else ""
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

results <- collect_results(results_dir)
if (!nrow(results)) stop("No result rows found in ", results_dir)
if (length(method_filter)) {
  results <- results[results$method %in% method_filter, , drop = FALSE]
  if (!nrow(results)) {
    stop("No rows found after METHOD_FILTER=", paste(method_filter, collapse = ","), call. = FALSE)
  }
}

write.csv(
  results,
  file.path(table_dir, paste0(run_label, "_", output_prefix, "completed_results.csv")),
  row.names = FALSE
)

group_cols <- intersect(
  c("block_size_mode", "loading_strength", "cross_loading_prob", "H_true", "G_true"),
  names(results)
)
setting_key <- do.call(interaction, c(
  as.data.frame(lapply(results[, group_cols, drop = FALSE], as.character), stringsAsFactors = FALSE),
  list(drop = TRUE, sep = " | ")
))

metrics <- list(
  factor_score_rmse = "factor score RMSE",
  lambda_rmse = "loading RMSE",
  alpha_rmse = "intercept RMSE",
  marginal_mu_rmse = "mixture mean RMSE",
  marginal_var_rmse = "mixture variance RMSE",
  marginal_weight_rmse = "mixture weight RMSE",
  seconds = "seconds"
)

cell_summary <- summarize_by_cell(results, metrics)
write.csv(
  cell_summary,
  file.path(table_dir, paste0(run_label, "_", output_prefix, "cell_summary.csv")),
  row.names = FALSE
)

vs_gibbs <- make_vs_gibbs_table(cell_summary, metrics)
if (nrow(vs_gibbs)) {
  write.csv(
    vs_gibbs,
    file.path(table_dir, paste0(run_label, "_", output_prefix, "product_vs_gibbs_overlap.csv")),
    row.names = FALSE
  )
}

plot_matched_product_vs_viroli_laplace(
  results = results,
  metrics = metrics,
  out_file = file.path(plot_dir, paste0(output_prefix, "matched_product_vs_viroli_laplace_scatter.png")),
  table_file = file.path(table_dir, paste0(run_label, "_", output_prefix, "product_vs_viroli_laplace_matched_reps.csv"))
)

for (key in levels(setting_key)) {
  d <- results[setting_key == key, , drop = FALSE]
  if (!nrow(d)) next
  base <- d[1L, group_cols, drop = FALSE]
  setting_name <- paste(
    paste0("block=", base$block_size_mode),
    paste0("strength=", base$loading_strength),
    paste0("cp=", base$cross_loading_prob),
    paste0("H=", base$H_true),
    paste0("G=", base$G_true),
    sep = "_"
  )
  title_base <- paste(
    ifelse(base$block_size_mode %in% c("ifeval_like", "ifeval_min30"), "Unbalanced", "Balanced"),
    "blocks,",
    "strength", base$loading_strength,
    paste0("cross prob ", base$cross_loading_prob),
    paste0("H=", base$H_true),
    paste0("G=", base$G_true)
  )
  plot_across_p_recovery_panel(
    d,
    metrics = metrics,
    out_file = file.path(plot_dir, paste0(output_prefix, "across_p_recovery_panel_", safe_token(setting_name), ".png")),
    title = paste(title_base, "- recovery across p")
  )
  for (metric in names(metrics)) {
    if (!metric %in% names(d)) next
    plot_metric_lines(
      d,
      metric = metric,
      ylab = metrics[[metric]],
      out_file = file.path(plot_dir, paste0(output_prefix, "lines_", metric, "_", safe_token(setting_name), ".png")),
      title = paste(title_base, "-", metrics[[metric]])
    )
    plot_metric_boxplots(
      d,
      metric = metric,
      ylab = metrics[[metric]],
      out_file = file.path(plot_dir, paste0(output_prefix, "boxplot_", metric, "_", safe_token(setting_name), ".png")),
      title = paste(title_base, "-", metrics[[metric]])
    )
  }
}

cat("Wrote completed-results table to:\n")
cat(file.path(table_dir, paste0(run_label, "_", output_prefix, "completed_results.csv")), "\n")
cat("Wrote cell summary table to:\n")
cat(file.path(table_dir, paste0(run_label, "_", output_prefix, "cell_summary.csv")), "\n")
if (nrow(vs_gibbs)) {
  cat("Wrote Product MAP vs Gibbs overlap table to:\n")
  cat(file.path(table_dir, paste0(run_label, "_", output_prefix, "product_vs_gibbs_overlap.csv")), "\n")
} else {
  cat("Product MAP vs Gibbs overlap table not written yet; no overlapping Gibbs rows found.\n")
}
cat("Wrote plots to:\n")
cat(plot_dir, "\n")
