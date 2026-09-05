#!/usr/bin/env Rscript

# Progress plots for the final signal-support simulation grid.
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
  final <- file.path(results_dir, "comparison_results.csv")
  if (file.exists(final)) return(read.csv(final, check.names = FALSE))
  chunk_dir <- file.path(results_dir, "chunks")
  task_dirs <- unique(dirname(list.files(
    chunk_dir,
    pattern = "comparison_results(_checkpoint)?\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )))
  if (!length(task_dirs)) return(data.frame())
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

run_label <- get_env("RUN_LABEL", "signal_support_grid_adaptive_product_workers")
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
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

results <- collect_results(results_dir)
if (!nrow(results)) stop("No result rows found in ", results_dir)

write.csv(
  results,
  file.path(table_dir, paste0(run_label, "_completed_results.csv")),
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
    ifelse(base$block_size_mode == "ifeval_like", "Unbalanced", "Balanced"),
    "blocks,",
    "strength", base$loading_strength,
    paste0("cross prob ", base$cross_loading_prob),
    paste0("H=", base$H_true),
    paste0("G=", base$G_true)
  )
  for (metric in names(metrics)) {
    if (!metric %in% names(d)) next
    plot_metric_lines(
      d,
      metric = metric,
      ylab = metrics[[metric]],
      out_file = file.path(plot_dir, paste0("lines_", metric, "_", safe_token(setting_name), ".png")),
      title = paste(title_base, "-", metrics[[metric]])
    )
    plot_metric_boxplots(
      d,
      metric = metric,
      ylab = metrics[[metric]],
      out_file = file.path(plot_dir, paste0("boxplot_", metric, "_", safe_token(setting_name), ".png")),
      title = paste(title_base, "-", metrics[[metric]])
    )
  }
}

cat("Wrote completed-results table to:\n")
cat(file.path(table_dir, paste0(run_label, "_completed_results.csv")), "\n")
cat("Wrote plots to:\n")
cat(plot_dir, "\n")
