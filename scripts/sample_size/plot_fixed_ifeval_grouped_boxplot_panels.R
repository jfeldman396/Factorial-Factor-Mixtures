#!/usr/bin/env Rscript

# Faceted boxplot panels for fixed IFEval-like simulation comparisons.
# Produces one PNG per metric with p as columns and H/G settings as rows.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

get_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

split_csv <- function(x) {
  if (!nzchar(x)) return(character(0))
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
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
  if (dir.exists(chunk_dir)) {
    task_dirs <- unique(dirname(list.files(
      chunk_dir,
      pattern = "comparison_results(_checkpoint)?\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )))
    if (length(task_dirs)) {
      files <- vapply(task_dirs, function(d) {
        f <- file.path(d, "comparison_results.csv")
        c <- file.path(d, "comparison_results_checkpoint.csv")
        if (file.exists(f)) f else c
      }, character(1L))
      files <- files[file.exists(files)]
      return(rbind_fill(lapply(files, function(path) {
        d <- read.csv(path, check.names = FALSE)
        d$source_chunk <- basename(dirname(path))
        d
      })))
    }
  }
  final <- file.path(results_dir, "comparison_results.csv")
  if (file.exists(final)) return(read.csv(final, check.names = FALSE))
  data.frame()
}

safe_token <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "value")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x)) y else x

method_label <- function(x) {
  labels <- c(
    independent_marginal_mixture = "Product MAP",
    viroli_laplace_gibbs = "Gibbs",
    viroli_gaussian_gibbs = "Gibbs (Gaussian)"
  )
  out <- labels[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  unname(out)
}

all_equal_components <- function(g_config, value) {
  vapply(strsplit(as.character(g_config), "-", fixed = TRUE), function(parts) {
    length(parts) > 0L && all(parts == as.character(value))
  }, logical(1L))
}

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
if (!dir.exists(repo_root)) repo_root <- getwd()

run_label <- get_env("RUN_LABEL", "fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10_canonical_comparison")
results_dir <- get_env("RESULTS_DIR", file.path(repo_root, "results", "full", run_label))
results_file <- get_env("RESULTS_FILE", file.path(results_dir, "comparison_results.csv"))
plot_dir <- get_env(
  "PLOT_DIR",
  file.path(repo_root, "results", "selected_plots", "sample_size", run_label, "grouped_boxplots")
)
table_dir <- get_env("TABLE_DIR", file.path(repo_root, "results", "selected_tables", "sample_size"))
method_filter <- split_csv(get_env("METHOD_FILTER", "independent_marginal_mixture,viroli_laplace_gibbs"))
g_component_filter <- get_env("G_COMPONENT_FILTER", "3")
output_tag <- get_env("OUTPUT_TAG", paste0("G", g_component_filter, "_product_vs_viroli_laplace"))

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

results <- if (file.exists(results_file)) {
  read.csv(results_file, check.names = FALSE)
} else {
  collect_results(results_dir)
}
if (!nrow(results)) {
  stop("No results found. Checked ", results_file, " and chunk files under ", results_dir)
}
if (length(method_filter)) {
  results <- results[results$method %in% method_filter, , drop = FALSE]
}
if (nzchar(g_component_filter)) {
  g_col <- if ("G_config" %in% names(results)) "G_config" else "G_true"
  results <- results[all_equal_components(results[[g_col]], g_component_filter), , drop = FALSE]
}
if (!nrow(results)) stop("No rows remain after filters.")

results$method_label <- factor(
  method_label(results$method),
  levels = unique(method_label(method_filter))
)
results$p_label <- factor(paste0("p=", results$p), levels = paste0("p=", sort(unique(results$p))))
results$n_label <- factor(paste0("n=", results$n), levels = paste0("n=", sort(unique(results$n))))
results$facet_label <- paste0("H=", results$H_true, ", G=", results$G_config)
facet_levels <- unique(results$facet_label[order(results$H_true, results$G_config)])
results$facet_label <- factor(results$facet_label, levels = facet_levels)

metrics <- c(
  factor_score_rmse = "Factor score RMSE",
  probability_rmse = "Probability RMSE",
  lambda_rmse = "Loading RMSE",
  alpha_rmse = "Intercept RMSE",
  marginal_mu_rmse = "Mixture mean RMSE",
  marginal_var_rmse = "Mixture variance RMSE",
  marginal_weight_rmse = "Mixture weight RMSE",
  seconds = "Runtime seconds"
)
metrics <- metrics[names(metrics) %in% names(results)]

summary_rows <- aggregate(
  results[names(metrics)],
  by = list(
    method = results$method_label,
    H_true = results$H_true,
    G_config = results$G_config,
    n = results$n,
    p = results$p
  ),
  FUN = function(z) mean(z, na.rm = TRUE)
)
write.csv(
  summary_rows,
  file.path(table_dir, paste0(run_label, "_", safe_token(output_tag), "_grouped_boxplot_cell_means.csv")),
  row.names = FALSE
)

plot_metric <- function(metric, label) {
  d <- results[is.finite(results[[metric]]), , drop = FALSE]
  if (!nrow(d)) return(invisible(FALSE))
  method_colors <- c(
    "Product MAP" = "#2b6db6",
    "Gibbs" = "#cf2f34",
    "Gibbs (Gaussian)" = "#2f9b57"
  )
  title <- paste0(label, " grouped by sample size")
  subtitle <- "Fixed IFEval-like DGP; p=1500/2000 have Product MAP only when Gibbs was not run"
  p <- ggplot(d, aes(x = n_label, y = .data[[metric]], fill = method_label, color = method_label)) +
    geom_boxplot(
      position = position_dodge(width = 0.78),
      width = 0.58,
      alpha = 0.58,
      outlier.alpha = 0.42,
      outlier.size = 1.15,
      linewidth = 0.45
    ) +
    facet_grid(facet_label ~ p_label, scales = "free_y") +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = label, fill = NULL, color = NULL) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = element_text(hjust = 0, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey86", color = "grey40"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  out_file <- file.path(plot_dir, paste0(safe_token(output_tag), "_grouped_boxplot_", metric, ".png"))
  ggsave(out_file, p, width = 15.8, height = 9.8, dpi = 180, bg = "white")
  out_file
}

written <- vapply(names(metrics), function(metric) plot_metric(metric, metrics[[metric]]), character(1L))
writeLines(c("Wrote grouped boxplots:", written))
writeLines(c("Wrote table:", file.path(table_dir, paste0(run_label, "_", safe_token(output_tag), "_grouped_boxplot_cell_means.csv"))))
