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
g_component_filter <- split_csv(get_env("G_COMPONENT_FILTER", "3"))
p_filter <- split_csv(get_env("P_FILTER", ""))
default_g_tag <- if (length(g_component_filter)) {
  paste0("G", paste(g_component_filter, collapse = "_G"))
} else {
  "all_G"
}
output_tag <- get_env("OUTPUT_TAG", paste0(default_g_tag, "_product_vs_viroli_laplace"))
plot_subtitle <- get_env("PLOT_SUBTITLE", "")
plot_width <- as.numeric(get_env("PLOT_WIDTH", "15.8"))
plot_height <- as.numeric(get_env("PLOT_HEIGHT", "9.8"))
plot_dpi <- as.integer(get_env("PLOT_DPI", "300"))
write_pdf <- as.logical(get_env("WRITE_PDF", "TRUE"))

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
if (length(g_component_filter)) {
  g_col <- if ("G_config" %in% names(results)) "G_config" else "G_true"
  keep_g <- Reduce(`|`, lapply(
    g_component_filter,
    function(g) all_equal_components(results[[g_col]], g)
  ))
  results <- results[keep_g, , drop = FALSE]
}
if (length(p_filter)) {
  results <- results[as.character(results$p) %in% p_filter, , drop = FALSE]
}
if (!nrow(results)) stop("No rows remain after filters.")

results$method_label <- factor(
  method_label(results$method),
  levels = unique(method_label(method_filter))
)
results$p_label <- factor(paste0("p=", results$p), levels = paste0("p=", sort(unique(results$p))))
results$n_label <- factor(paste0("n=", results$n), levels = paste0("n=", sort(unique(results$n))))
compact_g_label <- function(g_config) {
  vapply(strsplit(as.character(g_config), "-", fixed = TRUE), function(parts) {
    if (length(parts) > 0L && length(unique(parts)) == 1L) {
      paste0(parts[[1L]], "^", length(parts))
    } else {
      paste(parts, collapse = "-")
    }
  }, character(1L))
}
results$facet_label <- paste0("H=", results$H_true, ", G=", compact_g_label(results$G_config))
facet_frame <- unique(results[c("H_true", "G_config", "facet_label")])
facet_g <- suppressWarnings(as.integer(sub("-.*$", "", facet_frame$G_config)))
facet_frame <- facet_frame[order(-facet_frame$H_true, facet_g), , drop = FALSE]
facet_levels <- facet_frame$facet_label
results$facet_label <- factor(results$facet_label, levels = facet_levels)

metrics <- c(
  factor_score_rmse = "Factor score RMSE",
  probability_rmse = "Probability RMSE",
  lambda_rmse = "Loading RMSE",
  alpha_rmse = "Intercept RMSE",
  marginal_mu_rmse = "Mixture mean RMSE",
  marginal_var_rmse = "Mixture variance RMSE",
  marginal_weight_rmse = "Mixture weight RMSE",
  component_profile_hamming_accuracy = "Component-profile Hamming accuracy",
  component_profile_exact_accuracy = "Exact component-profile accuracy",
  mean_component_ari = "Mean component ARI",
  min_component_ari = "Minimum component ARI",
  mean_true_component_probability = "Mean probability on true component",
  component_brier_score = "Component Brier score",
  component_log_loss = "Component log loss",
  mean_component_entropy = "Mean component entropy",
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
  compared_p <- sort(unique(d$p))
  subtitle <- if (nzchar(plot_subtitle)) {
    plot_subtitle
  } else if (all(compared_p %in% c(500, 1000))) {
    "Fixed IFEval-like DGP; matched Product MAP and Gibbs cells"
  } else {
    "Fixed IFEval-like DGP; p=1500/2000 have Product MAP only when Gibbs was not run"
  }
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
      panel.spacing = grid::unit(0.10, "lines"),
      axis.title.y = element_text(face = "bold"),
      legend.position = "bottom"
    )
  file_stem <- file.path(plot_dir, paste0(safe_token(output_tag), "_grouped_boxplot_", metric))
  png_file <- paste0(file_stem, ".png")
  pdf_file <- paste0(file_stem, ".pdf")
  ggsave(png_file, p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  if (isTRUE(write_pdf)) {
    ggsave(
      pdf_file,
      p,
      width = plot_width,
      height = plot_height,
      device = grDevices::pdf,
      bg = "white",
      useDingbats = FALSE
    )
  }
  paste(c(png_file, if (isTRUE(write_pdf)) pdf_file), collapse = "; ")
}

written <- vapply(names(metrics), function(metric) plot_metric(metric, metrics[[metric]]), character(1L))
writeLines(c("Wrote grouped boxplots:", written))
writeLines(c("Wrote table:", file.path(table_dir, paste0(run_label, "_", safe_token(output_tag), "_grouped_boxplot_cell_means.csv"))))
