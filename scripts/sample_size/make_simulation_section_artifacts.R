#!/usr/bin/env Rscript

# Build the paper-facing figures and LaTeX tables for the simulation section.
# The matched p=500/1000 comparisons come from the completed 25-rep subtype
# rerun. Product MAP p=1500/2000 rows are appended from the completed recovery
# study. Runtime summaries use the original end-to-end recovery run so subtype
# post-processing is not included in fit time.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1L]])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

plot_root <- file.path(
  repo_root, "results", "selected_plots", "sample_size",
  "paper_simulation_section"
)
table_root <- file.path(
  repo_root, "results", "selected_tables", "sample_size",
  "paper_simulation_section"
)
supp_root <- file.path(plot_root, "supplement")
dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
dir.create(table_root, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_root, recursive = TRUE, showWarnings = FALSE)

arms <- c("separated", "asymmetric_pi", "moderate_overlap")
arm_labels <- c(
  separated = "Separated symmetric",
  asymmetric_pi = "Asymmetric weights",
  moderate_overlap = "Moderate overlap"
)
product_method <- "independent_marginal_mixture"
gibbs_method <- "viroli_laplace_gibbs"

rbind_fill <- function(xs) {
  xs <- xs[vapply(xs, nrow, integer(1L)) > 0L]
  all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    for (nm in setdiff(all_names, names(x))) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, xs)
}

read_results <- function(path) {
  if (!file.exists(path)) stop("Missing results file: ", path)
  read.csv(path, check.names = FALSE)
}

arm_paths <- function(arm) {
  saved_root <- Sys.getenv(
    "SIMULATION_RESULTS_DIR",
    file.path(repo_root, "results", "saved", "simulation")
  )
  list(
    analysis = file.path(saved_root, arm, "analysis_results.csv"),
    recovery = file.path(saved_root, arm, "recovery_results.csv"),
    subtype = file.path(saved_root, arm, "subtype_results.csv")
  )
}

assemble_arm <- function(arm) {
  paths <- arm_paths(arm)
  out <- read_results(paths$analysis)
  out$mixture_setting <- arm
  out
}

results <- do.call(rbind, lapply(arms, assemble_arm))
profile_results <- do.call(rbind, lapply(arms, function(arm) {
  d <- read_results(arm_paths(arm)$subtype)
  d$mixture_setting <- arm
  d
}))

method_label <- c(
  independent_marginal_mixture = "Spectral-MAP",
  viroli_laplace_gibbs = "Gibbs"
)
method_colors <- c("Spectral-MAP" = "#2C6DB2", "Gibbs" = "#D43D44")
method_fills <- c("Spectral-MAP" = "#7FA6D2", "Gibbs" = "#E68185")

compact_g <- function(g_config) {
  vapply(strsplit(as.character(g_config), "-", fixed = TRUE), function(x) {
    if (length(unique(x)) == 1L) as.integer(x[[1L]]) else NA_integer_
  }, integer(1L))
}

decorate <- function(d) {
  d$method_label <- factor(
    unname(method_label[d$method]),
    levels = c("Spectral-MAP", "Gibbs")
  )
  d$G_number <- compact_g(d$G_config)
  d$design <- factor(
    paste0("H = ", d$H_true, ", G = ", d$G_number),
    levels = c("H = 5, G = 2", "H = 5, G = 3", "H = 10, G = 2", "H = 10, G = 3")
  )
  d$p_label <- factor(
    paste0("p = ", d$p),
    levels = paste0("p = ", c(500, 1000, 1500, 2000))
  )
  d$n_label <- factor(
    paste0("n = ", d$n),
    levels = paste0("n = ", c(100, 200, 400))
  )
  d$arm_label <- factor(
    unname(arm_labels[d$mixture_setting]),
    levels = unname(arm_labels)
  )
  d
}

results <- decorate(results)
profile_results <- decorate(profile_results)

paper_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1),
      strip.background = element_rect(fill = "grey92", color = "grey45", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
      panel.spacing = unit(0.10, "lines"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(7, 8, 6, 7)
    )
}

save_plot <- function(plot, stem, width, height, dpi = 400) {
  png_path <- paste0(stem, ".png")
  pdf_path <- paste0(stem, ".pdf")
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(
    pdf_path, plot, width = width, height = height,
    device = grDevices::pdf, bg = "white", useDingbats = FALSE
  )
  c(png_path, pdf_path)
}

metric_labels <- c(
  factor_score_rmse = "Factor score RMSE",
  probability_rmse = "Response-probability RMSE",
  lambda_rmse = "Loading RMSE",
  alpha_rmse = "Intercept RMSE",
  marginal_mu_rmse = "Mixture mean RMSE",
  marginal_var_rmse = "Mixture variance RMSE",
  marginal_weight_rmse = "Mixture weight RMSE",
  component_profile_hamming_accuracy = "Hamming accuracy",
  component_profile_exact_accuracy = "Exact-profile accuracy",
  mean_component_ari = "Mean component ARI",
  seconds = "Runtime (seconds)"
)

boxplot_metric <- function(d, metric, y_label, include_all_p = TRUE) {
  keep_p <- if (include_all_p) c(500, 1000, 1500, 2000) else c(500, 1000)
  d <- d[d$p %in% keep_p & is.finite(d[[metric]]), , drop = FALSE]
  ggplot(
    d,
    aes(x = n_label, y = .data[[metric]], fill = method_label, color = method_label)
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.76), width = 0.56,
      alpha = 0.78, linewidth = 0.42,
      outlier.alpha = 0.42, outlier.size = 0.95
    ) +
    facet_grid(design ~ p_label, scales = "free_y") +
    scale_fill_manual(values = method_fills, drop = FALSE) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    labs(x = NULL, y = y_label, fill = NULL, color = NULL) +
    paper_theme(13.5) +
    theme(
      legend.key.width = unit(1.05, "lines"),
      legend.text = element_text(size = 12.2),
      axis.text = element_text(size = 11.7),
      axis.title.y = element_text(size = 13.2),
      strip.text = element_text(size = 12.2)
    )
}

core_metrics <- c(
  "factor_score_rmse", "probability_rmse", "lambda_rmse", "alpha_rmse",
  "marginal_mu_rmse", "marginal_var_rmse", "marginal_weight_rmse"
)
profile_metrics <- c(
  "component_profile_hamming_accuracy",
  "component_profile_exact_accuracy",
  "mean_component_ari"
)

written_plots <- character(0)
for (arm in arms) {
  arm_dir <- file.path(supp_root, arm)
  dir.create(arm_dir, recursive = TRUE, showWarnings = FALSE)
  d <- results[results$mixture_setting == arm, , drop = FALSE]
  dp <- profile_results[profile_results$mixture_setting == arm, , drop = FALSE]

  for (metric in core_metrics) {
    p <- boxplot_metric(d, metric, metric_labels[[metric]], include_all_p = TRUE)
    files <- save_plot(
      p,
      file.path(arm_dir, paste0("fig_", arm, "_", metric)),
      width = 11.5, height = 8.8
    )
    written_plots <- c(written_plots, files)
  }
  for (metric in profile_metrics) {
    p <- boxplot_metric(dp, metric, metric_labels[[metric]], include_all_p = FALSE)
    files <- save_plot(
      p,
      file.path(arm_dir, paste0("fig_", arm, "_", metric)),
      width = 9.0, height = 8.8
    )
    written_plots <- c(written_plots, files)
  }
}

# Main-text factor recovery figure: all p for Spectral-MAP and matched Gibbs
# comparisons where Gibbs was run.
main_factor <- boxplot_metric(
  results[results$mixture_setting == "separated", , drop = FALSE],
  "factor_score_rmse", "Factor score RMSE", include_all_p = TRUE
)
written_plots <- c(
  written_plots,
  save_plot(
    main_factor,
    file.path(plot_root, "fig_factor_recovery_separated"),
    width = 11.5, height = 8.8
  )
)

cell_summary <- function(d, metric) {
  keys <- c("method_label", "design", "n", "p")
  mean_d <- aggregate(d[[metric]], d[keys], mean, na.rm = TRUE)
  names(mean_d)[ncol(mean_d)] <- "estimate"
  se_d <- aggregate(d[[metric]], d[keys], function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  })
  names(se_d)[ncol(se_d)] <- "se"
  merge(mean_d, se_d, by = keys, all = TRUE)
}

extract_legend <- function(plot) {
  g <- ggplotGrob(plot + theme(legend.position = "bottom"))
  idx <- which(vapply(g$grobs, function(x) x$name, character(1L)) == "guide-box")
  if (!length(idx)) return(nullGrob())
  g$grobs[[idx[[1L]]]]
}

trend_plot <- function(d, metric, label, show_y = TRUE) {
  s <- cell_summary(d, metric)
  s$p_label <- factor(paste0("p = ", s$p), levels = c("p = 500", "p = 1000"))
  p <- ggplot(
    s,
    aes(
      x = n, y = estimate, color = method_label,
      linetype = p_label,
      group = interaction(method_label, p_label)
    )
  ) +
    geom_errorbar(
      aes(ymin = estimate - 2 * se, ymax = estimate + 2 * se),
      width = 8, linewidth = 0.38, alpha = 0.75
    ) +
    geom_line(linewidth = 0.72) +
    geom_point(size = 1.9) +
    facet_grid(design ~ ., scales = "free_y") +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_x_continuous(breaks = c(100, 200, 400)) +
    labs(
      title = label, x = "Sample size (n)",
      y = if (show_y) label else NULL,
      color = NULL, linetype = "Observed variables"
    ) +
    paper_theme(12.6) +
    theme(
      axis.text.x = element_text(angle = 0),
      legend.position = "bottom",
      legend.text = element_text(size = 10.8),
      axis.text = element_text(size = 10.8),
      axis.title = element_text(size = 11.5),
      strip.text.y = element_text(size = 10.7),
      plot.title = element_text(size = 13.0)
    )
  p
}

main_profile <- profile_results[
  profile_results$mixture_setting == "separated" &
    profile_results$p %in% c(500, 1000),
  , drop = FALSE
]
main_metrics <- c(
  marginal_mu_rmse = "Mixture mean RMSE",
  marginal_var_rmse = "Mixture variance RMSE",
  marginal_weight_rmse = "Mixture weight RMSE",
  component_profile_hamming_accuracy = "Hamming accuracy"
)
main_panels <- lapply(seq_along(main_metrics), function(i) {
  metric <- names(main_metrics)[[i]]
  trend_plot(main_profile, metric, main_metrics[[i]], show_y = TRUE)
})
legend <- extract_legend(main_panels[[1L]])
main_panels <- lapply(main_panels, function(p) p + theme(legend.position = "none"))
panel_row <- arrangeGrob(grobs = main_panels, ncol = 4)
main_composite <- arrangeGrob(panel_row, legend, ncol = 1, heights = c(12, 0.75))
for (ext in c("png", "pdf")) {
  out <- file.path(plot_root, paste0("fig_mixture_subtype_recovery_separated.", ext))
  if (ext == "png") {
    ggsave(out, main_composite, width = 11.8, height = 9.2, dpi = 400, bg = "white")
  } else {
    ggsave(
      out, main_composite, width = 11.8, height = 9.2,
      device = grDevices::pdf, bg = "white", useDingbats = FALSE
    )
  }
  written_plots <- c(written_plots, out)
}

# Generate the DGP figures in the paper bundle, then give the two figures used
# by the manuscript stable names.
dgp_dir <- plot_root
dgp_table_dir <- file.path(repo_root, "results", "saved", "simulation", "dgp")
dgp_script <- file.path(script_dir, "plot_paper_dgp_figures.R")
old_figure_dir <- Sys.getenv("FIGURE_DIR", unset = NA_character_)
old_table_dir <- Sys.getenv("TABLE_DIR", unset = NA_character_)
Sys.setenv(FIGURE_DIR = dgp_dir, TABLE_DIR = dgp_table_dir)
sys.source(dgp_script, envir = new.env(parent = globalenv()))
if (is.na(old_figure_dir)) Sys.unsetenv("FIGURE_DIR") else Sys.setenv(FIGURE_DIR = old_figure_dir)
if (is.na(old_table_dir)) Sys.unsetenv("TABLE_DIR") else Sys.setenv(TABLE_DIR = old_table_dir)
dgp_map <- c(
  dgp_loadings_p500 = "fig_data_generating_loadings",
  dgp_factor_marginals_G3 = "fig_data_generating_mixtures"
)
for (source_stem in names(dgp_map)) {
  for (ext in c("png", "pdf")) {
    source <- file.path(dgp_dir, paste0(source_stem, ".", ext))
    target <- file.path(plot_root, paste0(dgp_map[[source_stem]], ".", ext))
    if (!file.copy(source, target, overwrite = TRUE)) stop("Could not copy ", source)
    written_plots <- c(written_plots, target)
  }
}
unlink(
  list.files(dgp_dir, pattern = "^dgp_.*\\.(png|pdf)$", full.names = TRUE),
  force = TRUE
)

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
mean_sd <- function(x, digits = 3) {
  x <- x[is.finite(x)]
  paste0(fmt(mean(x), digits), " (", fmt(stats::sd(x), digits), ")")
}

# Main runtime table. Use original fit time from the completed separated run.
runtime_file <- arm_paths("separated")$recovery
runtime <- read_results(runtime_file)
runtime <- runtime[
  runtime$p %in% c(500, 1000) &
    runtime$method %in% c(product_method, gibbs_method),
  , drop = FALSE
]
runtime$G_number <- compact_g(runtime$G_config)
runtime_summary <- aggregate(
  seconds ~ method + H_true + G_number,
  data = runtime,
  FUN = mean
)
runtime_wide <- reshape(
  runtime_summary,
  idvar = c("H_true", "G_number"),
  timevar = "method",
  direction = "wide"
)
runtime_wide$speedup <-
  runtime_wide[[paste0("seconds.", gibbs_method)]] /
  runtime_wide[[paste0("seconds.", product_method)]]
ess <- aggregate(
  gibbs_median_parameter_ess ~ H_true + G_number,
  data = runtime[runtime$method == gibbs_method, , drop = FALSE],
  FUN = mean,
  na.action = na.pass
)
runtime_wide <- merge(runtime_wide, ess, by = c("H_true", "G_number"), all.x = TRUE)
runtime_wide <- runtime_wide[order(runtime_wide$H_true, runtime_wide$G_number), ]
write.csv(
  runtime_wide,
  file.path(table_root, "runtime_summary_by_H_G.csv"),
  row.names = FALSE
)

runtime_tex <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\caption{End-to-end computation time under the separated symmetric mixture configuration. Times are means in seconds, pooled over $n\\in\\{100,200,400\\}$ and $p\\in\\{500,1000\\}$, with 25 replications per design cell. Speed-up is the ratio of Gibbs time to Spectral-MAP time. ESS is the mean across replications of the median effective sample size over P-IFA parameters among the 1,000 retained Gibbs draws.}",
  "\\label{tab:runtime_comparison}",
  "\\begin{tabular}{ccrrrr}",
  "\\toprule",
  "$H$ & $G_h$ & Spectral-MAP & Gibbs & Speed-up & Gibbs ESS \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(runtime_wide))) {
  runtime_tex <- c(
    runtime_tex,
    sprintf(
      "%d & %d & %s & %s & \\textbf{%s$\\times$} & %s \\\\",
      runtime_wide$H_true[[i]], runtime_wide$G_number[[i]],
      fmt(runtime_wide[[paste0("seconds.", product_method)]][[i]], 1),
      fmt(runtime_wide[[paste0("seconds.", gibbs_method)]][[i]], 1),
      fmt(runtime_wide$speedup[[i]], 1),
      fmt(runtime_wide$gibbs_median_parameter_ess[[i]], 0)
    )
  )
}
runtime_tex <- c(runtime_tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(runtime_tex, file.path(table_root, "table_runtime_main.tex"))

# Pooled recovery table for the main text/supplement. All matched comparisons
# use the same 25 replications at p=500/1000.
summary_metrics <- c(
  "factor_score_rmse", "probability_rmse", "lambda_rmse", "alpha_rmse",
  "marginal_mu_rmse", "marginal_var_rmse", "marginal_weight_rmse",
  "component_profile_hamming_accuracy"
)
summary_rows <- list()
for (arm in arms) {
  for (method in c(product_method, gibbs_method)) {
    d <- profile_results[
      profile_results$mixture_setting == arm & profile_results$method == method,
      , drop = FALSE
    ]
    row <- data.frame(
      mixture_setting = arm_labels[[arm]],
      method = method_label[[method]],
      stringsAsFactors = FALSE
    )
    for (metric in summary_metrics) row[[metric]] <- mean_sd(d[[metric]])
    summary_rows[[length(summary_rows) + 1L]] <- row
  }
}
recovery_summary <- do.call(rbind, summary_rows)
write.csv(
  recovery_summary,
  file.path(table_root, "recovery_summary_by_mixture_setting.csv"),
  row.names = FALSE
)

recovery_tex <- c(
  "\\begin{table*}[t]",
  "\\centering",
  "\\caption{Recovery across the three factor-mixture configurations. Entries are mean RMSE (Monte Carlo standard deviation) or mean Hamming accuracy (standard deviation), pooled over the matched $p\\in\\{500,1000\\}$ design cells. Each design cell contains 25 replications. For RMSE, lower is better; for Hamming accuracy, higher is better.}",
  "\\label{tab:recovery_by_mixture_setting}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llcccccccc}",
  "\\toprule",
  "Configuration & Method & Factor & Probability & Loading & Intercept & Mean & Variance & Weight & Hamming \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(recovery_summary))) {
  r <- recovery_summary[i, ]
  recovery_tex <- c(
    recovery_tex,
    paste0(
      r$mixture_setting, " & ", r$method, " & ",
      paste(unlist(r[summary_metrics], use.names = FALSE), collapse = " & "),
      " \\\\"
    )
  )
  if (i %% 2L == 0L && i < nrow(recovery_summary)) {
    recovery_tex <- c(recovery_tex, "\\addlinespace")
  }
}
recovery_tex <- c(
  recovery_tex,
  "\\bottomrule", "\\end{tabular}%", "}", "\\end{table*}"
)
writeLines(recovery_tex, file.path(table_root, "table_recovery_supplement.tex"))

# Compare the Product MAP plug-in classifier with two defensible Gibbs
# summaries. The allocation-mode estimator takes the modal aligned sampled
# label for each subject and factor. The posterior-mean plug-in estimator first
# averages the aligned canonical draws and then classifies the posterior-mean
# factor score using the posterior-mean marginal mixture parameters.
subtype_rows <- list()
for (arm in arms) {
  product_d <- profile_results[
    profile_results$mixture_setting == arm &
      profile_results$method == product_method &
      profile_results$p %in% c(500, 1000),
    , drop = FALSE
  ]
  gibbs_d <- profile_results[
    profile_results$mixture_setting == arm &
      profile_results$method == gibbs_method &
      profile_results$p %in% c(500, 1000),
    , drop = FALSE
  ]
  subtype_rows[[length(subtype_rows) + 1L]] <- data.frame(
    mixture_setting = arm_labels[[arm]],
    estimator = c(
      "Spectral-MAP plug-in",
      "Gibbs allocation mode",
      "Gibbs posterior-mean plug-in"
    ),
    hamming_accuracy = c(
      mean_sd(product_d$component_profile_hamming_accuracy),
      mean_sd(gibbs_d$component_profile_hamming_accuracy),
      mean_sd(gibbs_d$posterior_mean_component_profile_hamming_accuracy)
    ),
    exact_profile_accuracy = c(
      mean_sd(product_d$component_profile_exact_accuracy),
      mean_sd(gibbs_d$component_profile_exact_accuracy),
      mean_sd(gibbs_d$posterior_mean_component_profile_exact_accuracy)
    ),
    stringsAsFactors = FALSE
  )
}
subtype_summary <- do.call(rbind, subtype_rows)
write.csv(
  subtype_summary,
  file.path(table_root, "subtype_estimator_summary.csv"),
  row.names = FALSE
)

subtype_tex <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\caption{Subtype recovery for Spectral-MAP and two Gibbs posterior summaries. Entries are means with Monte Carlo standard deviations in parentheses, pooled over the matched $p\\in\\{500,1000\\}$ design cells, with 25 replications per cell. The Gibbs allocation-mode estimator uses the modal aligned sampled component label for each subject and factor; the Gibbs posterior-mean plug-in estimator classifies the aligned posterior-mean score using the aligned posterior-mean mixture parameters. Higher values are better.}",
  "\\label{tab:subtype_estimators}",
  "\\begin{tabular}{llcc}",
  "\\toprule",
  "Configuration & Subtype estimator & Hamming & Exact profile \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(subtype_summary))) {
  r <- subtype_summary[i, ]
  config <- if ((i - 1L) %% 3L == 0L) r$mixture_setting else ""
  subtype_tex <- c(
    subtype_tex,
    paste0(
      config, " & ", r$estimator, " & ", r$hamming_accuracy,
      " & ", r$exact_profile_accuracy, " \\\\"
    )
  )
  if (i %% 3L == 0L && i < nrow(subtype_summary)) {
    subtype_tex <- c(subtype_tex, "\\addlinespace")
  }
}
subtype_tex <- c(
  subtype_tex,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"
)
writeLines(
  subtype_tex,
  file.path(table_root, "table_subtype_estimators.tex")
)

# Detailed runtime table for reproducibility and supplemental reporting.
runtime_cell <- aggregate(
  cbind(seconds, gibbs_median_parameter_ess) ~ method + n + p + H_true + G_number,
  data = runtime,
  FUN = function(x) mean(x, na.rm = TRUE)
)
write.csv(
  runtime_cell,
  file.path(table_root, "runtime_summary_by_n_p_H_G.csv"),
  row.names = FALSE
)

# Bundle the pre-existing, result-checked rotation tables with the section files.
file.copy(
  file.path(repo_root, "docs", "rotation_ablation_main_table.tex"),
  file.path(table_root, "table_rotation_ablation_main.tex"),
  overwrite = TRUE
)
file.copy(
  file.path(repo_root, "docs", "rotation_ablation_supplement_table.tex"),
  file.path(table_root, "table_rotation_ablation_supplement.tex"),
  overwrite = TRUE
)

include_tex <- c(
  "% Requires: \\usepackage{graphicx,booktabs}",
  "% Paths are relative to the repository root.",
  "\\begin{figure}[t]",
  "  \\centering",
  "  \\includegraphics[width=\\linewidth]{results/selected_plots/sample_size/paper_simulation_section/fig_data_generating_loadings.pdf}",
  "  \\includegraphics[width=\\linewidth]{results/selected_plots/sample_size/paper_simulation_section/fig_data_generating_mixtures.pdf}",
  "  \\caption{Data-generating loading structures and latent factor mixture distributions. Top: loading matrices for $p=500$ and $H\\in\\{5,10\\}$. Bottom: representative three-component factor marginals under separated symmetric, asymmetric-weight, and moderate-overlap configurations.}",
  "  \\label{fig:loadings_clusters}",
  "\\end{figure}",
  "",
  "\\begin{figure*}[t]",
  "  \\centering",
  "  \\includegraphics[width=\\textwidth]{results/selected_plots/sample_size/paper_simulation_section/fig_factor_recovery_separated.pdf}",
  "  \\caption{Factor-score recovery under the separated symmetric mixture configuration. Boxplots summarize 25 replications. Gibbs was run for $p\\in\\{500,1000\\}$; the larger-$p$ columns report Spectral-MAP.}",
  "  \\label{fig:factor_recovery_separated}",
  "\\end{figure*}",
  "",
  "\\begin{figure*}[t]",
  "  \\centering",
  "  \\includegraphics[width=\\textwidth]{results/selected_plots/sample_size/paper_simulation_section/fig_mixture_subtype_recovery_separated.pdf}",
  "  \\caption{Mixture-parameter and subtype recovery under the separated symmetric configuration. Points are means and bars show $\\pm2$ Monte Carlo standard errors over 25 replications. Line type distinguishes $p=500$ and $p=1000$. For Gibbs, Hamming accuracy uses the posterior mode of the aligned sampled allocations.}",
  "  \\label{fig:mixture_subtype_separated}",
  "\\end{figure*}",
  "",
  "\\input{results/selected_tables/sample_size/paper_simulation_section/table_subtype_estimators.tex}",
  "",
  "\\input{results/selected_tables/sample_size/paper_simulation_section/table_runtime_main.tex}",
  "\\input{results/selected_tables/sample_size/paper_simulation_section/table_rotation_ablation_main.tex}"
)
writeLines(include_tex, file.path(table_root, "simulation_section_includes.tex"))

manifest <- data.frame(
  artifact = c(
    "fig_data_generating_loadings.pdf",
    "fig_data_generating_mixtures.pdf",
    "fig_factor_recovery_separated.pdf",
    "fig_mixture_subtype_recovery_separated.pdf",
    "table_runtime_main.tex",
    "table_recovery_supplement.tex",
    "table_subtype_estimators.tex",
    "table_rotation_ablation_main.tex",
    "table_rotation_ablation_supplement.tex",
    "simulation_section_includes.tex"
  ),
  role = c(
    "Main DGP loading figure",
    "Main DGP mixture figure",
    "Main factor recovery figure",
    "Main mixture/subtype recovery figure",
    "Main runtime table",
    "Supplementary recovery table",
    "Subtype-estimator comparison table",
    "Main rotation ablation table",
    "Supplementary rotation ablation table",
    "Copy-paste LaTeX include block"
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(table_root, "artifact_manifest.csv"), row.names = FALSE)

cat("Wrote ", length(written_plots), " plot files.\n", sep = "")
cat("Figures: ", plot_root, "\n", sep = "")
cat("Tables:  ", table_root, "\n", sep = "")
