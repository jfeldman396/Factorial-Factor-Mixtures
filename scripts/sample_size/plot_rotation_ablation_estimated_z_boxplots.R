#!/usr/bin/env Rscript

# Paper-ready n/p/H/G panels for estimated probit-signal recovery.

options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.", call. = FALSE)
}

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

run_prefix <- "rotation_ablation_ifeval_u1_2_cp0_05_sep2_penalty3_5_"
arms <- c("separated", "asymmetric_pi", "moderate_overlap")
arm_titles <- c(
  separated = "Separated mixtures",
  asymmetric_pi = "Asymmetric mixture weights",
  moderate_overlap = "Moderate component overlap"
)
diagnostic_root <- file.path(repo_root, "results", "diagnostics")
plot_dir <- file.path(
  repo_root, "results", "selected_plots", "sample_size",
  "rotation_ablation_three_arms", "estimated_z_boxplots"
)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

read_arm <- function(arm) {
  path <- file.path(
    diagnostic_root,
    paste0(run_prefix, arm),
    "rotation_fastica_results.csv"
  )
  if (!file.exists(path)) stop("Missing result file: ", path, call. = FALSE)
  out <- read.csv(path, check.names = FALSE)
  out$mixture_scenario <- arm
  out
}

make_panel_plot <- function(data, arm) {
  data <- data[data$mixture_scenario == arm, , drop = FALSE]
  dataset_cols <- c(
    "mixture_scenario", "n", "p", "H", "G_config", "rep",
    "estimated_vs_oracle_signal_rmse"
  )
  data <- unique(data[dataset_cols])
  if (!nrow(data)) stop("No estimated-Z signal diagnostics for arm: ", arm, call. = FALSE)

  data$n_label <- factor(data$n, levels = sort(unique(data$n)), labels = paste0("n=", sort(unique(data$n))))
  data$p_label <- factor(data$p, levels = sort(unique(data$p)), labels = paste0("p=", sort(unique(data$p))))
  data$G <- as.integer(sub("-.*$", "", data$G_config))
  hg_levels <- c("H=5, G=2", "H=5, G=3", "H=10, G=2", "H=10, G=3")
  data$HG <- factor(paste0("H=", data$H, ", G=", data$G), levels = hg_levels)

  ggplot2::ggplot(data, ggplot2::aes(x = n_label, y = estimated_vs_oracle_signal_rmse)) +
    ggplot2::geom_boxplot(
      width = 0.58,
      linewidth = 0.55,
      color = "#245B8A",
      fill = "#79A7D3",
      alpha = 0.82,
      outlier.color = "#245B8A",
      outlier.size = 1.35,
      outlier.alpha = 0.65
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(HG), cols = ggplot2::vars(p_label), scales = "free_y") +
    ggplot2::labs(
      title = paste0("Estimated Probit Signal Recovery: ", arm_titles[[arm]]),
      x = "Sample size",
      y = "Estimated-vs-oracle signal RMSE"
    ) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 12)),
      axis.title = ggplot2::element_text(size = 17),
      axis.text = ggplot2::element_text(size = 12, color = "#222222"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      strip.text.x = ggplot2::element_text(size = 14, face = "bold"),
      strip.text.y = ggplot2::element_text(size = 13, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#E6E6E6", color = "#555555", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#E8E8E8", linewidth = 0.4),
      panel.spacing = grid::unit(0.55, "lines"),
      plot.margin = ggplot2::margin(12, 15, 8, 12)
    )
}

results <- do.call(rbind, lapply(arms, read_arm))

for (arm in arms) {
  figure <- make_panel_plot(results, arm)
  stem <- paste0("rotation_ablation_estimated_z_signal_rmse_", arm)
  ggplot2::ggsave(
    filename = file.path(plot_dir, paste0(stem, ".png")),
    plot = figure,
    width = 13.5,
    height = 10,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(plot_dir, paste0(stem, ".pdf")),
    plot = figure,
    width = 13.5,
    height = 10,
    units = "in",
    device = grDevices::pdf,
    bg = "white"
  )
}

cat("Wrote estimated-Z rotation-ablation boxplots to:\n", plot_dir, "\n")
