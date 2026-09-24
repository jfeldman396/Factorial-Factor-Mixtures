#!/usr/bin/env Rscript

# Validate and summarize the completed three-arm rotation ablation.
#
# The full per-replication files remain in results/diagnostics (gitignored).
# This script writes compact scientific summaries and paper-ready PNG/PDF
# figures to tracked results/selected_* directories.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

run_prefix <- "rotation_ablation_ifeval_u1_2_cp0_05_sep2_penalty3_5_"
arms <- c("separated", "asymmetric_pi", "moderate_overlap")
arm_labels <- c(
  separated = "Separated",
  asymmetric_pi = "Asymmetric weights",
  moderate_overlap = "Moderate overlap"
)
signal_levels <- c("estimated_Z", "oracle_Z")
signal_labels <- c(estimated_Z = "Estimated Z", oracle_Z = "Oracle Z")
method_levels <- c("fastica_only", "mixture_fastica_start", "mixture_identity_start")
method_labels <- c(
  fastica_only = "FastICA only",
  mixture_fastica_start = "FastICA + mixture",
  mixture_identity_start = "Identity + mixture"
)
method_colors <- c(
  fastica_only = "#C43C4A",
  mixture_fastica_start = "#2E7D5B",
  mixture_identity_start = "#2F6DAE"
)
method_pch <- c(fastica_only = 17, mixture_fastica_start = 16, mixture_identity_start = 15)
stage_levels <- c("pretrain", "refined")

expected <- list(
  n = c(100L, 200L, 400L),
  p = c(500L, 1000L, 1500L, 2000L),
  H = c(5L, 10L),
  rep = 1:25,
  signal_source = signal_levels,
  rotation_method = method_levels,
  stage = stage_levels
)
expected_g <- c(
  paste(rep(2L, 5L), collapse = "-"),
  paste(rep(3L, 5L), collapse = "-"),
  paste(rep(2L, 10L), collapse = "-"),
  paste(rep(3L, 10L), collapse = "-")
)

diagnostic_root <- file.path(repo_root, "results", "diagnostics")
plot_dir <- file.path(
  repo_root, "results", "selected_plots", "sample_size",
  "rotation_ablation_three_arms"
)
table_dir <- file.path(
  repo_root, "results", "selected_tables", "sample_size",
  "rotation_ablation_three_arms"
)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

mean_finite <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

sd_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x)
}

se_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

aggregate_metric <- function(data, metric, groups) {
  pieces <- split(data[[metric]], interaction(data[groups], drop = TRUE, lex.order = TRUE))
  keys <- unique(data[groups])
  key_id <- interaction(keys, drop = TRUE, lex.order = TRUE)
  stats <- lapply(pieces, function(x) {
    c(
      mean = mean_finite(x),
      median = median(x[is.finite(x)], na.rm = TRUE),
      sd = sd_finite(x),
      se = se_finite(x),
      n_replications = sum(is.finite(x))
    )
  })
  stat_df <- as.data.frame(do.call(rbind, stats))
  stat_df$.group <- names(pieces)
  keys$.group <- as.character(key_id)
  out <- merge(keys, stat_df, by = ".group", all.y = TRUE, sort = FALSE)
  out$.group <- NULL
  names(out)[(ncol(out) - 4L):ncol(out)] <- paste0(metric, "_", names(out)[(ncol(out) - 4L):ncol(out)])
  out
}

read_arm <- function(arm) {
  path <- file.path(
    diagnostic_root,
    paste0(run_prefix, arm),
    "rotation_fastica_results.csv"
  )
  if (!file.exists(path)) stop("Missing rotation-ablation result: ", path, call. = FALSE)
  data <- read.csv(path, check.names = FALSE)
  data$mixture_scenario <- arm
  data
}

read_eigengap_arm <- function(arm) {
  path <- file.path(
    diagnostic_root,
    paste0(run_prefix, arm),
    "rotation_fastica_eigengap_dataset_results.csv"
  )
  if (!file.exists(path)) stop("Missing eigengap result: ", path, call. = FALSE)
  data <- read.csv(path, check.names = FALSE)
  data$mixture_scenario <- arm
  data
}

results <- do.call(rbind, lapply(arms, read_arm))
eigengap <- do.call(rbind, lapply(arms, read_eigengap_arm))
rownames(results) <- NULL
rownames(eigengap) <- NULL

key_cols <- c(
  "mixture_scenario", "n", "p", "H", "G_config", "rep",
  "signal_source", "rotation_method", "stage"
)
dataset_cols <- c("mixture_scenario", "n", "p", "H", "G_config", "rep")
expected_datasets_per_arm <- length(expected$n) * length(expected$p) *
  length(expected$H) * 2L * length(expected$rep)
expected_rows_per_arm <- expected_datasets_per_arm * length(signal_levels) *
  length(method_levels) * length(stage_levels)

coverage <- lapply(arms, function(arm) {
  data <- results[results$mixture_scenario == arm, , drop = FALSE]
  dataset_keys <- do.call(paste, c(data[dataset_cols], sep = "|"))
  scientific_keys <- do.call(paste, c(data[key_cols], sep = "|"))
  core_metrics <- c(
    "factor_score_rmse", "lambda_rmse", "alpha_rmse", "probability_rmse",
    "marginal_mu_rmse", "marginal_var_rmse", "marginal_weight_rmse",
    "end_to_end_seconds"
  )
  data.frame(
    mixture_scenario = arm,
    observed_rows = nrow(data),
    expected_rows = expected_rows_per_arm,
    unique_datasets = length(unique(dataset_keys)),
    expected_datasets = expected_datasets_per_arm,
    duplicated_scientific_keys = sum(duplicated(scientific_keys)),
    missing_core_metrics = sum(!is.finite(as.matrix(data[core_metrics]))),
    n_values_match = identical(sort(unique(data$n)), expected$n),
    p_values_match = identical(sort(unique(data$p)), expected$p),
    H_values_match = identical(sort(unique(data$H)), expected$H),
    G_values_match = identical(sort(unique(data$G_config)), sort(expected_g)),
    rep_values_match = identical(sort(unique(data$rep)), expected$rep),
    signal_values_match = identical(sort(unique(data$signal_source)), sort(signal_levels)),
    method_values_match = identical(sort(unique(data$rotation_method)), sort(method_levels)),
    stage_values_match = identical(sort(unique(data$stage)), sort(stage_levels)),
    stringsAsFactors = FALSE
  )
})
coverage <- do.call(rbind, coverage)

if (any(coverage$observed_rows != coverage$expected_rows) ||
    any(coverage$unique_datasets != coverage$expected_datasets) ||
    any(coverage$duplicated_scientific_keys != 0L) ||
    any(coverage$missing_core_metrics != 0L) ||
    any(!as.matrix(coverage[grep("_match$", names(coverage))]))) {
  write.csv(coverage, file.path(table_dir, "rotation_ablation_coverage_validation.csv"), row.names = FALSE)
  stop("Rotation-ablation coverage validation failed; see the coverage table.", call. = FALSE)
}

write.csv(
  coverage,
  file.path(table_dir, "rotation_ablation_coverage_validation.csv"),
  row.names = FALSE
)

metrics <- c(
  "factor_score_rmse", "mean_factor_abs_cor", "lambda_rmse", "alpha_rmse",
  "probability_rmse", "marginal_mu_rmse", "marginal_var_rmse",
  "marginal_weight_rmse", "end_to_end_seconds", "rotation_seconds",
  "refinement_seconds"
)
summary_groups <- c(
  "mixture_scenario", "n", "p", "H", "G_config", "signal_source",
  "rotation_method", "stage"
)
summary_parts <- lapply(metrics, function(metric) aggregate_metric(results, metric, summary_groups))
cell_summary <- Reduce(
  function(left, right) merge(left, right, by = summary_groups, all = TRUE, sort = FALSE),
  summary_parts
)
cell_summary <- cell_summary[do.call(order, cell_summary[summary_groups]), , drop = FALSE]
write.csv(
  cell_summary,
  file.path(table_dir, "rotation_ablation_cell_summary.csv"),
  row.names = FALSE
)

signal_error_metrics <- c(
  "estimated_vs_oracle_signal_rmse", "estimated_vs_oracle_alpha_rmse",
  "estimated_vs_oracle_linear_predictor_rmse", "estimated_vs_oracle_probability_rmse",
  "estimated_vs_oracle_subspace_distance", "estimated_vs_full_oracle_z_rmse",
  "oracle_lowrank_vs_full_oracle_z_rmse"
)
dataset_once <- unique(results[c(dataset_cols, signal_error_metrics)])
signal_groups <- c("mixture_scenario", "n", "p", "H", "G_config")
signal_parts <- lapply(signal_error_metrics, function(metric) {
  aggregate_metric(dataset_once, metric, signal_groups)
})
signal_summary <- Reduce(
  function(left, right) merge(left, right, by = signal_groups, all = TRUE, sort = FALSE),
  signal_parts
)
signal_summary <- signal_summary[do.call(order, signal_summary[signal_groups]), , drop = FALSE]
write.csv(
  signal_summary,
  file.path(table_dir, "rotation_ablation_estimated_z_error_summary.csv"),
  row.names = FALSE
)

rank_metrics <- c(
  "estimated_z_largest_eigengap_selects_true_H",
  "estimated_z_largest_relative_eigengap_selects_true_H",
  "estimated_z_largest_sv_ratio_selects_true_H",
  "estimated_z_eigengap_at_true_H",
  "estimated_z_relative_eigengap_at_true_H",
  "estimated_z_sv_ratio_at_true_H",
  "eigengap_diagnostic_seconds"
)
rank_groups <- c("mixture_scenario", "n", "p", "H", "G_config")
rank_parts <- lapply(rank_metrics, function(metric) aggregate_metric(eigengap, metric, rank_groups))
rank_summary <- Reduce(
  function(left, right) merge(left, right, by = rank_groups, all = TRUE, sort = FALSE),
  rank_parts
)
rank_summary <- rank_summary[do.call(order, rank_summary[rank_groups]), , drop = FALSE]
write.csv(
  rank_summary,
  file.path(table_dir, "rotation_ablation_rank_selection_summary.csv"),
  row.names = FALSE
)

pair_cols <- c(
  "mixture_scenario", "n", "p", "H", "G_config", "rep",
  "signal_source", "stage"
)
contrast_metrics <- c(
  "factor_score_rmse", "lambda_rmse", "alpha_rmse", "probability_rmse",
  "marginal_mu_rmse", "marginal_var_rmse", "marginal_weight_rmse",
  "end_to_end_seconds"
)
fastica <- results[results$rotation_method == "fastica_only", c(pair_cols, contrast_metrics)]
names(fastica)[match(contrast_metrics, names(fastica))] <- paste0(contrast_metrics, "_fastica_only")
contrast_rows <- lapply(c("mixture_fastica_start", "mixture_identity_start"), function(method) {
  other <- results[results$rotation_method == method, c(pair_cols, contrast_metrics)]
  names(other)[match(contrast_metrics, names(other))] <- paste0(contrast_metrics, "_comparison")
  paired <- merge(fastica, other, by = pair_cols, all = FALSE)
  paired$comparison_method <- method
  for (metric in contrast_metrics) {
    paired[[paste0(metric, "_difference_vs_fastica")]] <-
      paired[[paste0(metric, "_comparison")]] - paired[[paste0(metric, "_fastica_only")]]
  }
  paired
})
paired_contrasts <- do.call(rbind, contrast_rows)
contrast_groups <- c("mixture_scenario", "signal_source", "stage", "comparison_method", "H", "G_config")
contrast_diff_metrics <- paste0(contrast_metrics, "_difference_vs_fastica")
contrast_parts <- lapply(contrast_diff_metrics, function(metric) {
  aggregate_metric(paired_contrasts, metric, contrast_groups)
})
contrast_summary <- Reduce(
  function(left, right) merge(left, right, by = contrast_groups, all = TRUE, sort = FALSE),
  contrast_parts
)
write.csv(
  contrast_summary,
  file.path(table_dir, "rotation_ablation_method_contrasts_vs_fastica.csv"),
  row.names = FALSE
)

best_method_groups <- c("mixture_scenario", "signal_source", "stage", "H", "G_config")
best_method <- aggregate(
  results$factor_score_rmse,
  results[c(best_method_groups, "rotation_method")],
  mean_finite
)
names(best_method)[ncol(best_method)] <- "mean_factor_score_rmse"
best_method <- do.call(rbind, lapply(split(best_method, interaction(best_method[best_method_groups], drop = TRUE)), function(data) {
  data$best_factor_method <- data$rotation_method == data$rotation_method[which.min(data$mean_factor_score_rmse)]
  data
}))
write.csv(
  best_method,
  file.path(table_dir, "rotation_ablation_best_factor_method.csv"),
  row.names = FALSE
)

open_device <- function(path, width, height, res = 180) {
  if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
    pdf(path, width = width / res, height = height / res, useDingbats = FALSE)
  } else {
    png(path, width = width, height = height, res = res)
  }
}

with_devices <- function(stem, width, height, draw) {
  for (extension in c("png", "pdf")) {
    path <- file.path(plot_dir, paste0(stem, ".", extension))
    open_device(path, width, height)
    draw()
    dev.off()
  }
}

plot_factor_prepost <- function() {
  paired <- reshape(
    cell_summary[c(
      "mixture_scenario", "n", "p", "H", "G_config", "signal_source",
      "rotation_method", "stage", "factor_score_rmse_mean"
    )],
    idvar = c("mixture_scenario", "n", "p", "H", "G_config", "signal_source", "rotation_method"),
    timevar = "stage",
    direction = "wide"
  )
  x_name <- "factor_score_rmse_mean.pretrain"
  y_name <- "factor_score_rmse_mean.refined"
  lim <- range(c(paired[[x_name]], paired[[y_name]]), finite = TRUE)
  pad <- 0.04 * diff(lim)
  lim <- lim + c(-pad, pad)
  op <- par(mfrow = c(3, 2), mar = c(4.5, 4.6, 2.8, 1), oma = c(0, 0, 3.8, 0))
  on.exit(par(op), add = TRUE)
  for (arm in arms) {
    for (signal in signal_levels) {
      data <- paired[paired$mixture_scenario == arm & paired$signal_source == signal, ]
      plot(
        data[[x_name]], data[[y_name]],
        xlim = lim, ylim = lim, asp = 1,
        xlab = "Pre-refinement factor RMSE",
        ylab = "Post-refinement factor RMSE",
        main = paste(arm_labels[[arm]], signal_labels[[signal]], sep = ": "),
        pch = method_pch[data$rotation_method],
        col = paste0(method_colors[data$rotation_method], "A6"),
        cex = 0.75
      )
      abline(0, 1, lty = 2, col = "#666666", lwd = 1.5)
      grid(col = "#E7E7E7")
      if (arm == arms[1L] && signal == signal_levels[1L]) {
        legend(
          "topleft", legend = unname(method_labels), col = method_colors,
          pch = method_pch, bty = "n", cex = 0.8
        )
      }
    }
  }
  mtext("Factor recovery before and after refinement", outer = TRUE, font = 2, cex = 1.35)
}

plot_estimated_z_error <- function() {
  op <- par(mfrow = c(3, 4), mar = c(4.2, 4.4, 2.7, 1), oma = c(0, 0, 3.8, 0))
  on.exit(par(op), add = TRUE)
  n_colors <- c(`100` = "#C43C4A", `200` = "#2F6DAE", `400` = "#2E7D5B")
  for (arm in arms) {
    for (H_value in expected$H) {
      for (g_value in c(2L, 3L)) {
        g_config <- paste(rep(g_value, H_value), collapse = "-")
        data <- signal_summary[
          signal_summary$mixture_scenario == arm &
            signal_summary$H == H_value & signal_summary$G_config == g_config,
          , drop = FALSE
        ]
        y <- data$estimated_vs_oracle_signal_rmse_mean
        ylim <- range(y, finite = TRUE)
        ylim <- ylim + c(-0.08, 0.12) * diff(ylim + c(-1e-8, 1e-8))
        plot(
          NA, xlim = range(expected$p), ylim = ylim,
          xlab = "Number of items (p)", ylab = "Estimated-Z signal RMSE",
          main = sprintf("%s: H=%d, G=%d", arm_labels[[arm]], H_value, g_value)
        )
        grid(col = "#E7E7E7")
        for (n_value in expected$n) {
          line <- data[data$n == n_value, ]
          line <- line[order(line$p), ]
          lines(line$p, line$estimated_vs_oracle_signal_rmse_mean,
                col = n_colors[as.character(n_value)], lwd = 2)
          points(line$p, line$estimated_vs_oracle_signal_rmse_mean,
                 col = n_colors[as.character(n_value)], pch = 16, cex = 0.8)
        }
        if (arm == arms[1L] && H_value == expected$H[1L] && g_value == 2L) {
          legend("topright", legend = paste0("n=", expected$n), col = n_colors,
                 lty = 1, pch = 16, bty = "n", cex = 0.72)
        }
      }
    }
  }
  mtext("Estimated-Z signal error across sample size and item count", outer = TRUE, font = 2, cex = 1.35)
}

plot_rank_accuracy <- function() {
  criteria <- c(
    estimated_z_largest_eigengap_selects_true_H_mean = "Raw eigengap",
    estimated_z_largest_relative_eigengap_selects_true_H_mean = "Relative eigengap",
    estimated_z_largest_sv_ratio_selects_true_H_mean = "Singular-value ratio"
  )
  colors <- c("#2F6DAE", "#2E7D5B", "#C43C4A")
  op <- par(mfrow = c(3, 2), mar = c(4.2, 4.5, 2.8, 1), oma = c(0, 0, 3.8, 0))
  on.exit(par(op), add = TRUE)
  for (arm in arms) {
    for (H_value in expected$H) {
      data <- rank_summary[rank_summary$mixture_scenario == arm & rank_summary$H == H_value, ]
      agg <- aggregate(data[names(criteria)], data[c("p")], mean_finite)
      plot(
        NA, xlim = range(expected$p), ylim = c(0, 1),
        xlab = "Number of items (p)", ylab = "Correct-rank selection rate",
        main = sprintf("%s: H=%d", arm_labels[[arm]], H_value)
      )
      grid(col = "#E7E7E7")
      for (j in seq_along(criteria)) {
        lines(agg$p, agg[[names(criteria)[j]]], col = colors[j], lwd = 2)
        points(agg$p, agg[[names(criteria)[j]]], col = colors[j], pch = 14 + j, cex = 0.85)
      }
      if (arm == arms[1L] && H_value == expected$H[1L]) {
        legend("bottomright", legend = unname(criteria), col = colors,
               lty = 1, pch = 15:17, bty = "n", cex = 0.78)
      }
    }
  }
  mtext("Estimated-Z rank selection", outer = TRUE, font = 2, cex = 1.35)
}

plot_runtime <- function() {
  data <- cell_summary[
    cell_summary$stage == "refined" & cell_summary$signal_source == "estimated_Z",
    , drop = FALSE
  ]
  op <- par(mfrow = c(3, 4), mar = c(4.2, 4.5, 2.7, 1), oma = c(0, 0, 3.8, 0))
  on.exit(par(op), add = TRUE)
  for (arm in arms) {
    for (H_value in expected$H) {
      for (g_value in c(2L, 3L)) {
        g_config <- paste(rep(g_value, H_value), collapse = "-")
        panel <- data[
          data$mixture_scenario == arm & data$H == H_value & data$G_config == g_config,
          , drop = FALSE
        ]
        agg <- aggregate(
          panel$end_to_end_seconds_mean,
          panel[c("p", "rotation_method")],
          mean_finite
        )
        names(agg)[ncol(agg)] <- "seconds"
        ylim <- range(agg$seconds, finite = TRUE)
        ylim <- ylim + c(-0.06, 0.12) * diff(ylim + c(-1e-8, 1e-8))
        plot(
          NA, xlim = range(expected$p), ylim = ylim,
          xlab = "Number of items (p)", ylab = "End-to-end seconds",
          main = sprintf("%s: H=%d, G=%d", arm_labels[[arm]], H_value, g_value)
        )
        grid(col = "#E7E7E7")
        for (method in method_levels) {
          line <- agg[agg$rotation_method == method, ]
          line <- line[order(line$p), ]
          lines(line$p, line$seconds, col = method_colors[[method]], lwd = 2)
          points(line$p, line$seconds, col = method_colors[[method]],
                 pch = method_pch[[method]], cex = 0.8)
        }
        if (arm == arms[1L] && H_value == expected$H[1L] && g_value == 2L) {
          legend("topleft", legend = unname(method_labels), col = method_colors,
                 lty = 1, pch = method_pch, bty = "n", cex = 0.68)
        }
      }
    }
  }
  mtext("Estimated-Z end-to-end runtime", outer = TRUE, font = 2, cex = 1.35)
}

with_devices("rotation_ablation_factor_pretrain_vs_refined", 2100, 2600, plot_factor_prepost)
with_devices("rotation_ablation_estimated_z_error", 3000, 2300, plot_estimated_z_error)
with_devices("rotation_ablation_rank_selection_accuracy", 2100, 2500, plot_rank_accuracy)
with_devices("rotation_ablation_runtime", 3000, 2300, plot_runtime)

overall <- aggregate(
  cbind(factor_score_rmse, lambda_rmse, probability_rmse, end_to_end_seconds) ~
    mixture_scenario + signal_source + stage + rotation_method,
  results,
  mean_finite
)
overall <- overall[do.call(order, overall[c("mixture_scenario", "signal_source", "stage", "rotation_method")]), ]
write.csv(
  overall,
  file.path(table_dir, "rotation_ablation_overall_method_summary.csv"),
  row.names = FALSE
)

cat("Validated rows:", nrow(results), "\n")
cat("Validated data sets:", nrow(eigengap), "\n")
cat("Figures:", plot_dir, "\n")
cat("Tables:", table_dir, "\n")
