#!/usr/bin/env Rscript

# Build compact selected outputs comparing serial and parallel runtime for the
# proposed product-mixture estimator and the joint-mixture Gibbs baseline.
#
# Inputs are the focused timing runs:
#   - serial,
#   - 4 workers,
#   - 18 workers.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(repo_root)

input_files <- c(
  serial = "results/full/parallel_runtime_serial_focused_H3_G331_p500_MAP_intercepts/comparison_results_checkpoint.csv",
  workers4 = "results/full/parallel_runtime_parallel_focused_H3_G331_p500_MAP_intercepts/comparison_results_checkpoint.csv",
  workers18 = "results/full/parallel_runtime_parallel18_focused_H3_G331_p500_MAP_intercepts/comparison_results_checkpoint.csv"
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop("Missing focused timing result files: ", paste(missing_files, collapse = ", "))
}

results <- do.call(rbind, lapply(names(input_files), function(mode) {
  d <- read.csv(input_files[[mode]])
  d$parallel_mode <- mode
  d
}))

results$method_label <- ifelse(
  results$method == "independent_marginal_mixture",
  "Product mixture",
  "Joint Gibbs"
)

group_vars <- c(
  "parallel_mode", "method", "method_label", "n", "p",
  "H_true", "G_true", "loading_design"
)

summary <- aggregate(
  cbind(
    seconds,
    mean_factor_abs_cor,
    lambda_rmse,
    joint_mu_rmse,
    joint_var_rmse,
    joint_weight_rmse,
    flat_parameter_corr
  ) ~ parallel_mode + method + method_label + n + p + H_true + G_true + loading_design,
  results,
  mean,
  na.rm = TRUE
)

counts <- aggregate(
  rep ~ parallel_mode + method + method_label + n + p + H_true + G_true + loading_design,
  results,
  length
)
names(counts)[names(counts) == "rep"] <- "n_reps"

summary <- merge(summary, counts, by = group_vars)

serial_seconds <- subset(
  summary,
  parallel_mode == "serial",
  select = c("method", "n", "p", "H_true", "G_true", "loading_design", "seconds")
)
names(serial_seconds)[names(serial_seconds) == "seconds"] <- "seconds_serial"

speedup <- merge(
  summary,
  serial_seconds,
  by = c("method", "n", "p", "H_true", "G_true", "loading_design")
)
speedup$speedup_vs_serial <- speedup$seconds_serial / speedup$seconds

table_dir <- file.path("results", "selected_tables", "sample_size")
plot_dir <- file.path("results", "selected_plots", "sample_size")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  summary,
  file.path(table_dir, "parallel_runtime_method_comparison_summary.csv"),
  row.names = FALSE
)
write.csv(
  speedup,
  file.path(table_dir, "parallel_runtime_method_speedup_summary.csv"),
  row.names = FALSE
)

plot_data <- subset(speedup, parallel_mode %in% c("serial", "workers4", "workers18"))
plot_data$parallel_mode <- factor(
  plot_data$parallel_mode,
  levels = c("serial", "workers4", "workers18"),
  labels = c("serial", "4 workers", "18 workers")
)

method_levels <- c("Product mixture", "Joint Gibbs")
method_cols <- c("Product mixture" = "#2f6db3", "Joint Gibbs" = "#b23a48")
method_pch <- c("Product mixture" = 19, "Joint Gibbs" = 17)
n_values <- sort(unique(plot_data$n))

png(
  file.path(plot_dir, "parallel_runtime_method_comparison.png"),
  width = 1600,
  height = 900,
  res = 150
)
old_par <- par(mfrow = c(1, length(n_values)), mar = c(5, 5, 4, 1) + 0.1, oma = c(0, 0, 2, 0))
on.exit(par(old_par), add = TRUE)

for (n_value in n_values) {
  d_n <- subset(plot_data, n == n_value)
  plot(
    NA,
    xlim = c(0.8, 3.2),
    ylim = c(0, max(d_n$seconds) * 1.18),
    xaxt = "n",
    xlab = "parallel setting",
    ylab = "seconds per fit",
    main = paste0("n=", n_value, ", H=3, G=(3,3,1)")
  )
  axis(1, at = 1:3, labels = levels(plot_data$parallel_mode))

  for (method_name in method_levels) {
    d_m <- d_n[d_n$method_label == method_name, ]
    d_m <- d_m[order(d_m$parallel_mode), ]
    x <- seq_len(nrow(d_m))
    serial_y <- d_m$seconds[d_m$parallel_mode == "serial"]
    if (length(serial_y) == 1L && is.finite(serial_y)) {
      abline(
        h = serial_y,
        col = method_cols[method_name],
        lty = 3,
        lwd = 1.5
      )
      text(
        x = 3.12,
        y = serial_y,
        labels = paste0("serial ", sprintf("%.1f", serial_y)),
        pos = 3,
        cex = 0.68,
        col = method_cols[method_name]
      )
    }
    lines(
      x,
      d_m$seconds,
      type = "b",
      lwd = 2,
      pch = method_pch[method_name],
      col = method_cols[method_name]
    )
    text(
      x,
      d_m$seconds,
      labels = sprintf("%.1f", d_m$seconds),
      pos = 3,
      cex = 0.75,
      col = method_cols[method_name]
    )
  }

  if (n_value == n_values[1L]) {
    legend(
      "topleft",
      legend = c(method_levels, "serial baseline"),
      col = c(method_cols[method_levels], "gray35"),
      pch = c(method_pch[method_levels], NA),
      lty = c(1, 1, 3),
      lwd = c(2, 2, 1.5),
      bty = "n"
    )
  }
}

mtext(
  "Serial and parallel runtime: product mixture vs joint Gibbs",
  outer = TRUE,
  cex = 1.15,
  font = 2
)
dev.off()

print(
  summary[order(summary$n, summary$method_label, summary$parallel_mode),
    c(
      "parallel_mode", "method_label", "n", "n_reps", "seconds",
      "mean_factor_abs_cor", "lambda_rmse", "flat_parameter_corr"
    )
  ],
  row.names = FALSE
)
print(
  speedup[order(speedup$n, speedup$method_label, speedup$parallel_mode),
    c("parallel_mode", "method_label", "n", "seconds", "seconds_serial", "speedup_vs_serial")
  ],
  row.names = FALSE
)

cat(
  "Wrote selected timing comparison to:\n",
  normalizePath(file.path(table_dir, "parallel_runtime_method_comparison_summary.csv")), "\n",
  normalizePath(file.path(table_dir, "parallel_runtime_method_speedup_summary.csv")), "\n",
  normalizePath(file.path(plot_dir, "parallel_runtime_method_comparison.png")), "\n"
)
