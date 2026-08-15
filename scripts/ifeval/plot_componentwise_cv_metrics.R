#!/usr/bin/env Rscript

# Visualize held-out CV metrics for models with factor-specific component
# counts, e.g. G_config = "3,3,1,2".  The script is intentionally independent
# of the fitting code: it only reads the checkpointed fold-score CSV.

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else getwd()
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(repo_root, "results", "ifeval_columnwise_G_cv_atmost1_gaussian_20260815")
)
score_path <- Sys.getenv(
  "SCORE_PATH",
  file.path(out_dir, "ifeval_rank_lambda_cv_fold_scores.csv")
)
top_n <- as.integer(Sys.getenv("TOP_N", "40"))

if (!file.exists(score_path)) {
  stop("Cannot find CV fold-score file: ", score_path)
}

scores <- read.csv(score_path, check.names = FALSE)
scores <- scores[scores$method == "independent_mixture_probit", , drop = FALSE]
scores <- scores[!is.na(scores$G_config) & nzchar(scores$G_config), , drop = FALSE]

key_cols <- c("method", "H", "G", "G_config", "lambda_l1_penalty")
split_key <- interaction(scores[, key_cols], drop = TRUE)
summary_list <- lapply(split(scores, split_key), function(d) {
  data.frame(
    method = d$method[1L],
    H = d$H[1L],
    G = d$G[1L],
    G_config = d$G_config[1L],
    lambda_l1_penalty = d$lambda_l1_penalty[1L],
    n_folds = nrow(d),
    mean_heldout_loglik_per_response = mean(d$heldout_loglik_per_response),
    se_heldout_loglik_per_response = if (nrow(d) > 1L) sd(d$heldout_loglik_per_response) / sqrt(nrow(d)) else NA_real_,
    mean_heldout_bic = mean(d$heldout_bic),
    mean_training_bic = mean(d$training_bic),
    mean_fit_seconds = mean(d$fit_seconds),
    stringsAsFactors = FALSE
  )
})
summary_scores <- do.call(rbind, summary_list)
summary_scores <- summary_scores[summary_scores$n_folds == 3L, , drop = FALSE]
summary_scores <- summary_scores[order(summary_scores$mean_heldout_loglik_per_response, decreasing = TRUE), ]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  summary_scores,
  file.path(out_dir, "componentwise_cv_summary_complete.csv"),
  row.names = FALSE
)

parse_G <- function(g_config, H) {
  vals <- as.integer(strsplit(gsub("[[:space:]]+", "", g_config), ",", fixed = TRUE)[[1L]])
  if (length(vals) == 1L && H > 1L) vals <- rep(vals, H)
  vals
}

component_matrix <- function(d) {
  max_H <- max(d$H)
  mat <- matrix(NA_integer_, nrow = nrow(d), ncol = max_H)
  for (i in seq_len(nrow(d))) {
    vals <- parse_G(d$G_config[i], d$H[i])
    mat[i, seq_along(vals)] <- vals
  }
  colnames(mat) <- paste0("F", seq_len(max_H))
  mat
}

component_cols <- c("1" = "#9CA3AF", "2" = "#2F6DAE", "3" = "#C43C4A")
metric_cols <- colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))

plot_top_component_configs <- function(d, top_n, path) {
  d <- head(d, min(top_n, nrow(d)))
  d <- d[order(d$mean_heldout_loglik_per_response), , drop = FALSE]
  mat <- component_matrix(d)
  row_labels <- paste0(
    "H=", d$H,
    ", G=(", d$G_config, ")",
    ", lambda=", d$lambda_l1_penalty
  )

  png(path, width = 1900, height = max(900, 34 * nrow(d) + 220), res = 180)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.35, 1))

  par(mar = c(4.5, 13, 4, 1))
  image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat),
    col = component_cols,
    breaks = c(0.5, 1.5, 2.5, 3.5),
    axes = FALSE,
    xlab = "factor coordinate",
    ylab = "",
    main = paste0("Top ", nrow(d), " component-wise CV fits")
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
  axis(2, at = seq_len(nrow(mat)), labels = row_labels, las = 2, cex.axis = 0.62)
  abline(v = seq(0.5, ncol(mat) + 0.5, by = 1), col = "white")
  abline(h = seq(0.5, nrow(mat) + 0.5, by = 1), col = "white")
  legend(
    "topleft",
    legend = paste("G =", names(component_cols)),
    fill = component_cols,
    horiz = TRUE,
    bty = "n",
    cex = 0.85
  )

  par(mar = c(4.5, 4.8, 4, 1))
  y <- seq_len(nrow(d))
  x <- d$mean_heldout_loglik_per_response
  se <- d$se_heldout_loglik_per_response
  xr <- range(c(x - 2 * se, x + 2 * se), na.rm = TRUE)
  plot(
    x,
    y,
    pch = 19,
    xlim = xr,
    yaxt = "n",
    xlab = "held-out log likelihood per response",
    ylab = "",
    main = "CV metric"
  )
  segments(x - 2 * se, y, x + 2 * se, y, col = "#4B5563")
  abline(v = max(x, na.rm = TRUE), col = "#B91C1C", lty = 2)
  grid(nx = NULL, ny = NA, col = "#E5E7EB")
}

plot_H_lambda_config_heatmaps <- function(d, path, max_rows_per_panel = 35L) {
  H_vals <- sort(unique(d$H))
  png(path, width = 2200, height = 650 * length(H_vals), res = 180)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  par(mfrow = c(length(H_vals), 1), mar = c(4.2, 12, 3.2, 4))

  global_range <- range(d$mean_heldout_loglik_per_response, na.rm = TRUE)
  cols <- metric_cols(80)

  for (H_now in H_vals) {
    dd <- d[d$H == H_now, , drop = FALSE]
    best_by_config <- aggregate(
      mean_heldout_loglik_per_response ~ G_config,
      dd,
      max
    )
    best_by_config <- best_by_config[order(best_by_config$mean_heldout_loglik_per_response, decreasing = TRUE), ]
    keep_configs <- head(best_by_config$G_config, min(max_rows_per_panel, nrow(best_by_config)))
    dd <- dd[dd$G_config %in% keep_configs, , drop = FALSE]

    lambdas <- sort(unique(dd$lambda_l1_penalty))
    configs <- rev(keep_configs)
    z <- matrix(NA_real_, nrow = length(lambdas), ncol = length(configs))
    for (i in seq_len(nrow(dd))) {
      z[match(dd$lambda_l1_penalty[i], lambdas), match(dd$G_config[i], configs)] <-
        dd$mean_heldout_loglik_per_response[i]
    }

    image(
      x = seq_along(lambdas),
      y = seq_along(configs),
      z = z,
      col = cols,
      zlim = global_range,
      axes = FALSE,
      xlab = "lambda",
      ylab = "component-wise G config",
      main = paste0("Held-out log likelihood by component-wise clusters, H=", H_now)
    )
    axis(1, at = seq_along(lambdas), labels = lambdas)
    axis(2, at = seq_along(configs), labels = configs, las = 2, cex.axis = 0.62)
    box()

    for (ix in seq_along(lambdas)) {
      for (iy in seq_along(configs)) {
        if (is.finite(z[ix, iy])) {
          text(ix, iy, labels = sprintf("%.3f", z[ix, iy]), cex = 0.48)
        }
      }
    }
  }
}

plot_coordinate_component_effects <- function(d, path) {
  d <- d[d$H >= 2L, , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(d))) {
    vals <- parse_G(d$G_config[i], d$H[i])
    rows[[i]] <- data.frame(
      H = d$H[i],
      factor = seq_along(vals),
      component_count = vals,
      lambda_l1_penalty = d$lambda_l1_penalty[i],
      heldout_ll = rep(d$mean_heldout_loglik_per_response[i], length(vals))
    )
  }
  long <- do.call(rbind, rows)
  H_vals <- sort(unique(long$H))

  png(path, width = 1900, height = 480 * length(H_vals), res = 180)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  par(mfrow = c(length(H_vals), 1), mar = c(4.2, 4.2, 3, 1))

  for (H_now in H_vals) {
    dd <- long[long$H == H_now, , drop = FALSE]
    agg <- aggregate(
      heldout_ll ~ factor + component_count,
      dd,
      function(x) c(mean = mean(x), q25 = quantile(x, 0.25), q75 = quantile(x, 0.75))
    )
    y_mean <- agg$heldout_ll[, 1L]
    y_q25 <- agg$heldout_ll[, 2L]
    y_q75 <- agg$heldout_ll[, 3L]
    x <- agg$factor + (agg$component_count - 2) * 0.18
    ylim <- range(c(y_q25, y_q75), na.rm = TRUE)
    plot(
      NA,
      xlim = c(0.5, H_now + 0.5),
      ylim = ylim,
      xaxt = "n",
      xlab = "factor coordinate",
      ylab = "held-out log likelihood per response",
      main = paste0("Marginal CV performance by component count, H=", H_now)
    )
    axis(1, at = seq_len(H_now), labels = paste0("F", seq_len(H_now)))
    grid(nx = NA, ny = NULL, col = "#E5E7EB")
    for (g in sort(unique(agg$component_count))) {
      idx <- agg$component_count == g
      segments(x[idx], y_q25[idx], x[idx], y_q75[idx], col = component_cols[as.character(g)], lwd = 2)
      points(x[idx], y_mean[idx], pch = 19, col = component_cols[as.character(g)], cex = 1.1)
    }
    legend(
      "bottomright",
      legend = paste("G =", names(component_cols)),
      col = component_cols,
      pch = 19,
      bty = "n",
      horiz = TRUE
    )
  }
}

plot_top_component_configs(
  summary_scores,
  top_n,
  file.path(out_dir, "componentwise_cv_top_configs.png")
)

best_slice <- summary_scores[which.max(summary_scores$mean_heldout_loglik_per_response), ]
slice_scores <- summary_scores[
  summary_scores$H == best_slice$H &
    summary_scores$lambda_l1_penalty == best_slice$lambda_l1_penalty,
  ,
  drop = FALSE
]
slice_scores <- slice_scores[order(slice_scores$mean_heldout_loglik_per_response, decreasing = TRUE), ]
plot_top_component_configs(
  slice_scores,
  nrow(slice_scores),
  file.path(
    out_dir,
    paste0(
      "componentwise_cv_all_configs_H",
      best_slice$H,
      "_lambda",
      best_slice$lambda_l1_penalty,
      ".png"
    )
  )
)

plot_H_lambda_config_heatmaps(
  summary_scores,
  file.path(out_dir, "componentwise_cv_metric_heatmaps_by_H_lambda.png")
)

plot_coordinate_component_effects(
  summary_scores,
  file.path(out_dir, "componentwise_cv_coordinate_component_effects.png")
)

cat("Wrote component-wise CV summary and plots to:\n", normalizePath(out_dir), "\n", sep = "")
