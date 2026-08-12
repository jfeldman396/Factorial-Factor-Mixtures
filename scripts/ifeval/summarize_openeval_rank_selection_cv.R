#!/usr/bin/env Rscript

# Summarize OpenEval random-cell holdout rank selection for:
#   1. independent-mixture probit factor model, G = 2;
#   2. independent-mixture probit factor model, G = 3;
#   3. ordinary Gaussian probit factor model.
#
# The script reports three rank diagnostics:
#
#   held-out predictive log likelihood:
#     mean log p(X_heldout | fit on training cells)
#
#   held-out BIC:
#     -2 log p(X_heldout | fit on training cells) + df log(N_heldout)
#
#   training BIC:
#     -2 log p(X_training | fit on training cells) + df log(N_training)
#
# The held-out BIC is best read as penalized held-out deviance; it is useful
# because it puts the held-out predictive score on a complexity-penalized scale.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
selected_table_dir <- file.path(repo_root, "results", "selected_tables", "ifeval")

matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  unset = file.path(repo_root, "data", "ifeval", "openeval_ifeval_only_binary_matrix.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(repo_root, "results", "reproduced_openeval_ifeval_rank_selection_comparison")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_binary_matrix_dims <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(X) <- "numeric"
  X <- X[, colSums(is.na(X)) == 0L, drop = FALSE]
  X <- X[, colSums(X == 1) > 0L & colSums(X == 0) > 0L, drop = FALSE]
  c(n = nrow(X), p = ncol(X), N = length(X))
}

method_component_count <- function(method) {
  if (grepl("G2", method)) return(2L)
  if (grepl("G3", method)) return(3L)
  NA_integer_
}

effective_df <- function(method, H, n, p) {
  # Count item intercepts, subject factor scores, and loading entries.  For the
  # independent-mixture model, also count each factor's marginal mixture
  # parameters: (G - 1) weights, G means, and G variances.
  H <- as.integer(H)
  df <- p + n * H + p * H
  G <- method_component_count(method)
  if (!is.na(G)) df <- df + H * (3L * G - 1L)
  df
}

read_method_summary <- function(path, method) {
  if (!file.exists(path)) return(NULL)
  d <- read.csv(path, check.names = FALSE)
  if ("method" %in% names(d)) d <- d[d$method == method, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$method <- method
  keep <- c(
    "method",
    "H",
    "mean_heldout_loglik_per_response",
    "sd_heldout_loglik_per_response",
    "se_heldout_loglik_per_response",
    "best_H",
    "selected_by_one_se",
    "one_se_threshold",
    "mean_fit_seconds"
  )
  missing <- setdiff(keep, names(d))
  for (col in missing) d[[col]] <- NA
  d[, keep, drop = FALSE]
}

fold_metrics_from_scores <- function(path, method, n, p, N_total) {
  if (!file.exists(path)) return(NULL)
  folds <- read.csv(path, check.names = FALSE)
  needed <- c(
    "H",
    "fold",
    "heldout_loglik_per_response",
    "train_loglik_per_response",
    "n_heldout"
  )
  if (!all(needed %in% names(folds))) return(NULL)

  out <- do.call(rbind, lapply(seq_len(nrow(folds)), function(i) {
    row <- folds[i, ]
    H <- as.integer(row$H)
    n_heldout <- as.numeric(row$n_heldout)
    n_train <- N_total - n_heldout
    df <- effective_df(method, H, n = n, p = p)
    heldout_ll <- as.numeric(row$heldout_loglik_per_response) * n_heldout
    train_ll <- as.numeric(row$train_loglik_per_response) * n_train
    data.frame(
      method = method,
      H = H,
      fold = as.integer(row$fold),
      effective_df = df,
      n_heldout = n_heldout,
      n_train = n_train,
      heldout_loglik = heldout_ll,
      train_loglik = train_ll,
      heldout_loglik_per_response = as.numeric(row$heldout_loglik_per_response),
      train_loglik_per_response = as.numeric(row$train_loglik_per_response),
      heldout_bic = -2 * heldout_ll + df * log(n_heldout),
      training_bic = -2 * train_ll + df * log(n_train),
      fit_seconds = if ("fit_seconds" %in% names(row)) as.numeric(row$fit_seconds) else NA_real_
    )
  }))
  out
}

summarize_folds_by_H <- function(folds) {
  if (is.null(folds) || nrow(folds) == 0L) return(NULL)
  do.call(rbind, lapply(split(folds, list(folds$method, folds$H), drop = TRUE), function(d) {
    data.frame(
      method = d$method[1L],
      H = d$H[1L],
      effective_df = d$effective_df[1L],
      mean_heldout_loglik_per_response = mean(d$heldout_loglik_per_response),
      sd_heldout_loglik_per_response = stats::sd(d$heldout_loglik_per_response),
      se_heldout_loglik_per_response = stats::sd(d$heldout_loglik_per_response) / sqrt(nrow(d)),
      mean_train_loglik_per_response = mean(d$train_loglik_per_response),
      sd_train_loglik_per_response = stats::sd(d$train_loglik_per_response),
      se_train_loglik_per_response = stats::sd(d$train_loglik_per_response) / sqrt(nrow(d)),
      mean_heldout_bic = mean(d$heldout_bic),
      sd_heldout_bic = stats::sd(d$heldout_bic),
      se_heldout_bic = stats::sd(d$heldout_bic) / sqrt(nrow(d)),
      mean_training_bic = mean(d$training_bic),
      sd_training_bic = stats::sd(d$training_bic),
      se_training_bic = stats::sd(d$training_bic) / sqrt(nrow(d)),
      mean_fit_seconds = mean(d$fit_seconds, na.rm = TRUE)
    )
  }))
}

summary_schema <- c(
  "method",
  "H",
  "effective_df",
  "mean_heldout_loglik_per_response",
  "sd_heldout_loglik_per_response",
  "se_heldout_loglik_per_response",
  "mean_train_loglik_per_response",
  "sd_train_loglik_per_response",
  "se_train_loglik_per_response",
  "mean_heldout_bic",
  "sd_heldout_bic",
  "se_heldout_bic",
  "mean_training_bic",
  "sd_training_bic",
  "se_training_bic",
  "best_H",
  "selected_by_one_se",
  "one_se_threshold",
  "mean_fit_seconds"
)

ensure_summary_schema <- function(d) {
  if (is.null(d)) return(NULL)
  missing <- setdiff(summary_schema, names(d))
  for (col in missing) d[[col]] <- NA
  d[, summary_schema, drop = FALSE]
}

merge_summary_fallback <- function(fold_summary, summary_fallback) {
  if (is.null(summary_fallback)) return(ensure_summary_schema(fold_summary))
  if (is.null(fold_summary)) {
    return(ensure_summary_schema(summary_fallback))
  }
  ensure_summary_schema(fold_summary)
}

add_selection_flags <- function(d) {
  d$best_H <- FALSE
  d$selected_by_one_se <- FALSE
  d$selected_by_heldout_bic <- FALSE
  d$selected_by_training_bic <- FALSE

  for (method in unique(d$method)) {
    idx <- which(d$method == method)
    dm <- d[idx, , drop = FALSE]
    if (any(is.finite(dm$mean_heldout_loglik_per_response))) {
      best_local <- which.max(dm$mean_heldout_loglik_per_response)
      best_global <- idx[best_local]
      threshold <- dm$mean_heldout_loglik_per_response[best_local] -
        dm$se_heldout_loglik_per_response[best_local]
      one_se_candidates <- idx[dm$mean_heldout_loglik_per_response >= threshold]
      d$best_H[best_global] <- TRUE
      d$one_se_threshold[idx] <- threshold
      if (length(one_se_candidates)) {
        d$selected_by_one_se[one_se_candidates[which.min(d$H[one_se_candidates])]] <- TRUE
      }
    }
    if (any(is.finite(dm$mean_heldout_bic))) {
      d$selected_by_heldout_bic[idx[which.min(dm$mean_heldout_bic)]] <- TRUE
    }
    if (any(is.finite(dm$mean_training_bic))) {
      d$selected_by_training_bic[idx[which.min(dm$mean_training_bic)]] <- TRUE
    }
  }
  d
}

method_specs <- data.frame(
  method = c("mixture_probit_G2", "mixture_probit_G3", "ordinary_probit_factor"),
  summary_env = c("G2_H_SUMMARY", "G3_H_SUMMARY", "ORDINARY_H_SUMMARY"),
  fold_env = c("G2_FOLD_SCORES", "G3_FOLD_SCORES", "ORDINARY_FOLD_SCORES"),
  default_summary = c(
    file.path(selected_table_dir, "cv_inputs", "mixture_G2_H_summary.csv"),
    file.path(selected_table_dir, "cv_inputs", "mixture_G3_H_summary.csv"),
    file.path(selected_table_dir, "cv_inputs", "ordinary_probit_H_summary.csv")
  ),
  default_folds = c(
    file.path(selected_table_dir, "cv_inputs", "mixture_G2_fold_scores.csv"),
    file.path(selected_table_dir, "cv_inputs", "mixture_G3_fold_scores.csv"),
    file.path(selected_table_dir, "cv_inputs", "ordinary_probit_factor_fold_scores.csv")
  ),
  stringsAsFactors = FALSE
)

dims <- read_binary_matrix_dims(matrix_path)
n_models <- unname(dims[["n"]])
n_items <- unname(dims[["p"]])
N_total <- unname(dims[["N"]])

fold_rows <- list()
summary_rows <- list()

for (i in seq_len(nrow(method_specs))) {
  spec <- method_specs[i, ]
  summary_path <- Sys.getenv(spec$summary_env, unset = spec$default_summary)
  fold_path <- Sys.getenv(spec$fold_env, unset = spec$default_folds)
  fallback <- read_method_summary(summary_path, spec$method)
  folds <- fold_metrics_from_scores(
    fold_path,
    method = spec$method,
    n = n_models,
    p = n_items,
    N_total = N_total
  )
  fold_rows[[length(fold_rows) + 1L]] <- folds
  summary_rows[[length(summary_rows) + 1L]] <- merge_summary_fallback(
    summarize_folds_by_H(folds),
    fallback
  )
}

fold_metrics <- do.call(rbind, fold_rows[!vapply(fold_rows, is.null, logical(1))])
if (!is.null(fold_metrics) && nrow(fold_metrics) > 0L) {
  write.csv(fold_metrics, file.path(out_dir, "openeval_ifeval_rank_selection_fold_metrics.csv"), row.names = FALSE)
}

summaries <- do.call(rbind, summary_rows[!vapply(summary_rows, is.null, logical(1))])
summaries <- summaries[order(summaries$method, summaries$H), ]
summaries <- add_selection_flags(summaries)
write.csv(summaries, file.path(out_dir, "openeval_ifeval_rank_selection_all_H.csv"), row.names = FALSE)
write.csv(
  summaries[, c(
    "method",
    "H",
    "effective_df",
    "mean_heldout_loglik_per_response",
    "se_heldout_loglik_per_response",
    "mean_heldout_bic",
    "se_heldout_bic",
    "mean_training_bic",
    "se_training_bic",
    "best_H",
    "selected_by_one_se",
    "selected_by_heldout_bic",
    "selected_by_training_bic"
  )],
  file.path(out_dir, "openeval_ifeval_rank_selection_loglik_bic_by_H.csv"),
  row.names = FALSE
)

selected <- do.call(rbind, lapply(split(summaries, summaries$method), function(d) {
  best <- d[d$best_H, ][1L, ]
  one_se <- d[d$selected_by_one_se, ][1L, ]
  heldout_bic <- d[d$selected_by_heldout_bic, ][1L, ]
  training_bic <- d[d$selected_by_training_bic, ][1L, ]
  data.frame(
    method = d$method[1L],
    best_heldout_loglik_H = best$H,
    one_se_H = one_se$H,
    heldout_bic_H = heldout_bic$H,
    training_bic_H = training_bic$H,
    best_mean_heldout_loglik_per_response = best$mean_heldout_loglik_per_response,
    min_mean_heldout_bic = heldout_bic$mean_heldout_bic,
    min_mean_training_bic = training_bic$mean_training_bic,
    mean_fit_seconds_at_best_loglik_H = best$mean_fit_seconds
  )
}))
selected <- selected[order(selected$method), ]
write.csv(selected, file.path(out_dir, "openeval_ifeval_selected_rank_by_method.csv"), row.names = FALSE)

cols <- c(
  mixture_probit_G2 = "#2B6CB0",
  mixture_probit_G3 = "#C05621",
  ordinary_probit_factor = "#4A5568"
)

plot_line_metric <- function(file_name, y_col, se_col, ylab, main, lower_is_better = FALSE) {
  ok <- is.finite(summaries[[y_col]])
  if (!any(ok)) return(invisible(FALSE))
  png(file.path(out_dir, file_name), width = 1500, height = 950, res = 160)
  y_low <- summaries[[y_col]] - if (se_col %in% names(summaries)) summaries[[se_col]] else 0
  y_high <- summaries[[y_col]] + if (se_col %in% names(summaries)) summaries[[se_col]] else 0
  plot(
    range(summaries$H),
    range(c(y_low[ok], y_high[ok]), na.rm = TRUE),
    type = "n",
    xlab = "number of factors H",
    ylab = ylab,
    main = main
  )
  for (m in names(cols)) {
    d <- summaries[summaries$method == m & is.finite(summaries[[y_col]]), ]
    if (!nrow(d)) next
    lines(d$H, d[[y_col]], type = "b", pch = 19, col = cols[m])
    if (se_col %in% names(d) && any(is.finite(d[[se_col]]))) {
      arrows(
        d$H,
        d[[y_col]] - d[[se_col]],
        d$H,
        d[[y_col]] + d[[se_col]],
        angle = 90,
        code = 3,
        length = 0.035,
        col = cols[m]
      )
    }
  }
  legend(if (lower_is_better) "topright" else "bottomleft",
         legend = names(cols), col = cols, pch = 19, lty = 1, bty = "n")
  dev.off()
  invisible(TRUE)
}

plot_line_metric(
  "openeval_ifeval_heldout_loglik_by_H_method.png",
  "mean_heldout_loglik_per_response",
  "se_heldout_loglik_per_response",
  "held-out log likelihood per response",
  "OpenEval rank selection by held-out predictive log score"
)
plot_line_metric(
  "openeval_ifeval_heldout_bic_by_H_method.png",
  "mean_heldout_bic",
  "se_heldout_bic",
  "held-out BIC",
  "OpenEval rank selection by held-out BIC",
  lower_is_better = TRUE
)
plot_line_metric(
  "openeval_ifeval_training_bic_by_H_method.png",
  "mean_training_bic",
  "se_training_bic",
  "training BIC",
  "OpenEval rank selection by training BIC",
  lower_is_better = TRUE
)
plot_line_metric(
  "openeval_ifeval_runtime_by_H_method.png",
  "mean_fit_seconds",
  "missing_se_column",
  "mean fit seconds per fold",
  "OpenEval CV runtime by H"
)

accuracy_by_fold_path <- Sys.getenv(
  "H3_ACCURACY_BY_FOLD",
  unset = file.path(selected_table_dir, "openeval_ifeval_H3_heldout_accuracy_by_fold.csv")
)
if (file.exists(accuracy_by_fold_path)) {
  acc_fold <- read.csv(accuracy_by_fold_path, check.names = FALSE)
  write.csv(
    acc_fold,
    file.path(out_dir, "openeval_ifeval_H3_heldout_accuracy_by_fold.csv"),
    row.names = FALSE
  )

  acc_summary <- do.call(rbind, lapply(split(acc_fold, acc_fold$method), function(d) {
    data.frame(
      method = d$method[1L],
      mean_accuracy = mean(d$acc),
      se_accuracy = stats::sd(d$acc) / sqrt(nrow(d)),
      mean_brier = mean(d$brier),
      se_brier = stats::sd(d$brier) / sqrt(nrow(d)),
      mean_loglik = mean(d$loglik),
      se_loglik = stats::sd(d$loglik) / sqrt(nrow(d))
    )
  }))
  acc_summary <- acc_summary[order(acc_summary$method), ]
  write.csv(
    acc_summary,
    file.path(out_dir, "openeval_ifeval_H3_heldout_accuracy_summary.csv"),
    row.names = FALSE
  )
}

cat("\nSelected ranks by method:\n")
print(selected)
cat("\nAll-H summaries saved in: ", normalizePath(out_dir), "\n", sep = "")
