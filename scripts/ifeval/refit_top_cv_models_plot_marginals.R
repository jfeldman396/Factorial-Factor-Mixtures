#!/usr/bin/env Rscript

# Refit the top component-wise IFEval CV models on the full analysis matrix and
# plot the fitted marginal factor distributions with mixture components overlaid.
#
# This script is intentionally separate from the long CV sweep.  It reads the
# checkpointed fold-score CSV, chooses the best complete CV settings, refits
# those settings on the full binary matrix, and writes marginal-mixture plots.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

score_path <- Sys.getenv(
  "SCORE_PATH",
  file.path(
    repo_root,
    "results", "full", "ifeval_threshold_sensitivity_2thirds_20260825",
    "threshold_0p5", "componentwise_cv", "ifeval_rank_lambda_cv_fold_scores.csv"
  )
)
matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  file.path(repo_root, "data", "ifeval_threshold_0p5", "openeval_ifeval_only_binary_matrix.csv")
)
item_metadata_path <- Sys.getenv(
  "ITEM_METADATA_PATH",
  file.path(repo_root, "data", "ifeval_threshold_0p5", "openeval_item_metadata.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(
    repo_root,
    "results", "full", "ifeval_threshold_sensitivity_2thirds_20260825",
    "threshold_0p5", "top3_full_data_refits"
  )
)
top_n <- as.integer(Sys.getenv("TOP_N", "3"))
workers <- as.integer(Sys.getenv("WORKERS", "4"))
pretrain_aug_iter <- as.integer(Sys.getenv("PRETRAIN_AUG_ITER", "20"))
refine_iter <- as.integer(Sys.getenv("REFINE_ITER", "20"))
mixture_max_iter <- as.integer(Sys.getenv("MIXTURE_MAX_ITER", "200"))

if (!file.exists(score_path)) stop("Cannot find SCORE_PATH: ", score_path)
if (!file.exists(matrix_path)) stop("Cannot find MATRIX_PATH: ", matrix_path)
if (!file.exists(item_metadata_path)) stop("Cannot find ITEM_METADATA_PATH: ", item_metadata_path)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(score_path, stringsAsFactors = FALSE, check.names = FALSE)
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
    se_heldout_loglik_per_response = if (nrow(d) > 1L) {
      stats::sd(d$heldout_loglik_per_response) / sqrt(nrow(d))
    } else {
      NA_real_
    },
    mean_fit_seconds = mean(d$fit_seconds),
    stringsAsFactors = FALSE
  )
})
cv_summary <- do.call(rbind, summary_list)
cv_summary <- cv_summary[cv_summary$n_folds == 3L, , drop = FALSE]
cv_summary <- cv_summary[order(cv_summary$mean_heldout_loglik_per_response, decreasing = TRUE), ]
top_models <- head(cv_summary, min(top_n, nrow(cv_summary)))

write.csv(
  top_models,
  file.path(out_dir, "top_cv_models_refit_targets.csv"),
  row.names = FALSE
)

safe_label <- function(G_config, lambda) {
  paste0("H", top_models$H[1L])
  gsub("[^A-Za-z0-9]+", "-", paste0("G", G_config, "_lambda", lambda))
}

cols <- c("#2F6DAE", "#C43C4A", "#2F855A", "#805AD5", "#B7791F")

run_fit_script_with_env <- function(env_values) {
  env_names <- names(env_values)
  old_env <- Sys.getenv(env_names, unset = NA_character_)
  names(old_env) <- env_names
  on.exit({
    for (nm in env_names) {
      if (is.na(old_env[[nm]])) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, as.list(stats::setNames(old_env[[nm]], nm)))
      }
    }
  }, add = TRUE)

  do.call(Sys.setenv, as.list(env_values))
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(script_dir)
  system2("Rscript", "fit_interpret_ifeval_mixture.R")
}

plot_marginals_for_fit <- function(fit_path, out_file, title_prefix = "") {
  fit <- readRDS(fit_path)
  H <- fit$H

  png(out_file, width = 1800, height = max(1000, 520 * ceiling(H / 2)), res = 180)
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  n_col <- min(2L, H)
  n_row <- ceiling(H / n_col)
  par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  for (h in seq_len(H)) {
    y <- fit$F_hat[, h]
    mf <- fit$mixture_fits[[h]]
    ord <- order(mf$mu)
    pi_ord <- mf$pi[ord]
    mu_ord <- mf$mu[ord]
    var_ord <- mf$var[ord]

    score_sd <- stats::sd(y)
    if (!is.finite(score_sd) || score_sd <= 0) score_sd <- 1
    xg <- seq(min(y) - 0.5 * score_sd, max(y) + 0.5 * score_sd, length.out = 700)
    comp <- sapply(seq_along(pi_ord), function(g) {
      pi_ord[g] * stats::dnorm(xg, mean = mu_ord[g], sd = sqrt(var_ord[g]))
    })
    if (is.null(dim(comp))) comp <- matrix(comp, ncol = 1L)
    mix <- rowSums(comp)

    hist(
      y,
      breaks = 18,
      freq = FALSE,
      col = "#F3F4F6",
      border = "#D1D5DB",
      main = paste0("F", h, " marginal mixture"),
      xlab = paste0("F", h, " score"),
      ylab = "density"
    )
    lines(xg, mix, col = "#111827", lwd = 2.5, lty = 2)
    for (g in seq_along(pi_ord)) {
      lines(xg, comp[, g], col = cols[g], lwd = 2)
      abline(v = mu_ord[g], col = cols[g], lwd = 1.5, lty = 3)
    }
    legend(
      "topright",
      legend = c("mixture density", paste0("cluster ", seq_along(pi_ord))),
      col = c("#111827", cols[seq_along(pi_ord)]),
      lty = c(2, rep(1, length(pi_ord))),
      lwd = c(2.5, rep(2, length(pi_ord))),
      bty = "n",
      cex = 0.72
    )
  }

  mtext(title_prefix, outer = TRUE, cex = 1.15, font = 2)
  invisible(out_file)
}

fit_records <- vector("list", nrow(top_models))
for (i in seq_len(nrow(top_models))) {
  H_i <- as.integer(top_models$H[i])
  G_i <- top_models$G_config[i]
  lambda_i <- top_models$lambda_l1_penalty[i]
  label_i <- gsub("[^A-Za-z0-9]+", "-", paste0("rank", i, "_H", H_i, "_G", G_i, "_lambda", lambda_i))
  model_dir <- file.path(out_dir, label_i)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  message("Refitting top model ", i, ": H=", H_i, ", G=(", G_i, "), lambda=", lambda_i)
  env <- c(
    MATRIX_PATH = matrix_path,
    ITEM_METADATA_PATH = item_metadata_path,
    OUT_DIR = model_dir,
    H_FIXED = as.character(H_i),
    G_FIXED = G_i,
    WORKERS = as.character(workers),
    PRETRAIN_AUG_ITER = as.character(pretrain_aug_iter),
    PRETRAIN_Z_UPDATE = "sample",
    REFINE_ITER = as.character(refine_iter),
    MIXTURE_MAX_ITER = as.character(mixture_max_iter),
    REQUIRE_MIXTURE_CONVERGENCE = "TRUE",
    REFINEMENT_LAMBDA_L1_PENALTY = as.character(lambda_i)
  )
  status <- run_fit_script_with_env(env)
  if (!identical(status, 0L)) {
    stop("Full-data refit failed for model ", i, " with status ", status)
  }

  fit_files <- list.files(model_dir, pattern = "_fit[.]rds$", full.names = TRUE)
  if (length(fit_files) != 1L) {
    stop("Expected one *_fit.rds in ", model_dir, "; found ", length(fit_files))
  }

  plot_file <- file.path(model_dir, paste0(label_i, "_factor_marginal_mixtures.png"))
  plot_marginals_for_fit(
    fit_files,
    plot_file,
    title_prefix = paste0(
      "IFEval full-data refit: rank ", i,
      ", H=", H_i,
      ", G=(", G_i, ")",
      ", lambda=", lambda_i,
      ", CV ll/resp=", sprintf("%.6f", top_models$mean_heldout_loglik_per_response[i])
    )
  )

  fit_records[[i]] <- data.frame(
    rank = i,
    H = H_i,
    G_config = G_i,
    lambda_l1_penalty = lambda_i,
    mean_heldout_loglik_per_response = top_models$mean_heldout_loglik_per_response[i],
    se_heldout_loglik_per_response = top_models$se_heldout_loglik_per_response[i],
    model_dir = model_dir,
    fit_path = fit_files,
    marginal_plot = plot_file,
    stringsAsFactors = FALSE
  )
}

fit_index <- do.call(rbind, fit_records)
write.csv(fit_index, file.path(out_dir, "top_cv_models_full_refit_index.csv"), row.names = FALSE)

message("Refit index written to: ", file.path(out_dir, "top_cv_models_full_refit_index.csv"))
message("Output directory: ", normalizePath(out_dir, mustWork = FALSE))
