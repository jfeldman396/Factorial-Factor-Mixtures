#!/usr/bin/env Rscript

# Tune the entrywise L1 sparsity penalty on Lambda for the selected IFEval
# independent-mixture probit factor model.
#
# The rank and mixture complexity are held fixed at H=3, G=3.  We pretrain once,
# then run refinement across a grid of lambda_l1 penalties.  The tuning summary
# reports the likelihood/sparsity tradeoff and selects a penalty by a BIC-style
# criterion using the number of nonzero loading entries as the varying degrees
# of freedom.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

source(file.path(script_dir, "math500_intercept_imfm_fit.R"))

parse_num_grid <- function(x, default) {
  if (is.na(x) || !nzchar(x)) return(default)
  as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))
}

read_binary_matrix <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- raw[[1L]]
  X <- X[, colSums(is.na(X)) == 0L, drop = FALSE]
  X <- X[, colSums(X == 1) > 0L & colSums(X == 0) > 0L, drop = FALSE]
  X
}

matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  file.path(bundle_root, "data", "openeval_ifeval_only_binary_matrix.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(bundle_root, "results", "reproduced_openeval_ifeval_lambda_sparsity_tuning")
)
H_fixed <- as.integer(Sys.getenv("H_FIXED", "3"))
G_fixed <- as.integer(Sys.getenv("G_FIXED", "3"))
workers <- as.integer(Sys.getenv("WORKERS", "8"))
penalty_grid <- parse_num_grid(Sys.getenv("LAMBDA_L1_GRID", unset = ""), c(0, 0.5, 1, 2, 3, 4, 6, 8, 10, 12))
zero_threshold <- as.numeric(Sys.getenv("ZERO_THRESHOLD", unset = "1e-6"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

X <- read_binary_matrix(matrix_path)
n <- nrow(X)
p <- ncol(X)
N <- length(X)

set.seed(20260724L)
message("Pretraining once for lambda sparsity tuning.")
pre <- fit_binary_probit_pretraining_intercept(
  X = X,
  H = H_fixed,
  G_fixed = G_fixed,
  n_aug_iter = as.integer(Sys.getenv("PRETRAIN_AUG_ITER", "8")),
  z_update = "expectation",
  n_random_starts = 1L,
  max_outer = as.integer(Sys.getenv("PRETRAIN_MAX_OUTER", "4")),
  n_mix_starts = 3L,
  mixture_update = "map",
  mu_prior_kappa = 0.05,
  var_prior_shape = 4,
  var_prior_scale = 0.35,
  weight_prior_alpha = 1.2,
  loading_penalty = as.numeric(Sys.getenv("PRETRAIN_LOADING_PENALTY", "0.05")),
  objective_tolerance = 5e-4,
  objective_tolerance_scale = "per_response",
  min_aug_iter = as.integer(Sys.getenv("PRETRAIN_MIN_AUG_ITER", "5")),
  parallel = TRUE,
  workers = workers,
  seed = 20260724L,
  verbose = FALSE
)

rows <- vector("list", length(penalty_grid))

for (idx in seq_along(penalty_grid)) {
  pen <- penalty_grid[idx]
  message("Lambda sparsity tuning: penalty=", pen, " (", idx, "/", length(penalty_grid), ")")
  t0 <- proc.time()[["elapsed"]]
  ref <- fit_binary_probit_refinement_intercept(
    X = X,
    pretrain_fit = pre,
    n_refine_iter = as.integer(Sys.getenv("REFINE_ITER", "8")),
    maxit_per_subject = as.integer(Sys.getenv("MAXIT_PER_SUBJECT", "60")),
    n_mix_starts = 3L,
    min_mixture_var = 0.05,
    mixture_update = "map",
    mu_prior_kappa = 0.05,
    var_prior_shape = 4,
    var_prior_scale = 0.35,
    weight_prior_alpha = 1.2,
    mixture_prior_weight = 0.2,
    lambda_l1_penalty = pen,
    objective_tolerance = 2e-4,
    objective_tolerance_scale = "relative_total",
    min_refine_iter = as.integer(Sys.getenv("MIN_REFINE_ITER", "4")),
    keep_best_binary_iterate = TRUE,
    parallel = TRUE,
    workers = workers,
    verbose = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  ll <- binary_probit_loglik_alpha(ref$X, ref$F_hat, ref$Lambda_hat, ref$alpha_hat)
  abs_L <- abs(ref$Lambda_hat)
  nnz <- sum(abs_L > zero_threshold)
  mixture_df <- H_fixed * (3L * G_fixed - 1L)
  effective_df <- p + n * H_fixed + nnz + mixture_df
  aic <- -2 * ll + 2 * effective_df
  bic <- -2 * ll + log(N) * effective_df

  rows[[idx]] <- data.frame(
    lambda_l1_penalty = pen,
    binary_loglik = ll,
    binary_loglik_per_response = ll / N,
    effective_df = effective_df,
    nnz_lambda = nnz,
    lambda_density = nnz / length(ref$Lambda_hat),
    n_abs_gt_0_1 = sum(abs_L > 0.1),
    n_abs_gt_0_5 = sum(abs_L > 0.5),
    n_abs_gt_1_0 = sum(abs_L > 1.0),
    mean_abs_lambda = mean(abs_L),
    median_abs_lambda = median(abs_L),
    aic = aic,
    bic = bic,
    refinement_iterations = ref$joint_refinement$n_completed,
    selected_refinement_iteration = ref$joint_refinement$selected_iteration,
    elapsed_sec = elapsed
  )

  saveRDS(
    ref,
    file.path(out_dir, sprintf("ifeval_H%d_G%d_lambda_l1_%s_fit.rds", H_fixed, G_fixed, gsub("[.]", "p", as.character(pen))))
  )
}

summary <- do.call(rbind, rows)
summary$selected_by_bic <- summary$bic == min(summary$bic)
summary$selected_by_aic <- summary$aic == min(summary$aic)
summary <- summary[order(summary$lambda_l1_penalty), ]
write.csv(summary, file.path(out_dir, "ifeval_lambda_sparsity_tuning_summary.csv"), row.names = FALSE)

selected <- summary[which.min(summary$bic), , drop = FALSE]
write.csv(selected, file.path(out_dir, "ifeval_selected_lambda_sparsity_penalty.csv"), row.names = FALSE)
writeLines(as.character(selected$lambda_l1_penalty[1L]), file.path(out_dir, "selected_lambda_penalty.txt"))

png(file.path(out_dir, "ifeval_lambda_sparsity_tuning_path.png"), width = 1700, height = 950, res = 160)
op <- par(mfrow = c(2, 2), mar = c(5, 5, 3, 1), oma = c(0, 0, 2, 0))
plot(summary$lambda_l1_penalty, summary$binary_loglik_per_response,
     type = "b", pch = 19, col = "#2B6CB0",
     xlab = "entrywise Lambda L1 penalty", ylab = "binary loglik per response",
     main = "Fit")
abline(v = selected$lambda_l1_penalty, lty = 2, col = "#BC4749")
plot(summary$lambda_l1_penalty, summary$bic,
     type = "b", pch = 19, col = "#4A5568",
     xlab = "entrywise Lambda L1 penalty", ylab = "BIC-style criterion",
     main = "Selection criterion")
abline(v = selected$lambda_l1_penalty, lty = 2, col = "#BC4749")
plot(summary$lambda_l1_penalty, summary$lambda_density,
     type = "b", pch = 19, col = "#C05621",
     xlab = "entrywise Lambda L1 penalty", ylab = "nonzero loading fraction",
     main = "Sparsity")
abline(v = selected$lambda_l1_penalty, lty = 2, col = "#BC4749")
plot(summary$lambda_l1_penalty, summary$n_abs_gt_0_5,
     type = "b", pch = 19, col = "#805AD5",
     xlab = "entrywise Lambda L1 penalty", ylab = "number |lambda| > 0.5",
     main = "Large loadings")
abline(v = selected$lambda_l1_penalty, lty = 2, col = "#BC4749")
mtext("IFEval Lambda sparsity tuning path; dashed line = BIC-selected penalty", outer = TRUE, font = 2)
par(op)
dev.off()

cat("\nLambda sparsity tuning summary:\n")
print(summary)
cat("\nSelected penalty by BIC-style criterion: ", selected$lambda_l1_penalty, "\n", sep = "")
cat("Outputs saved in: ", normalizePath(out_dir), "\n", sep = "")
