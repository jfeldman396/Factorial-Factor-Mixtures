#!/usr/bin/env Rscript

# Fit the H=3, G=3 independent-mixture probit model to an IFEval-like
# simulated dataset and compare the fitted parameters to the known DGP.

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

get_env <- function(name, default, FUN = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  FUN(value)
}

sim_dir <- get_env(
  "SIM_DIR",
  file.path(bundle_root, "results", "ifeval_like_H3_G3_simulated_dgp_n400"),
  as.character
)
out_dir <- get_env(
  "OUT_DIR",
  file.path(bundle_root, "results", "ifeval_like_H3_G3_recovery_n400"),
  as.character
)
workers <- get_env("WORKERS", 2L, as.integer)
seed <- get_env("SEED", 20260811L, as.integer)
lambda_l1_penalty <- get_env("REFINEMENT_LAMBDA_L1_PENALTY", 10, as.numeric)
pretrain_iter <- get_env("PRETRAIN_AUG_ITER", 8L, as.integer)
refine_iter <- get_env("REFINE_ITER", 8L, as.integer)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_matrix_with_id <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  x <- as.matrix(raw[, -1L, drop = FALSE])
  rownames(x) <- raw[[1L]]
  storage.mode(x) <- "numeric"
  x
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L || sd(x) <= 1e-12 || sd(y) <= 1e-12) return(NA_real_)
  suppressWarnings(cor(x, y))
}

permutation_matrix <- function(idx) {
  out <- matrix(0, length(idx), length(idx))
  for (j in seq_along(idx)) out[idx[j], j] <- 1
  out
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], all_permutations(x[-i]))
  }))
}

best_signed_permutation <- function(L_true, L_est) {
  H <- ncol(L_true)
  perms <- all_permutations(seq_len(H))
  best <- NULL
  for (rr in seq_len(nrow(perms))) {
    perm <- perms[rr, ]
    for (mask in seq_len(2^H) - 1L) {
      signs <- ifelse(as.integer(intToBits(mask))[seq_len(H)] == 1L, -1, 1)
      L_aligned <- sweep(L_est[, perm, drop = FALSE], 2L, signs, "*")
      score <- sum((L_true - L_aligned)^2)
      if (is.null(best) || score < best$score) {
        best <- list(perm = perm, signs = signs, score = score, Lambda = L_aligned)
      }
    }
  }
  best
}

align_matrix <- function(M, align) {
  sweep(M[, align$perm, drop = FALSE], 2L, align$signs, "*")
}

extract_ordered_mixture <- function(mixture_fits, align, true_mixture) {
  rows <- list()
  idx <- 0L
  H <- length(mixture_fits)
  for (h in seq_len(H)) {
    est_fit <- mixture_fits[[align$perm[h]]]
    est_mu <- align$signs[h] * est_fit$mu
    est_var <- est_fit$var
    est_pi <- est_fit$pi / sum(est_fit$pi)
    est_ord <- order(est_mu)
    true_h <- true_mixture[true_mixture$factor == paste0("F", h), ]
    true_h <- true_h[order(true_h$mean), ]
    for (g in seq_len(nrow(true_h))) {
      idx <- idx + 1L
      eg <- est_ord[g]
      rows[[idx]] <- data.frame(
        factor = paste0("F", h),
        component = g,
        true_weight = true_h$weight[g],
        est_weight = est_pi[eg],
        true_mean = true_h$mean[g],
        est_mean = est_mu[eg],
        true_variance = true_h$variance[g],
        est_variance = est_var[eg]
      )
    }
  }
  do.call(rbind, rows)
}

binary_matrix_path <- file.path(sim_dir, "ifeval_like_H3_G3_binary_matrix.csv")
factor_path <- file.path(sim_dir, "ifeval_like_H3_G3_true_factors.csv")
item_path <- file.path(sim_dir, "ifeval_like_H3_G3_dgp_item_parameters.csv")
mixture_path <- file.path(sim_dir, "ifeval_like_H3_G3_dgp_mixture_parameters.csv")

X <- read_matrix_with_id(binary_matrix_path)
F_true <- read_matrix_with_id(factor_path)
item_params <- read.csv(item_path, check.names = FALSE)
mixture_true <- read.csv(mixture_path, check.names = FALSE)

Lambda_true <- as.matrix(item_params[, c("F1", "F2", "F3")])
rownames(Lambda_true) <- item_params$item_id
alpha_true <- item_params$alpha
names(alpha_true) <- item_params$item_id

message("Fitting simulated IFEval-like data from: ", sim_dir)
message("n=", nrow(X), ", p=", ncol(X), ", H=3, G=3, lambda_l1=", lambda_l1_penalty)

set.seed(seed)
t0 <- proc.time()[["elapsed"]]
pre <- fit_binary_probit_pretraining_intercept(
  X = X,
  H = 3L,
  G_fixed = 3L,
  n_aug_iter = pretrain_iter,
  z_update = "expectation",
  n_random_starts = 1L,
  max_outer = 4L,
  n_mix_starts = 3L,
  mixture_update = "map",
  mu_prior_kappa = 0.05,
  var_prior_shape = 4,
  var_prior_scale = 0.35,
  weight_prior_alpha = 1.2,
  loading_penalty = 0.05,
  objective_tolerance = 5e-4,
  objective_tolerance_scale = "per_response",
  min_aug_iter = min(5L, pretrain_iter),
  parallel = workers > 1L,
  workers = workers,
  seed = seed,
  verbose = FALSE
)

ref <- fit_binary_probit_refinement_intercept(
  X = X,
  pretrain_fit = pre,
  n_refine_iter = refine_iter,
  maxit_per_subject = 60L,
  n_mix_starts = 3L,
  min_mixture_var = 0.05,
  mixture_update = "map",
  mu_prior_kappa = 0.05,
  var_prior_shape = 4,
  var_prior_scale = 0.35,
  weight_prior_alpha = 1.2,
  mixture_prior_weight = 0.2,
  lambda_l1_penalty = lambda_l1_penalty,
  objective_tolerance = 2e-4,
  objective_tolerance_scale = "relative_total",
  min_refine_iter = min(4L, refine_iter),
  keep_best_binary_iterate = TRUE,
  parallel = workers > 1L,
  workers = workers,
  verbose = FALSE
)
elapsed <- proc.time()[["elapsed"]] - t0

fit <- orient_factors_by_accuracy(ref)
rownames(fit$F_hat) <- rownames(X)
rownames(fit$Lambda_hat) <- colnames(X)
names(fit$alpha_hat) <- colnames(X)
align <- best_signed_permutation(Lambda_true[colnames(X), , drop = FALSE], fit$Lambda_hat[colnames(X), , drop = FALSE])
Lambda_aligned <- align$Lambda
F_aligned <- align_matrix(fit$F_hat[rownames(F_true), , drop = FALSE], align)
alpha_est <- fit$alpha_hat[colnames(X)]

eta_true <- sweep(F_true %*% t(Lambda_true[colnames(X), , drop = FALSE]), 2L, alpha_true[colnames(X)], "+")
eta_est <- sweep(F_aligned %*% t(Lambda_aligned[colnames(X), , drop = FALSE]), 2L, alpha_est[colnames(X)], "+")
prob_true <- pnorm(eta_true)
prob_est <- pnorm(eta_est)

mix_cmp <- extract_ordered_mixture(fit$mixture_fits, align, mixture_true)
mix_cmp$mean_error <- mix_cmp$est_mean - mix_cmp$true_mean
mix_cmp$variance_error <- mix_cmp$est_variance - mix_cmp$true_variance
mix_cmp$weight_error <- mix_cmp$est_weight - mix_cmp$true_weight

factor_summary <- data.frame(
  factor = paste0("F", seq_len(3L)),
  estimated_column = align$perm,
  sign = align$signs,
  factor_cor = vapply(seq_len(3L), function(h) safe_cor(F_true[, h], F_aligned[, h]), numeric(1)),
  lambda_cor = vapply(seq_len(3L), function(h) safe_cor(Lambda_true[, h], Lambda_aligned[, h]), numeric(1)),
  lambda_rmse = sqrt(colMeans((Lambda_true - Lambda_aligned)^2)),
  true_lambda_l2 = sqrt(colSums(Lambda_true^2)),
  est_lambda_l2 = sqrt(colSums(Lambda_aligned^2))
)

summary <- data.frame(
  sim_dir = sim_dir,
  n = nrow(X),
  p = ncol(X),
  H = 3L,
  G = 3L,
  lambda_l1_penalty = lambda_l1_penalty,
  elapsed_sec = elapsed,
  pretraining_converged = isTRUE(pre$pretraining$converged),
  pretraining_completed_iter = pre$pretraining$n_completed,
  refinement_converged = isTRUE(fit$joint_refinement$converged),
  refinement_completed_iter = fit$joint_refinement$n_completed,
  selected_refinement_iteration = fit$joint_refinement$selected_iteration,
  alpha_cor = safe_cor(alpha_true[colnames(X)], alpha_est[colnames(X)]),
  alpha_rmse = sqrt(mean((alpha_true[colnames(X)] - alpha_est[colnames(X)])^2)),
  lambda_cor_flat = safe_cor(as.numeric(Lambda_true), as.numeric(Lambda_aligned)),
  lambda_rmse_flat = sqrt(mean((Lambda_true - Lambda_aligned)^2)),
  factor_mean_abs_cor = mean(abs(factor_summary$factor_cor), na.rm = TRUE),
  factor_min_abs_cor = min(abs(factor_summary$factor_cor), na.rm = TRUE),
  eta_cor_flat = safe_cor(as.numeric(eta_true), as.numeric(eta_est)),
  eta_rmse_flat = sqrt(mean((eta_true - eta_est)^2)),
  prob_cor_flat = safe_cor(as.numeric(prob_true), as.numeric(prob_est)),
  prob_rmse_flat = sqrt(mean((prob_true - prob_est)^2)),
  mixture_weight_cor = safe_cor(mix_cmp$true_weight, mix_cmp$est_weight),
  mixture_mean_cor = safe_cor(mix_cmp$true_mean, mix_cmp$est_mean),
  mixture_variance_cor = safe_cor(mix_cmp$true_variance, mix_cmp$est_variance),
  mixture_weight_rmse = sqrt(mean((mix_cmp$true_weight - mix_cmp$est_weight)^2)),
  mixture_mean_rmse = sqrt(mean((mix_cmp$true_mean - mix_cmp$est_mean)^2)),
  mixture_variance_rmse = sqrt(mean((mix_cmp$true_variance - mix_cmp$est_variance)^2))
)

write.csv(summary, file.path(out_dir, "ifeval_like_H3_G3_recovery_summary.csv"), row.names = FALSE)
write.csv(factor_summary, file.path(out_dir, "ifeval_like_H3_G3_factor_alignment_summary.csv"), row.names = FALSE)
write.csv(mix_cmp, file.path(out_dir, "ifeval_like_H3_G3_mixture_parameter_recovery.csv"), row.names = FALSE)
write.csv(
  data.frame(item_id = colnames(X), alpha_true = alpha_true[colnames(X)], alpha_est = alpha_est[colnames(X)]),
  file.path(out_dir, "ifeval_like_H3_G3_alpha_recovery.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(item_id = rownames(Lambda_true), Lambda_true, Lambda_aligned, check.names = FALSE),
  file.path(out_dir, "ifeval_like_H3_G3_lambda_recovery_aligned.csv"),
  row.names = FALSE
)
saveRDS(fit, file.path(out_dir, "ifeval_like_H3_G3_fit.rds"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  mix_long <- rbind(
    data.frame(mix_cmp[, c("factor", "component")], source = "truth", variance = mix_cmp$true_variance),
    data.frame(mix_cmp[, c("factor", "component")], source = "estimate", variance = mix_cmp$est_variance)
  )
  p_mix <- ggplot(mix_long, aes(x = factor(component), y = variance, fill = source)) +
    geom_col(position = "dodge", width = 0.72) +
    facet_wrap(~factor) +
    labs(title = "IFEval-like simulated DGP: mixture variance recovery", x = "ordered component", y = "variance", fill = NULL) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, "ifeval_like_H3_G3_mixture_variance_recovery.png"), p_mix, width = 7.5, height = 4.5, dpi = 200)

  lam_df <- data.frame(lambda_true = as.numeric(Lambda_true), lambda_est = as.numeric(Lambda_aligned))
  p_lam <- ggplot(lam_df, aes(x = lambda_true, y = lambda_est)) +
    geom_point(alpha = 0.45, size = 1) +
    geom_abline(slope = 1, intercept = 0, color = "#b2182b") +
    coord_equal() +
    labs(title = "Aligned Lambda recovery", x = "true lambda", y = "estimated lambda") +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, "ifeval_like_H3_G3_lambda_recovery_scatter.png"), p_lam, width = 5.5, height = 5, dpi = 200)

  prob_df <- data.frame(prob_true = as.numeric(prob_true), prob_est = as.numeric(prob_est))
  p_prob <- ggplot(prob_df[sample.int(nrow(prob_df), min(10000L, nrow(prob_df))), ], aes(x = prob_true, y = prob_est)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "#b2182b") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "Probability recovery", x = "true probability", y = "estimated probability") +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, "ifeval_like_H3_G3_probability_recovery_scatter.png"), p_prob, width = 5.5, height = 5, dpi = 200)
}

cat("\nRecovery summary:\n")
print(summary, row.names = FALSE)
cat("\nFactor alignment:\n")
print(factor_summary, row.names = FALSE)
cat("\nMixture parameter recovery:\n")
print(mix_cmp, row.names = FALSE)
cat("\nOutputs written to:", normalizePath(out_dir), "\n")
