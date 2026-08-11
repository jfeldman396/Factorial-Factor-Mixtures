#!/usr/bin/env Rscript

# Add IFEval-like observation features to the original H=3, G=3 simulation:
# heterogeneous item intercepts, less balanced loading patterns, and the fitted
# IFEval sparse Lambda.  The latent mixture is kept fixed by default so the
# observation-side features can be compared directly.

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

split_csv <- function(x) {
  y <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  y[nzchar(y)]
}

parse_int_csv <- function(x) as.integer(split_csv(x))
parse_chr_csv <- function(x) split_csv(x)

sim_n <- get_env("SIM_N", 400L, as.integer)
sim_p <- get_env("SIM_P", 500L, as.integer)
rep_values <- get_env("REP_VALUES", 1:2, parse_int_csv)
variants <- get_env(
  "VARIANTS",
  c(
    "balanced_strong_zero_alpha",
    "balanced_strong_ifeval_alpha",
    "balanced_moderate_ifeval_alpha",
    "ifeval_thresholded_alpha",
    "ifeval_full_alpha"
  ),
  parse_chr_csv
)
out_dir <- get_env(
  "OUT_DIR",
  file.path(bundle_root, "results", "ifeval_features_added_to_original_g3_sim"),
  as.character
)
seed_base <- get_env("SEED", 20260812L, as.integer)
workers <- get_env("WORKERS", 1L, as.integer)
lambda_l1_penalty <- get_env("REFINEMENT_LAMBDA_L1_PENALTY", 10, as.numeric)
pretrain_iter <- get_env("PRETRAIN_AUG_ITER", 15L, as.integer)
refine_iter <- get_env("REFINE_ITER", 12L, as.integer)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_path <- file.path(
  bundle_root,
  "results",
  "reproduced_openeval_ifeval_lambda_sparsity_tuning",
  "ifeval_H3_G3_lambda_l1_10_fit.rds"
)
thresholded_lambda_path <- file.path(
  bundle_root,
  "results",
  "loadings_crossloadings",
  "ifeval_sparse_lambda_matrix_threshold_0p5.csv"
)
ifeval_fit <- readRDS(fit_path)
ifeval_alpha <- as.numeric(ifeval_fit$alpha_hat)
names(ifeval_alpha) <- names(ifeval_fit$alpha_hat)
ifeval_full_lambda <- as.matrix(ifeval_fit$Lambda_hat)
thresholded <- read.csv(thresholded_lambda_path, check.names = FALSE)
ifeval_thresholded_lambda <- as.matrix(thresholded[, c("F1", "F2", "F3")])
rownames(ifeval_thresholded_lambda) <- thresholded$item_id
ifeval_thresholded_alpha <- ifeval_alpha[thresholded$item_id]

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L || sd(x) <= 1e-12 || sd(y) <= 1e-12) return(NA_real_)
  suppressWarnings(cor(x, y))
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], all_permutations(x[-i]))))
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
      if (is.null(best) || score < best$score) best <- list(perm = perm, signs = signs, Lambda = L_aligned, score = score)
    }
  }
  best
}

align_matrix <- function(M, align) sweep(M[, align$perm, drop = FALSE], 2L, align$signs, "*")

balanced_block_sizes <- function(p, H) {
  block_sizes <- rep(floor(p / H), H)
  remainder <- p - sum(block_sizes)
  if (remainder > 0L) block_sizes[seq_len(remainder)] <- block_sizes[seq_len(remainder)] + 1L
  block_sizes
}

make_balanced_loading <- function(p, H, primary_range, rare_cross = TRUE) {
  block_sizes <- balanced_block_sizes(p, H)
  block_id <- rep(seq_len(H), times = block_sizes)
  Lambda <- matrix(0, p, H)
  for (j in seq_len(p)) {
    h <- block_id[j]
    Lambda[j, h] <- runif(1, primary_range[1L], primary_range[2L])
    if (rare_cross) {
      for (k in setdiff(seq_len(H), h)) {
        if (runif(1) < 0.035) Lambda[j, k] <- runif(1, 0.12, 0.28)
      }
    }
  }
  rownames(Lambda) <- sprintf("item_%03d", seq_len(p))
  colnames(Lambda) <- paste0("F", seq_len(H))
  list(Lambda = Lambda, block_id = block_id)
}

make_original_g3_params <- function(H) {
  replicate(
    H,
    list(pi = c(0.30, 0.40, 0.30), mu = c(-1, 0, 1), sd = c(0.25, 0.50, 0.75)),
    simplify = FALSE
  )
}

simulate_factors <- function(n, H) {
  params <- make_original_g3_params(H)
  F <- matrix(NA_real_, n, H)
  component <- matrix(NA_integer_, n, H)
  standardized_params <- vector("list", H)
  for (h in seq_len(H)) {
    draw <- sample_standardized_mixture(n, params[[h]])
    F[, h] <- draw$x
    component[, h] <- draw$component
    standardized_params[[h]] <- draw$parameters
  }
  colnames(F) <- paste0("F", seq_len(H))
  colnames(component) <- paste0("F", seq_len(H))
  list(F = F, component = component, mixture_params = standardized_params)
}

make_observation_model <- function(variant, p) {
  if (variant == "balanced_strong_zero_alpha") {
    loading <- make_balanced_loading(p, 3L, c(1.65, 2.35))
    alpha <- rep(0, p)
  } else if (variant == "balanced_strong_ifeval_alpha") {
    loading <- make_balanced_loading(p, 3L, c(1.65, 2.35))
    alpha <- sample(ifeval_alpha, p, replace = p > length(ifeval_alpha))
  } else if (variant == "balanced_moderate_ifeval_alpha") {
    loading <- make_balanced_loading(p, 3L, c(0.75, 1.25))
    alpha <- sample(ifeval_alpha, p, replace = p > length(ifeval_alpha))
  } else if (variant == "ifeval_thresholded_alpha") {
    Lambda <- ifeval_thresholded_lambda
    alpha <- ifeval_thresholded_alpha
    loading <- list(Lambda = Lambda, block_id = max.col(abs(Lambda), ties.method = "first"))
  } else if (variant == "ifeval_full_alpha") {
    Lambda <- ifeval_full_lambda
    alpha <- ifeval_alpha[rownames(Lambda)]
    loading <- list(Lambda = Lambda, block_id = max.col(abs(Lambda), ties.method = "first"))
  } else {
    stop("Unknown variant: ", variant)
  }
  Lambda <- loading$Lambda
  Lambda <- Lambda[seq_len(min(nrow(Lambda), p)), , drop = FALSE]
  if (nrow(Lambda) < p) stop("Requested p exceeds available Lambda rows for variant ", variant)
  alpha <- as.numeric(alpha[seq_len(p)])
  names(alpha) <- rownames(Lambda)
  list(alpha = alpha, Lambda = Lambda, block_id = loading$block_id[seq_len(p)])
}

extract_mixture_recovery <- function(true_params, mixture_fits, align) {
  rows <- list()
  idx <- 0L
  for (h in seq_along(true_params)) {
    fit_h <- mixture_fits[[align$perm[h]]]
    est_mu <- align$signs[h] * fit_h$mu
    est_var <- fit_h$var
    est_weight <- fit_h$pi / sum(fit_h$pi)
    est_order <- order(est_mu)
    true_order <- order(true_params[[h]]$mu)
    for (g in seq_along(true_order)) {
      idx <- idx + 1L
      tg <- true_order[g]
      eg <- est_order[g]
      rows[[idx]] <- data.frame(
        factor = paste0("F", h),
        component = g,
        true_weight = true_params[[h]]$pi[tg],
        est_weight = est_weight[eg],
        true_mean = true_params[[h]]$mu[tg],
        est_mean = est_mu[eg],
        true_variance = true_params[[h]]$sd[tg]^2,
        est_variance = est_var[eg]
      )
    }
  }
  do.call(rbind, rows)
}

fit_one <- function(variant, rep_id) {
  seed <- seed_base + 10000L * match(variant, variants) + rep_id
  set.seed(seed)
  obs <- make_observation_model(variant, sim_p)
  factors <- simulate_factors(sim_n, 3L)
  eta <- sweep(factors$F %*% t(obs$Lambda), 2L, obs$alpha, "+")
  prob <- pnorm(eta)
  X <- 1L * (eta + matrix(rnorm(sim_n * sim_p), sim_n, sim_p) > 0)
  colnames(X) <- rownames(obs$Lambda)
  rownames(X) <- sprintf("sim_%03d", seq_len(sim_n))
  rownames(factors$F) <- rownames(X)

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
    seed = seed + 101L,
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

  align <- best_signed_permutation(obs$Lambda, fit$Lambda_hat[colnames(X), , drop = FALSE])
  Lambda_aligned <- align$Lambda
  F_aligned <- align_matrix(fit$F_hat, align)
  mix <- extract_mixture_recovery(factors$mixture_params, fit$mixture_fits, align)
  eta_est <- sweep(F_aligned %*% t(Lambda_aligned), 2L, fit$alpha_hat[colnames(X)], "+")

  item_active_counts <- rowSums(abs(obs$Lambda) >= 0.5)
  metric <- data.frame(
    variant = variant,
    rep = rep_id,
    n = sim_n,
    p = sim_p,
    lambda_l1_penalty = lambda_l1_penalty,
    alpha_sd = sd(obs$alpha),
    alpha_min = min(obs$alpha),
    alpha_median = median(obs$alpha),
    alpha_max = max(obs$alpha),
    mean_abs_lambda = mean(abs(obs$Lambda)),
    max_abs_lambda = max(abs(obs$Lambda)),
    n_lambda_abs_gt_0_5 = sum(abs(obs$Lambda) >= 0.5),
    n_items_no_abs_gt_0_5 = sum(item_active_counts == 0L),
    n_items_multi_abs_gt_0_5 = sum(item_active_counts >= 2L),
    prob_extreme = mean(prob < 0.05 | prob > 0.95),
    row_extreme = mean(rowMeans(X) <= 0.05 | rowMeans(X) >= 0.95),
    item_extreme = mean(colMeans(X) <= 0.05 | colMeans(X) >= 0.95),
    alpha_cor = safe_cor(obs$alpha, fit$alpha_hat[colnames(X)]),
    lambda_cor = safe_cor(as.numeric(obs$Lambda), as.numeric(Lambda_aligned)),
    factor_mean_abs_cor = mean(abs(vapply(seq_len(3L), function(h) safe_cor(factors$F[, h], F_aligned[, h]), numeric(1)))),
    eta_cor = safe_cor(as.numeric(eta), as.numeric(eta_est)),
    prob_cor = safe_cor(as.numeric(prob), as.numeric(pnorm(eta_est))),
    mixture_weight_cor = safe_cor(mix$true_weight, mix$est_weight),
    mixture_mean_cor = safe_cor(mix$true_mean, mix$est_mean),
    mixture_variance_cor = safe_cor(mix$true_variance, mix$est_variance),
    mixture_variance_rmse = sqrt(mean((mix$true_variance - mix$est_variance)^2)),
    pretraining_converged = isTRUE(pre$pretraining$converged),
    pretraining_completed_iter = pre$pretraining$n_completed,
    refinement_converged = isTRUE(fit$joint_refinement$converged),
    refinement_completed_iter = fit$joint_refinement$n_completed,
    elapsed_sec = elapsed
  )
  mix$variant <- variant
  mix$rep <- rep_id
  list(metric = metric, mixture = mix)
}

metrics_file <- file.path(out_dir, "feature_grid_metrics_checkpoint.csv")
mixture_file <- file.path(out_dir, "feature_grid_mixture_checkpoint.csv")
metrics_all <- if (file.exists(metrics_file)) read.csv(metrics_file) else data.frame()
mixture_all <- if (file.exists(mixture_file)) read.csv(mixture_file) else data.frame()
completed <- if (nrow(metrics_all)) paste(metrics_all$variant, metrics_all$rep, sep = "|") else character(0)

cat("Writing outputs to:", out_dir, "\n")
cat("Variants:", paste(variants, collapse = ", "), "\n")
cat("Rep values:", paste(rep_values, collapse = ", "), "\n")

for (variant in variants) {
  for (rep_id in rep_values) {
    key <- paste(variant, rep_id, sep = "|")
    if (key %in% completed) {
      cat("Skipping completed", key, "\n")
      next
    }
    cat(sprintf("\nRunning variant=%s rep=%d\n", variant, rep_id))
    out <- fit_one(variant, rep_id)
    metrics_all <- rbind(metrics_all, out$metric)
    mixture_all <- rbind(mixture_all, out$mixture)
    write.csv(metrics_all, metrics_file, row.names = FALSE)
    write.csv(mixture_all, mixture_file, row.names = FALSE)
    print(out$metric[, c("variant", "rep", "prob_extreme", "lambda_cor", "factor_mean_abs_cor", "mixture_variance_cor", "mixture_variance_rmse", "elapsed_sec")], row.names = FALSE)
  }
}

metric_cols <- c(
  "alpha_sd",
  "mean_abs_lambda",
  "max_abs_lambda",
  "n_lambda_abs_gt_0_5",
  "n_items_no_abs_gt_0_5",
  "n_items_multi_abs_gt_0_5",
  "prob_extreme",
  "row_extreme",
  "item_extreme",
  "alpha_cor",
  "lambda_cor",
  "factor_mean_abs_cor",
  "eta_cor",
  "prob_cor",
  "mixture_weight_cor",
  "mixture_mean_cor",
  "mixture_variance_cor",
  "mixture_variance_rmse",
  "pretraining_converged",
  "refinement_converged",
  "elapsed_sec"
)
summary <- do.call(rbind, lapply(split(metrics_all, metrics_all$variant), function(d) {
  vals <- vapply(metric_cols, function(nm) mean(d[[nm]], na.rm = TRUE), numeric(1))
  vals[!is.finite(vals)] <- NA_real_
  data.frame(
    variant = unique(d$variant)[1L],
    as.data.frame(as.list(vals), check.names = FALSE),
    n_reps = length(unique(d$rep)),
    check.names = FALSE
  )
}))
rownames(summary) <- NULL
write.csv(summary, file.path(out_dir, "feature_grid_summary.csv"), row.names = FALSE)
write.csv(metrics_all, file.path(out_dir, "feature_grid_metrics.csv"), row.names = FALSE)
write.csv(mixture_all, file.path(out_dir, "feature_grid_mixture.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  summary$variant <- factor(summary$variant, levels = variants)
  p1 <- ggplot(summary, aes(x = variant, y = mixture_variance_cor, fill = variant)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "grey70") +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(title = "Variance recovery after adding IFEval-like observation features", x = NULL, y = "mean variance correlation") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
  ggsave(file.path(out_dir, "feature_grid_variance_correlation.png"), p1, width = 8, height = 4.8, dpi = 200)

  p2 <- ggplot(summary, aes(x = prob_extreme, y = mixture_variance_cor, label = variant)) +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_point(size = 2.5, color = "#2B6CB0") +
    geom_text(nudge_y = 0.06, size = 3) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(title = "Variance recovery vs probability saturation", x = "Pr(prob < .05 or > .95)", y = "mean variance correlation") +
    theme_minimal(base_size = 11)
  ggsave(file.path(out_dir, "feature_grid_variance_vs_saturation.png"), p2, width = 8, height = 5, dpi = 200)
}

cat("\nFinal summary:\n")
print(summary[order(summary$variant), ], row.names = FALSE)
cat("\nOutputs written to:", normalizePath(out_dir), "\n")
