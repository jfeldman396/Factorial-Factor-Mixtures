#!/usr/bin/env Rscript

# Short H=3, G=3 diagnostic: balanced moderate primary blocks, but with
# substantially denser signed cross-loadings.

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

parse_int_csv <- function(x) as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))

sim_n <- get_env("SIM_N", 400L, as.integer)
sim_p <- get_env("SIM_P", 500L, as.integer)
H <- get_env("SIM_H", 3L, as.integer)
G <- get_env("SIM_G", 3L, as.integer)
rep_values <- get_env("REP_VALUES", 1:2, parse_int_csv)
seed_base <- get_env("SEED", 20260820L, as.integer)
workers <- get_env("WORKERS", 1L, as.integer)
lambda_l1_penalty <- get_env("REFINEMENT_LAMBDA_L1_PENALTY", 10, as.numeric)
pretrain_iter <- get_env("PRETRAIN_AUG_ITER", 12L, as.integer)
refine_iter <- get_env("REFINE_ITER", 10L, as.integer)
cross_prob <- get_env("CROSS_PROB", 0.25, as.numeric)
cross_min <- get_env("CROSS_MIN", 0.20, as.numeric)
cross_max <- get_env("CROSS_MAX", 0.60, as.numeric)

out_dir <- get_env(
  "OUT_DIR",
  file.path(bundle_root, "results", sprintf("signed_crossloading_short_test_H%d_G%d", H, G)),
  as.character
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit_path <- file.path(
  bundle_root,
  "results",
  "reproduced_openeval_ifeval_lambda_sparsity_tuning",
  "ifeval_H3_G3_lambda_l1_10_fit.rds"
)
ifeval_fit <- readRDS(fit_path)
ifeval_alpha <- as.numeric(ifeval_fit$alpha_hat)

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
  H_local <- ncol(L_true)
  perms <- all_permutations(seq_len(H_local))
  best <- NULL
  for (rr in seq_len(nrow(perms))) {
    perm <- perms[rr, ]
    for (mask in seq_len(2^H_local) - 1L) {
      signs <- ifelse(as.integer(intToBits(mask))[seq_len(H_local)] == 1L, -1, 1)
      L_aligned <- sweep(L_est[, perm, drop = FALSE], 2L, signs, "*")
      score <- sum((L_true - L_aligned)^2)
      if (is.null(best) || score < best$score) {
        best <- list(perm = perm, signs = signs, Lambda = L_aligned, score = score)
      }
    }
  }
  best
}

align_matrix <- function(M, align) sweep(M[, align$perm, drop = FALSE], 2L, align$signs, "*")

balanced_block_sizes <- function(p, H_local) {
  block_sizes <- rep(floor(p / H_local), H_local)
  remainder <- p - sum(block_sizes)
  if (remainder > 0L) block_sizes[seq_len(remainder)] <- block_sizes[seq_len(remainder)] + 1L
  block_sizes
}

make_signed_cross_loading <- function(p, H_local) {
  block_sizes <- balanced_block_sizes(p, H_local)
  block_id <- rep(seq_len(H_local), times = block_sizes)
  Lambda <- matrix(0, p, H_local)
  for (j in seq_len(p)) {
    h <- block_id[j]
    Lambda[j, h] <- runif(1, 0.75, 1.25)
    for (k in setdiff(seq_len(H_local), h)) {
      if (runif(1) < cross_prob) {
        Lambda[j, k] <- sample(c(-1, 1), 1L) * runif(1, cross_min, cross_max)
      }
    }
  }
  rownames(Lambda) <- sprintf("item_%03d", seq_len(p))
  colnames(Lambda) <- paste0("F", seq_len(H_local))
  list(Lambda = Lambda, block_id = block_id)
}

make_original_g3_params <- function(H_local) {
  replicate(
    H_local,
    list(pi = c(0.30, 0.40, 0.30), mu = c(-1, 0, 1), sd = c(0.25, 0.50, 0.75)),
    simplify = FALSE
  )
}

simulate_factors <- function(n, H_local) {
  params <- make_original_g3_params(H_local)
  Fmat <- matrix(NA_real_, n, H_local)
  component <- matrix(NA_integer_, n, H_local)
  standardized_params <- vector("list", H_local)
  for (h in seq_len(H_local)) {
    draw <- sample_standardized_mixture(n, params[[h]])
    Fmat[, h] <- draw$x
    component[, h] <- draw$component
    standardized_params[[h]] <- draw$parameters
  }
  colnames(Fmat) <- paste0("F", seq_len(H_local))
  colnames(component) <- paste0("F", seq_len(H_local))
  list(F = Fmat, component = component, mixture_params = standardized_params)
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

write_lambda_heatmap <- function(Lambda, file, title, zlim = NULL) {
  if (is.null(zlim)) zlim <- c(-max(abs(Lambda)), max(abs(Lambda)))
  block_sizes <- balanced_block_sizes(nrow(Lambda), ncol(Lambda))
  block_ends <- cumsum(block_sizes)
  tick_items <- sort(unique(c(1L, block_ends[-length(block_ends)] + 1L, block_ends)))
  png(file, width = 1200, height = 1800, res = 180)
  par(mar = c(5, 5, 4, 6), xaxs = "i", yaxs = "i")
  pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(256)
  image(
    x = seq_len(ncol(Lambda)),
    y = seq_len(nrow(Lambda)),
    z = t(Lambda[nrow(Lambda):1, , drop = FALSE]),
    col = pal,
    zlim = zlim,
    axes = FALSE,
    xlab = "factor",
    ylab = "item",
    main = title
  )
  axis(1, at = seq_len(ncol(Lambda)), labels = colnames(Lambda))
  axis(2, at = nrow(Lambda) - tick_items + 1, labels = tick_items, las = 1)
  if (length(block_ends) > 1L) {
    abline(h = nrow(Lambda) - block_ends[-length(block_ends)] + 1.5, col = "#444444", lwd = 1)
  }
  box()
  usr <- par("usr")
  legend_x <- usr[2] + 0.28
  legend_y <- seq(usr[3], usr[4], length.out = 256)
  rect(legend_x, legend_y[-256], legend_x + 0.18, legend_y[-1], col = pal, border = NA, xpd = NA)
  axis(4, at = seq(usr[3], usr[4], length.out = 5), labels = sprintf("%.2f", seq(zlim[1], zlim[2], length.out = 5)), las = 1, xpd = NA)
  mtext("lambda", side = 4, line = 4)
  invisible(dev.off())
}

fit_one <- function(rep_id) {
  seed <- seed_base + rep_id
  set.seed(seed)
  loading <- make_signed_cross_loading(sim_p, H)
  Lambda <- loading$Lambda
  alpha <- sample(ifeval_alpha, sim_p, replace = sim_p > length(ifeval_alpha))
  names(alpha) <- rownames(Lambda)
  factors <- simulate_factors(sim_n, H)
  eta <- sweep(factors$F %*% t(Lambda), 2L, alpha, "+")
  prob <- pnorm(eta)
  X <- 1L * (eta + matrix(rnorm(sim_n * sim_p), sim_n, sim_p) > 0)
  colnames(X) <- rownames(Lambda)
  rownames(X) <- sprintf("sim_%03d", seq_len(sim_n))
  rownames(factors$F) <- rownames(X)

  t0 <- proc.time()[["elapsed"]]
  pre <- fit_binary_probit_pretraining_intercept(
    X = X,
    H = H,
    G_fixed = G,
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

  align <- best_signed_permutation(Lambda, fit$Lambda_hat[colnames(X), , drop = FALSE])
  Lambda_aligned <- align$Lambda
  colnames(Lambda_aligned) <- colnames(Lambda)
  F_aligned <- align_matrix(fit$F_hat, align)
  colnames(F_aligned) <- colnames(factors$F)
  mix <- extract_mixture_recovery(factors$mixture_params, fit$mixture_fits, align)
  eta_est <- sweep(F_aligned %*% t(Lambda_aligned), 2L, fit$alpha_hat[colnames(X)], "+")

  item_active_counts <- rowSums(abs(Lambda) >= 0.5)
  metric <- data.frame(
    variant = sprintf("balanced_moderate_signed_cross_dense_H%d_G%d", H, G),
    rep = rep_id,
    n = sim_n,
    p = sim_p,
    H = H,
    G = G,
    cross_prob = cross_prob,
    cross_min = cross_min,
    cross_max = cross_max,
    lambda_l1_penalty = lambda_l1_penalty,
    alpha_sd = sd(alpha),
    mean_abs_lambda = mean(abs(Lambda)),
    max_abs_lambda = max(abs(Lambda)),
    n_lambda_abs_gt_0_5 = sum(abs(Lambda) >= 0.5),
    n_nonprimary_abs_gt_0_1 = sum(abs(Lambda[cbind(rep(seq_len(sim_p), each = H), rep(seq_len(H), times = sim_p))]) > 0.1) - sim_p,
    n_items_no_abs_gt_0_5 = sum(item_active_counts == 0L),
    n_items_multi_abs_gt_0_5 = sum(item_active_counts >= 2L),
    prob_extreme = mean(prob < 0.05 | prob > 0.95),
    row_extreme = mean(rowMeans(X) <= 0.05 | rowMeans(X) >= 0.95),
    item_extreme = mean(colMeans(X) <= 0.05 | colMeans(X) >= 0.95),
    alpha_cor = safe_cor(alpha, fit$alpha_hat[colnames(X)]),
    lambda_cor = safe_cor(as.numeric(Lambda), as.numeric(Lambda_aligned)),
    factor_mean_abs_cor = mean(abs(vapply(seq_len(H), function(h) safe_cor(factors$F[, h], F_aligned[, h]), numeric(1)))),
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
  mix$variant <- metric$variant
  mix$rep <- rep_id

  if (rep_id == rep_values[1L]) {
    write.csv(data.frame(item = rownames(Lambda), block = loading$block_id, Lambda, check.names = FALSE), file.path(out_dir, "rep1_true_lambda.csv"), row.names = FALSE)
    write.csv(data.frame(item = rownames(Lambda_aligned), Lambda_aligned, check.names = FALSE), file.path(out_dir, "rep1_estimated_lambda_aligned.csv"), row.names = FALSE)
    lim <- c(-1.25, 1.25)
    write_lambda_heatmap(Lambda, file.path(out_dir, "rep1_true_lambda_heatmap.png"), "Signed dense-cross DGP Lambda, rep 1", lim)
    write_lambda_heatmap(Lambda_aligned, file.path(out_dir, "rep1_estimated_lambda_aligned_heatmap.png"), "Aligned estimated Lambda, rep 1", lim)
  }

  list(metric = metric, mixture = mix)
}

cat("Writing outputs to:", out_dir, "\n")
cat("Settings: n=", sim_n, ", p=", sim_p, ", H=", H, ", G=", G,
    ", cross_prob=", cross_prob, ", cross_range=[", cross_min, ", ", cross_max,
    "], lambda_l1=", lambda_l1_penalty, "\n", sep = "")

metrics_all <- data.frame()
mixture_all <- data.frame()
for (rep_id in rep_values) {
  cat("\nRunning signed-cross rep=", rep_id, "\n", sep = "")
  out <- fit_one(rep_id)
  metrics_all <- rbind(metrics_all, out$metric)
  mixture_all <- rbind(mixture_all, out$mixture)
  write.csv(metrics_all, file.path(out_dir, "signed_cross_metrics_checkpoint.csv"), row.names = FALSE)
  write.csv(mixture_all, file.path(out_dir, "signed_cross_mixture_checkpoint.csv"), row.names = FALSE)
  print(out$metric[, c("rep", "prob_extreme", "n_items_multi_abs_gt_0_5", "lambda_cor", "factor_mean_abs_cor", "mixture_variance_cor", "mixture_variance_rmse", "elapsed_sec")], row.names = FALSE)
}

metric_cols <- c(
  "alpha_sd",
  "mean_abs_lambda",
  "max_abs_lambda",
  "n_lambda_abs_gt_0_5",
  "n_nonprimary_abs_gt_0_1",
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
summary <- as.data.frame(as.list(vapply(metric_cols, function(nm) mean(metrics_all[[nm]], na.rm = TRUE), numeric(1))), check.names = FALSE)
summary$variant <- sprintf("balanced_moderate_signed_cross_dense_H%d_G%d", H, G)
summary$n_reps <- length(unique(metrics_all$rep))
summary <- summary[, c("variant", setdiff(names(summary), "variant"))]
write.csv(metrics_all, file.path(out_dir, "signed_cross_metrics.csv"), row.names = FALSE)
write.csv(mixture_all, file.path(out_dir, "signed_cross_mixture.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "signed_cross_summary.csv"), row.names = FALSE)

cat("\nFinal summary:\n")
print(summary, row.names = FALSE)
cat("\nOutputs written to:", normalizePath(out_dir), "\n")
