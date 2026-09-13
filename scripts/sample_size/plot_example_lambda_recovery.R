#!/usr/bin/env Rscript

# Example loading-matrix recovery diagnostics.
#
# This script recreates one fixed-IFEval-like simulation cell, fits Product MAP
# and Viroli-Laplace Gibbs, aligns the fitted loading matrices to the true
# columns/signs, and plots true/estimated/error heatmaps.  It is intentionally
# separate from the long-running simulation driver so that visual diagnostics do
# not change the main result files.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
file_arg <- gsub("~\\+~", " ", file_arg)
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "../.."))

source(file.path(repo_root, "R", "binary_probit_pretraining.R"))
source(file.path(repo_root, "R", "binary_probit_refinement.R"))
source(file.path(repo_root, "R", "probit_ifa_em_svd_pretraining.R"))
source(file.path(repo_root, "R", "riemannian_rotation.R"))
source(file.path(repo_root, "R", "viroli_probit_independent_gibbs.R"))
source(file.path(repo_root, "R", "sample_size_dgp.R"))

get_env <- function(name, default, FUN = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  FUN(value)
}

parse_int_csv <- function(x) as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1L]]))

normalize_G_counts_local <- function(G, H) {
  G <- as.integer(G)
  if (length(G) == 1L) G <- rep(G, H)
  if (length(G) != H) stop("G must have length 1 or H.", call. = FALSE)
  G
}

format_G_config <- function(G) paste(as.integer(G), collapse = "-")

make_viroli_smoke_mixture_params <- function(H, G, sep) {
  G <- normalize_G_counts_local(G, H)
  lapply(seq_len(H), function(h) {
    Gh <- G[h]
    if (Gh == 1L) {
      list(pi = 1, mu = 0, sd = 1)
    } else if (Gh == 2L) {
      list(pi = c(0.50, 0.50), mu = sep * c(-1, 1), sd = c(0.55, 0.85))
    } else if (Gh == 3L) {
      list(pi = c(0.30, 0.40, 0.30), mu = sep * c(-1.35, 0, 1.35), sd = c(0.45, 0.65, 0.45))
    } else {
      stop("This diagnostic supports G_h in {1,2,3}.", call. = FALSE)
    }
  })
}

make_item_intercepts <- function(p, H, block_id, seed) {
  make_sample_size_item_intercepts(
    p = p,
    H = H,
    block_id = block_id,
    seed = seed,
    mode = "ifeval_like",
    intercept_sd = 0.45,
    intercept_block_span = 1.6,
    intercept_clip = 1.75
  )
}

simulate_fixed_ifeval_cell <- function(n, p, H, G, rep_value, seed_base = 20260731L,
                                       row_idx = NULL, dgp_p_max = 2000L,
                                       sep = 1,
                                       primary_loading_range = c(2, 3),
                                       cross_loading_range = c(2, 3),
                                       cross_loading_prob = 0.05) {
  G <- normalize_G_counts_local(G, H)
  if (is.null(row_idx)) {
    n_grid <- c(100L, 200L, 400L)
    p_grid <- c(500L, 1000L, 1500L, 2000L)
    np <- expand.grid(n = n_grid, p = p_grid, KEEP.OUT.ATTRS = FALSE)
    np <- np[order(np$p, np$n), , drop = FALSE]
    row_idx <- match(paste(n, p), paste(np$n, np$p))
    if (is.na(row_idx)) stop("Could not infer row_idx for n=", n, ", p=", p, call. = FALSE)
  }

  loading_design <- "balanced_moderate_dense_signed_cross"
  loading_designs <- normalize_sample_size_loading_design(loading_design)
  loading_index <- 1L
  block_index <- match("ifeval_min30", c("balanced", "ifeval_like", "moderate_ifeval_like", "ifeval_min30"))
  strength_index <- match("strong", c("default", "weak", "strong"))
  cross_prob_offset <- as.integer(round(10000 * cross_loading_prob))

  seed <- seed_base + 100000L * H + 50000L * sum(G * seq_along(G)) +
    10000L * rep_value + 1000L * row_idx
  loading_parameter_seed <- seed_base + 100000L * H +
    5000L * loading_index + 1000L * block_index +
    100L * strength_index + cross_prob_offset + as.integer(round(1000 * sep))
  mixture_parameter_seed <- seed_base + 100000L * H +
    50000L * sum(G * seq_along(G)) + as.integer(round(1000 * sep))

  mixture_params <- make_viroli_smoke_mixture_params(H, G, sep)

  set.seed(loading_parameter_seed)
  master_loading <- make_sample_size_loadings(
    design = loading_designs,
    p = dgp_p_max,
    H = H,
    block_size_mode = "ifeval_min30",
    loading_sign_mode = "block",
    primary_loading_range = primary_loading_range,
    cross_loading_range = cross_loading_range,
    cross_loading_prob = cross_loading_prob,
    cross_sign_mode = "random"
  )
  loading_out <- subset_sample_size_loading_output(
    master_loading,
    p = p,
    H = H,
    block_size_mode = "ifeval_min30"
  )
  master_alpha <- make_item_intercepts(
    p = dgp_p_max,
    H = H,
    block_id = master_loading$block_id,
    seed = loading_parameter_seed + 3571L
  )
  alpha <- master_alpha[loading_out$master_rows]

  set.seed(seed)
  F <- matrix(NA_real_, n, H)
  component <- matrix(NA_integer_, n, H)
  standardized_params <- vector("list", H)
  for (h in seq_len(H)) {
    draw_h <- sample_standardized_mixture(n, mixture_params[[h]])
    F[, h] <- draw_h$x
    component[, h] <- draw_h$component
    standardized_params[[h]] <- draw_h$parameters
  }
  Z_latent <- sweep(F %*% t(loading_out$Lambda), 2L, alpha, "+") +
    matrix(rnorm(n * p), n, p)
  X_binary <- 1L * (Z_latent > 0)

  list(
    X_binary = X_binary,
    F = F,
    Lambda = loading_out$Lambda,
    alpha = alpha,
    component = component,
    mixture_params = standardized_params,
    block_id = loading_out$block_id,
    block_sizes = loading_out$block_sizes,
    seed = seed,
    row_idx = row_idx,
    loading_parameter_seed = loading_parameter_seed,
    mixture_parameter_seed = mixture_parameter_seed
  )
}

align_loadings_to_truth <- function(Lambda_true, Lambda_est_raw, F_true = NULL, F_est = NULL) {
  H_true <- ncol(Lambda_true)
  H_est <- ncol(Lambda_est_raw)
  if (H_true > H_est) stop("Cannot align fewer estimated columns than true columns.")
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' is required for loading-column alignment.", call. = FALSE)
  }
  dist_mat <- matrix(NA_real_, H_true, H_est)
  sign_mat <- matrix(1, H_true, H_est)
  for (h_true in seq_len(H_true)) {
    for (h_est in seq_len(H_est)) {
      pos_dist <- sum((Lambda_true[, h_true] - Lambda_est_raw[, h_est])^2)
      neg_dist <- sum((Lambda_true[, h_true] + Lambda_est_raw[, h_est])^2)
      if (pos_dist <= neg_dist) {
        dist_mat[h_true, h_est] <- pos_dist
        sign_mat[h_true, h_est] <- 1
      } else {
        dist_mat[h_true, h_est] <- neg_dist
        sign_mat[h_true, h_est] <- -1
      }
    }
  }
  est_index <- as.integer(clue::solve_LSAP(dist_mat))
  signs <- sign_mat[cbind(seq_len(H_true), est_index)]
  matched_abs_cor <- rep(NA_real_, H_true)
  if (!is.null(F_true) && !is.null(F_est)) {
    C <- suppressWarnings(cor(F_true, F_est))
    matched_abs_cor <- abs(C[cbind(seq_len(H_true), est_index)])
  }
  list(
    est_index = est_index,
    signs = signs,
    matched_abs_cor = matched_abs_cor,
    mean_abs_cor = mean(matched_abs_cor, na.rm = TRUE)
  )
}

align_lambda_to_truth <- function(Lambda_true, Lambda_est_raw, align) {
  out <- matrix(0, nrow(Lambda_true), ncol(Lambda_true))
  out[] <- sweep(Lambda_est_raw[, align$est_index, drop = FALSE], 2L, align$signs, "*")
  dimnames(out) <- dimnames(Lambda_true)
  out
}

align_factor_scores_to_truth <- function(F_est_raw, align) {
  sweep(F_est_raw[, align$est_index, drop = FALSE], 2L, align$signs, "*")
}

fit_var <- function(fit) {
  if (!is.null(fit$var)) fit$var else fit$sd^2
}

sort_mixture_fit_by_mean <- function(fit) {
  ord <- order(fit$mu)
  fit$pi <- fit$pi[ord]
  fit$mu <- fit$mu[ord]
  if (!is.null(fit$var)) fit$var <- fit$var[ord]
  if (!is.null(fit$sd)) fit$sd <- fit$sd[ord]
  fit
}

align_mixture_fits_to_truth <- function(mixture_fits, align) {
  H <- length(align$est_index)
  out <- vector("list", H)
  for (h in seq_len(H)) {
    fit_h <- mixture_fits[[align$est_index[h]]]
    out[[h]] <- list(
      pi = fit_h$pi,
      mu = align$signs[h] * fit_h$mu,
      var = fit_var(fit_h),
      G = length(fit_h$pi)
    )
  }
  out
}

make_mixture_parameter_table <- function(true_mixture_params, mixture_fits, method) {
  rows <- vector("list", sum(vapply(true_mixture_params, function(z) length(z$pi), integer(1))))
  idx <- 0L
  for (h in seq_along(true_mixture_params)) {
    true_h <- sort_mixture_fit_by_mean(list(
      pi = true_mixture_params[[h]]$pi,
      mu = true_mixture_params[[h]]$mu,
      var = true_mixture_params[[h]]$sd^2
    ))
    fit_h <- sort_mixture_fit_by_mean(mixture_fits[[h]])
    for (g in seq_along(true_h$pi)) {
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        method = method,
        factor = h,
        component = g,
        true_weight = true_h$pi[g],
        est_weight = fit_h$pi[g],
        true_mu = true_h$mu[g],
        est_mu = fit_h$mu[g],
        true_var = true_h$var[g],
        est_var = fit_h$var[g],
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out$weight_error <- out$est_weight - out$true_weight
  out$mu_error <- out$est_mu - out$true_mu
  out$var_error <- out$est_var - out$true_var
  out
}

summarize_mixture_parameter_table <- function(table) {
  do.call(rbind, lapply(split(table, table$method), function(x) {
    data.frame(
      method = x$method[1L],
      marginal_mu_rmse = sqrt(mean(x$mu_error^2)),
      marginal_var_rmse = sqrt(mean(x$var_error^2)),
      marginal_weight_rmse = sqrt(mean(x$weight_error^2)),
      stringsAsFactors = FALSE
    )
  }))
}

plot_factor_score_panel <- function(F_true, product_F, viroli_F, out_file, title) {
  H <- ncol(F_true)
  methods <- list("Product MAP" = product_F, "Viroli Laplace" = viroli_F)
  xylim <- range(c(F_true, product_F, viroli_F), finite = TRUE)
  pad <- 0.06 * diff(xylim)
  if (!is.finite(pad) || pad <= 0) pad <- 0.25
  xylim <- xylim + c(-pad, pad)

  png(out_file, width = 2600, height = 1050, res = 180)
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    invisible(dev.off())
  }, add = TRUE)
  layout(matrix(seq_len(length(methods) * H), nrow = length(methods), byrow = TRUE))
  par(mar = c(4.1, 4.3, 3.0, 1.0), oma = c(0, 0, 3.0, 0))
  for (method_name in names(methods)) {
    F_hat <- methods[[method_name]]
    for (h in seq_len(H)) {
      r <- suppressWarnings(cor(F_true[, h], F_hat[, h]))
      plot(
        F_true[, h],
        F_hat[, h],
        pch = 16,
        cex = 0.55,
        col = rgb(33, 102, 172, 110, maxColorValue = 255),
        xlim = xylim,
        ylim = xylim,
        xlab = paste0("true F", h),
        ylab = paste0(method_name, " F", h),
        main = sprintf("%s F%d, cor=%.3f", method_name, h, r)
      )
      abline(0, 1, lty = 2, col = "gray35")
      abline(h = 0, v = 0, lty = 3, col = "gray75")
      box()
    }
  }
  mtext(title, outer = TRUE, cex = 1.25, font = 2)
}

plot_factor_mixture_panel <- function(F_true, product_F, viroli_F,
                                      true_mixture_params, product_mixtures,
                                      viroli_mixtures, out_file, title) {
  H <- ncol(F_true)
  score_sets <- list(True = F_true, "Product MAP" = product_F, "Viroli Laplace" = viroli_F)
  mixture_sets <- list(
    True = lapply(true_mixture_params, function(z) list(pi = z$pi, mu = z$mu, var = z$sd^2)),
    "Product MAP" = product_mixtures,
    "Viroli Laplace" = viroli_mixtures
  )
  component_cols <- c("#2f70bf", "#cc333f", "#2e8b57", "#9a6fb0", "#d9a21b")

  png(out_file, width = 2800, height = 1600, res = 180)
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    invisible(dev.off())
  }, add = TRUE)
  layout(matrix(seq_len(length(score_sets) * H), nrow = length(score_sets), byrow = TRUE))
  par(mar = c(4.0, 4.2, 3.0, 1.0), oma = c(0, 0, 3.2, 0))
  for (method_name in names(score_sets)) {
    F_mat <- score_sets[[method_name]]
    mix_list <- mixture_sets[[method_name]]
    for (h in seq_len(H)) {
      fit_h <- sort_mixture_fit_by_mean(mix_list[[h]])
      var_h <- fit_var(fit_h)
      x_all <- c(
        F_true[, h],
        product_F[, h],
        viroli_F[, h],
        fit_h$mu - 4 * sqrt(var_h),
        fit_h$mu + 4 * sqrt(var_h)
      )
      xlim <- range(x_all, finite = TRUE)
      grid <- seq(xlim[1L], xlim[2L], length.out = 300L)
      comp_density <- vapply(seq_along(fit_h$pi), function(g) {
        fit_h$pi[g] * dnorm(grid, mean = fit_h$mu[g], sd = sqrt(var_h[g]))
      }, numeric(length(grid)))
      total_density <- rowSums(comp_density)
      hist_out <- hist(F_mat[, h], breaks = 14L, plot = FALSE)
      ylim <- c(0, 1.08 * max(hist_out$density, total_density, finite = TRUE))
      plot(
        hist_out,
        freq = FALSE,
        col = "#eef2f7",
        border = "#cad2dc",
        xlim = xlim,
        ylim = ylim,
        xlab = paste0("F", h, " score"),
        ylab = if (h == 1L) "density" else "",
        main = sprintf("%s F%d", method_name, h)
      )
      lines(grid, total_density, lwd = 2.0, col = "black")
      for (g in seq_along(fit_h$pi)) {
        col_g <- component_cols[((g - 1L) %% length(component_cols)) + 1L]
        lines(grid, comp_density[, g], lwd = 1.7, col = col_g)
        segments(fit_h$mu[g], 0, fit_h$mu[g], 0.10 * ylim[2L], col = col_g, lwd = 2)
      }
      rug(F_mat[, h], col = rgb(0, 0, 0, 55, maxColorValue = 255), ticksize = 0.03)
      box()
    }
  }
  mtext(title, outer = TRUE, cex = 1.25, font = 2)
}

plot_lambda_panel <- function(true_Lambda, product_Lambda, viroli_Lambda, block_id,
                              out_file, title) {
  row_order <- unlist(lapply(seq_len(ncol(true_Lambda)), function(h) {
    rows <- which(block_id == h)
    rows[order(abs(true_Lambda[rows, h]), decreasing = TRUE)]
  }), use.names = FALSE)
  true_s <- true_Lambda[row_order, , drop = FALSE]
  product_s <- product_Lambda[row_order, , drop = FALSE]
  viroli_s <- viroli_Lambda[row_order, , drop = FALSE]
  product_err <- product_s - true_s
  viroli_err <- viroli_s - true_s
  err_diff <- abs(product_err) - abs(viroli_err)
  block_sizes <- tabulate(block_id[row_order], nbins = ncol(true_Lambda))
  block_cuts <- cumsum(block_sizes)

  zlim <- c(-1, 1) * max(abs(c(true_s, product_s, viroli_s)), na.rm = TRUE)
  elim <- c(-1, 1) * max(abs(c(product_err, viroli_err)), na.rm = TRUE)
  dlim <- c(-1, 1) * max(abs(err_diff), na.rm = TRUE)
  pal <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)
  diff_pal <- colorRampPalette(c("#2f70bf", "#f7f7f7", "#c9343a"))(101)

  draw_heatmap <- function(M, main, zlim_use, palette_use, ylab = "") {
    image(
      x = seq_len(ncol(M)),
      y = seq_len(nrow(M)),
      z = t(M[nrow(M):1, , drop = FALSE]),
      col = palette_use,
      zlim = zlim_use,
      axes = FALSE,
      xlab = "factor",
      ylab = ylab,
      main = main,
      useRaster = TRUE
    )
    axis(1, at = seq_len(ncol(M)), labels = paste0("F", seq_len(ncol(M))), las = 1)
    axis(2, at = pretty(seq_len(nrow(M)), 4), labels = rev(pretty(seq_len(nrow(M)), 4)), las = 1)
    box()
    y_cuts <- nrow(M) - block_cuts + 0.5
    abline(h = y_cuts[-length(y_cuts)], lty = 3, col = "gray35")
  }

  png(out_file, width = 2600, height = 1700, res = 180)
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    invisible(dev.off())
  }, add = TRUE)
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7, 7), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.28))
  par(mar = c(4.3, 4.6, 3.0, 1.0), oma = c(0, 0, 3.2, 0))
  draw_heatmap(true_s, "true Lambda", zlim, pal, "ordered items")
  draw_heatmap(product_s, "Product MAP aligned", zlim, pal)
  draw_heatmap(viroli_s, "Viroli Laplace aligned", zlim, pal)
  draw_heatmap(product_err, "Product MAP error", elim, pal, "ordered items")
  draw_heatmap(viroli_err, "Viroli Laplace error", elim, pal)
  draw_heatmap(err_diff, "|Product error| - |Viroli error|", dlim, diff_pal)
  par(mar = c(0.5, 2, 0.5, 2))
  plot.new()
  legend(
    "center",
    legend = c(
      "Rows sorted by true primary block, then primary magnitude",
      "Red/blue heatmaps: signed loading value or signed error",
      "Bottom-right: red means Product MAP has larger absolute error; blue means Viroli has larger absolute error"
    ),
    bty = "n",
    cex = 0.95
  )
  mtext(title, outer = TRUE, cex = 1.25, font = 2)
}

evaluate_lambda_fit <- function(sim, fit_name, F_hat, Lambda_hat, alpha_hat, seconds) {
  align <- align_loadings_to_truth(sim$Lambda, Lambda_hat, sim$F, F_hat)
  Lambda_aligned <- align_lambda_to_truth(sim$Lambda, Lambda_hat, align)
  F_aligned <- sweep(F_hat[, align$est_index, drop = FALSE], 2L, align$signs, "*")
  eta_true <- sweep(sim$F %*% t(sim$Lambda), 2L, sim$alpha, "+")
  eta_hat <- sweep(F_aligned %*% t(Lambda_aligned), 2L, alpha_hat, "+")
  data.frame(
    method = fit_name,
    lambda_rmse = sqrt(mean((sim$Lambda - Lambda_aligned)^2)),
    lambda_corr = suppressWarnings(cor(as.vector(sim$Lambda), as.vector(Lambda_aligned))),
    mean_factor_abs_cor = align$mean_abs_cor,
    factor_score_rmse = sqrt(mean((scale(sim$F) - scale(F_aligned))^2)),
    probability_rmse = sqrt(mean((pnorm(eta_true) - pnorm(eta_hat))^2)),
    seconds = seconds,
    stringsAsFactors = FALSE
  )
}

n_value <- get_env("N_VALUE", 100L, as.integer)
p_value <- get_env("P_VALUE", 500L, as.integer)
H_value <- get_env("H_TRUE", 5L, as.integer)
G_value <- get_env("G_TRUE", 3L, as.integer)
rep_value <- get_env("REP_VALUE", 1L, as.integer)
lasso_penalty <- get_env("LASSO_PENALTY", 10, as.numeric)
sep_value <- get_env("SEPARATIONS", 2, as.numeric)
primary_loading_range <- get_env("PRIMARY_LOADING_RANGE", c(2, 3), parse_int_csv)
cross_loading_range <- get_env("CROSS_LOADING_RANGE", c(2, 3), parse_int_csv)
cross_loading_prob <- get_env("CROSS_LOADING_PROB", 0.05, as.numeric)
product_workers <- get_env("PRODUCT_WORKERS", 18L, as.integer)
viroli_workers <- get_env("VIROLI_WORKERS", 4L, as.integer)
viroli_iter <- get_env("VIROLI_ITER", 2000L, as.integer)
viroli_burn <- get_env("VIROLI_BURN", 1000L, as.integer)
out_dir <- get_env(
  "OUT_DIR",
  file.path(repo_root, "results", "selected_plots", "sample_size",
            "fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10", "lambda_recovery_examples"),
  as.character
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

G_config <- rep(G_value, H_value)
G_label <- format_G_config(G_config)
tag <- sprintf("lambda_recovery_n%d_p%d_H%d_G%d_rep%d_sep%s_penalty%s",
               n_value, p_value, H_value, G_value, rep_value,
               gsub("[.]", "p", format(sep_value, trim = TRUE)),
               gsub("[.]", "p", format(lasso_penalty, trim = TRUE)))

cat(sprintf("Simulating fixed cell n=%d, p=%d, H=%d, G=%s, rep=%d, sep=%s\n",
            n_value, p_value, H_value, G_label, rep_value, format(sep_value, trim = TRUE)))
sim <- simulate_fixed_ifeval_cell(
  n = n_value,
  p = p_value,
  H = H_value,
  G = G_config,
  rep_value = rep_value,
  sep = sep_value,
  primary_loading_range = primary_loading_range,
  cross_loading_range = cross_loading_range,
  cross_loading_prob = cross_loading_prob
)

cat("Fitting Product MAP...\n")
t0 <- proc.time()[["elapsed"]]
product <- fit_binary_probit_em_svd_pretrain_then_refine(
  X = sim$X_binary,
  H = H_value,
  G_fixed = G_config,
  em_max_iter = 50L,
  em_tol_loglik = 1e-5,
  em_tol_L = 1e-4,
  em_init_method = "both",
  em_init_z = "expectation",
  rotation_random_starts = 1L,
  rotation_loading_l1_penalty = lasso_penalty,
  rotation_max_outer = 20L,
  n_mix_starts = 3L,
  mixture_max_iter = 100L,
  mixture_update = "map",
  mu_prior_mean = 0,
  mu_prior_kappa = 0.05,
  var_prior_shape = 3,
  var_prior_scale = 2,
  weight_prior_alpha = 1,
  refine_mu_prior_mean = 0,
  refine_mu_prior_kappa = 0.05,
  refine_var_prior_shape = 3,
  refine_var_prior_scale = 2,
  refine_weight_prior_alpha = 1,
  min_mixture_var = 0.05,
  grid_size = 21L,
  rotation_optimizer = "riemannian",
  riemannian_rotation_steps = 10L,
  rotation_objective_tolerance = 1e-4,
  rotation_min_outer = 2L,
  require_mixture_convergence_for_rotation_stop = TRUE,
  loading_penalty = lasso_penalty,
  n_refine_iter = 50L,
  factor_update = "marginal",
  lambda_l1_penalty = lasso_penalty,
  lasso_backend = "glmnet",
  glmnet_standardize = FALSE,
  objective_tolerance = 1e-3,
  min_refine_iter = 3L,
  return_best_refinement_iteration = TRUE,
  refinement_selection_objective = "posterior_objective",
  require_mixture_convergence_for_stop = TRUE,
  mixture_refit = "em",
  enforce_monotone_refinement = TRUE,
  monotone_tolerance = 1e-8,
  parallel = TRUE,
  workers = product_workers,
  seed = sim$seed + 11L,
  verbose = FALSE
)
product_seconds <- proc.time()[["elapsed"]] - t0

cat("Fitting Viroli-Laplace Gibbs...\n")
viroli <- fit_viroli_probit_independent_gibbs(
  X = sim$X_binary,
  H = H_value,
  G = G_config,
  n_iter = viroli_iter,
  burn = viroli_burn,
  thin = 1L,
  tau_lambda = 1.5,
  tau_intercept = 5,
  lambda_l1_penalty = lasso_penalty,
  alpha_dirichlet = 1,
  min_scale = 1e-4,
  normalize_each_draw = TRUE,
  parallel = TRUE,
  workers = viroli_workers,
  compute_parameter_ess = FALSE,
  seed = sim$seed + 37L,
  verbose = FALSE
)

product_align <- align_loadings_to_truth(
  sim$Lambda,
  product$refine_fit$Lambda_hat,
  sim$F,
  product$refine_fit$F_hat
)
viroli_align <- align_loadings_to_truth(sim$Lambda, viroli$Lambda_hat, sim$F, viroli$F_hat)
product_lambda <- align_lambda_to_truth(sim$Lambda, product$refine_fit$Lambda_hat, product_align)
viroli_lambda <- align_lambda_to_truth(sim$Lambda, viroli$Lambda_hat, viroli_align)
product_F <- align_factor_scores_to_truth(product$refine_fit$F_hat, product_align)
viroli_F <- align_factor_scores_to_truth(viroli$F_hat, viroli_align)
product_mixture_fits <- align_mixture_fits_to_truth(product$refine_fit$mixture_fits, product_align)
viroli_mixture_fits <- align_mixture_fits_to_truth(viroli$mixture_fits, viroli_align)

mixture_parameter_table <- rbind(
  make_mixture_parameter_table(sim$mixture_params, product_mixture_fits, "Product MAP"),
  make_mixture_parameter_table(sim$mixture_params, viroli_mixture_fits, "Viroli Laplace")
)
mixture_parameter_summary <- summarize_mixture_parameter_table(mixture_parameter_table)

metrics <- rbind(
  evaluate_lambda_fit(
    sim,
    "Product MAP",
    product$refine_fit$F_hat,
    product$refine_fit$Lambda_hat,
    product$refine_fit$alpha_hat,
    product_seconds
  ),
  evaluate_lambda_fit(
    sim,
    "Viroli Laplace",
    viroli$F_hat,
    viroli$Lambda_hat,
    viroli$alpha_hat,
    viroli$seconds
  )
)
metrics <- merge(metrics, mixture_parameter_summary, by = "method", all.x = TRUE, sort = FALSE)
metrics$n <- n_value
metrics$p <- p_value
metrics$H <- H_value
metrics$G <- G_label
metrics$rep <- rep_value
metrics$lasso_penalty <- lasso_penalty
metrics$product_pretraining_converged <- isTRUE(product$pretrain_fit$pretraining_converged)
metrics$product_rotation_converged <- isTRUE(product$pretrain_fit$rotation_converged)
metrics$product_refinement_converged <- isTRUE(product$refine_fit$joint_refinement$converged)

write.csv(metrics, file.path(out_dir, paste0(tag, "_metrics.csv")), row.names = FALSE)
write.csv(mixture_parameter_table, file.path(out_dir, paste0(tag, "_mixture_parameters.csv")), row.names = FALSE)
write.csv(sim$Lambda, file.path(out_dir, paste0(tag, "_true_lambda.csv")), row.names = FALSE)
write.csv(product_lambda, file.path(out_dir, paste0(tag, "_product_map_lambda_aligned.csv")), row.names = FALSE)
write.csv(viroli_lambda, file.path(out_dir, paste0(tag, "_viroli_laplace_lambda_aligned.csv")), row.names = FALSE)
write.csv(product_F, file.path(out_dir, paste0(tag, "_product_map_factor_scores_aligned.csv")), row.names = FALSE)
write.csv(viroli_F, file.path(out_dir, paste0(tag, "_viroli_laplace_factor_scores_aligned.csv")), row.names = FALSE)

plot_lambda_panel(
  true_Lambda = sim$Lambda,
  product_Lambda = product_lambda,
  viroli_Lambda = viroli_lambda,
  block_id = sim$block_id,
  out_file = file.path(out_dir, paste0(tag, "_heatmap_panel.png")),
  title = sprintf(
    "Example Lambda recovery: n=%d, p=%d, H=%d, G=%s, rep=%d, sep=%s, penalty=%s",
    n_value, p_value, H_value, G_label, rep_value, format(sep_value, trim = TRUE),
    format(lasso_penalty, trim = TRUE)
  )
)
plot_factor_score_panel(
  F_true = sim$F,
  product_F = product_F,
  viroli_F = viroli_F,
  out_file = file.path(out_dir, paste0(tag, "_factor_scores.png")),
  title = sprintf(
    "Example factor-score recovery: n=%d, p=%d, H=%d, G=%s, rep=%d, sep=%s, penalty=%s",
    n_value, p_value, H_value, G_label, rep_value, format(sep_value, trim = TRUE),
    format(lasso_penalty, trim = TRUE)
  )
)
plot_factor_mixture_panel(
  F_true = sim$F,
  product_F = product_F,
  viroli_F = viroli_F,
  true_mixture_params = sim$mixture_params,
  product_mixtures = product_mixture_fits,
  viroli_mixtures = viroli_mixture_fits,
  out_file = file.path(out_dir, paste0(tag, "_factor_mixtures.png")),
  title = sprintf(
    "Example fitted factor marginals: n=%d, p=%d, H=%d, G=%s, rep=%d, sep=%s, penalty=%s",
    n_value, p_value, H_value, G_label, rep_value, format(sep_value, trim = TRUE),
    format(lasso_penalty, trim = TRUE)
  )
)

cat("Wrote:\n")
cat(file.path(out_dir, paste0(tag, "_heatmap_panel.png")), "\n")
cat(file.path(out_dir, paste0(tag, "_factor_scores.png")), "\n")
cat(file.path(out_dir, paste0(tag, "_factor_mixtures.png")), "\n")
cat(file.path(out_dir, paste0(tag, "_metrics.csv")), "\n")
cat(file.path(out_dir, paste0(tag, "_mixture_parameters.csv")), "\n")
print(metrics)
