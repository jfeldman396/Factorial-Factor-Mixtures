#!/usr/bin/env Rscript

# Compare product-mixture rotation starts to direct FastICA rotation.
#
# This diagnostic isolates the stage-one signal and rotation stages.  For each
# simulated data set we construct two low-rank signals:
#
#   1. estimated_Z: the EM-SVD probit low-rank signal estimated from X.
#   2. oracle_Z: the latent Gaussian Z from the simulator, centered and projected
#      to rank H before taking the same SVD scores/loadings.
#
# From each signal we create three pretraining objects:
#
#   1. fastica_only: use FastICA scores directly, then fit marginal mixtures.
#   2. mixture_fastica_start: initialize at FastICA, then optimize the profiled
#      marginal-mixture rotation criterion.
#   3. mixture_identity_start: initialize at the unrotated SVD scores, then
#      optimize the same marginal-mixture rotation criterion.
#
# All pretraining objects are passed through the identical MAP refinement
# routine.  The output has pretraining and refined metrics so we can see whether
# the bottleneck is binary signal estimation, rotation, or the refinement stage.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

source(file.path(repo_root, "R", "binary_probit_pretraining.R"))
source(file.path(repo_root, "R", "binary_probit_refinement.R"))
source(file.path(repo_root, "R", "canonical_factor_normalization.R"))
source(file.path(repo_root, "R", "sample_size_dgp.R"))
source(file.path(repo_root, "R", "probit_ifa_em_svd_pretraining.R"))
source(file.path(repo_root, "R", "riemannian_rotation.R"))

get_env <- function(name, default, coercer = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  coercer(value)
}

split_csv <- function(x) {
  out <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

parse_ints <- function(x) as.integer(split_csv(x))

parse_numeric_pair <- function(x) {
  out <- as.numeric(split_csv(x))
  if (length(out) != 2L || any(!is.finite(out))) {
    stop("Expected a comma-separated numeric pair.", call. = FALSE)
  }
  sort(out)
}

parse_optional_numeric_pair <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parse_numeric_pair(value)
}

parse_optional_numeric <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  out <- as.numeric(value)
  if (length(out) != 1L || !is.finite(out)) {
    stop(name, " must be a single finite number.", call. = FALSE)
  }
  out
}

parse_numeric_vector <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  out <- if (nzchar(value)) as.numeric(split_csv(value)) else as.numeric(default)
  if (!length(out) || any(!is.finite(out))) {
    stop(name, " must be a comma-separated vector of finite numbers.", call. = FALSE)
  }
  out
}

parse_penalty_by_n <- function(x) {
  x <- trimws(as.character(x))
  if (!nzchar(x)) return(data.frame(n = integer(), penalty = numeric()))
  pieces <- split_csv(x)
  out <- lapply(pieces, function(piece) {
    kv <- trimws(strsplit(piece, ":", fixed = TRUE)[[1L]])
    if (length(kv) != 2L || !nzchar(kv[1L]) || !nzchar(kv[2L])) {
      stop("PENALTY_BY_N entries must look like n:penalty, e.g. 100:5,400:8.", call. = FALSE)
    }
    data.frame(n = as.integer(kv[1L]), penalty = as.numeric(kv[2L]))
  })
  out <- do.call(rbind, out)
  if (any(!is.finite(out$n)) || any(!is.finite(out$penalty))) {
    stop("PENALTY_BY_N contains a non-finite n or penalty.", call. = FALSE)
  }
  out
}

safe_cor <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok]))
}

rbind_fill <- function(x) {
  x <- x[!vapply(x, is.null, logical(1L))]
  if (!length(x)) return(data.frame())
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) d[[nm]] <- NA
    d[, all_names, drop = FALSE]
  })
  do.call(rbind, x)
}

normalize_G_counts <- function(G, H) normalize_G_fixed(G, H)

g_vector <- function(type, H) {
  type <- match.arg(type, c("all2", "all3"))
  if (type == "all2") rep(2L, H) else rep(3L, H)
}

g_label <- function(G) paste(as.integer(G), collapse = "-")

dataset_key <- function(block_size_mode, loading_design, n, p, H, G_config, rep) {
  paste(block_size_mode, loading_design, n, p, H, G_config, rep, sep = "|")
}

dataset_key_from_results <- function(results) {
  block_size_mode <- if ("block_size_mode" %in% names(results)) {
    results$block_size_mode
  } else {
    rep("ifeval_min30", nrow(results))
  }
  loading_design <- if ("loading_design" %in% names(results)) {
    results$loading_design
  } else {
    rep("balanced_moderate_dense_signed_cross", nrow(results))
  }
  dataset_key(
    block_size_mode,
    loading_design,
    results$n,
    results$p,
    results$H,
    results$G_config,
    results$rep
  )
}

make_viroli_smoke_mixture_params_local <- function(H, G, sep) {
  G <- normalize_G_counts(G, H)
  lapply(seq_len(H), function(h) {
    Gh <- G[h]
    if (Gh == 1L) {
      list(pi = 1, mu = 0, sd = 1)
    } else if (Gh == 2L) {
      if (length(viroli_smoke_g2_pi) != 2L || length(viroli_smoke_g2_sd) != 2L) {
        stop("The G=2 mixture weights and standard deviations must have length 2.", call. = FALSE)
      }
      list(
        pi = viroli_smoke_g2_pi / sum(viroli_smoke_g2_pi),
        mu = sep * viroli_smoke_g2_mu_multiplier * c(-1, 1),
        sd = viroli_smoke_g2_sd
      )
    } else if (Gh == 3L) {
      if (length(viroli_smoke_g3_pi) != 3L || length(viroli_smoke_g3_sd) != 3L) {
        stop("The G=3 mixture weights and standard deviations must have length 3.", call. = FALSE)
      }
      list(
        pi = viroli_smoke_g3_pi / sum(viroli_smoke_g3_pi),
        mu = sep * viroli_smoke_g3_mu_multiplier * c(-1, 0, 1),
        sd = viroli_smoke_g3_sd
      )
    } else {
      stop("This diagnostic supports G_h in {1, 2, 3}.", call. = FALSE)
    }
  })
}

estimated_z_eigengap_diagnostics <- function(lowrank, true_H, max_candidate_rank) {
  L <- as.matrix(lowrank$L)
  max_available <- min(nrow(L), ncol(L)) - 1L
  max_candidate_rank <- min(as.integer(max_candidate_rank), max_available)
  if (max_candidate_rank < 1L) stop("No candidate ranks are available for eigengap selection.")

  singular_values <- svd(L, nu = 0L, nv = 0L)$d
  required <- max_candidate_rank + 1L
  singular_values <- c(singular_values, rep(0, max(0L, required - length(singular_values))))
  singular_values <- singular_values[seq_len(required)]
  eigenvalues <- singular_values^2 / (nrow(L) * ncol(L))
  raw_gaps <- eigenvalues[seq_len(max_candidate_rank)] - eigenvalues[seq.int(2L, required)]
  relative_gaps <- raw_gaps / pmax(eigenvalues[seq_len(max_candidate_rank)], .Machine$double.eps)
  ratio_gaps <- singular_values[seq_len(max_candidate_rank)] /
    pmax(singular_values[seq.int(2L, required)], .Machine$double.eps)

  raw_rank <- which.max(raw_gaps)
  relative_rank <- which.max(relative_gaps)
  ratio_rank <- which.max(ratio_gaps)
  true_gap <- if (true_H <= max_candidate_rank) raw_gaps[true_H] else NA_real_
  true_relative_gap <- if (true_H <= max_candidate_rank) relative_gaps[true_H] else NA_real_
  true_ratio <- if (true_H <= max_candidate_rank) ratio_gaps[true_H] else NA_real_

  summary <- data.frame(
    eigengap_fit_rank = if (!is.null(lowrank$requested_rank)) lowrank$requested_rank else required,
    eigengap_max_candidate_rank = max_candidate_rank,
    estimated_z_h_hat_largest_eigengap = raw_rank,
    estimated_z_h_hat_largest_relative_eigengap = relative_rank,
    estimated_z_h_hat_largest_sv_ratio = ratio_rank,
    estimated_z_largest_eigengap = raw_gaps[raw_rank],
    estimated_z_largest_relative_eigengap = relative_gaps[relative_rank],
    estimated_z_largest_sv_ratio = ratio_gaps[ratio_rank],
    estimated_z_eigengap_at_true_H = true_gap,
    estimated_z_relative_eigengap_at_true_H = true_relative_gap,
    estimated_z_sv_ratio_at_true_H = true_ratio,
    estimated_z_largest_eigengap_selects_true_H = as.integer(raw_rank == true_H),
    estimated_z_largest_relative_eigengap_selects_true_H = as.integer(relative_rank == true_H),
    estimated_z_largest_sv_ratio_selects_true_H = as.integer(ratio_rank == true_H),
    stringsAsFactors = FALSE
  )
  spectrum <- as.list(setNames(eigenvalues, paste0("estimated_z_eigenvalue_", seq_along(eigenvalues))))
  gaps <- as.list(setNames(raw_gaps, paste0("estimated_z_eigengap_", seq_along(raw_gaps))))
  relative <- as.list(setNames(relative_gaps, paste0("estimated_z_relative_eigengap_", seq_along(relative_gaps))))
  ratios <- as.list(setNames(ratio_gaps, paste0("estimated_z_sv_ratio_", seq_along(ratio_gaps))))
  cbind(summary, as.data.frame(c(spectrum, gaps, relative, ratios), check.names = FALSE))
}

simulate_rotation_dgp <- function(
    n,
    p,
    H,
    G,
    sep,
    rep,
    p_master,
    base_seed,
    loading_seed,
    mixture_seed,
    loading_design,
    block_size_mode,
    intercept_mode) {
  G <- normalize_G_counts(G, H)

  set.seed(mixture_seed + 1000003L * H + 1009L * sum(G))
  mixture_params <- make_viroli_smoke_mixture_params_local(H, G, sep = sep)

  loading_design <- normalize_sample_size_loading_design(loading_design)
  design_id <- match(loading_design, sort(unique(normalize_sample_size_loading_design(names(sample_size_loading_design_aliases())))))
  if (is.na(design_id)) design_id <- 1L
  block_id <- match(block_size_mode, c("balanced", "ifeval_like", "moderate_ifeval_like", "ifeval_min30"))
  if (is.na(block_id)) stop("Unsupported block_size_mode: ", block_size_mode, call. = FALSE)

  set.seed(loading_seed + 1000003L * H + 1009L * design_id + 101L * block_id)
  master_loading <- make_sample_size_loadings(
    design = loading_design,
    p = p_master,
    H = H,
    block_size_mode = block_size_mode,
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
    block_size_mode = block_size_mode
  )

  alpha_master <- make_sample_size_item_intercepts(
    p = p_master,
    H = H,
    block_id = master_loading$block_id,
    seed = loading_seed + 3571L + H,
    mode = intercept_mode,
    intercept_sd = intercept_sd,
    intercept_block_span = intercept_block_span,
    intercept_clip = intercept_clip
  )
  alpha <- alpha_master[loading_out$master_rows]

  data_seed <- base_seed + 104729L * rep + 1009L * n + 17L * p +
    7919L * H + 101L * sum(G) + 10007L * design_id + 100003L * block_id
  set.seed(data_seed)
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
    Z_latent = Z_latent,
    F = F,
    Lambda = loading_out$Lambda,
    alpha = alpha,
    component = component,
    mixture_params = standardized_params,
    block_id = loading_out$block_id,
    block_sizes = loading_out$block_sizes,
    loading_design = loading_design,
    loading_design_label = sample_size_loading_design_label(loading_design),
    block_size_mode = block_size_mode,
    data_seed = data_seed,
    loading_seed = loading_seed,
    mixture_seed = mixture_seed,
    prevalence = colMeans(X_binary)
  )
}

align_by_loadings <- function(Lambda_true, Lambda_est, F_true, F_est) {
  H <- ncol(Lambda_true)
  if (ncol(Lambda_est) < H) stop("Estimated loading matrix has too few columns.")
  dist_mat <- matrix(NA_real_, H, ncol(Lambda_est))
  sign_mat <- matrix(1, H, ncol(Lambda_est))
  for (h in seq_len(H)) {
    for (k in seq_len(ncol(Lambda_est))) {
      pos <- sum((Lambda_true[, h] - Lambda_est[, k])^2)
      neg <- sum((Lambda_true[, h] + Lambda_est[, k])^2)
      if (pos <= neg) {
        dist_mat[h, k] <- pos
        sign_mat[h, k] <- 1
      } else {
        dist_mat[h, k] <- neg
        sign_mat[h, k] <- -1
      }
    }
  }
  est_index <- as.integer(clue::solve_LSAP(dist_mat))
  signs <- sign_mat[cbind(seq_len(H), est_index)]
  C <- suppressWarnings(cor(F_true, F_est))
  matched_abs_cor <- abs(C[cbind(seq_len(H), est_index)])
  list(
    est_index = est_index,
    signs = signs,
    matched_abs_cor = matched_abs_cor,
    mean_abs_cor = mean(matched_abs_cor, na.rm = TRUE)
  )
}

align_lambda <- function(Lambda_true, Lambda_est, align) {
  H <- ncol(Lambda_true)
  out <- matrix(0, nrow(Lambda_true), H)
  out[, seq_len(H)] <- sweep(Lambda_est[, align$est_index, drop = FALSE], 2L, align$signs, "*")
  out
}

align_mixture_fits <- function(mixture_fits, align, H) {
  out <- vector("list", H)
  for (h in seq_len(H)) {
    fit_h <- mixture_fits[[align$est_index[h]]]
    var_h <- if (!is.null(fit_h$var)) fit_h$var else fit_h$sd^2
    mu_h <- align$signs[h] * fit_h$mu
    ord <- order(mu_h)
    out[[h]] <- list(pi = fit_h$pi[ord], mu = mu_h[ord], var = var_h[ord])
  }
  out
}

marginal_mixture_summary <- function(true_mixture_params, mixture_fits, G) {
  H <- length(true_mixture_params)
  G <- normalize_G_counts(G, H)
  rows <- vector("list", sum(G))
  idx <- 0L
  for (h in seq_len(H)) {
    true_h <- true_mixture_params[[h]]
    fit_h <- mixture_fits[[h]]
    ord <- order(fit_h$mu)
    fit_pi <- fit_h$pi[ord]
    fit_mu <- fit_h$mu[ord]
    fit_var <- fit_h$var[ord]
    for (g in seq_len(G[h])) {
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        factor = h,
        component = g,
        true_weight = true_h$pi[g],
        est_weight = fit_pi[g],
        true_mu = true_h$mu[g],
        est_mu = fit_mu[g],
        true_var = true_h$sd[g]^2,
        est_var = fit_var[g],
        stringsAsFactors = FALSE
      )
    }
  }
  tab <- do.call(rbind, rows)
  data.frame(
    marginal_mu_rmse = sqrt(mean((tab$est_mu - tab$true_mu)^2)),
    marginal_var_rmse = sqrt(mean((tab$est_var - tab$true_var)^2)),
    marginal_weight_rmse = sqrt(mean((tab$est_weight - tab$true_weight)^2)),
    stringsAsFactors = FALSE
  )
}

standardize_scores <- function(F) {
  out <- scale(F)
  out[!is.finite(out)] <- 0
  out
}

evaluate_estimate <- function(stage, rotation_method, sim, H, G, F_hat, Lambda_hat, alpha_hat, mixture_fits, seconds) {
  align <- align_by_loadings(sim$Lambda, Lambda_hat, sim$F, F_hat)
  F_aligned <- sweep(F_hat[, align$est_index, drop = FALSE], 2L, align$signs, "*")
  Lambda_aligned <- align_lambda(sim$Lambda, Lambda_hat, align)
  mixture_aligned <- align_mixture_fits(mixture_fits, align, H)
  mix <- marginal_mixture_summary(sim$mixture_params, mixture_aligned, G)
  eta_true <- sweep(sim$F %*% t(sim$Lambda), 2L, sim$alpha, "+")
  eta_hat <- sweep(F_aligned %*% t(Lambda_aligned), 2L, alpha_hat, "+")

  cbind(
    data.frame(
      stage = stage,
      rotation_method = rotation_method,
      mean_factor_abs_cor = align$mean_abs_cor,
      factor_score_rmse = sqrt(mean((standardize_scores(sim$F) - standardize_scores(F_aligned))^2)),
      factor_score_raw_rmse = sqrt(mean((sim$F - F_aligned)^2)),
      lambda_rmse = sqrt(mean((sim$Lambda - Lambda_aligned)^2)),
      lambda_corr = safe_cor(as.numeric(sim$Lambda), as.numeric(Lambda_aligned)),
      alpha_rmse = sqrt(mean((sim$alpha - alpha_hat)^2)),
      alpha_corr = safe_cor(sim$alpha, alpha_hat),
      linear_predictor_rmse = sqrt(mean((eta_true - eta_hat)^2)),
      probability_rmse = sqrt(mean((pnorm(eta_true) - pnorm(eta_hat))^2)),
      seconds = seconds,
      stringsAsFactors = FALSE
    ),
    mix
  )
}

orthonormal_left_basis <- function(M, H) {
  dec <- svd(as.matrix(M), nu = H, nv = 0)
  dec$u[, seq_len(H), drop = FALSE]
}

subspace_distance <- function(A, B, H) {
  U <- orthonormal_left_basis(A, H)
  V <- orthonormal_left_basis(B, H)
  sqrt(sum((tcrossprod(U) - tcrossprod(V))^2) / (2 * H))
}

make_oracle_z_lowrank <- function(X, Z_latent, H) {
  t0 <- proc.time()[["elapsed"]]
  projection <- rank_H_centered_projection(Z_latent, H)
  elapsed <- proc.time()[["elapsed"]] - t0
  list(
    alpha = projection$alpha,
    L = projection$L,
    probit_loglik = probit_ifa_loglik_lowrank(X, projection$alpha, projection$L),
    history = data.frame(
      iteration = 1L,
      probit_loglik = probit_ifa_loglik_lowrank(X, projection$alpha, projection$L),
      relative_loglik_change = NA_real_,
      relative_L_change = NA_real_,
      projection_update = "oracle_Z_rank_H_projection",
      stochastic_svd_draws = NA_integer_,
      stochastic_svd_average = NA_character_,
      e_step_seconds = NA_real_,
      projection_seconds = elapsed,
      objective_seconds = NA_real_,
      iteration_seconds = elapsed
    ),
    converged = TRUE,
    n_completed = 1L,
    start_name = "oracle_Z",
    selected_start = "oracle_Z",
    init_method = "oracle_Z",
    init_z = "oracle_Z",
    projection_update = "oracle_Z_rank_H_projection",
    signal_seconds = elapsed
  )
}

build_pretrain_object <- function(
    X,
    H,
    G,
    lowrank,
    spectral,
    rotation_method,
    R,
    F_hat,
    Lambda_hat,
    mixture_fits,
    rotation_history,
    rotation_seconds,
    rotation_converged,
    rotation_completed_outer,
    selection_score) {
  responsibilities <- lapply(seq_len(H), function(h) {
    mixture_responsibilities(F_hat[, h], mixture_fits[[h]])
  })
  class_map <- sapply(responsibilities, max.col, ties.method = "first")
  if (H == 1L) class_map <- matrix(class_map, ncol = 1L)
  colnames(class_map) <- paste0("factor_", seq_len(H))
  fitted <- sweep(F_hat %*% t(Lambda_hat), 2L, lowrank$alpha, "+")
  residual <- probit_ifa_truncated_mean(X, lowrank$alpha, lowrank$L) - fitted

  list(
    model = paste0("rotation_diagnostic_", rotation_method),
    X = X,
    H = H,
    alpha_hat = lowrank$alpha,
    L_hat = lowrank$L,
    S = spectral$S,
    svd_fit = spectral,
    lowrank_fit = lowrank,
    R = R,
    F_hat = F_hat,
    Lambda_hat = Lambda_hat,
    Lambda_ls = Lambda_hat,
    mixture_fits = mixture_fits,
    G_hat = vapply(mixture_fits, function(z) length(z$pi), integer(1L)),
    G_fixed = normalize_G_counts(G, H),
    mixture_update = mixture_update,
    rotation_optimizer = rotation_method,
    responsibilities = responsibilities,
    class_map = class_map,
    profile_id = binary_profile_id(class_map),
    fitted = fitted,
    residual = residual,
    history = data.frame(
      stage = c("lowrank_em_svd", rotation_method),
      iteration = c(lowrank$n_completed, rotation_completed_outer),
      probit_loglik = c(lowrank$probit_loglik, lowrank$probit_loglik),
      mixture_loglik = c(NA_real_, mixture_loglik_total(F_hat, mixture_fits)),
      objective = c(lowrank$probit_loglik, lowrank$probit_loglik + selection_score)
    ),
    em_history = lowrank$history,
    pretraining_converged = isTRUE(lowrank$converged),
    pretraining_completed_iter = lowrank$n_completed,
    selected_pretraining_iteration = lowrank$n_completed,
    selected_em_start = lowrank$selected_start,
    rotation_completed_outer = rotation_completed_outer,
    rotation_converged = isTRUE(rotation_converged),
    rotation_history = rotation_history,
    rotation_seconds = rotation_seconds,
    z_update = "none_em_svd",
    fix_psi_identity = TRUE,
    estimate_intercept = TRUE
  )
}

fit_mixture_criterion_pretrain <- function(
    X,
    H,
    G,
    lowrank,
    spectral,
    rotation_method,
    seed,
    n_random_starts_use = 0L,
    n_ica_starts_use = 0L,
    include_identity_start_use = TRUE) {
  t0 <- proc.time()[["elapsed"]]
  rotation <- rotate_em_svd_scores_with_mixtures(
    S = spectral$S,
    G_fixed = G,
    loading_basis = spectral$B,
    rotation_loading_l1_penalty = rotation_loading_l1_penalty,
    n_random_starts = n_random_starts_use,
    n_ica_starts = n_ica_starts_use,
    include_identity_start = include_identity_start_use,
    ica_functions = fastica_functions,
    ica_max_iter = fastica_max_iter,
    ica_tol = fastica_tol,
    max_outer = rotation_max_outer,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha,
    rotation_objective_tolerance = rotation_objective_tolerance,
    rotation_min_outer = rotation_min_outer,
    require_mixture_convergence_for_rotation_stop = TRUE,
    rotation_optimizer = "riemannian",
    riemannian_rotation_steps = riemannian_rotation_steps,
    seed = seed,
    parallel = parallel_enabled,
    workers = workers,
    verbose = FALSE
  )
  Lambda_hat <- spectral$B %*% rotation$R
  if (pretrain_loading_penalty > 0) Lambda_hat <- soft_threshold(Lambda_hat, pretrain_loading_penalty)
  colnames(Lambda_hat) <- paste0("factor_", seq_len(H))
  build_pretrain_object(
    X = X,
    H = H,
    G = G,
    lowrank = lowrank,
    spectral = spectral,
    rotation_method = rotation_method,
    R = rotation$R,
    F_hat = rotation$F_hat,
    Lambda_hat = Lambda_hat,
    mixture_fits = rotation$fits,
    rotation_history = rotation$rotation_history,
    rotation_seconds = proc.time()[["elapsed"]] - t0,
    rotation_converged = rotation$rotation_converged,
    rotation_completed_outer = rotation$rotation_completed_outer,
    selection_score = rotation$selection_score
  )
}

fit_fastica_direct_pretrain <- function(
    X,
    H,
    G,
    lowrank,
    spectral,
    seed,
    rotation_method = "fastica_only") {
  t0 <- proc.time()[["elapsed"]]
  S <- spectral$S
  starts <- fastica_rotation_starts(
    S = S,
    n_starts = fastica_starts,
    ica_functions = fastica_functions,
    max_iter = fastica_max_iter,
    tol = fastica_tol,
    seed = seed,
    verbose = FALSE
  )
  if (!length(starts)) starts <- list(fastica_unavailable_identity = diag(H))

  fits_by_start <- lapply(seq_along(starts), function(s) {
    R <- starts[[s]]
    F_hat <- S %*% R
    fits <- fit_column_mixtures_fixed_G(
      F = F_hat,
      G_fixed = G,
      n_starts = n_mix_starts,
      max_iter = mixture_max_iter,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      parallel = parallel_enabled,
      workers = workers
    )
    score <- rotation_selection_score(
      F = F_hat,
      fits = fits,
      R = R,
      loading_basis = spectral$B,
      loading_l1_penalty = rotation_loading_l1_penalty
    )
    list(R = R, F_hat = F_hat, fits = fits, score = score, start_name = names(starts)[s])
  })

  scores <- vapply(fits_by_start, `[[`, numeric(1L), "score")
  best <- if (identical(fastica_selection, "best_mixture")) which.max(scores) else 1L
  best_fit <- fits_by_start[[best]]
  Lambda_hat <- spectral$B %*% best_fit$R
  if (pretrain_loading_penalty > 0) Lambda_hat <- soft_threshold(Lambda_hat, pretrain_loading_penalty)
  colnames(Lambda_hat) <- paste0("factor_", seq_len(H))
  rotation_history <- data.frame(
    outer_iteration = 0L,
    mixture_loglik = mixture_loglik_total(best_fit$F_hat, best_fit$fits),
    rotation_selection_score = best_fit$score,
    loading_l1 = rotation_loading_l1(spectral$B, best_fit$R),
    start_name = best_fit$start_name,
    fastica_selection = fastica_selection,
    all_mixtures_converged = all(vapply(best_fit$fits, function(z) isTRUE(z$converged), logical(1L))),
    stringsAsFactors = FALSE
  )

  build_pretrain_object(
    X = X,
    H = H,
    G = G,
    lowrank = lowrank,
    spectral = spectral,
    rotation_method = rotation_method,
    R = best_fit$R,
    F_hat = best_fit$F_hat,
    Lambda_hat = Lambda_hat,
    mixture_fits = best_fit$fits,
    rotation_history = rotation_history,
    rotation_seconds = proc.time()[["elapsed"]] - t0,
    rotation_converged = TRUE,
    rotation_completed_outer = 0L,
    selection_score = best_fit$score
  )
}

refine_pretrain <- function(X, pretrain, seed) {
  t0 <- proc.time()[["elapsed"]]
  set.seed(seed)
  fit <- fit_binary_probit_refinement(
    X = X,
    pretrain_fit = pretrain,
    n_refine_iter = refine_iter,
    maxit_per_subject = maxit_per_subject,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mixture_refit = "em",
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha,
    factor_update = "marginal",
    min_mixture_var = min_mixture_var,
    lambda_l1_penalty = lambda_l1_penalty,
    lasso_backend = lasso_backend,
    glmnet_standardize = FALSE,
    objective_tolerance = refinement_objective_tolerance,
    min_refine_iter = refine_min_iter,
    stopping_objective = "posterior_objective",
    return_best_refinement_iteration = TRUE,
    refinement_selection_objective = "posterior_objective",
    require_mixture_convergence_for_stop = TRUE,
    enforce_monotone_refinement = TRUE,
    normalize_factor_scale = TRUE,
    normalize_factor_location = TRUE,
    parallel = parallel_enabled,
    workers = workers,
    verbose = FALSE
  )
  fit <- canonical_normalize_refined_fit(
    fit,
    min_scale = min_mixture_var,
    sign_rule = "none",
    order_rule = "none"
  )
  fit$fit_seconds <- proc.time()[["elapsed"]] - t0
  fit
}

metadata_pair <- function(x) {
  if (is.null(x)) return(c(NA_real_, NA_real_))
  as.numeric(x)
}

fit_one_dataset <- function(block_size_mode, loading_design, n, p, H, G, rep) {
  cat(sprintf(
    "[%s] block=%s loading=%s n=%d p=%d H=%d G=%s rep=%d\n",
    format(Sys.time(), "%H:%M:%S"),
    block_size_mode,
    sample_size_loading_design_label(loading_design),
    n, p, H, g_label(G), rep
  ))
  sim <- simulate_rotation_dgp(
    n = n,
    p = p,
    H = H,
    G = G,
    sep = separation,
    rep = rep,
    p_master = max(p_values),
    base_seed = base_seed,
    loading_seed = loading_seed,
    mixture_seed = mixture_seed,
    loading_design = loading_design,
    block_size_mode = block_size_mode,
    intercept_mode = intercept_mode
  )

  true_signal <- sim$F %*% t(sim$Lambda)

  lowrank_t0 <- proc.time()[["elapsed"]]
  lowrank <- fit_lowrank_probit_em_svd(
    X = sim$X_binary,
    H = H,
    max_iter = em_max_iter,
    tol_loglik = em_tol_loglik,
    tol_L = em_tol_L,
    tol_subspace = em_tol_subspace,
    init_method = em_init_method,
    init_z = em_init_z,
    projection_update = "expectation",
    seed = base_seed + 10000019L * rep + 101L * H + p,
    parallel = parallel_enabled,
    workers = workers,
    verbose = FALSE
  )
  lowrank_seconds <- proc.time()[["elapsed"]] - lowrank_t0
  lowrank$signal_seconds <- lowrank_seconds

  eigengap_t0 <- proc.time()[["elapsed"]]
  eigengap_diagnostics <- data.frame(
    eigengap_fit_rank = NA_integer_,
    eigengap_max_candidate_rank = NA_integer_,
    estimated_z_h_hat_largest_eigengap = NA_integer_,
    estimated_z_h_hat_largest_relative_eigengap = NA_integer_,
    estimated_z_h_hat_largest_sv_ratio = NA_integer_,
    estimated_z_largest_eigengap = NA_real_,
    estimated_z_largest_relative_eigengap = NA_real_,
    estimated_z_largest_sv_ratio = NA_real_,
    estimated_z_eigengap_at_true_H = NA_real_,
    estimated_z_relative_eigengap_at_true_H = NA_real_,
    estimated_z_sv_ratio_at_true_H = NA_real_,
    estimated_z_largest_eigengap_selects_true_H = NA_integer_,
    estimated_z_largest_relative_eigengap_selects_true_H = NA_integer_,
    estimated_z_largest_sv_ratio_selects_true_H = NA_integer_,
    stringsAsFactors = FALSE
  )
  if (isTRUE(run_eigengap_diagnostic)) {
    rank_fit <- min(eigengap_fit_rank, min(n, p) - 1L)
    candidate_max <- min(eigengap_max_candidate_rank, rank_fit - 1L)
    eigengap_fit <- fit_lowrank_probit_em_svd(
      X = sim$X_binary,
      H = rank_fit,
      max_iter = em_max_iter,
      tol_loglik = em_tol_loglik,
      tol_L = em_tol_L,
      tol_subspace = em_tol_subspace,
      init_method = em_init_method,
      init_z = em_init_z,
      projection_update = "expectation",
      seed = base_seed + 10000019L * rep + 101L * H + p + 900001L,
      parallel = parallel_enabled,
      workers = workers,
      verbose = FALSE
    )
    eigengap_fit$requested_rank <- rank_fit
    eigengap_diagnostics <- estimated_z_eigengap_diagnostics(
      lowrank = eigengap_fit,
      true_H = H,
      max_candidate_rank = candidate_max
    )
  }
  eigengap_seconds <- proc.time()[["elapsed"]] - eigengap_t0

  oracle_lowrank <- make_oracle_z_lowrank(sim$X_binary, sim$Z_latent, H)

  eta_lowrank <- sweep(lowrank$L, 2L, lowrank$alpha, "+")
  eta_oracle_lowrank <- sweep(oracle_lowrank$L, 2L, oracle_lowrank$alpha, "+")
  estimated_oracle_gap <- list(
    estimated_vs_oracle_signal_rmse = sqrt(mean((lowrank$L - oracle_lowrank$L)^2)),
    estimated_vs_oracle_alpha_rmse = sqrt(mean((lowrank$alpha - oracle_lowrank$alpha)^2)),
    estimated_vs_oracle_linear_predictor_rmse = sqrt(mean((eta_lowrank - eta_oracle_lowrank)^2)),
    estimated_vs_oracle_probability_rmse = sqrt(mean((pnorm(eta_lowrank) - pnorm(eta_oracle_lowrank))^2)),
    estimated_vs_oracle_subspace_distance = subspace_distance(lowrank$L, oracle_lowrank$L, H),
    estimated_vs_full_oracle_z_rmse = sqrt(mean((eta_lowrank - sim$Z_latent)^2)),
    oracle_lowrank_vs_full_oracle_z_rmse = sqrt(mean((eta_oracle_lowrank - sim$Z_latent)^2))
  )
  primary_range_meta <- metadata_pair(primary_loading_range)
  cross_range_meta <- metadata_pair(cross_loading_range)
  cross_prob_meta <- if (is.null(cross_loading_prob)) NA_real_ else cross_loading_prob

  signals <- list(
    estimated_Z = list(
      signal_source = "estimated_Z",
      lowrank = lowrank,
      signal_seconds = lowrank_seconds,
      signal_converged = isTRUE(lowrank$converged),
      signal_completed_iter = lowrank$n_completed,
      stage1_signal_rmse = sqrt(mean((true_signal - lowrank$L)^2)),
      stage1_subspace_distance = subspace_distance(true_signal, lowrank$L, H)
    ),
    oracle_Z = list(
      signal_source = "oracle_Z",
      lowrank = oracle_lowrank,
      signal_seconds = oracle_lowrank$signal_seconds,
      signal_converged = TRUE,
      signal_completed_iter = 1L,
      stage1_signal_rmse = sqrt(mean((true_signal - oracle_lowrank$L)^2)),
      stage1_subspace_distance = subspace_distance(true_signal, oracle_lowrank$L, H)
    )
  )

  rotation_arms <- list(
    fastica_only = list(kind = "fastica"),
    mixture_fastica_start = list(
      kind = "mixture",
      n_random_starts = 0L,
      n_ica_starts = max(1L, fastica_starts),
      include_identity_start = FALSE
    ),
    mixture_identity_start = list(
      kind = "mixture",
      n_random_starts = 0L,
      n_ica_starts = 0L,
      include_identity_start = TRUE
    )
  )
  rotation_arms <- rotation_arms[rotation_methods]

  out <- list()
  idx <- 0L
  for (signal_name in names(signals)) {
    signal <- signals[[signal_name]]
    spectral <- spectral_scores_from_lowrank_signal(signal$lowrank$L, H)
    for (method in names(rotation_arms)) {
      arm <- rotation_arms[[method]]
      method_id <- match(method, names(rotation_arms))
      signal_id <- match(signal_name, names(signals))
      pretrain <- if (identical(arm$kind, "fastica")) {
        fit_fastica_direct_pretrain(
          X = sim$X_binary,
          H = H,
          G = G,
          lowrank = signal$lowrank,
          spectral = spectral,
          seed = base_seed + 3000017L * rep + 1009L * signal_id + p + H,
          rotation_method = method
        )
      } else {
        fit_mixture_criterion_pretrain(
          X = sim$X_binary,
          H = H,
          G = G,
          lowrank = signal$lowrank,
          spectral = spectral,
          rotation_method = method,
          seed = base_seed + 2000003L * rep + 1009L * signal_id + 997L * method_id + p + H,
          n_random_starts_use = arm$n_random_starts,
          n_ica_starts_use = arm$n_ica_starts,
          include_identity_start_use = arm$include_identity_start
        )
      }

      pre_seconds <- signal$signal_seconds + pretrain$rotation_seconds
      pre_eval <- evaluate_estimate(
        stage = "pretrain",
        rotation_method = method,
        sim = sim,
        H = H,
        G = G,
        F_hat = pretrain$F_hat,
        Lambda_hat = pretrain$Lambda_hat,
        alpha_hat = pretrain$alpha_hat,
        mixture_fits = pretrain$mixture_fits,
        seconds = pre_seconds
      )

      refined <- refine_pretrain(
        X = sim$X_binary,
        pretrain = pretrain,
        seed = base_seed + 4000037L * rep + 1009L * signal_id + 997L * method_id + p + H
      )
      refined_seconds <- signal$signal_seconds + pretrain$rotation_seconds + refined$fit_seconds
      ref_eval <- evaluate_estimate(
        stage = "refined",
        rotation_method = method,
        sim = sim,
        H = H,
        G = G,
        F_hat = refined$F_hat,
        Lambda_hat = refined$Lambda_hat,
        alpha_hat = refined$alpha_hat,
        mixture_fits = refined$mixture_fits,
        seconds = refined_seconds
      )

      common_extra <- data.frame(
        n = n,
        p = p,
        H = H,
        G_config = g_label(G),
        rep = rep,
        block_size_mode = sim$block_size_mode,
        loading_design = sim$loading_design,
        loading_design_label = sim$loading_design_label,
        data_seed = sim$data_seed,
        signal_source = signal$signal_source,
        mixture_scenario = mixture_scenario,
        separation = separation,
        viroli_smoke_g2_pi = paste(viroli_smoke_g2_pi, collapse = "-"),
        viroli_smoke_g2_mu_multiplier = viroli_smoke_g2_mu_multiplier,
        viroli_smoke_g2_sd = paste(viroli_smoke_g2_sd, collapse = "-"),
        viroli_smoke_g3_pi = paste(viroli_smoke_g3_pi, collapse = "-"),
        viroli_smoke_g3_mu_multiplier = viroli_smoke_g3_mu_multiplier,
        viroli_smoke_g3_sd = paste(viroli_smoke_g3_sd, collapse = "-"),
        primary_loading_min = primary_range_meta[1L],
        primary_loading_max = primary_range_meta[2L],
        cross_loading_min = cross_range_meta[1L],
        cross_loading_max = cross_range_meta[2L],
        cross_loading_prob = cross_prob_meta,
        pretrain_loading_penalty = pretrain_loading_penalty,
        rotation_loading_l1_penalty = rotation_loading_l1_penalty,
        lambda_l1_penalty = lambda_l1_penalty,
        workers = workers,
        parallel = parallel_enabled,
        fastica_selection = fastica_selection,
        min_block_size = min(sim$block_sizes),
        median_block_size = median(sim$block_sizes),
        max_block_size = max(sim$block_sizes),
        stage1_signal_rmse = signal$stage1_signal_rmse,
        stage1_subspace_distance = signal$stage1_subspace_distance,
        estimated_vs_oracle_signal_rmse = estimated_oracle_gap$estimated_vs_oracle_signal_rmse,
        estimated_vs_oracle_alpha_rmse = estimated_oracle_gap$estimated_vs_oracle_alpha_rmse,
        estimated_vs_oracle_linear_predictor_rmse = estimated_oracle_gap$estimated_vs_oracle_linear_predictor_rmse,
        estimated_vs_oracle_probability_rmse = estimated_oracle_gap$estimated_vs_oracle_probability_rmse,
        estimated_vs_oracle_subspace_distance = estimated_oracle_gap$estimated_vs_oracle_subspace_distance,
        estimated_vs_full_oracle_z_rmse = estimated_oracle_gap$estimated_vs_full_oracle_z_rmse,
        oracle_lowrank_vs_full_oracle_z_rmse = estimated_oracle_gap$oracle_lowrank_vs_full_oracle_z_rmse,
        signal_seconds = signal$signal_seconds,
        lowrank_seconds = signal$signal_seconds,
        lowrank_converged = signal$signal_converged,
        lowrank_completed_iter = signal$signal_completed_iter,
        eigengap_diagnostic_seconds = eigengap_seconds,
        rotation_converged = isTRUE(pretrain$rotation_converged),
        rotation_completed_outer = pretrain$rotation_completed_outer,
        rotation_seconds = pretrain$rotation_seconds,
        stringsAsFactors = FALSE
      )
      common_extra <- cbind(common_extra, eigengap_diagnostics)
      pre_extra <- cbind(
        common_extra,
        data.frame(
          refinement_seconds = 0,
          end_to_end_seconds = pre_seconds,
          refinement_converged = NA,
          refinement_completed_iter = NA_integer_
        )
      )
      ref_extra <- cbind(
        common_extra,
        data.frame(
          refinement_seconds = refined$fit_seconds,
          end_to_end_seconds = refined_seconds,
          refinement_converged = isTRUE(refined$converged),
          refinement_completed_iter = if (!is.null(refined$completed_iter)) refined$completed_iter else NA_integer_
        )
      )
      idx <- idx + 1L
      out[[idx]] <- cbind(pre_extra, pre_eval)
      idx <- idx + 1L
      out[[idx]] <- cbind(ref_extra, ref_eval)
    }
  }
  rbind_fill(out)
}

ensure_design_columns <- function(results) {
  if (!"block_size_mode" %in% names(results)) {
    results$block_size_mode <- "ifeval_min30"
  }
  if (!"loading_design" %in% names(results)) {
    results$loading_design <- "balanced_moderate_dense_signed_cross"
  }
  if (!"loading_design_label" %in% names(results)) {
    results$loading_design_label <- sample_size_loading_design_label(results$loading_design)
  }
  results
}

summarize_results <- function(results) {
  results <- ensure_design_columns(results)
  metrics <- c(
    "factor_score_rmse", "lambda_rmse", "alpha_rmse", "probability_rmse",
    "marginal_mu_rmse", "marginal_var_rmse", "marginal_weight_rmse",
    "estimated_vs_oracle_signal_rmse", "estimated_vs_oracle_alpha_rmse",
    "estimated_vs_oracle_linear_predictor_rmse", "estimated_vs_oracle_probability_rmse",
    "estimated_vs_oracle_subspace_distance", "estimated_vs_full_oracle_z_rmse",
    "oracle_lowrank_vs_full_oracle_z_rmse",
    "seconds", "end_to_end_seconds", "signal_seconds", "rotation_seconds",
    "refinement_seconds"
  )
  split_vars <- c(
    "stage", "signal_source", "rotation_method", "block_size_mode",
    "loading_design", "loading_design_label", "n", "p", "H", "G_config"
  )
  pieces <- list()
  idx <- 0L
  for (metric in metrics) {
    if (!metric %in% names(results)) next
    agg <- aggregate(
      results[[metric]],
      results[split_vars],
      function(x) c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE),
                    sd = sd(x, na.rm = TRUE), n = sum(is.finite(x)))
    )
    vals <- if (is.matrix(agg$x)) {
      as.data.frame(agg$x)
    } else {
      do.call(data.frame, agg$x)
    }
    names(vals) <- paste0(metric, "_", c("mean", "median", "sd", "n"))
    idx <- idx + 1L
    pieces[[idx]] <- cbind(agg[split_vars], vals)
  }
  if (!length(pieces)) return(data.frame())
  Reduce(function(a, b) merge(a, b, by = split_vars, all = TRUE), pieces)
}

method_colors <- c(
  fastica_only = "#c4313b",
  mixture_fastica_start = "#2b8a57",
  mixture_identity_start = "#2b6db6"
)
method_labels <- c(
  fastica_only = "FastICA only",
  mixture_fastica_start = "FastICA start + mixture rotation",
  mixture_identity_start = "identity start + mixture rotation"
)
method_short_labels <- c(
  fastica_only = "FastICA",
  mixture_fastica_start = "Mix+ICA",
  mixture_identity_start = "Mix+Id"
)
method_pch <- c(fastica_only = 17, mixture_fastica_start = 16, mixture_identity_start = 15)

plot_metric_boxplots <- function(results, out_file) {
  results <- ensure_design_columns(results)
  refined <- results[results$stage == "refined", , drop = FALSE]
  metrics <- c(
    factor_score_rmse = "factor scores",
    lambda_rmse = "loadings",
    alpha_rmse = "intercepts",
    probability_rmse = "probabilities",
    marginal_mu_rmse = "mix means",
    marginal_var_rmse = "mix variances"
  )
  png(out_file, width = 1900, height = 1200, res = 170)
  op <- par(mfrow = c(2, 3), mar = c(4.7, 4.6, 3, 1), oma = c(0, 0, 3, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  refined$plot_group <- paste(refined$signal_source, refined$rotation_method, sep = ": ")
  ordered_groups <- as.vector(outer(c("estimated_Z", "oracle_Z"), names(method_colors), paste, sep = ": "))
  ordered_groups <- ordered_groups[ordered_groups %in% unique(refined$plot_group)]
  group_methods <- sub("^.*: ", "", ordered_groups)
  group_sources <- sub(": .*$", "", ordered_groups)
  group_labels <- paste0(ifelse(group_sources == "estimated_Z", "est", "oracle"), "\n", method_short_labels[group_methods])
  for (metric in names(metrics)) {
    vals <- split(refined[[metric]], refined$plot_group)[ordered_groups]
    boxplot(
      vals,
      col = method_colors[group_methods],
      border = "#333333",
      ylab = if (metric == "probability_rmse") "RMSE" else "RMSE",
      names = group_labels,
      las = 1,
      cex.axis = 0.72,
      main = metrics[[metric]]
    )
    if (length(vals) == 6L) abline(v = 3.5, col = "#777777", lty = 3)
    stripchart(vals, vertical = TRUE, method = "jitter", add = TRUE, pch = 16, cex = 0.6, col = "#00000055")
    grid(nx = NA, ny = NULL, col = "#eeeeee")
  }
  mtext("Rotation diagnostic: refined recovery across all cells", outer = TRUE, font = 2, cex = 1.2)
}

plot_runtime_boxplot <- function(results, out_file) {
  results <- ensure_design_columns(results)
  refined <- results[results$stage == "refined", , drop = FALSE]
  refined$plot_group <- paste(refined$signal_source, refined$rotation_method, sep = ": ")
  ordered_groups <- as.vector(outer(c("estimated_Z", "oracle_Z"), names(method_colors), paste, sep = ": "))
  ordered_groups <- ordered_groups[ordered_groups %in% unique(refined$plot_group)]
  group_methods <- sub("^.*: ", "", ordered_groups)
  group_sources <- sub(": .*$", "", ordered_groups)
  group_labels <- paste0(ifelse(group_sources == "estimated_Z", "est", "oracle"), "\n", method_short_labels[group_methods])

  png(out_file, width = 1200, height = 800, res = 170)
  op <- par(mar = c(4.8, 5, 3.2, 1), oma = c(0, 0, 0, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  vals <- split(refined$end_to_end_seconds, refined$plot_group)[ordered_groups]
  boxplot(
    vals,
    col = method_colors[group_methods],
    border = "#333333",
    names = group_labels,
    ylab = "end-to-end seconds",
    main = "Rotation diagnostic: end-to-end runtime",
    cex.axis = 0.82
  )
  if (length(vals) == 6L) abline(v = 3.5, col = "#777777", lty = 3)
  stripchart(vals, vertical = TRUE, method = "jitter", add = TRUE, pch = 16, cex = 0.75, col = "#00000055")
  grid(nx = NA, ny = NULL, col = "#eeeeee")
}

plot_cell_lines <- function(results, out_file) {
  results <- ensure_design_columns(results)
  refined <- results[results$stage == "refined", , drop = FALSE]
  metrics <- c(factor_score_rmse = "factor score RMSE", lambda_rmse = "loading RMSE")
  cells <- unique(refined[c("signal_source", "block_size_mode", "loading_design_label", "H", "G_config")])
  cells <- cells[order(cells$signal_source, cells$block_size_mode, cells$loading_design_label,
                       cells$H, cells$G_config), , drop = FALSE]
  png(out_file, width = 2000, height = max(1300, 350 * nrow(cells)), res = 170)
  op <- par(mfrow = c(nrow(cells), length(metrics)), mar = c(4.8, 4.8, 3, 1), oma = c(0, 0, 3, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    dd <- refined[
      refined$signal_source == cell$signal_source &
        refined$block_size_mode == cell$block_size_mode &
        refined$loading_design_label == cell$loading_design_label &
        refined$H == cell$H &
        refined$G_config == cell$G_config,
      ,
      drop = FALSE
    ]
    for (metric in names(metrics)) {
      agg <- aggregate(dd[[metric]], dd[c("rotation_method", "n", "p")], mean, na.rm = TRUE)
      names(agg)[names(agg) == "x"] <- "value"
      yr <- range(agg$value, finite = TRUE)
      plot(
        NA,
        xlim = range(p_values),
        ylim = yr + c(-0.06, 0.10) * diff(yr + c(-1e-6, 1e-6)),
        xlab = "p",
        ylab = "RMSE",
        main = sprintf(
          "%s | %s | %s | H=%d, G=%s: %s",
          cell$signal_source,
          cell$block_size_mode,
          cell$loading_design_label,
          cell$H,
          cell$G_config,
          metrics[[metric]]
        )
      )
      grid(col = "#eeeeee")
      for (method in unique(agg$rotation_method)) {
        for (n_val in sort(unique(agg$n))) {
          aa <- agg[agg$rotation_method == method & agg$n == n_val, , drop = FALSE]
          aa <- aa[order(aa$p), ]
          lty <- if (n_val == min(n_values)) 2 else 1
          points(aa$p, aa$value, pch = method_pch[[method]],
                 col = method_colors[[method]], cex = if (n_val == min(n_values)) 0.9 else 1.1)
          lines(aa$p, aa$value, col = method_colors[[method]], lwd = 2, lty = lty)
        }
      }
      if (i == 1L && metric == names(metrics)[1L]) {
        legend(
          "topright",
          legend = c(unname(method_short_labels), paste0("n=", min(n_values)), paste0("n=", max(n_values))),
          col = c(method_colors, "#333333", "#333333"),
          lty = c(rep(1, length(method_colors)), 2, 1),
          pch = c(method_pch, NA, NA),
          bty = "n",
          cex = 0.66
        )
      }
    }
  }
  mtext("Rotation diagnostic: refined recovery by p", outer = TRUE, font = 2, cex = 1.2)
}

plot_pretrain_refine_delta <- function(results, out_file) {
  results <- ensure_design_columns(results)
  metrics <- c(factor_score_rmse = "factor score RMSE", lambda_rmse = "loading RMSE",
               marginal_mu_rmse = "mixture mean RMSE", probability_rmse = "probability RMSE")
  png(out_file, width = 1700, height = 1100, res = 170)
  op <- par(mfrow = c(2, 2), mar = c(4.8, 4.8, 3, 1), oma = c(4.2, 0, 3, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  for (metric in names(metrics)) {
    agg <- aggregate(results[[metric]], results[c("stage", "signal_source", "rotation_method", "H")], mean, na.rm = TRUE)
    names(agg)[names(agg) == "x"] <- "value"
    yr <- range(agg$value, finite = TRUE)
    plot(
      NA,
      xlim = c(0.8, 2.2),
      ylim = yr + c(-0.08, 0.08) * diff(yr + c(-1e-6, 1e-6)),
      xaxt = "n",
      xlab = "",
      ylab = "RMSE",
      main = metrics[[metric]]
    )
    axis(1, at = c(1, 2), labels = c("pretrain", "refined"))
    grid(nx = NA, col = "#eeeeee")
    for (method in unique(agg$rotation_method)) {
      for (source in unique(agg$signal_source)) {
        for (H_val in sort(unique(agg$H))) {
          aa <- agg[
            agg$rotation_method == method &
              agg$signal_source == source &
              agg$H == H_val,
            ,
            drop = FALSE
          ]
          aa$x <- match(aa$stage, c("pretrain", "refined"))
          aa <- aa[order(aa$x), ]
          lty <- if (source == "estimated_Z") 1 else 2
          lines(aa$x, aa$value, col = method_colors[[method]], lwd = if (H_val == max(h_values)) 2.6 else 1.8,
                lty = lty)
          points(aa$x, aa$value, col = method_colors[[method]], pch = method_pch[[method]])
        }
      }
    }
  }
  legend(
    "bottom",
    inset = -0.02,
    xpd = NA,
    horiz = TRUE,
    legend = c(unname(method_short_labels), "estimated Z", "oracle Z"),
    col = c(method_colors, "#333333", "#333333"),
    lty = c(rep(1, length(method_colors)), 1, 2),
    pch = c(method_pch, NA, NA),
    bty = "n",
    cex = 0.8
  )
  mtext("Rotation diagnostic: pretraining versus refinement", outer = TRUE, font = 2, cex = 1.2)
}

write_pretrain_refine_factor_pairs <- function(results, out_file) {
  results <- ensure_design_columns(results)
  key_cols <- c(
    "signal_source", "rotation_method", "block_size_mode", "loading_design",
    "loading_design_label", "n", "p", "H", "G_config", "rep"
  )
  metric_cols <- c(
    "factor_score_rmse",
    "factor_score_raw_rmse",
    "mean_factor_abs_cor",
    "lambda_rmse",
    "probability_rmse"
  )
  keep_cols <- c(key_cols, metric_cols)
  pre <- results[results$stage == "pretrain", keep_cols, drop = FALSE]
  ref <- results[results$stage == "refined", keep_cols, drop = FALSE]
  names(pre)[match(metric_cols, names(pre))] <- paste0("pretrain_", metric_cols)
  names(ref)[match(metric_cols, names(ref))] <- paste0("refined_", metric_cols)
  paired <- merge(pre, ref, by = key_cols, all = FALSE)
  paired$delta_factor_score_rmse <- paired$refined_factor_score_rmse - paired$pretrain_factor_score_rmse
  paired$delta_factor_score_raw_rmse <- paired$refined_factor_score_raw_rmse - paired$pretrain_factor_score_raw_rmse
  paired$delta_mean_factor_abs_cor <- paired$refined_mean_factor_abs_cor - paired$pretrain_mean_factor_abs_cor
  paired$delta_lambda_rmse <- paired$refined_lambda_rmse - paired$pretrain_lambda_rmse
  paired$delta_probability_rmse <- paired$refined_probability_rmse - paired$pretrain_probability_rmse
  paired <- paired[order(paired$signal_source, paired$rotation_method, paired$H,
                         paired$G_config, paired$n, paired$p, paired$rep), , drop = FALSE]
  write.csv(paired, out_file, row.names = FALSE)
  invisible(paired)
}

write_eigengap_outputs <- function(results, out_dir) {
  if (!("estimated_z_h_hat_largest_eigengap" %in% names(results))) return(invisible(NULL))
  id_cols <- intersect(
    c("mixture_scenario", "block_size_mode", "loading_design", "n", "p", "H", "G_config", "rep"),
    names(results)
  )
  diagnostic_cols <- grep("^(eigengap_|estimated_z_)", names(results), value = TRUE)
  dataset_rows <- unique(results[c(id_cols, diagnostic_cols)])
  dataset_file <- file.path(out_dir, "rotation_fastica_eigengap_dataset_results.csv")
  write.csv(dataset_rows, dataset_file, row.names = FALSE)

  accuracy_cols <- intersect(
    c(
      "estimated_z_largest_eigengap_selects_true_H",
      "estimated_z_largest_relative_eigengap_selects_true_H",
      "estimated_z_largest_sv_ratio_selects_true_H",
      "estimated_z_eigengap_at_true_H",
      "estimated_z_relative_eigengap_at_true_H",
      "estimated_z_sv_ratio_at_true_H",
      "eigengap_diagnostic_seconds"
    ),
    names(dataset_rows)
  )
  group_cols <- intersect(c("mixture_scenario", "n", "p", "H", "G_config"), names(dataset_rows))
  summary <- aggregate(
    dataset_rows[accuracy_cols],
    dataset_rows[group_cols],
    mean,
    na.rm = TRUE
  )
  counts <- aggregate(
    rep(1L, nrow(dataset_rows)),
    dataset_rows[group_cols],
    length
  )
  names(counts)[ncol(counts)] <- "n_replications"
  summary <- merge(summary, counts, by = group_cols, all.x = TRUE, sort = FALSE)
  summary_file <- file.path(out_dir, "rotation_fastica_eigengap_setting_summary.csv")
  write.csv(summary, summary_file, row.names = FALSE)

  cat("  ", dataset_file, "\n", sep = "")
  cat("  ", summary_file, "\n", sep = "")
  invisible(list(dataset = dataset_rows, summary = summary))
}

write_results_and_plots <- function(results, results_file, summary_file, out_dir) {
  results <- ensure_design_columns(results)
  summary <- summarize_results(results)
  write.csv(results, results_file, row.names = FALSE)
  write.csv(summary, summary_file, row.names = FALSE)

  paired_file <- file.path(out_dir, "rotation_fastica_pretrain_refine_factor_pairs.csv")
  refined_boxplot <- file.path(out_dir, "rotation_fastica_refined_boxplots.png")
  by_p_plot <- file.path(out_dir, "rotation_fastica_factor_loading_by_p.png")
  pretrain_plot <- file.path(out_dir, "rotation_fastica_pretrain_vs_refine.png")
  runtime_plot <- file.path(out_dir, "rotation_fastica_end_to_end_runtime.png")

  write_pretrain_refine_factor_pairs(results, paired_file)
  plot_metric_boxplots(results, refined_boxplot)
  plot_cell_lines(results, by_p_plot)
  plot_pretrain_refine_delta(results, pretrain_plot)
  plot_runtime_boxplot(results, runtime_plot)

  cat("\nWrote:\n")
  cat("  ", results_file, "\n", sep = "")
  cat("  ", summary_file, "\n", sep = "")
  cat("  ", paired_file, "\n", sep = "")
  cat("  ", refined_boxplot, "\n", sep = "")
  cat("  ", by_p_plot, "\n", sep = "")
  cat("  ", pretrain_plot, "\n", sep = "")
  cat("  ", runtime_plot, "\n", sep = "")
  write_eigengap_outputs(results, out_dir)

  refined <- results[results$stage == "refined", , drop = FALSE]
  if (nrow(refined)) {
    compact <- aggregate(
      cbind(factor_score_rmse, lambda_rmse, marginal_mu_rmse, probability_rmse, end_to_end_seconds) ~
        signal_source + rotation_method,
      refined,
      median,
      na.rm = TRUE
    )
    cat("\nRefined median metrics across diagnostic grid:\n")
    print(compact, row.names = FALSE)
  }

  invisible(summary)
}

n_values <- get_env("N_VALUES", c(100L, 200L), parse_ints)
p_values <- get_env("P_VALUES", c(500L, 1000L), parse_ints)
h_values <- get_env("H_VALUES", c(5L, 10L), parse_ints)
g_types <- split_csv(get_env("G_TYPES", "all2,all3"))
rep_values <- get_env("REP_VALUES", 1:2, parse_ints)

base_seed <- get_env("SEED", 20260731L, as.integer)
loading_seed <- get_env("LOADING_SEED", 20772531L, as.integer)
mixture_seed <- get_env("MIXTURE_SEED", 22262731L, as.integer)
separation <- get_env("SEPARATION", 2, as.numeric)
mixture_scenario <- get_env("MIXTURE_SCENARIO", "separated", as.character)
viroli_smoke_g2_pi <- parse_numeric_vector("VIROLI_SMOKE_G2_PI", c(0.50, 0.50))
viroli_smoke_g2_mu_multiplier <- get_env("VIROLI_SMOKE_G2_MU_MULTIPLIER", 1.35, as.numeric)
viroli_smoke_g2_sd <- parse_numeric_vector("VIROLI_SMOKE_G2_SD", c(0.45, 0.45))
viroli_smoke_g3_pi <- parse_numeric_vector("VIROLI_SMOKE_G3_PI", c(0.30, 0.40, 0.30))
viroli_smoke_g3_mu_multiplier <- get_env("VIROLI_SMOKE_G3_MU_MULTIPLIER", 1.35, as.numeric)
viroli_smoke_g3_sd <- parse_numeric_vector("VIROLI_SMOKE_G3_SD", c(0.45, 0.65, 0.45))
primary_loading_range <- parse_optional_numeric_pair("PRIMARY_LOADING_RANGE", NULL)
cross_loading_range <- parse_optional_numeric_pair("CROSS_LOADING_RANGE", NULL)
cross_loading_prob <- parse_optional_numeric("CROSS_LOADING_PROB", NULL)
loading_designs <- get_env(
  "LOADING_DESIGNS",
  "balanced_moderate_dense_signed_cross",
  split_csv
)
loading_designs <- normalize_sample_size_loading_design(loading_designs)
block_size_modes <- split_csv(get_env("BLOCK_SIZE_MODES", "ifeval_min30"))
block_size_modes <- match.arg(
  block_size_modes,
  choices = c("balanced", "ifeval_like", "moderate_ifeval_like", "ifeval_min30"),
  several.ok = TRUE
)
intercept_sd <- get_env("INTERCEPT_SD", 0.45, as.numeric)
intercept_block_span <- get_env("INTERCEPT_BLOCK_SPAN", 1.6, as.numeric)
intercept_clip <- get_env("INTERCEPT_CLIP", 1.75, as.numeric)
intercept_mode <- get_env("INTERCEPT_MODE", "ifeval_like", as.character)
intercept_mode <- match.arg(intercept_mode, c("none", "ifeval_like", "random", "block", "viroli_smoke"))

em_max_iter <- get_env("EM_SVD_ITER", 35L, as.integer)
em_tol_loglik <- get_env("EM_SVD_TOL_LOGLIK", 1e-5, as.numeric)
em_tol_L <- get_env("EM_SVD_TOL_L", 1e-4, as.numeric)
em_tol_subspace <- get_env("EM_SVD_TOL_SUBSPACE", 2e-3, as.numeric)
em_init_method <- get_env("EM_SVD_INIT", "both")
em_init_z <- get_env("EM_SVD_INIT_Z", "expectation")
run_eigengap_diagnostic <- get_env("RUN_EIGENGAP_DIAGNOSTIC", TRUE, as.logical)
eigengap_fit_rank <- get_env("EIGENGAP_FIT_RANK", max(h_values) + 5L, as.integer)
eigengap_max_candidate_rank <- get_env(
  "EIGENGAP_MAX_CANDIDATE_RANK",
  eigengap_fit_rank - 1L,
  as.integer
)
if (eigengap_fit_rank <= max(h_values)) {
  stop("EIGENGAP_FIT_RANK must exceed the largest true H in the grid.", call. = FALSE)
}

rotation_random_starts <- get_env("ROTATION_RANDOM_STARTS", 1L, as.integer)
rotation_max_outer <- get_env("ROTATION_ITER", 12L, as.integer)
rotation_min_outer <- get_env("ROTATION_MIN_ITER", 2L, as.integer)
riemannian_rotation_steps <- get_env("RIEMANNIAN_ROTATION_STEPS", 8L, as.integer)
rotation_objective_tolerance <- get_env("ROTATION_OBJECTIVE_TOLERANCE", 1e-4, as.numeric)
fastica_starts <- get_env("FASTICA_STARTS", 1L, as.integer)
fastica_functions <- split_csv(get_env("FASTICA_FUNCTIONS", "logcosh"))
fastica_max_iter <- get_env("FASTICA_MAX_ITER", 400L, as.integer)
fastica_tol <- get_env("FASTICA_TOL", 1e-5, as.numeric)
fastica_selection <- get_env("FASTICA_SELECTION", "first")
fastica_selection <- match.arg(fastica_selection, c("first", "best_mixture"))
rotation_methods <- split_csv(get_env(
  "ROTATION_METHODS",
  "fastica_only,mixture_fastica_start,mixture_identity_start"
))
rotation_methods <- match.arg(
  rotation_methods,
  choices = c("fastica_only", "mixture_fastica_start", "mixture_identity_start"),
  several.ok = TRUE
)

n_mix_starts <- get_env("N_MIX_STARTS", 3L, as.integer)
mixture_max_iter <- get_env("MIXTURE_MAX_ITER", 80L, as.integer)
mixture_update <- get_env("MIXTURE_UPDATE", "map")
mu_prior_mean <- get_env("MU_PRIOR_MEAN", 0, as.numeric)
mu_prior_kappa <- get_env("MU_PRIOR_KAPPA", 0.05, as.numeric)
var_prior_shape <- get_env("VAR_PRIOR_SHAPE", 3, as.numeric)
var_prior_scale <- get_env("VAR_PRIOR_SCALE", 2, as.numeric)
weight_prior_alpha <- get_env("WEIGHT_PRIOR_ALPHA", 1, as.numeric)
min_mixture_var <- get_env("MIN_MIXTURE_VAR", 0.05, as.numeric)

base_pretrain_loading_penalty <- get_env("PRETRAIN_LOADING_PENALTY", 5, as.numeric)
base_rotation_loading_l1_penalty <- get_env("ROTATION_LOADING_L1_PENALTY", 5, as.numeric)
base_lambda_l1_penalty <- get_env("LAMBDA_L1_PENALTY", 5, as.numeric)
penalty_by_n <- parse_penalty_by_n(get_env("PENALTY_BY_N", ""))
pretrain_loading_penalty <- base_pretrain_loading_penalty
rotation_loading_l1_penalty <- base_rotation_loading_l1_penalty
lambda_l1_penalty <- base_lambda_l1_penalty
lasso_backend <- get_env("LASSO_BACKEND", "glmnet")
refine_iter <- get_env("REFINE_ITER", 25L, as.integer)
refine_min_iter <- get_env("REFINE_MIN_ITER", 3L, as.integer)
refinement_objective_tolerance <- get_env("REFINE_OBJECTIVE_TOLERANCE", 1e-3, as.numeric)
maxit_per_subject <- get_env("MAXIT_PER_SUBJECT", 50L, as.integer)

parallel_enabled <- get_env("PARALLEL", TRUE, as.logical)
workers <- get_env("WORKERS", 4L, as.integer)
run_label <- get_env("RUN_LABEL", "rotation_vs_fastica")
plot_only <- get_env("PLOT_ONLY", FALSE, as.logical)
resume_existing <- get_env("RESUME_EXISTING", FALSE, as.logical)
out_dir <- get_env(
  "OUT_DIR",
  file.path(repo_root, "results", "diagnostics", run_label)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set_penalties_for_n <- function(n) {
  if (nrow(penalty_by_n)) {
    hit <- match(as.integer(n), penalty_by_n$n)
    if (!is.na(hit)) {
      pretrain_loading_penalty <<- penalty_by_n$penalty[hit]
      rotation_loading_l1_penalty <<- penalty_by_n$penalty[hit]
      lambda_l1_penalty <<- penalty_by_n$penalty[hit]
      return(invisible(TRUE))
    }
  }
  pretrain_loading_penalty <<- base_pretrain_loading_penalty
  rotation_loading_l1_penalty <<- base_rotation_loading_l1_penalty
  lambda_l1_penalty <<- base_lambda_l1_penalty
  invisible(FALSE)
}

load_resume_rows <- function(results_file, expected_rows_per_dataset) {
  if (!file.exists(results_file)) {
    return(list(rows = list(), completed_keys = character(), kept = 0L, dropped = 0L))
  }
  existing <- read.csv(results_file)
  if (!nrow(existing)) {
    return(list(rows = list(), completed_keys = character(), kept = 0L, dropped = 0L))
  }
  existing <- ensure_design_columns(existing)
  id_cols <- intersect(
    c(
      "block_size_mode", "loading_design", "n", "p", "H", "G_config", "rep",
      "signal_source", "rotation_method", "stage"
    ),
    names(existing)
  )
  if (length(id_cols) == 10L) {
    existing <- existing[!duplicated(existing[id_cols], fromLast = TRUE), , drop = FALSE]
  }
  keys <- dataset_key_from_results(existing)
  counts <- table(keys)
  completed_keys <- names(counts[counts >= expected_rows_per_dataset])
  keep <- keys %in% completed_keys
  kept <- sum(keep)
  dropped <- sum(!keep)
  rows <- if (kept) list(existing[keep, , drop = FALSE]) else list()
  list(rows = rows, completed_keys = completed_keys, kept = kept, dropped = dropped)
}

cat("Rotation diagnostic output:", out_dir, "\n")
cat("Signal sources: estimated_Z, oracle_Z\n")
cat("Rotation arms:", paste(rotation_methods, collapse = ","), "\n")
cat("Workers:", workers, " parallel:", parallel_enabled, "\n")
if (nrow(penalty_by_n)) {
  cat("Penalty by n:", paste(paste0(penalty_by_n$n, ":", penalty_by_n$penalty), collapse = ","), "\n")
} else {
  cat("Penalties: pretrain=", base_pretrain_loading_penalty,
      " rotation=", base_rotation_loading_l1_penalty,
      " refinement=", base_lambda_l1_penalty, "\n", sep = "")
}

results_file <- file.path(out_dir, "rotation_fastica_results.csv")
summary_file <- file.path(out_dir, "rotation_fastica_summary.csv")

if (isTRUE(plot_only)) {
  if (!file.exists(results_file)) {
    stop("PLOT_ONLY=TRUE but no rotation_fastica_results.csv exists in ", out_dir, call. = FALSE)
  }
  results <- read.csv(results_file)
  cat("Plot-only grid from saved CSV: n=", paste(sort(unique(results$n)), collapse = ","),
      " p=", paste(sort(unique(results$p)), collapse = ","),
      " H=", paste(sort(unique(results$H)), collapse = ","),
      " G=", paste(sort(unique(results$G_config)), collapse = ","),
      " reps=", paste(sort(unique(results$rep)), collapse = ","), "\n", sep = "")
  write_results_and_plots(results, results_file, summary_file, out_dir)
  quit(save = "no", status = 0L)
}

cat("Grid: block=", paste(block_size_modes, collapse = ","),
    " loading=", paste(sample_size_loading_design_label(loading_designs), collapse = ","),
    " n=", paste(n_values, collapse = ","), " p=", paste(p_values, collapse = ","),
    " H=", paste(h_values, collapse = ","), " G=", paste(g_types, collapse = ","),
    " reps=", paste(rep_values, collapse = ","), "\n", sep = "")

expected_rows_per_dataset <- 2L * length(rotation_methods) * 2L
resume_state <- if (isTRUE(resume_existing)) {
  load_resume_rows(results_file, expected_rows_per_dataset)
} else {
  list(rows = list(), completed_keys = character(), kept = 0L, dropped = 0L)
}
if (isTRUE(resume_existing)) {
  cat(
    "Resume mode: kept ", resume_state$kept,
    " rows from ", length(resume_state$completed_keys),
    " completed dataset cells; dropped ", resume_state$dropped,
    " incomplete/duplicate rows.\n",
    sep = ""
  )
}

all_rows <- resume_state$rows
completed_keys <- resume_state$completed_keys
idx <- length(all_rows)
total_dataset_cells <- length(block_size_modes) * length(loading_designs) *
  length(h_values) * length(g_types) * length(n_values) * length(p_values) *
  length(rep_values)
seen_dataset_cells <- 0L
for (block_size_mode in block_size_modes) {
  for (loading_design in loading_designs) {
    for (H in h_values) {
      for (g_type in g_types) {
        G <- g_vector(g_type, H)
        for (n in n_values) {
          set_penalties_for_n(n)
          for (p in p_values) {
            for (rep in rep_values) {
              seen_dataset_cells <- seen_dataset_cells + 1L
              key <- dataset_key(block_size_mode, loading_design, n, p, H, g_label(G), rep)
              if (key %in% completed_keys) {
                if (seen_dataset_cells %% 25L == 0L || seen_dataset_cells == total_dataset_cells) {
                  cat(sprintf(
                    "[%s] skipped completed cells through %d/%d\n",
                    format(Sys.time(), "%H:%M:%S"),
                    seen_dataset_cells,
                    total_dataset_cells
                  ))
                }
                next
              }
              cat(sprintf(
                "Active penalties: pretrain=%.3f rotation=%.3f refinement=%.3f; cell %d/%d\n",
                pretrain_loading_penalty,
                rotation_loading_l1_penalty,
                lambda_l1_penalty,
                seen_dataset_cells,
                total_dataset_cells
              ))
              idx <- idx + 1L
              all_rows[[idx]] <- fit_one_dataset(
                block_size_mode = block_size_mode,
                loading_design = loading_design,
                n = n,
                p = p,
                H = H,
                G = G,
                rep = rep
              )
              partial <- rbind_fill(all_rows)
              write.csv(partial, results_file, row.names = FALSE)
              completed_keys <- c(completed_keys, key)
            }
          }
        }
      }
    }
  }
}

results <- rbind_fill(all_rows)
write_results_and_plots(results, results_file, summary_file, out_dir)
