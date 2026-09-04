#!/usr/bin/env Rscript

# Deterministic EM-SVD pretraining for probit independent factor analysis.
#
# This module implements the EM-SVD low-rank probit signal pretraining step
# summarized in writeup/final_simulation_design/final_simulation_design_algorithms.pdf.
# It is intentionally kept separate from R/binary_probit_pretraining.R so the
# older sampled-Z pretraining remains available as a legacy ablation.
#
# Required helpers from R/binary_probit_pretraining.R:
#   initialize_binary_intercepts()
#   estimate_mixture_ica_unknown_G()
#   mixture_responsibilities()
#   binary_profile_id()
#   normalize_G_fixed()
#   project_to_orthogonal()
#
# Required helper from R/binary_probit_refinement.R for the wrapper:
#   fit_binary_probit_refinement()

probit_ifa_loglik_lowrank <- function(X, alpha, L) {
  X <- as.matrix(X)
  eta <- sweep(L, 2L, alpha, "+")
  s <- 2 * X - 1
  sum(pnorm(s * eta, log.p = TRUE))
}

probit_ifa_inverse_mills <- function(x, max_value = 1e6) {
  # Compute phi(x) / Phi(x) on the log scale.  This is the only delicate
  # numerical operation in the deterministic E-step.
  out <- exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
  out[!is.finite(out)] <- max_value
  pmin(out, max_value)
}

probit_ifa_truncated_mean <- function(X, alpha, L) {
  # E[Z_ij | X_ij, alpha_j + L_ij] under the probit augmentation.
  X <- as.matrix(X)
  eta <- sweep(L, 2L, alpha, "+")
  s <- 2 * X - 1
  eta + s * probit_ifa_inverse_mills(s * eta)
}

probit_ifa_truncated_sample <- function(X, alpha, L, seed = NULL) {
  # Draw Z_ij | X_ij, alpha_j + L_ij from the probit augmentation law.  This is
  # the stochastic counterpart of probit_ifa_truncated_mean().
  if (!exists("rtruncnorm_binary_vec", mode = "function", inherits = TRUE)) {
    stop("rtruncnorm_binary_vec() must be sourced before probit_ifa_truncated_sample().")
  }
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  eta <- sweep(as.matrix(L), 2L, as.numeric(alpha), "+")
  Z <- eta
  observed <- !is.na(X)

  one <- observed & X == 1
  zero <- observed & X == 0
  if (any(one)) Z[one] <- rtruncnorm_binary_vec(eta[one], 1, 0, Inf)
  if (any(zero)) Z[zero] <- rtruncnorm_binary_vec(eta[zero], 1, -Inf, 0)
  if (any(!observed)) Z[!observed] <- rnorm(sum(!observed), eta[!observed], 1)
  Z
}

rank_H_centered_projection <- function(W, H) {
  # Solve min_{alpha,L: rank(L)<=H, colMeans(L)=0} ||W - 1 alpha' - L||_F^2.
  alpha <- colMeans(W)
  Wc <- sweep(W, 2L, alpha, "-")
  dec <- svd(Wc, nu = H, nv = H)
  dH <- dec$d[seq_len(H)]
  U <- dec$u[, seq_len(H), drop = FALSE]
  V <- dec$v[, seq_len(H), drop = FALSE]
  L <- U %*% diag(dH, H, H) %*% t(V)
  list(alpha = as.numeric(alpha), L = L, svd = dec)
}

stochastic_rank_H_centered_projection <- function(
    X,
    alpha,
    L,
    H,
    n_draws = 5L,
    average = c("projected_signal", "sample_mean"),
    seed = NULL,
    parallel = FALSE,
    workers = NULL) {
  # Monte Carlo version of the EM-SVD projection.  There are two useful
  # variants:
  #
  #   sample_mean:
  #       Draw Z^(b), average the sampled latent matrices, then take one rank-H
  #       centered SVD.  With many draws this approaches the usual deterministic
  #       E-step projection P_H(E[Z | X]).
  #
  #   projected_signal:
  #       Draw Z^(b), take a rank-H centered SVD of each draw, average the
  #       resulting low-rank fitted signals, then project that average back to
  #       rank H.  This estimates E[P_H(Z) | X], which keeps augmentation
  #       variability in the low-rank step instead of collapsing immediately to
  #       the conditional mean.
  average <- match.arg(average)
  n_draws <- max(1L, as.integer(n_draws))
  seeds <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, n_draws)
  } else {
    as.integer(seed + 7919L * seq_len(n_draws))
  }

  draws <- parallel_lapply(
    seq_len(n_draws),
    function(b) {
      Zb <- probit_ifa_truncated_sample(X, alpha, L, seed = seeds[b])
      if (average == "sample_mean") return(Zb)
      proj <- rank_H_centered_projection(Zb, H)
      sweep(proj$L, 2L, proj$alpha, "+")
    },
    parallel = parallel,
    workers = workers
  )
  W_bar <- Reduce("+", draws) / n_draws
  projection <- rank_H_centered_projection(W_bar, H)
  projection$n_stochastic_draws <- n_draws
  projection$stochastic_average <- average
  projection
}

fit_lowrank_probit_em_svd_one_start <- function(
    X,
    H,
    max_iter = 50L,
    tol_loglik = 1e-5,
    tol_L = NULL,
    projection_update = c("expectation", "sample_once", "stochastic_average"),
    stochastic_svd_draws = 5L,
    stochastic_svd_average = c("projected_signal", "sample_mean"),
    alpha_init = NULL,
    L_init = NULL,
    seed = 1L,
    parallel = FALSE,
    workers = NULL,
    verbose = FALSE) {
  projection_update <- match.arg(projection_update)
  stochastic_svd_average <- match.arg(stochastic_svd_average)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  H <- as.integer(H)
  if (H < 1L || H > min(n, p)) stop("H must be between 1 and min(n, p).")

  alpha <- if (is.null(alpha_init)) initialize_binary_intercepts(X) else as.numeric(alpha_init)
  if (length(alpha) != p) stop("alpha_init has wrong length.")

  L <- if (is.null(L_init)) matrix(0, n, p) else as.matrix(L_init)
  if (!identical(dim(L), dim(X))) stop("L_init must have the same dimensions as X.")
  L <- sweep(L, 2L, colMeans(L), "-")

  history <- vector("list", max_iter)
  converged <- FALSE
  old_loglik <- probit_ifa_loglik_lowrank(X, alpha, L)

  for (iter in seq_len(max_iter)) {
    iter_start <- Sys.time()
    e_step_start <- Sys.time()
    if (projection_update == "expectation") {
      W <- probit_ifa_truncated_mean(X, alpha, L)
    } else {
      W <- NULL
    }
    e_step_seconds <- as.numeric(difftime(Sys.time(), e_step_start, units = "secs"))
    projection_start <- Sys.time()
    projection <- if (projection_update == "expectation") {
      rank_H_centered_projection(W, H)
    } else if (projection_update == "sample_once") {
      stochastic_rank_H_centered_projection(
        X = X,
        alpha = alpha,
        L = L,
        H = H,
        n_draws = 1L,
        average = "sample_mean",
        seed = seed + 104729L * iter,
        parallel = FALSE,
        workers = workers
      )
    } else {
      stochastic_rank_H_centered_projection(
        X = X,
        alpha = alpha,
        L = L,
        H = H,
        n_draws = stochastic_svd_draws,
        average = stochastic_svd_average,
        seed = seed + 104729L * iter,
        parallel = parallel,
        workers = workers
      )
    }
    alpha_new <- projection$alpha
    L_new <- projection$L
    projection_seconds <- as.numeric(difftime(Sys.time(), projection_start, units = "secs"))
    objective_start <- Sys.time()
    loglik <- probit_ifa_loglik_lowrank(X, alpha_new, L_new)
    objective_seconds <- as.numeric(difftime(Sys.time(), objective_start, units = "secs"))

    rel_loglik <- abs(loglik - old_loglik) / (1 + abs(old_loglik))
    rel_L <- sqrt(sum((L_new - L)^2)) / (1 + sqrt(sum(L^2)))
    history[[iter]] <- data.frame(
      iteration = iter,
      probit_loglik = loglik,
      relative_loglik_change = rel_loglik,
      relative_L_change = rel_L,
      projection_update = projection_update,
      stochastic_svd_draws = if (projection_update == "stochastic_average") stochastic_svd_draws else NA_integer_,
      stochastic_svd_average = if (projection_update == "stochastic_average") stochastic_svd_average else NA_character_,
      e_step_seconds = e_step_seconds,
      projection_seconds = projection_seconds,
      objective_seconds = objective_seconds,
      iteration_seconds = as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
    )

    if (isTRUE(verbose)) {
      message(
        "  EM-SVD iter ", iter,
        ": probit loglik=", round(loglik, 3),
        ", projection=", projection_update,
        ", rel ll=", signif(rel_loglik, 3),
        ", rel L=", signif(rel_L, 3)
      )
    }

    alpha <- alpha_new
    L <- L_new
    old_loglik <- loglik

    L_done <- is.null(tol_L) || rel_L <= tol_L
    if (rel_loglik <= tol_loglik && L_done) {
      converged <- TRUE
      break
    }
  }

  history <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  list(
    alpha = alpha,
    L = L,
    probit_loglik = old_loglik,
    history = history,
    converged = converged,
    n_completed = nrow(history)
  )
}

random_centered_rank_H_matrix <- function(n, p, H, scale = 0.05) {
  A <- matrix(rnorm(n * H), n, H)
  B <- matrix(rnorm(p * H), p, H)
  L <- scale * A %*% t(B) / sqrt(H)
  sweep(L, 2L, colMeans(L), "-")
}

augmented_z_svd_lowrank_start <- function(
    X,
    H,
    alpha = NULL,
    seed = 1L,
    z_start = c("sample", "expectation", "stochastic_average"),
    stochastic_svd_draws = 5L,
    stochastic_svd_average = c("projected_signal", "sample_mean"),
    parallel = FALSE,
    workers = NULL) {
  # Viroli-inspired low-rank start.  Begin with itemwise probit intercepts,
  # build an augmented latent matrix under that intercept-only model, and use
  # its rank-H centered SVD as the initial probit signal L.  With
  # z_start="stochastic_average", average several sampled SVD projections before
  # returning the start.
  z_start <- match.arg(z_start)
  stochastic_svd_average <- match.arg(stochastic_svd_average)
  X <- as.matrix(X)
  if (is.null(alpha)) alpha <- initialize_binary_intercepts(X)
  alpha <- as.numeric(alpha)
  if (z_start == "sample") {
    Z0 <- initialize_binary_Z(X, seed = seed, alpha = alpha)
  } else {
    if (z_start == "stochastic_average") {
      projection <- stochastic_rank_H_centered_projection(
        X = X,
        alpha = alpha,
        L = matrix(0, nrow(X), ncol(X)),
        H = H,
        n_draws = stochastic_svd_draws,
        average = stochastic_svd_average,
        seed = seed,
        parallel = parallel,
        workers = workers
      )
      return(list(alpha = projection$alpha, L = projection$L))
    }
    Z0 <- probit_ifa_truncated_mean(
      X = X,
      alpha = alpha,
      L = matrix(0, nrow(X), ncol(X))
    )
  }
  projection <- rank_H_centered_projection(Z0, H)
  list(alpha = projection$alpha, L = projection$L)
}

fit_lowrank_probit_em_svd <- function(
    X,
    H,
    max_iter = 50L,
    tol_loglik = 1e-5,
    tol_L = NULL,
    init_method = c("intercept_only", "viroli_svd", "both"),
    init_z = c("sample", "expectation", "stochastic_average"),
    projection_update = c("expectation", "sample_once", "stochastic_average"),
    stochastic_svd_draws = 5L,
    stochastic_svd_average = c("projected_signal", "sample_mean"),
    n_random_starts = 0L,
    random_start_scale = 0.05,
    seed = 1L,
    parallel = FALSE,
    workers = NULL,
    verbose = FALSE) {
  X <- as.matrix(X)
  init_method <- match.arg(init_method)
  init_z <- match.arg(init_z)
  projection_update <- match.arg(projection_update)
  stochastic_svd_average <- match.arg(stochastic_svd_average)
  set.seed(seed)

  starts <- list()
  if (init_method %in% c("intercept_only", "both")) {
    starts[[length(starts) + 1L]] <- list(
      name = "intercept_only",
      alpha = NULL,
      L = matrix(0, nrow(X), ncol(X))
    )
  }
  if (init_method %in% c("viroli_svd", "both")) {
    augmented_start <- augmented_z_svd_lowrank_start(
      X = X,
      H = H,
      seed = seed + 7919L,
      z_start = init_z,
      stochastic_svd_draws = stochastic_svd_draws,
      stochastic_svd_average = stochastic_svd_average,
      parallel = parallel,
      workers = workers
    )
    starts[[length(starts) + 1L]] <- list(
      name = paste0("viroli_", init_z, "_svd"),
      alpha = augmented_start$alpha,
      L = augmented_start$L
    )
  }
  if (n_random_starts > 0L) {
    for (s in seq_len(n_random_starts)) {
      starts[[length(starts) + 1L]] <- list(
        name = paste0("random_", s),
        alpha = NULL,
        L = random_centered_rank_H_matrix(nrow(X), ncol(X), H, scale = random_start_scale)
      )
    }
  }

  fits <- lapply(seq_along(starts), function(s) {
    if (isTRUE(verbose)) message("EM-SVD start: ", starts[[s]]$name)
    fit <- fit_lowrank_probit_em_svd_one_start(
      X = X,
      H = H,
      max_iter = max_iter,
      tol_loglik = tol_loglik,
      tol_L = tol_L,
      projection_update = projection_update,
      stochastic_svd_draws = stochastic_svd_draws,
      stochastic_svd_average = stochastic_svd_average,
      alpha_init = starts[[s]]$alpha,
      L_init = starts[[s]]$L,
      seed = seed + 1009L * s,
      parallel = parallel,
      workers = workers,
      verbose = verbose
    )
    fit$start_name <- starts[[s]]$name
    fit
  })

  scores <- vapply(fits, function(z) z$probit_loglik, numeric(1))
  best <- which.max(scores)
  out <- fits[[best]]
  out$all_start_loglik <- scores
  names(out$all_start_loglik) <- vapply(starts, `[[`, character(1), "name")
  out$selected_start <- starts[[best]]$name
  out$init_method <- init_method
  out$init_z <- init_z
  out$projection_update <- projection_update
  out$stochastic_svd_draws <- stochastic_svd_draws
  out$stochastic_svd_average <- stochastic_svd_average
  out
}

spectral_scores_from_lowrank_signal <- function(L, H) {
  L <- as.matrix(L)
  n <- nrow(L)
  dec <- svd(L, nu = H, nv = H)
  U <- dec$u[, seq_len(H), drop = FALSE]
  V <- dec$v[, seq_len(H), drop = FALSE]
  dH <- dec$d[seq_len(H)]
  S <- sqrt(n) * U
  B <- V %*% diag(dH / sqrt(n), H, H)
  colnames(S) <- paste0("factor_", seq_len(H))
  colnames(B) <- paste0("factor_", seq_len(H))
  list(S = S, B = B, svd = dec)
}

ensure_riemannian_rotation_available <- function() {
  if (exists("estimate_mixture_ica_riemannian", mode = "function", inherits = TRUE)) {
    return(invisible(TRUE))
  }

  candidates <- unique(c(
    file.path(getwd(), "R", "riemannian_rotation.R"),
    file.path(dirname(getwd()), "R", "riemannian_rotation.R"),
    file.path(dirname(dirname(getwd())), "R", "riemannian_rotation.R")
  ))
  for (path in candidates) {
    if (file.exists(path)) {
      source(path)
      if (exists("estimate_mixture_ica_riemannian", mode = "function", inherits = TRUE)) {
        return(invisible(TRUE))
      }
    }
  }

  stop(
    "rotation_optimizer = 'riemannian' requires R/riemannian_rotation.R; ",
    "source that file or run from the repository root."
  )
}

rotate_em_svd_scores_with_mixtures <- function(
    S,
    G_fixed,
    loading_basis = NULL,
    rotation_loading_l1_penalty = 0,
    n_random_starts = 3L,
    n_ica_starts = 0L,
    ica_functions = c("logcosh", "exp"),
    ica_max_iter = 200L,
    ica_tol = 1e-4,
    max_outer = 5L,
    n_mix_starts = 5L,
    mixture_max_iter = 50L,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 1.5,
    weight_prior_alpha = 1,
    grid_size = 31L,
    rotation_objective_tolerance = 1e-4,
    rotation_min_outer = 2L,
    require_mixture_convergence_for_rotation_stop = TRUE,
    rotation_optimizer = c("riemannian", "givens"),
    rotation_sweep = c("full", "hybrid", "multi_disjoint", "disjoint", "promising"),
    riemannian_rotation_steps = 10L,
    riemannian_eta0 = 1,
    riemannian_beta = 0.5,
    riemannian_min_eta = 1e-8,
    riemannian_grad_tol = 1e-6,
    riemannian_update = c("cayley", "expm"),
    seed = 1L,
    parallel = FALSE,
    workers = NULL,
    verbose = FALSE) {
  mixture_update <- match.arg(mixture_update)
  rotation_optimizer <- match.arg(rotation_optimizer)
  rotation_sweep <- match.arg(rotation_sweep)
  riemannian_update <- match.arg(riemannian_update)

  if (rotation_optimizer == "riemannian") {
    ensure_riemannian_rotation_available()
    return(estimate_mixture_ica_riemannian(
      S = S,
      G_fixed = G_fixed,
      rotation_loading_basis = loading_basis,
      rotation_loading_l1_penalty = rotation_loading_l1_penalty,
      n_random_starts = n_random_starts,
      n_ica_starts = n_ica_starts,
      ica_functions = ica_functions,
      ica_max_iter = ica_max_iter,
      ica_tol = ica_tol,
      outer_maxit = max_outer,
      outer_min_iter = rotation_min_outer,
      rotation_steps = riemannian_rotation_steps,
      eta0 = riemannian_eta0,
      beta = riemannian_beta,
      min_eta = riemannian_min_eta,
      grad_tol = riemannian_grad_tol,
      update = riemannian_update,
      n_mix_starts = n_mix_starts,
      mixture_max_iter = mixture_max_iter,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      rel_tol = rotation_objective_tolerance,
      require_mixture_convergence_for_rotation_stop = require_mixture_convergence_for_rotation_stop,
      seed = seed,
      parallel = parallel,
      workers = workers,
      verbose = verbose
    ))
  }

  estimate_mixture_ica_unknown_G(
    S = S,
    G_fixed = G_fixed,
    rotation_loading_basis = loading_basis,
    rotation_loading_l1_penalty = rotation_loading_l1_penalty,
    n_random_starts = n_random_starts,
    n_ica_starts = n_ica_starts,
    ica_functions = ica_functions,
    ica_max_iter = ica_max_iter,
    ica_tol = ica_tol,
    max_outer = max_outer,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha,
    grid_size = grid_size,
    rotation_objective_tolerance = rotation_objective_tolerance,
    rotation_min_outer = rotation_min_outer,
    require_mixture_convergence_for_rotation_stop = require_mixture_convergence_for_rotation_stop,
    rotation_sweep = rotation_sweep,
    seed = seed,
    parallel = parallel,
    workers = workers,
    verbose = verbose
  )
}

pretrain_probit_ifa_em_svd <- function(
    X,
    H,
    G_fixed,
    em_max_iter = 50L,
    em_tol_loglik = 1e-5,
    em_tol_L = NULL,
    em_init_method = c("intercept_only", "viroli_svd", "both"),
    em_init_z = c("sample", "expectation", "stochastic_average"),
    em_projection_update = c("expectation", "sample_once", "stochastic_average"),
    stochastic_svd_draws = 5L,
    stochastic_svd_average = c("projected_signal", "sample_mean"),
    em_random_starts = 0L,
    em_random_start_scale = 0.05,
    rotation_loading_l1_penalty = 0,
    rotation_random_starts = 3L,
    rotation_ica_starts = 0L,
    rotation_ica_functions = c("logcosh", "exp"),
    rotation_ica_max_iter = 200L,
    rotation_ica_tol = 1e-4,
    rotation_max_outer = 5L,
    n_mix_starts = 5L,
    mixture_max_iter = 50L,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 1.5,
    weight_prior_alpha = 1,
    grid_size = 31L,
    rotation_objective_tolerance = 1e-4,
    rotation_min_outer = 2L,
    require_mixture_convergence_for_rotation_stop = TRUE,
    rotation_optimizer = c("riemannian", "givens"),
    rotation_sweep = c("full", "hybrid", "multi_disjoint", "disjoint", "promising"),
    riemannian_rotation_steps = 10L,
    riemannian_eta0 = 1,
    riemannian_beta = 0.5,
    riemannian_min_eta = 1e-8,
    riemannian_grad_tol = 1e-6,
    riemannian_update = c("cayley", "expm"),
    loading_penalty = 0,
    seed = 1L,
    parallel = FALSE,
    workers = NULL,
    verbose = FALSE) {
  X <- as.matrix(X)
  H <- as.integer(H)
  G_fixed <- normalize_G_fixed(G_fixed, H)
  mixture_update <- match.arg(mixture_update)
  rotation_optimizer <- match.arg(rotation_optimizer)
  rotation_sweep <- match.arg(rotation_sweep)
  riemannian_update <- match.arg(riemannian_update)
  em_init_method <- match.arg(em_init_method)
  em_init_z <- match.arg(em_init_z)
  em_projection_update <- match.arg(em_projection_update)
  stochastic_svd_average <- match.arg(stochastic_svd_average)

  lowrank <- fit_lowrank_probit_em_svd(
    X = X,
    H = H,
    max_iter = em_max_iter,
    tol_loglik = em_tol_loglik,
    tol_L = em_tol_L,
    init_method = em_init_method,
    init_z = em_init_z,
    projection_update = em_projection_update,
    stochastic_svd_draws = stochastic_svd_draws,
    stochastic_svd_average = stochastic_svd_average,
    n_random_starts = em_random_starts,
    random_start_scale = em_random_start_scale,
    seed = seed,
    parallel = parallel,
    workers = workers,
    verbose = verbose
  )

  spectral <- spectral_scores_from_lowrank_signal(lowrank$L, H)
  rotation <- rotate_em_svd_scores_with_mixtures(
    S = spectral$S,
    G_fixed = G_fixed,
    loading_basis = spectral$B,
    rotation_loading_l1_penalty = rotation_loading_l1_penalty,
    n_random_starts = rotation_random_starts,
    n_ica_starts = rotation_ica_starts,
    ica_functions = rotation_ica_functions,
    ica_max_iter = rotation_ica_max_iter,
    ica_tol = rotation_ica_tol,
    max_outer = rotation_max_outer,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha,
    grid_size = grid_size,
    rotation_objective_tolerance = rotation_objective_tolerance,
    rotation_min_outer = rotation_min_outer,
    require_mixture_convergence_for_rotation_stop = require_mixture_convergence_for_rotation_stop,
    rotation_optimizer = rotation_optimizer,
    rotation_sweep = rotation_sweep,
    riemannian_rotation_steps = riemannian_rotation_steps,
    riemannian_eta0 = riemannian_eta0,
    riemannian_beta = riemannian_beta,
    riemannian_min_eta = riemannian_min_eta,
    riemannian_grad_tol = riemannian_grad_tol,
    riemannian_update = riemannian_update,
    seed = seed + 10000L,
    parallel = parallel,
    workers = workers,
    verbose = verbose
  )

  F_hat <- rotation$F_hat
  Lambda_hat <- spectral$B %*% rotation$R
  if (loading_penalty > 0) {
    Lambda_hat <- soft_threshold(Lambda_hat, loading_penalty)
  }
  rownames(Lambda_hat) <- colnames(X)
  colnames(Lambda_hat) <- paste0("factor_", seq_len(H))

  responsibilities <- lapply(seq_len(H), function(h) {
    mixture_responsibilities(F_hat[, h], rotation$fits[[h]])
  })
  class_map <- sapply(responsibilities, max.col, ties.method = "first")
  if (H == 1L) class_map <- matrix(class_map, ncol = 1L)
  colnames(class_map) <- paste0("factor_", seq_len(H))

  fitted <- sweep(F_hat %*% t(Lambda_hat), 2L, lowrank$alpha, "+")
  residual <- probit_ifa_truncated_mean(X, lowrank$alpha, lowrank$L) - fitted

  history <- data.frame(
    stage = c("lowrank_em_svd", "mixture_rotation"),
    iteration = c(lowrank$n_completed, rotation$rotation_completed_outer),
    probit_loglik = c(lowrank$probit_loglik, lowrank$probit_loglik),
    mixture_loglik = c(NA_real_, rotation$loglik),
    objective = c(lowrank$probit_loglik, lowrank$probit_loglik + rotation$loglik)
  )

  list(
    model = "probit_ifa_em_svd_spectral_mixture_pretraining",
    X = X,
    H = H,
    H_selection = NULL,
    H_selection_strategy = "fixed_supplied_H",
    S = spectral$S,
    svd_fit = spectral,
    lowrank_fit = lowrank,
    alpha_hat = lowrank$alpha,
    L_hat = lowrank$L,
    R = rotation$R,
    F_hat = F_hat,
    Lambda_hat = Lambda_hat,
    Lambda_ls = Lambda_hat,
    loading_update = "em_svd_lowrank_signal",
    loading_penalty = loading_penalty,
    mixture_fits = rotation$fits,
    G_hat = rotation$G_hat,
    G_fixed = G_fixed,
    pretrain_G_selection = "fixed",
    mixture_update = mixture_update,
    rotation_optimizer = rotation_optimizer,
    rotation_fit = rotation,
    responsibilities = responsibilities,
    class_map = class_map,
    profile_id = binary_profile_id(class_map),
    fitted = fitted,
    residual = residual,
    Psi_hat = diag(1, ncol(X)),
    Psi_hat_unconstrained = diag(pmax(colMeans(residual^2), 1e-4), ncol(X)),
    history = history,
    em_history = lowrank$history,
    pretraining_converged = isTRUE(lowrank$converged),
    pretraining_completed_iter = lowrank$n_completed,
    selected_pretraining_iteration = lowrank$n_completed,
    em_init_method = lowrank$init_method,
    em_init_z = lowrank$init_z,
    selected_em_start = lowrank$selected_start,
    all_em_start_loglik = lowrank$all_start_loglik,
    em_projection_update = lowrank$projection_update,
    stochastic_svd_draws = lowrank$stochastic_svd_draws,
    stochastic_svd_average = lowrank$stochastic_svd_average,
    rotation_completed_outer = rotation$rotation_completed_outer,
    rotation_converged = isTRUE(rotation$rotation_converged),
    rotation_history = rotation$rotation_history,
    rotation_step_history = rotation$rotation_step_history,
    z_update = "none_em_svd",
    fix_psi_identity = TRUE,
    estimate_intercept = TRUE,
    call = match.call()
  )
}

fit_binary_probit_em_svd_pretrain_then_refine <- function(
    X,
    H,
    G_fixed,
    n_refine_iter = 5L,
    n_mix_starts = 3L,
    mixture_max_iter = 20L,
    em_init_method = c("intercept_only", "viroli_svd", "both"),
    em_init_z = c("sample", "expectation", "stochastic_average"),
    em_projection_update = c("expectation", "sample_once", "stochastic_average"),
    stochastic_svd_draws = 5L,
    stochastic_svd_average = c("projected_signal", "sample_mean"),
    mixture_update = c("map", "mle"),
    mixture_refit = c("em", "fixed_responsibility_mstep", "posterior_moment_mstep"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 1.5,
    weight_prior_alpha = 1,
    refine_mu_prior_mean = NULL,
    refine_mu_prior_kappa = NULL,
    refine_var_prior_shape = NULL,
    refine_var_prior_scale = NULL,
    refine_weight_prior_alpha = NULL,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    min_mixture_var = 1e-3,
    lambda_l1_penalty = 0,
    lasso_backend = c("proximal", "glmnet"),
    glmnet_standardize = FALSE,
    objective_tolerance = 1e-5,
    min_refine_iter = 1L,
    enforce_monotone_refinement = TRUE,
    monotone_tolerance = 1e-8,
    return_best_refinement_iteration = FALSE,
    refinement_selection_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    require_mixture_convergence_for_stop = FALSE,
    seed = 1L,
    parallel = FALSE,
    workers = NULL,
    verbose = FALSE,
    ...) {
  mixture_update <- match.arg(mixture_update)
  mixture_refit <- match.arg(mixture_refit)
  factor_update <- match.arg(factor_update)
  refinement_selection_objective <- match.arg(refinement_selection_objective)
  lasso_backend <- match.arg(lasso_backend)
  em_init_method <- match.arg(em_init_method)
  em_init_z <- match.arg(em_init_z)
  em_projection_update <- match.arg(em_projection_update)
  stochastic_svd_average <- match.arg(stochastic_svd_average)
  if (is.null(refine_mu_prior_mean)) refine_mu_prior_mean <- mu_prior_mean
  if (is.null(refine_mu_prior_kappa)) refine_mu_prior_kappa <- mu_prior_kappa
  if (is.null(refine_var_prior_shape)) refine_var_prior_shape <- var_prior_shape
  if (is.null(refine_var_prior_scale)) refine_var_prior_scale <- var_prior_scale
  if (is.null(refine_weight_prior_alpha)) refine_weight_prior_alpha <- weight_prior_alpha

  pretrain_fit <- pretrain_probit_ifa_em_svd(
    X = X,
    H = H,
    G_fixed = G_fixed,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    em_init_method = em_init_method,
    em_init_z = em_init_z,
    em_projection_update = em_projection_update,
    stochastic_svd_draws = stochastic_svd_draws,
    stochastic_svd_average = stochastic_svd_average,
    mixture_update = mixture_update,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha,
    seed = seed,
    parallel = parallel,
    workers = workers,
    verbose = verbose,
    ...
  )

  refine_fit <- fit_binary_probit_refinement(
    X = X,
    pretrain_fit = pretrain_fit,
    n_refine_iter = n_refine_iter,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mixture_refit = mixture_refit,
    mu_prior_mean = refine_mu_prior_mean,
    mu_prior_kappa = refine_mu_prior_kappa,
    var_prior_shape = refine_var_prior_shape,
    var_prior_scale = refine_var_prior_scale,
    weight_prior_alpha = refine_weight_prior_alpha,
    factor_update = factor_update,
    min_mixture_var = min_mixture_var,
    lambda_l1_penalty = lambda_l1_penalty,
    lasso_backend = lasso_backend,
    glmnet_standardize = glmnet_standardize,
    objective_tolerance = objective_tolerance,
    min_refine_iter = min_refine_iter,
    stopping_objective = "posterior_objective",
    enforce_monotone_refinement = enforce_monotone_refinement,
    monotone_tolerance = monotone_tolerance,
    return_best_refinement_iteration = return_best_refinement_iteration,
    refinement_selection_objective = refinement_selection_objective,
    require_mixture_convergence_for_stop = require_mixture_convergence_for_stop,
    parallel = parallel,
    workers = workers,
    verbose = verbose,
    ...
  )

  list(pretrain_fit = pretrain_fit, refine_fit = refine_fit)
}
