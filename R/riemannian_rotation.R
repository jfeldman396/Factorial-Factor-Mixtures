#!/usr/bin/env Rscript

# Riemannian orthogonal rotation for the spectral-mixture pretraining step.
#
# This is an experimental alternative to cyclic pairwise Givens sweeps.  It
# keeps the same profiled product-of-marginals objective, but updates all
# rotation directions simultaneously on the orthogonal group.
#
# Expected prerequisites, sourced from R/binary_probit_pretraining.R:
#   fit_column_mixtures_fixed_G()
#   mixture_responsibilities()
#   mixture_loglik_total()
#   rotation_selection_score()
#   random_orthogonal()
#   project_to_orthogonal()
#   fastica_rotation_starts()
#   normalize_G_fixed()
#   parallel_lapply()
#   resolve_workers()

skew_part <- function(A) {
  0.5 * (A - t(A))
}

riemannian_mixture_score_matrix <- function(F, fits, sigma_min = 1e-6) {
  # Compute psi_ih = d log q_h(F_ih) / d F_ih for every rotated score entry.
  F <- as.matrix(F)
  n <- nrow(F)
  H <- ncol(F)
  Psi <- matrix(0, n, H)

  for (h in seq_len(H)) {
    fit_h <- fits[[h]]
    resp <- mixture_responsibilities(F[, h], fit_h)
    var_h <- pmax(as.numeric(fit_h$var), sigma_min^2)
    centered <- sweep(matrix(F[, h], nrow = n, ncol = length(var_h)), 2L, fit_h$mu, "-")
    Psi[, h] <- -rowSums(sweep(resp * centered, 2L, var_h, "/"))
  }

  Psi
}

riemannian_rotation_objective <- function(
    S,
    R,
    fits,
    loading_basis = NULL,
    loading_l1_penalty = 0) {
  F <- as.matrix(S) %*% as.matrix(R)
  rotation_selection_score(
    F = F,
    fits = fits,
    R = R,
    loading_basis = loading_basis,
    loading_l1_penalty = loading_l1_penalty
  )
}

riemannian_rotation_skew_gradient <- function(
    S,
    R,
    fits,
    loading_basis = NULL,
    loading_l1_penalty = 0,
    sigma_min = 1e-6) {
  # Euclidean gradient for the mixture part is S' Psi.  If requested, include
  # the subgradient of -lambda ||B R||_1, where B is the loading basis.
  S <- as.matrix(S)
  R <- as.matrix(R)
  F <- S %*% R
  Psi <- riemannian_mixture_score_matrix(F, fits, sigma_min = sigma_min)
  G <- crossprod(S, Psi)

  if (!is.null(loading_basis) && loading_l1_penalty > 0) {
    Lambda_current <- as.matrix(loading_basis) %*% R
    G <- G - loading_l1_penalty * crossprod(as.matrix(loading_basis), sign(Lambda_current))
  }

  skew_part(crossprod(R, G))
}

cayley_rotation_update <- function(R, A, eta) {
  H <- ncol(R)
  I <- diag(H)
  R %*% ((I + 0.5 * eta * A) %*% solve(I - 0.5 * eta * A))
}

expm_rotation_update <- function(R, A, eta) {
  if (!requireNamespace("expm", quietly = TRUE)) {
    return(cayley_rotation_update(R, A, eta))
  }
  R %*% expm::expm(eta * A)
}

riemannian_update_rotation <- function(R, A, eta, update = c("cayley", "expm")) {
  update <- match.arg(update)
  out <- if (update == "expm") expm_rotation_update(R, A, eta) else cayley_rotation_update(R, A, eta)
  project_to_orthogonal(out)
}

riemannian_rotation_block <- function(
    S,
    R,
    fits,
    n_steps = 10L,
    eta0 = 1,
    beta = 0.5,
    min_eta = 1e-8,
    grad_tol = 1e-6,
    update = c("cayley", "expm"),
    loading_basis = NULL,
    loading_l1_penalty = 0,
    monotone_tolerance = 1e-10,
    verbose = FALSE) {
  # Holding mixture parameters fixed, do several monotone manifold-gradient
  # ascent steps.  Mixtures are intentionally not refit during line search.
  update <- match.arg(update)
  S <- as.matrix(S)
  R <- as.matrix(R)
  n_steps <- max(1L, as.integer(n_steps))
  eta0 <- as.numeric(eta0)
  beta <- as.numeric(beta)
  if (!is.finite(eta0) || eta0 <= 0) stop("eta0 must be positive.")
  if (!is.finite(beta) || beta <= 0 || beta >= 1) stop("beta must be in (0, 1).")

  current_objective <- riemannian_rotation_objective(
    S, R, fits,
    loading_basis = loading_basis,
    loading_l1_penalty = loading_l1_penalty
  )
  history <- vector("list", n_steps)
  converged <- FALSE

  for (step in seq_len(n_steps)) {
    step_start <- Sys.time()
    A <- riemannian_rotation_skew_gradient(
      S = S,
      R = R,
      fits = fits,
      loading_basis = loading_basis,
      loading_l1_penalty = loading_l1_penalty
    )
    grad_norm <- sqrt(sum(A^2))
    accepted <- FALSE
    eta <- eta0
    candidate_objective <- current_objective
    candidate_R <- R
    n_backtrack <- 0L

    if (is.finite(grad_norm) && grad_norm > grad_tol) {
      while (eta >= min_eta) {
        trial_R <- riemannian_update_rotation(R, A, eta, update = update)
        trial_objective <- riemannian_rotation_objective(
          S, trial_R, fits,
          loading_basis = loading_basis,
          loading_l1_penalty = loading_l1_penalty
        )
        if (is.finite(trial_objective) &&
            trial_objective >= current_objective - monotone_tolerance) {
          accepted <- TRUE
          candidate_R <- trial_R
          candidate_objective <- trial_objective
          break
        }
        eta <- eta * beta
        n_backtrack <- n_backtrack + 1L
      }
    }

    history[[step]] <- data.frame(
      block_step = step,
      fixed_mixture_objective_before = current_objective,
      fixed_mixture_objective_after = candidate_objective,
      objective_improvement = candidate_objective - current_objective,
      grad_norm = grad_norm,
      eta = if (accepted) eta else NA_real_,
      n_backtrack = n_backtrack,
      accepted = accepted,
      orthogonality_error = sqrt(sum((crossprod(candidate_R) - diag(ncol(candidate_R)))^2)),
      step_seconds = as.numeric(difftime(Sys.time(), step_start, units = "secs")),
      stringsAsFactors = FALSE
    )

    if (isTRUE(verbose)) {
      message(
        "    Riemannian step ", step,
        ": objective=", round(candidate_objective, 3),
        ", gain=", signif(candidate_objective - current_objective, 3),
        ", ||A||=", signif(grad_norm, 3),
        ", accepted=", accepted
      )
    }

    if (!accepted || !is.finite(grad_norm) || grad_norm <= grad_tol) {
      converged <- is.finite(grad_norm) && grad_norm <= grad_tol
      break
    }

    R <- candidate_R
    current_objective <- candidate_objective
  }

  history <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  list(
    R = project_to_orthogonal(R),
    F = S %*% R,
    fixed_mixture_objective = current_objective,
    gradient_norm = if (nrow(history)) tail(history$grad_norm, 1L) else NA_real_,
    converged = converged,
    history = history
  )
}

estimate_mixture_ica_riemannian <- function(
    S,
    G_fixed,
    n_random_starts = 1L,
    n_ica_starts = 0L,
    ica_functions = c("logcosh", "exp"),
    ica_max_iter = 200L,
    ica_tol = 1e-4,
    outer_maxit = 20L,
    outer_min_iter = 2L,
    rotation_steps = 10L,
    eta0 = 1,
    beta = 0.5,
    min_eta = 1e-8,
    rel_tol = 1e-4,
    grad_tol = 1e-6,
    update = c("cayley", "expm"),
    n_mix_starts = 3L,
    mixture_max_iter = 50L,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    require_mixture_convergence_for_rotation_stop = TRUE,
    rotation_loading_basis = NULL,
    rotation_loading_l1_penalty = 0,
    seed = 1L,
    parallel = FALSE,
    workers = NULL,
    verbose = FALSE) {
  update <- match.arg(update)
  mixture_update <- match.arg(mixture_update)
  workers <- resolve_workers(workers)
  set.seed(seed)
  S <- as.matrix(S)
  H <- ncol(S)
  G_fixed <- normalize_G_fixed(G_fixed, H)

  if (!is.null(rotation_loading_basis)) {
    rotation_loading_basis <- as.matrix(rotation_loading_basis)
    if (ncol(rotation_loading_basis) != H) {
      stop("rotation_loading_basis must have one column per factor.")
    }
  }
  rotation_loading_l1_penalty <- max(0, as.numeric(rotation_loading_l1_penalty))

  starts <- list(identity = diag(H))
  ica_starts <- fastica_rotation_starts(
    S = S,
    n_starts = n_ica_starts,
    ica_functions = ica_functions,
    max_iter = ica_max_iter,
    tol = ica_tol,
    seed = seed + 4242L,
    verbose = verbose
  )
  if (length(ica_starts)) starts <- c(starts, ica_starts)
  for (s in seq_len(max(0L, as.integer(n_random_starts)))) {
    starts[[paste0("random_", s)]] <- random_orthogonal(H)
  }

  fit_one_start <- function(start_index, use_parallel_mixtures = FALSE) {
    set.seed(seed + 100L * start_index)
    start_name <- names(starts)[start_index]
    if (isTRUE(verbose)) message("  Riemannian rotation start: ", start_name)
    R <- starts[[start_name]]
    F <- S %*% R

    fits <- fit_column_mixtures_fixed_G(
      F = F,
      G_fixed = G_fixed,
      n_starts = n_mix_starts,
      max_iter = mixture_max_iter,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      parallel = use_parallel_mixtures,
      workers = workers
    )

    previous_profiled_score <- rotation_selection_score(
      F = F,
      fits = fits,
      R = R,
      loading_basis = rotation_loading_basis,
      loading_l1_penalty = rotation_loading_l1_penalty
    )
    previous_mixture_loglik <- mixture_loglik_total(F, fits)
    previous_grad_norm <- sqrt(sum(riemannian_rotation_skew_gradient(
      S, R, fits,
      loading_basis = rotation_loading_basis,
      loading_l1_penalty = rotation_loading_l1_penalty
    )^2))

    rotation_history <- vector("list", outer_maxit + 1L)
    rotation_history[[1L]] <- data.frame(
      outer_iteration = 0L,
      mixture_loglik = previous_mixture_loglik,
      rotation_selection_score = previous_profiled_score,
      loading_l1 = rotation_loading_l1(rotation_loading_basis, R),
      gradient_norm = previous_grad_norm,
      fixed_score_before_refit = NA_real_,
      relative_mixture_loglik_improvement = NA_real_,
      relative_rotation_score_improvement = NA_real_,
      all_mixtures_converged = all(vapply(fits, function(z) isTRUE(z$converged), logical(1))),
      accepted_rotation_steps = 0L,
      rotation_update_seconds = NA_real_,
      mixture_refit_seconds = NA_real_,
      objective_seconds = NA_real_,
      iteration_seconds = NA_real_,
      stringsAsFactors = FALSE
    )
    step_histories <- vector("list", outer_maxit)
    rotation_converged <- FALSE
    rotation_completed_outer <- 0L

    for (outer in seq_len(max(1L, as.integer(outer_maxit)))) {
      iter_start <- Sys.time()
      block_start <- Sys.time()
      block <- riemannian_rotation_block(
        S = S,
        R = R,
        fits = fits,
        n_steps = rotation_steps,
        eta0 = eta0,
        beta = beta,
        min_eta = min_eta,
        grad_tol = grad_tol,
        update = update,
        loading_basis = rotation_loading_basis,
        loading_l1_penalty = rotation_loading_l1_penalty,
        verbose = verbose
      )
      rotation_update_seconds <- as.numeric(difftime(Sys.time(), block_start, units = "secs"))
      R <- block$R
      F <- block$F
      step_hist <- block$history
      if (!is.null(step_hist) && nrow(step_hist)) {
        step_hist$outer_iteration <- outer
        step_hist$start_name <- start_name
      }
      step_histories[[outer]] <- step_hist

      fixed_score_before_refit <- block$fixed_mixture_objective
      mixture_refit_start <- Sys.time()
      fits <- fit_column_mixtures_fixed_G(
        F = F,
        G_fixed = G_fixed,
        n_starts = n_mix_starts,
        max_iter = mixture_max_iter,
        mixture_update = mixture_update,
        mu_prior_mean = mu_prior_mean,
        mu_prior_kappa = mu_prior_kappa,
        var_prior_shape = var_prior_shape,
        var_prior_scale = var_prior_scale,
        weight_prior_alpha = weight_prior_alpha,
        previous_fits = fits,
        parallel = use_parallel_mixtures,
        workers = workers
      )
      mixture_refit_seconds <- as.numeric(difftime(Sys.time(), mixture_refit_start, units = "secs"))

      objective_start <- Sys.time()
      current_mixture_loglik <- mixture_loglik_total(F, fits)
      current_profiled_score <- rotation_selection_score(
        F = F,
        fits = fits,
        R = R,
        loading_basis = rotation_loading_basis,
        loading_l1_penalty = rotation_loading_l1_penalty
      )
      current_grad_norm <- sqrt(sum(riemannian_rotation_skew_gradient(
        S, R, fits,
        loading_basis = rotation_loading_basis,
        loading_l1_penalty = rotation_loading_l1_penalty
      )^2))
      objective_seconds <- as.numeric(difftime(Sys.time(), objective_start, units = "secs"))

      relative_score_improvement <- (current_profiled_score - previous_profiled_score) /
        (1 + abs(previous_profiled_score))
      relative_mixture_improvement <- (current_mixture_loglik - previous_mixture_loglik) /
        (1 + abs(previous_mixture_loglik))
      all_mixtures_converged <- all(vapply(fits, function(z) isTRUE(z$converged), logical(1)))
      accepted_steps <- if (!is.null(step_hist) && nrow(step_hist)) sum(step_hist$accepted %in% TRUE) else 0L
      iteration_seconds <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

      rotation_history[[outer + 1L]] <- data.frame(
        outer_iteration = outer,
        mixture_loglik = current_mixture_loglik,
        rotation_selection_score = current_profiled_score,
        loading_l1 = rotation_loading_l1(rotation_loading_basis, R),
        gradient_norm = current_grad_norm,
        fixed_score_before_refit = fixed_score_before_refit,
        relative_mixture_loglik_improvement = relative_mixture_improvement,
        relative_rotation_score_improvement = relative_score_improvement,
        all_mixtures_converged = all_mixtures_converged,
        accepted_rotation_steps = accepted_steps,
        rotation_update_seconds = rotation_update_seconds,
        mixture_refit_seconds = mixture_refit_seconds,
        objective_seconds = objective_seconds,
        iteration_seconds = iteration_seconds,
        stringsAsFactors = FALSE
      )
      rotation_completed_outer <- outer

      if (isTRUE(verbose)) {
        message(
          "  Riemannian outer ", outer,
          ": score=", round(current_profiled_score, 3),
          ", rel gain=", signif(relative_score_improvement, 3),
          ", ||A||=", signif(current_grad_norm, 3),
          ", accepted steps=", accepted_steps
        )
      }

      if (outer >= outer_min_iter &&
          is.finite(relative_score_improvement) &&
          relative_score_improvement <= rel_tol &&
          (!isTRUE(require_mixture_convergence_for_rotation_stop) || all_mixtures_converged)) {
        rotation_converged <- TRUE
        break
      }
      if (is.finite(current_grad_norm) && current_grad_norm <= grad_tol &&
          (!isTRUE(require_mixture_convergence_for_rotation_stop) || all_mixtures_converged)) {
        rotation_converged <- TRUE
        break
      }

      previous_profiled_score <- current_profiled_score
      previous_mixture_loglik <- current_mixture_loglik
    }

    selection_score <- rotation_selection_score(
      F = F,
      fits = fits,
      R = R,
      loading_basis = rotation_loading_basis,
      loading_l1_penalty = rotation_loading_l1_penalty
    )
    loglik <- mixture_loglik_total(F, fits)

    list(
      F_hat = F,
      R = project_to_orthogonal(R),
      fits = fits,
      G_hat = vapply(fits, function(z) length(z$pi), integer(1)),
      loglik = loglik,
      selection_score = selection_score,
      loading_l1 = rotation_loading_l1(rotation_loading_basis, R),
      rotation_history = do.call(rbind, rotation_history[!vapply(rotation_history, is.null, logical(1))]),
      rotation_step_history = do.call(rbind, step_histories[!vapply(step_histories, is.null, logical(1))]),
      rotation_converged = rotation_converged,
      rotation_completed_outer = rotation_completed_outer,
      rotation_objective_tolerance = rel_tol,
      rotation_min_outer = outer_min_iter,
      require_mixture_convergence_for_rotation_stop = isTRUE(require_mixture_convergence_for_rotation_stop),
      G_selection = "fixed",
      mixture_update = mixture_update,
      rotation_sweep = "riemannian",
      riemannian_rotation_steps = rotation_steps,
      riemannian_eta0 = eta0,
      riemannian_beta = beta,
      riemannian_grad_tol = grad_tol,
      riemannian_update = update,
      rotation_loading_l1_penalty = rotation_loading_l1_penalty
    )
  }

  use_start_parallel <- isTRUE(parallel) && workers > 1L && length(starts) > 1L
  results <- parallel_lapply(
    seq_along(starts),
    function(s) fit_one_start(s, use_parallel_mixtures = !use_start_parallel && isTRUE(parallel)),
    parallel = use_start_parallel,
    workers = workers
  )
  names(results) <- names(starts)

  scores <- vapply(results, function(z) z$selection_score, numeric(1))
  best <- which.max(scores)
  out <- results[[best]]
  out$start_name <- names(results)[best]
  out$all_start_scores <- scores
  out$all_start_mixture_loglik <- vapply(results, function(z) z$loglik, numeric(1))
  out$all_start_loading_l1 <- vapply(results, function(z) z$loading_l1, numeric(1))
  out$all_start_loglik <- out$all_start_mixture_loglik
  out
}
