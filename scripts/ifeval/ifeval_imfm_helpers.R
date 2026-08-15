#!/usr/bin/env Rscript

# Shared IFEval helpers for the binary-probit independent-mixture factor model.
#
# Model:
#   X_ij | f_i, alpha_j, lambda_j ~ Bernoulli(Phi(alpha_j + lambda_j' f_i))
#   f_ih independently follows a one-dimensional Gaussian mixture.
#
# This file intentionally contains reusable fitting code only.  Analysis
# entry-points live in:
#   - cv_ifeval_rank_lambda_models.R
#   - fit_interpret_ifeval_mixture.R

script_dir <- local({
  frames <- sys.frames()
  files <- vapply(frames, function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  files <- files[!is.na(files)]
  if (length(files) > 0L) dirname(normalizePath(tail(files, 1L))) else getwd()
})

repo_root <- normalizePath(file.path(script_dir, "../.."))
source(file.path(repo_root, "R", "binary_probit_pretraining.R"))
source(file.path(repo_root, "R", "binary_probit_refinement.R"))

# ----------------------------------------------------------------------------
# Small utilities
# ----------------------------------------------------------------------------

clip01 <- function(p, n) {
  eps <- 0.5 / max(1L, n)
  pmin(pmax(p, eps), 1 - eps)
}

ordered_component_labels <- function(F_hat, mixture_fits) {
  # Convert arbitrary mixture labels into ordered low/mid/high labels by
  # sorting each factor's component means.
  labs <- matrix(NA_integer_, nrow(F_hat), ncol(F_hat))
  resp_list <- vector("list", ncol(F_hat))

  for (h in seq_len(ncol(F_hat))) {
    resp_h <- mixture_responsibilities(F_hat[, h], mixture_fits[[h]])
    raw <- max.col(resp_h, ties.method = "first")
    ord <- order(mixture_fits[[h]]$mu)
    map <- integer(length(ord))
    map[ord] <- seq_along(ord)
    labs[, h] <- map[raw]
    resp_list[[h]] <- resp_h[, ord, drop = FALSE]
  }

  colnames(labs) <- paste0("factor_", seq_len(ncol(F_hat)))
  list(class_map = labs, responsibilities = resp_list)
}

profile_id_from_class_map <- function(class_map) {
  apply(class_map, 1L, paste, collapse = "-")
}

orient_factors_by_accuracy <- function(fit) {
  # Factor signs are unidentified. For reporting, orient each factor so larger
  # values are positively associated with overall row accuracy.
  acc <- rowMeans(fit$X)

  for (h in seq_len(ncol(fit$F_hat))) {
    rho <- suppressWarnings(cor(fit$F_hat[, h], acc))
    if (is.finite(rho) && rho < 0) {
      fit$F_hat[, h] <- -fit$F_hat[, h]
      fit$Lambda_hat[, h] <- -fit$Lambda_hat[, h]
      fit$mixture_fits[[h]]$mu <- -fit$mixture_fits[[h]]$mu
    }
  }

  ord <- ordered_component_labels(fit$F_hat, fit$mixture_fits)
  fit$class_map <- ord$class_map
  fit$responsibilities <- ord$responsibilities
  fit$profile_id <- profile_id_from_class_map(fit$class_map)
  fit
}

# ----------------------------------------------------------------------------
# Intercept-aware binary-probit likelihood and Z augmentation
# ----------------------------------------------------------------------------

binary_probit_loglik_alpha <- function(X, F_hat, Lambda, alpha) {
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  alpha <- as.numeric(alpha)

  eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  sum(X * log(p1) + (1 - X) * log(p0))
}

sample_binary_Z_given_model_alpha <- function(X, F_hat, Lambda, alpha, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  mean_mat <- sweep(F_hat %*% t(Lambda), 2L, as.numeric(alpha), "+")
  Z <- matrix(NA_real_, nrow(X), ncol(X), dimnames = dimnames(X))

  for (j in seq_len(ncol(X))) {
    lower <- ifelse(X[, j] == 1, 0, -Inf)
    upper <- ifelse(X[, j] == 1, Inf, 0)
    Z[, j] <- rtruncnorm_binary_vec(mean_mat[, j], 1, lower, upper)
  }

  Z
}

expected_binary_Z_given_model_alpha <- function(X, F_hat, Lambda, alpha) {
  X <- as.matrix(X)
  mean_mat <- sweep(F_hat %*% t(Lambda), 2L, as.numeric(alpha), "+")
  Z_mean <- matrix(NA_real_, nrow(X), ncol(X), dimnames = dimnames(X))

  for (j in seq_len(ncol(X))) {
    lower <- ifelse(X[, j] == 1, 0, -Inf)
    upper <- ifelse(X[, j] == 1, Inf, 0)
    Z_mean[, j] <- truncnorm_binary_moments_vec(
      mean = mean_mat[, j],
      sd = 1,
      lower = lower,
      upper = upper
    )$mean
  }

  Z_mean
}

initialize_binary_Z_intercept <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  alpha <- qnorm(clip01(colMeans(X), n))
  F0 <- matrix(0, nrow(X), 1L)
  Lambda0 <- matrix(0, ncol(X), 1L)
  expected_binary_Z_given_model_alpha(X, F0, Lambda0, alpha)
}

# ----------------------------------------------------------------------------
# Intercept-aware loading updates
# ----------------------------------------------------------------------------

update_working_loadings_with_intercept <- function(
    Z,
    F_hat,
    loading_penalty = 0,
    min_residual_var = 1e-4) {
  Z <- as.matrix(Z)
  F_hat <- as.matrix(F_hat)
  n <- nrow(Z)
  p <- ncol(Z)

  # Since SVD is performed on column-centered Z, the rotated scores are nearly
  # mean zero.  The intercept is therefore the column mean of the current Z.
  alpha <- colMeans(Z)
  Z_centered <- sweep(Z, 2L, alpha, "-")

  Lambda_ls <- crossprod(Z_centered, F_hat) / n
  Lambda <- soft_threshold(Lambda_ls, loading_penalty)
  fitted <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  residual <- Z - fitted
  psi <- pmax(colMeans(residual^2), min_residual_var)

  rownames(Lambda) <- colnames(Z)
  colnames(Lambda) <- paste0("factor_", seq_len(ncol(F_hat)))

  list(
    alpha = alpha,
    Lambda = Lambda,
    Lambda_ls = Lambda_ls,
    loading_update = if (loading_penalty > 0) "intercept_lasso_soft_threshold" else "intercept_least_squares",
    loading_penalty = loading_penalty,
    Psi = diag(psi, p),
    fitted = fitted,
    residual = residual
  )
}

probit_item_negloglik_alpha <- function(y, F_hat, theta, lambda_l1_penalty = 0) {
  alpha <- theta[1L]
  beta <- theta[-1L]
  eta <- alpha + as.numeric(F_hat %*% beta)
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  -sum(y * log(p1) + (1 - y) * log(p0)) + lambda_l1_penalty * sum(abs(beta))
}

probit_item_negloglik_grad_alpha <- function(y, F_hat, theta) {
  alpha <- theta[1L]
  beta <- theta[-1L]
  eta <- alpha + as.numeric(F_hat %*% beta)
  phi <- dnorm(eta)
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  score_eta <- y * phi / p1 - (1 - y) * phi / p0
  -c(sum(score_eta), as.numeric(crossprod(F_hat, score_eta)))
}

fit_probit_lasso_item_alpha <- function(
    y,
    F_hat,
    theta_init,
    lambda_l1_penalty = 0,
    maxit = 300L,
    tol = 1e-6) {
  theta <- as.numeric(theta_init)
  if (length(theta) != ncol(F_hat) + 1L || any(!is.finite(theta))) {
    theta <- c(qnorm(clip01(mean(y), length(y))), rep(0, ncol(F_hat)))
  }

  old_obj <- probit_item_negloglik_alpha(y, F_hat, theta, lambda_l1_penalty)
  for (iter in seq_len(maxit)) {
    grad <- probit_item_negloglik_grad_alpha(y, F_hat, theta)
    step <- 1
    accepted <- FALSE

    for (bt in seq_len(35L)) {
      theta_trial <- theta - step * grad
      theta_trial[-1L] <- soft_threshold(theta_trial[-1L], step * lambda_l1_penalty)
      trial_obj <- probit_item_negloglik_alpha(y, F_hat, theta_trial, lambda_l1_penalty)
      if (is.finite(trial_obj) && trial_obj <= old_obj + 1e-10) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }

    if (!accepted) break
    max_change <- max(abs(theta_trial - theta))
    theta <- theta_trial
    old_obj <- trial_obj
    if (max_change <= tol * (1 + max(abs(theta)))) break
  }

  theta
}

update_binary_probit_loadings_glm_alpha <- function(
    X,
    F_hat,
    Lambda_init = NULL,
    alpha_init = NULL,
    lambda_l1_penalty = 0,
    lasso_maxit = 300L,
    lasso_tol = 1e-6,
    parallel = FALSE,
    workers = NULL) {
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  n <- nrow(X)
  p <- ncol(X)
  H <- ncol(F_hat)

  alpha <- qnorm(clip01(colMeans(X), n))
  Lambda <- matrix(0, p, H, dimnames = list(colnames(X), paste0("factor_", seq_len(H))))

  dat <- as.data.frame(F_hat)
  names(dat) <- paste0("factor_", seq_len(H))
  form <- as.formula(paste("y ~", paste(names(dat), collapse = " + ")))

  fit_one_item <- function(j) {
    y <- X[, j]
    empirical_alpha <- qnorm(clip01(mean(y), n))
    beta_init <- if (!is.null(Lambda_init)) Lambda_init[j, ] else rep(0, H)
    a_init <- if (!is.null(alpha_init) && is.finite(alpha_init[j])) alpha_init[j] else empirical_alpha

    # Constant columns are perfectly separated.  Keep an intercept-only fit.
    if (length(unique(y)) == 1L) {
      return(c(empirical_alpha, rep(0, H)))
    }

    if (lambda_l1_penalty > 0) {
      return(fit_probit_lasso_item_alpha(
        y = y,
        F_hat = F_hat,
        theta_init = c(a_init, beta_init),
        lambda_l1_penalty = lambda_l1_penalty,
        maxit = lasso_maxit,
        tol = lasso_tol
      ))
    }

    dat_j <- dat
    dat_j$y <- y
    fit_j <- tryCatch(
      suppressWarnings(glm(form, data = dat_j, family = binomial(link = "probit"))),
      error = function(e) NULL
    )

    if (!is.null(fit_j) && all(is.finite(coef(fit_j)))) {
      coef(fit_j)
    } else {
      c(a_init, beta_init)
    }
  }

  rows <- parallel_lapply(seq_len(p), fit_one_item, parallel = parallel, workers = workers)
  mat <- do.call(rbind, rows)
  alpha[] <- mat[, 1L]
  Lambda[,] <- mat[, -1L, drop = FALSE]

  list(alpha = alpha, Lambda = Lambda)
}

# ----------------------------------------------------------------------------
# H pruning with intercepts
# ----------------------------------------------------------------------------

prune_factors_by_probit_drop_alpha <- function(
    X,
    F_hat,
    Lambda,
    alpha,
    mixture_fits,
    loading_threshold = 0.08,
    min_active_loadings = 10L,
    min_loading_energy_share = 0.025,
    likelihood_penalty_multiplier = 0.25,
    min_keep = 1L,
    penalty_mode = "active_log_n") {
  penalty_mode <- match.arg(penalty_mode)
  H <- ncol(F_hat)
  n <- nrow(X)
  p <- ncol(X)
  n_response <- length(X)
  full_ll <- binary_probit_loglik_alpha(X, F_hat, Lambda, alpha)
  active_count <- colSums(abs(Lambda) > loading_threshold)
  loading_energy <- colSums(Lambda^2)
  loading_energy_share <- loading_energy / pmax(sum(loading_energy), 1e-12)
  mixture_param_count <- vapply(mixture_fits, function(fit) {
    G <- length(fit$pi)
    3L * G - 1L
  }, integer(1))
  ll_drop <- numeric(H)

  for (h in seq_len(H)) {
    keep_h <- setdiff(seq_len(H), h)
    if (length(keep_h) == 0L) {
      minus_ll <- binary_probit_loglik_alpha(
        X,
        matrix(0, nrow(X), 1L),
        matrix(0, ncol(X), 1L),
        alpha
      )
    } else {
      minus_ll <- binary_probit_loglik_alpha(
        X,
        F_hat[, keep_h, drop = FALSE],
        Lambda[, keep_h, drop = FALSE],
        alpha
      )
    }
    ll_drop[h] <- full_ll - minus_ll
  }

  effective_parameter_count <- active_count
  penalty_sample_size <- n
  penalty <- likelihood_penalty_multiplier *
    effective_parameter_count *
    log(max(penalty_sample_size, 2L))
  keep <- active_count >= min_active_loadings &
    loading_energy_share >= min_loading_energy_share &
    (2 * ll_drop) > penalty

  min_keep <- max(1L, min(as.integer(min_keep), H))
  if (sum(keep) < min_keep) {
    rank_score <- 2 * ll_drop - penalty
    keep[order(rank_score, loading_energy, decreasing = TRUE)[seq_len(min_keep)]] <- TRUE
  }

  list(
    kept = which(keep),
    dropped = which(!keep),
    diagnostics = data.frame(
      factor = seq_len(H),
      active_loadings = active_count,
      loading_energy = loading_energy,
      loading_energy_share = loading_energy_share,
      probit_loglik_drop = ll_drop,
      effective_parameter_count = effective_parameter_count,
      penalty_mode = penalty_mode,
      penalty_sample_size = penalty_sample_size,
      penalty = penalty,
      keep = keep
    )
  )
}

# ----------------------------------------------------------------------------
# Intercept-aware pretraining
# ----------------------------------------------------------------------------

fit_binary_probit_pretraining_intercept <- function(
    X,
    H = NULL,
    H_max = min(10L, nrow(as.matrix(X)) - 1L, ncol(as.matrix(X))),
    G_fixed = 3L,
    n_aug_iter = 8L,
    z_update = c("expectation", "sample"),
    n_random_starts = 1L,
    max_outer = 4L,
    n_mix_starts = 3L,
    mixture_max_iter = 20L,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.05,
    var_prior_shape = 4,
    var_prior_scale = 0.35,
    weight_prior_alpha = 1.2,
    loading_penalty = 0.02,
    H_prune_after_iter = 4L,
    H_prune_loading_threshold = 0.08,
    H_prune_min_active_loadings = 10L,
    H_prune_min_loading_energy_share = 0.025,
    H_prune_likelihood_penalty_multiplier = 0.25,
    H_prune_penalty_mode = "active_log_n",
    H_prune_min_keep = 1L,
    objective_tolerance = 2e-4,
    objective_tolerance_scale = c("relative_total", "per_response"),
    min_aug_iter = 5L,
    grid_size = 21L,
    parallel = FALSE,
    workers = NULL,
    seed = 20260717L,
    verbose = TRUE) {
  z_update <- match.arg(z_update)
  mixture_update <- match.arg(mixture_update)
  H_prune_penalty_mode <- match.arg(H_prune_penalty_mode)
  objective_tolerance_scale <- match.arg(objective_tolerance_scale)
  workers <- resolve_workers(workers)

  X <- as.matrix(X)
  if (!all(X %in% c(0, 1))) stop("X must be binary.")
  if (is.null(colnames(X))) colnames(X) <- paste0("problem_", seq_len(ncol(X)))

  set.seed(seed)
  Z <- initialize_binary_Z_intercept(X)
  H_selection <- NULL
  H_pruning <- NULL
  H_pruned <- !is.null(H)
  if (is.null(H)) {
    H <- as.integer(min(H_max, nrow(X) - 1L, ncol(X)))
    H_selection <- list(strategy = "overfit_prune_after_augmented_Z_updates", H_start = H)
    if (verbose) message("Starting overfit-prune pretraining at H = ", H, ".")
  }

  history <- vector("list", n_aug_iter)
  trace <- vector("list", n_aug_iter)
  current <- NULL
  converged <- FALSE

  for (iter in seq_len(n_aug_iter)) {
    if (verbose) message("Intercept pretraining iteration ", iter, " of ", n_aug_iter, ".")

    svd_out <- svd_scores_from_Z(Z, H = H, center_Z = TRUE)
    S <- svd_out$S
    G_iter <- if (length(G_fixed) == 1L) rep(as.integer(G_fixed), H) else as.integer(G_fixed)
    if (length(G_iter) != H) {
      stop("G_fixed must be scalar or a length-H component-count vector.")
    }

    rotation_out <- estimate_mixture_ica_unknown_G(
      S = S,
      G_fixed = G_iter,
      n_random_starts = n_random_starts,
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
      seed = seed + 1000L + iter,
      parallel = parallel,
      workers = workers,
      verbose = verbose
    )

    F_hat <- rotation_out$F_hat
    loading_out <- update_working_loadings_with_intercept(
      Z = Z,
      F_hat = F_hat,
      loading_penalty = loading_penalty
    )

    H_before_prune <- H
    prune_diagnostics <- NULL
    pruned_factors <- integer(0)

    if (!isTRUE(H_pruned) && iter >= as.integer(H_prune_after_iter)) {
      prune_out <- prune_factors_by_probit_drop_alpha(
        X = X,
        F_hat = F_hat,
        Lambda = loading_out$Lambda,
        alpha = loading_out$alpha,
        mixture_fits = rotation_out$fits,
        loading_threshold = H_prune_loading_threshold,
        min_active_loadings = H_prune_min_active_loadings,
        min_loading_energy_share = H_prune_min_loading_energy_share,
        likelihood_penalty_multiplier = H_prune_likelihood_penalty_multiplier,
        min_keep = H_prune_min_keep,
        penalty_mode = H_prune_penalty_mode
      )

      F_hat <- F_hat[, prune_out$kept, drop = FALSE]
      loading_out$Lambda <- loading_out$Lambda[, prune_out$kept, drop = FALSE]
      loading_out$Lambda_ls <- loading_out$Lambda_ls[, prune_out$kept, drop = FALSE]
      loading_out$fitted <- sweep(F_hat %*% t(loading_out$Lambda), 2L, loading_out$alpha, "+")
      loading_out$residual <- Z - loading_out$fitted
      rotation_out$fits <- rotation_out$fits[prune_out$kept]
      H <- ncol(F_hat)
      H_pruned <- TRUE
      pruned_factors <- prune_out$dropped
      prune_diagnostics <- prune_out$diagnostics
      H_pruning <- list(
        iteration = iter,
        H_before = H_before_prune,
        H_after = H,
        kept = prune_out$kept,
        dropped = prune_out$dropped,
        diagnostics = prune_out$diagnostics
      )
      if (verbose) message("Pruned H from ", H_before_prune, " to ", H, ".")
    }

    ord <- ordered_component_labels(F_hat, rotation_out$fits)
    class_map <- ord$class_map
    responsibilities <- ord$responsibilities
    G_hat <- vapply(rotation_out$fits, function(z) length(z$pi), integer(1))
    probit_loglik <- binary_probit_loglik_alpha(X, F_hat, loading_out$Lambda, loading_out$alpha)
    mixture_loglik <- mixture_prior_loglik(F_hat, rotation_out$fits)
    full_loglik <- probit_loglik + mixture_loglik
    reconstruction_mse <- mean(loading_out$residual^2)
    previous <- if (iter > 1L) history[[iter - 1L]]$probit_loglik else NA_real_
    abs_change <- if (iter > 1L && is.finite(previous)) abs(probit_loglik - previous) else NA_real_
    rel_change <- if (iter > 1L && is.finite(previous)) {
      abs_change / (1 + abs(previous))
    } else {
      NA_real_
    }
    per_response_change <- abs_change / length(X)
    stopping_change <- if (objective_tolerance_scale == "per_response") {
      per_response_change
    } else {
      rel_change
    }

    history[[iter]] <- data.frame(
      iteration = iter,
      H = H,
      G_hat = paste(G_hat, collapse = ","),
      probit_loglik = probit_loglik,
      probit_loglik_per_response = probit_loglik / length(X),
      mixture_loglik = mixture_loglik,
      full_loglik = full_loglik,
      reconstruction_mse = reconstruction_mse,
      absolute_probit_change = abs_change,
      relative_probit_change = rel_change,
      per_response_probit_change = per_response_change,
      stopping_change = stopping_change,
      stopping_scale = objective_tolerance_scale,
      H_before_prune = H_before_prune,
      H_after_prune = H,
      n_pruned_factors = length(pruned_factors),
      converged = FALSE
    )

    current <- list(
      model = "binary_probit_pretraining_intercept_independent_mixture_factor",
      X = X,
      Z = Z,
      H = H,
      F_hat = F_hat,
      Lambda_hat = loading_out$Lambda,
      alpha_hat = loading_out$alpha,
      mixture_fits = rotation_out$fits,
      G_hat = G_hat,
      class_map = class_map,
      responsibilities = responsibilities,
      profile_id = profile_id_from_class_map(class_map),
      svd_fit = svd_out,
      R = rotation_out$R,
      rotation_fit = rotation_out,
      fitted = loading_out$fitted,
      residual = loading_out$residual,
      Psi_hat = diag(1, ncol(X))
    )

    trace[[iter]] <- list(
      iteration = iter,
      F_hat = F_hat,
      Lambda_hat = loading_out$Lambda,
      alpha_hat = loading_out$alpha,
      mixture_fits = rotation_out$fits,
      class_map = class_map,
      G_hat = G_hat,
      probit_loglik = probit_loglik,
      mixture_loglik = mixture_loglik,
      prune_diagnostics = prune_diagnostics
    )

    if (verbose) {
      message(
        "  H = ", H,
        "; G = [", paste(G_hat, collapse = ", "),
        "]; probit loglik = ", round(probit_loglik, 2),
        "; per response = ", round(probit_loglik / length(X), 5),
        "; reconstruction MSE = ", round(reconstruction_mse, 4)
      )
    }

    if (!is.null(objective_tolerance) &&
        iter >= min_aug_iter &&
        is.finite(stopping_change) &&
        stopping_change <= objective_tolerance) {
      converged <- TRUE
      history[[iter]]$converged <- TRUE
      if (verbose) message("Stopping pretraining by probit loglik tolerance.")
      break
    }

    if (z_update == "sample") {
      Z <- sample_binary_Z_given_model_alpha(
        X = X,
        F_hat = F_hat,
        Lambda = loading_out$Lambda,
        alpha = loading_out$alpha,
        seed = seed + 2000L + iter
      )
    } else {
      Z <- expected_binary_Z_given_model_alpha(
        X = X,
        F_hat = F_hat,
        Lambda = loading_out$Lambda,
        alpha = loading_out$alpha
      )
    }
  }

  history_out <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  trace_out <- trace[!vapply(trace, is.null, logical(1))]

  current$pretraining <- list(
    n_aug_iter = n_aug_iter,
    n_completed = nrow(history_out),
    converged = converged,
    z_update = z_update,
    H_selection = H_selection,
    H_pruning = H_pruning,
    G_fixed = G_fixed,
    loading_penalty = loading_penalty,
    objective_tolerance = objective_tolerance,
    objective_tolerance_scale = objective_tolerance_scale,
    history = history_out,
    iteration_trace = trace_out,
    intercept = TRUE,
    parallel = list(enabled = isTRUE(parallel), workers = workers)
  )

  current
}

# ----------------------------------------------------------------------------
# Intercept-aware refinement
# ----------------------------------------------------------------------------

update_one_factor_score_alpha <- function(
    x_i,
    f_init,
    Lambda,
    alpha,
    mixture_fits,
    mixture_prior_weight = 1,
    maxit = 50L) {
  x_i <- as.numeric(x_i)
  f_init <- as.numeric(f_init)
  Lambda <- as.matrix(Lambda)
  alpha <- as.numeric(alpha)
  H <- length(f_init)

  objective <- function(f) {
    eta <- alpha + as.numeric(Lambda %*% f)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    ll <- sum(x_i * log(p1) + (1 - x_i) * log(p0))

    lp <- 0
    for (h in seq_len(H)) {
      lp <- lp + single_factor_log_prior_grad(f[h], mixture_fits[[h]])$log_density
    }

    -(ll + mixture_prior_weight * lp)
  }

  gradient <- function(f) {
    eta <- alpha + as.numeric(Lambda %*% f)
    phi <- dnorm(eta)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    score_eta <- x_i * phi / p1 - (1 - x_i) * phi / p0
    grad <- as.numeric(crossprod(Lambda, score_eta))

    for (h in seq_len(H)) {
      grad[h] <- grad[h] + mixture_prior_weight *
        single_factor_log_prior_grad(f[h], mixture_fits[[h]])$grad
    }

    -grad
  }

  opt <- tryCatch(
    optim(
      par = f_init,
      fn = objective,
      gr = gradient,
      method = "L-BFGS-B",
      control = list(maxit = maxit, factr = 1e7)
    ),
    error = function(e) NULL
  )

  if (is.null(opt) || !all(is.finite(opt$par))) f_init else opt$par
}

update_factor_scores_joint_map_alpha <- function(
    X,
    F_hat,
    Lambda,
    alpha,
    mixture_fits,
    mixture_prior_weight = 1,
    maxit_per_subject = 50L,
    parallel = FALSE,
    workers = NULL) {
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)

  rows <- parallel_lapply(seq_len(nrow(X)), function(i) {
    update_one_factor_score_alpha(
      x_i = X[i, ],
      f_init = F_hat[i, ],
      Lambda = Lambda,
      alpha = alpha,
      mixture_fits = mixture_fits,
      mixture_prior_weight = mixture_prior_weight,
      maxit = maxit_per_subject
    )
  }, parallel = parallel, workers = workers)

  out <- F_hat
  out[,] <- do.call(rbind, rows)
  out
}

normalize_factor_scale_alpha <- function(F_hat, Lambda, mixture_fits, target_scale = 1) {
  scaled <- normalize_refinement_factor_scale(
    F_hat = F_hat,
    Lambda = Lambda,
    mixture_fits = mixture_fits,
    target_scale = target_scale,
    scale_method = "sd"
  )
  scaled
}

fit_binary_probit_refinement_intercept <- function(
    X,
    pretrain_fit,
    n_refine_iter = 8L,
    maxit_per_subject = 60L,
    n_mix_starts = 3L,
    mixture_max_iter = 20L,
    min_mixture_var = 0.05,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.05,
    var_prior_shape = 4,
    var_prior_scale = 0.35,
    weight_prior_alpha = 1.2,
    mixture_prior_weight = 0.35,
    lambda_l1_penalty = 0.25,
    lasso_maxit = 300L,
    lasso_tol = 1e-6,
    objective_tolerance = 2e-4,
    objective_tolerance_scale = c("relative_total", "per_response"),
    min_refine_iter = 4L,
    require_mixture_convergence = FALSE,
    keep_best_binary_iterate = TRUE,
    parallel = FALSE,
    workers = NULL,
    verbose = TRUE) {
  mixture_update <- match.arg(mixture_update)
  objective_tolerance_scale <- match.arg(objective_tolerance_scale)
  require_mixture_convergence <- isTRUE(require_mixture_convergence)
  workers <- resolve_workers(workers)
  X <- as.matrix(X)

  F_hat <- pretrain_fit$F_hat
  Lambda <- pretrain_fit$Lambda_hat
  alpha <- pretrain_fit$alpha_hat
  mixture_fits <- pretrain_fit$mixture_fits

  if (verbose) message("Pre-estimating intercept probit loadings before refinement.")
  load_out <- update_binary_probit_loadings_glm_alpha(
    X = X,
    F_hat = F_hat,
    Lambda_init = Lambda,
    alpha_init = alpha,
    lambda_l1_penalty = lambda_l1_penalty,
    lasso_maxit = lasso_maxit,
    lasso_tol = lasso_tol,
    parallel = parallel,
    workers = workers
  )
  alpha <- load_out$alpha
  Lambda <- load_out$Lambda

  history <- vector("list", n_refine_iter + 1L)
  trace <- vector("list", n_refine_iter + 1L)
  current_binary <- binary_probit_loglik_alpha(X, F_hat, Lambda, alpha)
  current_mixture <- mixture_prior_loglik(F_hat, mixture_fits)
  history[[1L]] <- data.frame(
    iteration = 0L,
    binary_loglik = current_binary,
    binary_loglik_per_response = current_binary / length(X),
    mixture_loglik = current_mixture,
    joint_objective = current_binary + mixture_prior_weight * current_mixture,
    absolute_binary_change = NA_real_,
    relative_binary_change = NA_real_,
    per_response_binary_change = NA_real_,
    stopping_change = NA_real_,
    stopping_scale = objective_tolerance_scale,
    G_hat = paste(vapply(mixture_fits, function(z) length(z$pi), integer(1)), collapse = ","),
    all_mixtures_converged = all(vapply(mixture_fits, function(z) isTRUE(z$converged), logical(1))),
    converged = FALSE
  )
  trace[[1L]] <- list(
    iteration = 0L,
    F_hat = F_hat,
    Lambda_hat = Lambda,
    alpha_hat = alpha,
    mixture_fits = mixture_fits
  )

  converged <- FALSE
  n_completed <- 0L

  for (iter in seq_len(n_refine_iter)) {
    if (verbose) message("Intercept refinement iteration ", iter, " of ", n_refine_iter, ".")

    F_hat <- update_factor_scores_joint_map_alpha(
      X = X,
      F_hat = F_hat,
      Lambda = Lambda,
      alpha = alpha,
      mixture_fits = mixture_fits,
      mixture_prior_weight = mixture_prior_weight,
      maxit_per_subject = maxit_per_subject,
      parallel = parallel,
      workers = workers
    )

    scaled <- normalize_factor_scale_alpha(
      F_hat = F_hat,
      Lambda = Lambda,
      mixture_fits = mixture_fits,
      target_scale = 1
    )
    F_hat <- scaled$F_hat
    Lambda <- scaled$Lambda
    mixture_fits <- scaled$mixture_fits

    load_out <- update_binary_probit_loadings_glm_alpha(
      X = X,
      F_hat = F_hat,
      Lambda_init = Lambda,
      alpha_init = alpha,
      lambda_l1_penalty = lambda_l1_penalty,
      lasso_maxit = lasso_maxit,
      lasso_tol = lasso_tol,
      parallel = parallel,
      workers = workers
    )
    alpha <- load_out$alpha
    Lambda <- load_out$Lambda

    mixture_fits <- update_mixture_fits_refinement(
      F_hat = F_hat,
      mixture_fits = mixture_fits,
      n_starts = n_mix_starts,
      max_iter = mixture_max_iter,
      min_var = min_mixture_var,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      parallel = parallel,
      workers = workers
    )
    all_mixtures_converged <- all(vapply(mixture_fits, function(z) isTRUE(z$converged), logical(1)))

    current_binary <- binary_probit_loglik_alpha(X, F_hat, Lambda, alpha)
    current_mixture <- mixture_prior_loglik(F_hat, mixture_fits)
    previous_binary <- history[[iter]]$binary_loglik
    abs_change <- abs(current_binary - previous_binary)
    rel_change <- abs_change / (1 + abs(previous_binary))
    per_response_change <- abs_change / length(X)
    stopping_change <- if (objective_tolerance_scale == "per_response") {
      per_response_change
    } else {
      rel_change
    }

    history[[iter + 1L]] <- data.frame(
      iteration = iter,
      binary_loglik = current_binary,
      binary_loglik_per_response = current_binary / length(X),
      mixture_loglik = current_mixture,
      joint_objective = current_binary + mixture_prior_weight * current_mixture,
      absolute_binary_change = abs_change,
      relative_binary_change = rel_change,
      per_response_binary_change = per_response_change,
      stopping_change = stopping_change,
      stopping_scale = objective_tolerance_scale,
      G_hat = paste(vapply(mixture_fits, function(z) length(z$pi), integer(1)), collapse = ","),
      all_mixtures_converged = all_mixtures_converged,
      converged = FALSE
    )
    trace[[iter + 1L]] <- list(
      iteration = iter,
      F_hat = F_hat,
      Lambda_hat = Lambda,
      alpha_hat = alpha,
      mixture_fits = mixture_fits
    )

    n_completed <- iter
    if (verbose) {
      message(
        "  binary loglik = ", round(current_binary, 2),
        "; per response = ", round(current_binary / length(X), 5),
        "; mixture loglik = ", round(current_mixture, 2),
        "; stopping change = ", signif(stopping_change, 3),
        "; all mixtures converged = ", all_mixtures_converged
      )
    }

    if (!is.null(objective_tolerance) &&
        is.finite(objective_tolerance) &&
        iter >= min_refine_iter &&
        is.finite(stopping_change) &&
        stopping_change <= objective_tolerance &&
        (!require_mixture_convergence || isTRUE(all_mixtures_converged))) {
      converged <- TRUE
      history[[iter + 1L]]$converged <- TRUE
      if (verbose) message("Stopping refinement by binary loglik tolerance.")
      break
    }
  }

  history_out <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  trace_out <- trace[!vapply(trace, is.null, logical(1))]
  best_candidates <- seq_len(nrow(history_out))
  if (require_mixture_convergence && "all_mixtures_converged" %in% names(history_out)) {
    converged_candidates <- which(history_out$all_mixtures_converged %in% TRUE)
    if (length(converged_candidates) > 0L) best_candidates <- converged_candidates
  }
  best_trace_index <- best_candidates[which.max(history_out$binary_loglik[best_candidates])]
  best_iteration <- history_out$iteration[best_trace_index]

  if (isTRUE(keep_best_binary_iterate) && length(trace_out) >= best_trace_index) {
    best_state <- trace_out[[best_trace_index]]
    F_hat <- best_state$F_hat
    Lambda <- best_state$Lambda_hat
    alpha <- best_state$alpha_hat
    mixture_fits <- best_state$mixture_fits
    if (verbose && best_iteration != tail(history_out$iteration, 1L)) {
      message("Restoring refinement iterate ", best_iteration, " with best binary loglik.")
    }
  }

  ord <- ordered_component_labels(F_hat, mixture_fits)
  class_map <- ord$class_map
  responsibilities <- ord$responsibilities

  out <- pretrain_fit
  out$model <- "joint_refined_binary_probit_intercept_independent_mixture_factor"
  out$F_hat <- F_hat
  out$Lambda_hat <- Lambda
  out$alpha_hat <- alpha
  out$mixture_fits <- mixture_fits
  out$G_hat <- vapply(mixture_fits, function(z) length(z$pi), integer(1))
  out$class_map <- class_map
  out$responsibilities <- responsibilities
  out$profile_id <- profile_id_from_class_map(class_map)
  out$joint_refinement <- list(
    n_refine_iter = n_refine_iter,
    n_completed = n_completed,
    converged = converged,
    objective_tolerance = objective_tolerance,
    objective_tolerance_scale = objective_tolerance_scale,
    min_refine_iter = min_refine_iter,
    require_mixture_convergence = require_mixture_convergence,
    keep_best_binary_iterate = keep_best_binary_iterate,
    selected_iteration = best_iteration,
    selected_binary_loglik = history_out$binary_loglik[best_trace_index],
    mixture_prior_weight = mixture_prior_weight,
    lambda_l1_penalty = lambda_l1_penalty,
    min_mixture_var = min_mixture_var,
    history = history_out,
    refinement_trace = trace_out,
    intercept = TRUE,
    parallel = list(enabled = isTRUE(parallel), workers = workers)
  )

  out
}
