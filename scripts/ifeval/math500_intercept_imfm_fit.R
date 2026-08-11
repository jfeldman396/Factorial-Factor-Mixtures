#!/usr/bin/env Rscript

# ============================================================================
# Math500 binary-probit independent-mixture factor model with item intercepts
#
# This script fits the pretraining + refinement workflow to accuracy_math500.csv.
# Rows are LLM/evaluation settings, columns are Math500 problems, and entries are
# binary correctness indicators.
#
# Model:
#   X_ij | f_i, alpha_j, lambda_j ~ Bernoulli(Phi(alpha_j + lambda_j' f_i))
#   f_ih independently follows a 1D Gaussian mixture.
#
# The item intercept alpha_j captures baseline problem difficulty.  The latent
# factors f_i then represent model-specific deviations from these difficulties.
#
# Workflow:
#   1. Initialize latent probit Z from the intercept-only truncation model.
#   2. Repeat pretraining:
#        a. SVD on centered Z.
#        b. Rotate SVD scores toward independent Gaussian-mixture coordinates.
#        c. Estimate alpha and Lambda from Z ~ alpha + F Lambda'.
#        d. Optionally prune an overfit starting H by probit likelihood evidence.
#        e. Update Z by E[Z | X, F, alpha, Lambda].
#   3. Refine:
#        a. Update each model's factor vector by MAP under the probit likelihood
#           plus the independent mixture prior.
#        b. Update alpha and Lambda by itemwise probit regressions.
#        c. Refit each marginal mixture with fixed G.
#   4. Save summaries and plots for the LLM evaluations.
# ============================================================================

script_dir <- local({
  frames <- sys.frames()
  files <- vapply(frames, function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  files <- files[!is.na(files)]
  if (length(files) > 0L) dirname(normalizePath(tail(files, 1L))) else getwd()
})

source(file.path(script_dir, "binary_probit_pretraining_algorithm_commented.R"))
source(file.path(script_dir, "binary_probit_refinement_algorithm_commented.R"))

# ----------------------------------------------------------------------------
# Small utilities
# ----------------------------------------------------------------------------

clip01 <- function(p, n) {
  eps <- 0.5 / max(1L, n)
  pmin(pmax(p, eps), 1 - eps)
}

parse_prompt_type <- function(model_id) {
  ifelse(grepl("_zero_shot$", model_id), "zero_shot",
         ifelse(grepl("_one_shot$", model_id), "one_shot", "unknown"))
}

parse_base_model <- function(model_id) {
  sub("_(zero_shot|one_shot)$", "", model_id)
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
  # Factor signs are unidentified.  For reporting, orient each factor so larger
  # values are positively associated with overall Math500 accuracy.
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
    penalty_mode = c("active_log_n", "bic_effective")) {
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

  effective_parameter_count <- if (penalty_mode == "bic_effective") {
    n + active_count + mixture_param_count
  } else {
    active_count
  }
  penalty_sample_size <- if (penalty_mode == "bic_effective") n_response else n
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
    H_prune_penalty_mode = c("active_log_n", "bic_effective"),
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
    G_iter <- rep(as.integer(G_fixed), H)

    rotation_out <- estimate_mixture_ica_unknown_G(
      S = S,
      G_selection = "fixed",
      G_fixed = G_iter,
      n_random_starts = n_random_starts,
      max_outer = max_outer,
      n_mix_starts = n_mix_starts,
      mixture_penalty_multiplier = 1,
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
    keep_best_binary_iterate = TRUE,
    parallel = FALSE,
    workers = NULL,
    verbose = TRUE) {
  mixture_update <- match.arg(mixture_update)
  objective_tolerance_scale <- match.arg(objective_tolerance_scale)
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
      G_selection = "fixed",
      n_starts = n_mix_starts,
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
        "; stopping change = ", signif(stopping_change, 3)
      )
    }

    if (iter >= min_refine_iter &&
        is.finite(stopping_change) &&
        stopping_change <= objective_tolerance) {
      converged <- TRUE
      history[[iter + 1L]]$converged <- TRUE
      if (verbose) message("Stopping refinement by binary loglik tolerance.")
      break
    }
  }

  history_out <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  trace_out <- trace[!vapply(trace, is.null, logical(1))]
  best_trace_index <- which.max(history_out$binary_loglik)
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

# ----------------------------------------------------------------------------
# Summaries and plots
# ----------------------------------------------------------------------------

summarize_mixture_profiles_intercept <- function(fit) {
  do.call(rbind, lapply(seq_len(fit$H), function(h) {
    cls <- fit$class_map[, h]
    tab <- tabulate(cls, nbins = fit$G_hat[h])
    ord <- order(fit$mixture_fits[[h]]$mu)
    data.frame(
      factor = h,
      group = seq_len(fit$G_hat[h]),
      n = tab,
      prop = tab / length(cls),
      mean = fit$mixture_fits[[h]]$mu[ord],
      sd = sqrt(fit$mixture_fits[[h]]$var[ord]),
      weight = fit$mixture_fits[[h]]$pi[ord]
    )
  }))
}

save_model_outputs <- function(fit, model_id, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  X <- fit$X
  score_df <- data.frame(
    model_id = model_id,
    base_model = parse_base_model(model_id),
    prompt_type = parse_prompt_type(model_id),
    accuracy = rowMeans(X),
    fit$F_hat,
    fit$class_map,
    profile_id = fit$profile_id,
    check.names = FALSE
  )
  names(score_df)[5:(4 + fit$H)] <- paste0("factor_", seq_len(fit$H))
  names(score_df)[(5 + fit$H):(4 + 2 * fit$H)] <- paste0("group_factor_", seq_len(fit$H))

  load_df <- data.frame(
    problem = colnames(X),
    empirical_accuracy = colMeans(X),
    alpha = fit$alpha_hat,
    fit$Lambda_hat,
    check.names = FALSE
  )
  names(load_df)[4:ncol(load_df)] <- paste0("loading_factor_", seq_len(fit$H))

  write.csv(score_df, file.path(output_dir, "math500_model_scores_profiles.csv"), row.names = FALSE)
  write.csv(load_df, file.path(output_dir, "math500_problem_intercepts_loadings.csv"), row.names = FALSE)
  write.csv(summarize_mixture_profiles_intercept(fit),
            file.path(output_dir, "math500_factor_mixture_groups.csv"),
            row.names = FALSE)
  write.csv(fit$pretraining$history,
            file.path(output_dir, "math500_pretraining_history.csv"),
            row.names = FALSE)
  write.csv(fit$joint_refinement$history,
            file.path(output_dir, "math500_refinement_history.csv"),
            row.names = FALSE)
  saveRDS(fit, file.path(output_dir, "math500_intercept_imfm_fit.rds"))

  invisible(list(scores = score_df, loadings = load_df))
}

plot_factor_heatmap <- function(score_df, H, output_dir) {
  png(file.path(output_dir, "math500_refined_factor_heatmap.png"), width = 1200, height = 1500, res = 150)
  op <- par(mar = c(4, 7, 3, 2))
  on.exit(par(op), add = TRUE)
  ord <- order(score_df$accuracy, decreasing = TRUE)
  F_mat <- as.matrix(score_df[ord, paste0("factor_", seq_len(H)), drop = FALSE])
  pal <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)
  br <- seq(-max(abs(F_mat)), max(abs(F_mat)), length.out = 102)
  image(
    x = seq_len(H),
    y = seq_len(nrow(F_mat)),
    z = t(F_mat[nrow(F_mat):1L, , drop = FALSE]),
    col = pal,
    breaks = br,
    axes = FALSE,
    xlab = "latent factor",
    ylab = "LLM evaluations ordered by accuracy, best at top",
    main = "Refined Math500 latent factor scores"
  )
  axis(1, at = seq_len(H), labels = paste0("F", seq_len(H)))
  axis(2, labels = FALSE)
  box()
  dev.off()
}

plot_factor_scores_llm_columns <- function(score_df, H, output_dir) {
  # Transposed continuous factor-score heatmap: LLM/evaluation labels on
  # columns and latent factors on rows.  Columns are ordered by accuracy.
  png(
    file.path(output_dir, "math500_factor_scores_llm_columns.png"),
    width = 3400,
    height = 900,
    res = 180
  )
  op <- par(mar = c(12, 5, 4, 2))
  on.exit(par(op), add = TRUE)

  ord <- order(score_df$accuracy, decreasing = TRUE)
  F_mat <- t(as.matrix(score_df[ord, paste0("factor_", seq_len(H)), drop = FALSE]))
  colnames(F_mat) <- score_df$model_id[ord]
  rownames(F_mat) <- paste0("F", seq_len(H))

  max_abs <- max(abs(F_mat), na.rm = TRUE)
  pal <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)
  br <- seq(-max_abs, max_abs, length.out = 102)

  image(
    x = seq_len(ncol(F_mat)),
    y = seq_len(nrow(F_mat)),
    z = t(F_mat),
    col = pal,
    breaks = br,
    axes = FALSE,
    xlab = "",
    ylab = "latent factor",
    main = "Refined factor scores by LLM evaluation, ordered by accuracy"
  )
  axis(2, at = seq_len(H), labels = rownames(F_mat), las = 1)
  axis(
    1,
    at = seq_len(ncol(F_mat)),
    labels = colnames(F_mat),
    las = 2,
    cex.axis = 0.28,
    tick = FALSE
  )
  box()
  dev.off()
}

plot_group_heatmap <- function(score_df, H, output_dir) {
  png(file.path(output_dir, "math500_latent_group_heatmap.png"), width = 1200, height = 1500, res = 150)
  op <- par(mar = c(4, 7, 3, 2))
  on.exit(par(op), add = TRUE)
  ord <- order(score_df$accuracy, decreasing = TRUE)
  G_mat <- as.matrix(score_df[ord, paste0("group_factor_", seq_len(H)), drop = FALSE])
  image(
    x = seq_len(H),
    y = seq_len(nrow(G_mat)),
    z = t(G_mat[nrow(G_mat):1L, , drop = FALSE]),
    col = c("#3B6EA8", "#F2C14E", "#9E2F44"),
    breaks = c(0.5, 1.5, 2.5, 3.5),
    axes = FALSE,
    xlab = "latent factor",
    ylab = "LLM evaluations ordered by accuracy, best at top",
    main = "Ordered latent group labels by factor"
  )
  axis(1, at = seq_len(H), labels = paste0("F", seq_len(H)))
  axis(2, labels = FALSE)
  box()
  legend("topright", fill = c("#3B6EA8", "#F2C14E", "#9E2F44"),
         legend = c("low", "middle", "high"), cex = 0.8, bty = "n")
  dev.off()
}

plot_group_heatmap_llm_columns <- function(score_df, H, output_dir) {
  # Transposed latent group heatmap: LLM/evaluation labels on columns and
  # factors on rows.  This is more natural for scanning which latent profile
  # each model receives across factors.
  png(
    file.path(output_dir, "math500_latent_group_heatmap_llm_columns.png"),
    width = 3400,
    height = 900,
    res = 180
  )
  op <- par(mar = c(12, 5, 4, 2))
  on.exit(par(op), add = TRUE)

  ord <- order(score_df$accuracy, decreasing = TRUE)
  G_mat <- t(as.matrix(score_df[ord, paste0("group_factor_", seq_len(H)), drop = FALSE]))
  colnames(G_mat) <- score_df$model_id[ord]
  rownames(G_mat) <- paste0("F", seq_len(H))

  image(
    x = seq_len(ncol(G_mat)),
    y = seq_len(nrow(G_mat)),
    z = t(G_mat),
    col = c("#3B6EA8", "#F2C14E", "#9E2F44"),
    breaks = c(0.5, 1.5, 2.5, 3.5),
    axes = FALSE,
    xlab = "",
    ylab = "latent factor",
    main = "Latent factor groups by LLM evaluation, ordered by accuracy"
  )
  axis(2, at = seq_len(H), labels = rownames(G_mat), las = 1)
  axis(
    1,
    at = seq_len(ncol(G_mat)),
    labels = colnames(G_mat),
    las = 2,
    cex.axis = 0.28,
    tick = FALSE
  )
  box()
  legend(
    "topright",
    fill = c("#3B6EA8", "#F2C14E", "#9E2F44"),
    legend = c("low", "middle", "high"),
    cex = 0.8,
    bty = "n"
  )
  dev.off()
}

plot_accuracy_scatter <- function(score_df, H, output_dir) {
  if (H < 2L) return(invisible(NULL))
  png(file.path(output_dir, "math500_factor1_factor2_accuracy_scatter.png"), width = 1200, height = 900, res = 150)
  op <- par(mar = c(5, 5, 3, 6))
  on.exit(par(op), add = TRUE)
  acc <- score_df$accuracy
  pal <- colorRampPalette(c("#355C9A", "#F2C14E", "#9E2F44"))(100)
  idx <- pmax(1, pmin(100, as.integer(cut(acc, breaks = 100, labels = FALSE))))
  plot(
    score_df$factor_1,
    score_df$factor_2,
    pch = ifelse(score_df$prompt_type == "one_shot", 16, 17),
    col = pal[idx],
    xlab = "factor 1",
    ylab = "factor 2",
    main = "Math500 LLM evaluations in first two refined factors"
  )
  top_idx <- order(score_df$accuracy, decreasing = TRUE)[seq_len(min(10L, nrow(score_df)))]
  text(
    score_df$factor_1[top_idx],
    score_df$factor_2[top_idx],
    labels = seq_len(length(top_idx)),
    cex = 0.65,
    pos = 3
  )
  legend("right", inset = c(-0.2, 0), xpd = TRUE, bty = "n",
         legend = c("one_shot", "zero_shot"), pch = c(16, 17), col = "gray30")
  dev.off()
}

plot_factor_marginals <- function(fit, output_dir) {
  H <- fit$H
  nrow_plot <- ceiling(H / 2)
  png(file.path(output_dir, "math500_factor_marginals_with_mixture_centers.png"),
      width = 1500, height = max(900, 300 * nrow_plot), res = 150)
  op <- par(mfrow = c(nrow_plot, 2), mar = c(4, 4, 3, 1))
  on.exit(par(op), add = TRUE)
  for (h in seq_len(H)) {
    hist(fit$F_hat[, h], breaks = 18, col = "#D8DEE9", border = "white",
         main = paste0("factor ", h, " marginal"), xlab = "refined factor score")
    ord <- order(fit$mixture_fits[[h]]$mu)
    abline(v = fit$mixture_fits[[h]]$mu[ord], col = c("#3B6EA8", "#F2C14E", "#9E2F44"),
           lwd = 2, lty = 2)
  }
  dev.off()
}

plot_item_difficulty <- function(fit, output_dir) {
  png(file.path(output_dir, "math500_problem_difficulty_intercepts.png"), width = 1200, height = 800, res = 150)
  op <- par(mar = c(5, 5, 3, 1))
  on.exit(par(op), add = TRUE)
  hist(fit$alpha_hat, breaks = 30, col = "#D8DEE9", border = "white",
       xlab = "item intercept alpha", main = "Math500 problem difficulty intercepts")
  abline(v = 0, col = "gray35", lwd = 2)
  dev.off()
}

plot_top_loading_heatmap <- function(fit, output_dir, top_n = 80L) {
  Lambda <- fit$Lambda_hat
  keep <- order(rowSums(abs(Lambda)), decreasing = TRUE)[seq_len(min(top_n, nrow(Lambda)))]
  mat <- Lambda[keep, , drop = FALSE]
  png(file.path(output_dir, "math500_top_problem_loading_heatmap.png"), width = 1200, height = 1400, res = 150)
  op <- par(mar = c(4, 7, 3, 2))
  on.exit(par(op), add = TRUE)
  pal <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)
  br <- seq(-max(abs(mat)), max(abs(mat)), length.out = 102)
  image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1L, , drop = FALSE]),
    col = pal,
    breaks = br,
    axes = FALSE,
    xlab = "latent factor",
    ylab = "top-loading Math500 problems",
    main = "Largest refined problem loadings"
  )
  axis(1, at = seq_len(ncol(mat)), labels = paste0("F", seq_len(ncol(mat))))
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.5)
  box()
  dev.off()
}

make_all_plots <- function(fit, score_df, output_dir) {
  plot_factor_heatmap(score_df, fit$H, output_dir)
  plot_factor_scores_llm_columns(score_df, fit$H, output_dir)
  plot_group_heatmap(score_df, fit$H, output_dir)
  plot_group_heatmap_llm_columns(score_df, fit$H, output_dir)
  plot_accuracy_scatter(score_df, fit$H, output_dir)
  plot_factor_marginals(fit, output_dir)
  plot_item_difficulty(fit, output_dir)
  plot_top_loading_heatmap(fit, output_dir)
}

load_math500_problem_metadata <- function(
    url = "https://huggingface.co/datasets/HuggingFaceH4/MATH-500/raw/main/test.jsonl") {
  # Download the public MATH-500 JSONL metadata and align row k to problem_k.
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The jsonlite package is required to parse MATH-500 metadata.")
  }
  rows <- readLines(url, warn = FALSE)
  meta <- do.call(rbind, lapply(seq_along(rows), function(i) {
    z <- jsonlite::fromJSON(rows[[i]])
    data.frame(
      problem_index = i - 1L,
      problem = paste0("problem_", i - 1L),
      subject = z$subject,
      level = z$level,
      unique_id = z$unique_id,
      answer = z$answer,
      problem_text = z$problem,
      stringsAsFactors = FALSE
    )
  }))
  meta
}

save_problem_interpretation_tables <- function(
    output_dir,
    metadata = NULL,
    top_n = 25L) {
  # Join fitted problem loadings to MATH-500 text metadata.  This produces both
  # a full item table and compact top-positive/top-negative tables per factor.
  if (is.null(metadata)) metadata <- load_math500_problem_metadata()

  load_df <- read.csv(
    file.path(output_dir, "math500_problem_intercepts_loadings.csv"),
    check.names = FALSE
  )
  merged <- merge(load_df, metadata, by = "problem", all.x = TRUE, sort = FALSE)
  merged$problem_snippet <- gsub("[[:space:]]+", " ", merged$problem_text)
  merged$problem_snippet <- substr(merged$problem_snippet, 1L, 220L)

  write.csv(
    merged,
    file.path(output_dir, "math500_problem_loading_metadata_full.csv"),
    row.names = FALSE
  )

  loading_cols <- grep("^loading_factor_", names(merged), value = TRUE)
  top_tables <- list()
  subject_tables <- list()

  for (h in seq_along(loading_cols)) {
    col_h <- loading_cols[h]
    pos_idx <- order(merged[[col_h]], decreasing = TRUE)[seq_len(min(top_n, nrow(merged)))]
    neg_idx <- order(merged[[col_h]], decreasing = FALSE)[seq_len(min(top_n, nrow(merged)))]
    abs_idx <- order(abs(merged[[col_h]]), decreasing = TRUE)[seq_len(min(top_n, nrow(merged)))]

    one_top <- rbind(
      data.frame(direction = "positive", factor = h, merged[pos_idx, ], check.names = FALSE),
      data.frame(direction = "negative", factor = h, merged[neg_idx, ], check.names = FALSE),
      data.frame(direction = "absolute", factor = h, merged[abs_idx, ], check.names = FALSE)
    )
    top_tables[[h]] <- one_top

    subject_tables[[h]] <- do.call(rbind, lapply(c("positive", "negative", "absolute"), function(dir) {
      idx <- switch(dir, positive = pos_idx, negative = neg_idx, absolute = abs_idx)
      tab <- as.data.frame(table(merged$subject[idx]), stringsAsFactors = FALSE)
      names(tab) <- c("subject", "n")
      tab$factor <- h
      tab$direction <- dir
      tab$mean_level <- tapply(merged$level[idx], merged$subject[idx], mean)[tab$subject]
      tab[order(tab$n, decreasing = TRUE), ]
    }))
  }

  top_out <- do.call(rbind, top_tables)
  subject_out <- do.call(rbind, subject_tables)
  write.csv(
    top_out,
    file.path(output_dir, "math500_top_loading_problem_texts.csv"),
    row.names = FALSE
  )
  write.csv(
    subject_out,
    file.path(output_dir, "math500_top_loading_subject_summaries.csv"),
    row.names = FALSE
  )

  invisible(list(full = merged, top = top_out, subjects = subject_out))
}

# ----------------------------------------------------------------------------
# End-to-end Math500 run
# ----------------------------------------------------------------------------

run_math500_intercept_imfm <- function(
    data_path = file.path(script_dir, "accuracy_math500.csv"),
    output_dir = file.path(script_dir, "math500_intercept_imfm_results"),
    H_max = 10L,
    G_fixed = 3L,
    pretrain_loading_penalty = 0.02,
    pretrain_objective_tolerance = 2e-4,
    pretrain_objective_tolerance_scale = c("relative_total", "per_response"),
    H_prune_min_loading_energy_share = 0.025,
    H_prune_likelihood_penalty_multiplier = 0.25,
    H_prune_penalty_mode = c("active_log_n", "bic_effective"),
    refinement_lambda_l1_penalty = 0.25,
    refinement_mixture_prior_weight = 0.35,
    refinement_objective_tolerance = 2e-4,
    refinement_objective_tolerance_scale = c("relative_total", "per_response"),
    keep_best_binary_iterate = TRUE,
    parallel = TRUE,
    workers = NULL,
    seed = 20260717L,
    verbose = TRUE) {
  pretrain_objective_tolerance_scale <- match.arg(pretrain_objective_tolerance_scale)
  H_prune_penalty_mode <- match.arg(H_prune_penalty_mode)
  refinement_objective_tolerance_scale <- match.arg(refinement_objective_tolerance_scale)
  workers <- resolve_workers(workers)
  df <- read.csv(data_path, check.names = FALSE)
  model_id <- as.character(df[[1L]])
  X <- as.matrix(df[, -1L])
  storage.mode(X) <- "numeric"
  rownames(X) <- model_id
  colnames(X) <- paste0("problem_", colnames(X))

  pre <- fit_binary_probit_pretraining_intercept(
    X = X,
    H = NULL,
    H_max = H_max,
    G_fixed = G_fixed,
    n_aug_iter = 8L,
    z_update = "expectation",
    n_random_starts = 1L,
    max_outer = 4L,
    n_mix_starts = 3L,
    mixture_update = "map",
    mu_prior_kappa = 0.05,
    var_prior_shape = 4,
    var_prior_scale = 0.35,
    weight_prior_alpha = 1.2,
    loading_penalty = pretrain_loading_penalty,
    H_prune_after_iter = 4L,
    H_prune_loading_threshold = 0.08,
    H_prune_min_active_loadings = 10L,
    H_prune_min_loading_energy_share = H_prune_min_loading_energy_share,
    H_prune_likelihood_penalty_multiplier = H_prune_likelihood_penalty_multiplier,
    H_prune_penalty_mode = H_prune_penalty_mode,
    H_prune_min_keep = 1L,
    objective_tolerance = pretrain_objective_tolerance,
    objective_tolerance_scale = pretrain_objective_tolerance_scale,
    min_aug_iter = 5L,
    parallel = parallel,
    workers = workers,
    seed = seed,
    verbose = verbose
  )

  ref <- fit_binary_probit_refinement_intercept(
    X = X,
    pretrain_fit = pre,
    n_refine_iter = 8L,
    maxit_per_subject = 60L,
    n_mix_starts = 3L,
    min_mixture_var = 0.05,
    mixture_update = "map",
    mu_prior_kappa = 0.05,
    var_prior_shape = 4,
    var_prior_scale = 0.35,
    weight_prior_alpha = 1.2,
    mixture_prior_weight = refinement_mixture_prior_weight,
    lambda_l1_penalty = refinement_lambda_l1_penalty,
    objective_tolerance = refinement_objective_tolerance,
    objective_tolerance_scale = refinement_objective_tolerance_scale,
    min_refine_iter = 4L,
    keep_best_binary_iterate = keep_best_binary_iterate,
    parallel = parallel,
    workers = workers,
    verbose = verbose
  )
  ref$H <- ncol(ref$F_hat)
  ref <- orient_factors_by_accuracy(ref)

  saved <- save_model_outputs(ref, model_id, output_dir)
  make_all_plots(ref, saved$scores, output_dir)

  top_models <- saved$scores[order(saved$scores$accuracy, decreasing = TRUE), ][seq_len(15), ]
  write.csv(top_models, file.path(output_dir, "math500_top_15_models.csv"), row.names = FALSE)

  summary <- data.frame(
    n_models = nrow(X),
    n_problems = ncol(X),
    selected_H = ref$H,
    G_fixed_per_factor = G_fixed,
    pretrain_loading_penalty = pretrain_loading_penalty,
    pretrain_objective_tolerance = pretrain_objective_tolerance,
    pretrain_objective_tolerance_scale = pretrain_objective_tolerance_scale,
    H_prune_min_loading_energy_share = H_prune_min_loading_energy_share,
    H_prune_likelihood_penalty_multiplier = H_prune_likelihood_penalty_multiplier,
    H_prune_penalty_mode = H_prune_penalty_mode,
    refinement_lambda_l1_penalty = refinement_lambda_l1_penalty,
    refinement_mixture_prior_weight = refinement_mixture_prior_weight,
    refinement_objective_tolerance = refinement_objective_tolerance,
    refinement_objective_tolerance_scale = refinement_objective_tolerance_scale,
    refinement_selected_iteration = ref$joint_refinement$selected_iteration,
    pretraining_iterations = ref$pretraining$n_completed,
    refinement_iterations = ref$joint_refinement$n_completed,
    pretraining_final_probit_loglik = tail(ref$pretraining$history$probit_loglik, 1L),
    refinement_final_binary_loglik = tail(ref$joint_refinement$history$binary_loglik, 1L),
    refinement_selected_binary_loglik = ref$joint_refinement$selected_binary_loglik,
    mean_accuracy = mean(X),
    min_model_accuracy = min(rowMeans(X)),
    max_model_accuracy = max(rowMeans(X))
  )
  write.csv(summary, file.path(output_dir, "math500_fit_summary.csv"), row.names = FALSE)

  list(fit = ref, scores = saved$scores, summary = summary, output_dir = output_dir)
}

cmd_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
cmd_file <- if (length(cmd_file_arg) > 0L) sub("^--file=", "", cmd_file_arg[1L]) else NA_character_
this_file <- file.path(script_dir, "math500_intercept_imfm_fit.R")
is_script_execution <- !is.na(cmd_file) &&
  identical(normalizePath(cmd_file, mustWork = FALSE), normalizePath(this_file, mustWork = FALSE))
if (identical(environment(), globalenv()) && !interactive() && is_script_execution) {
  out <- run_math500_intercept_imfm()
  print(out$summary)
  cat("\nOutputs saved in:\n", out$output_dir, "\n", sep = "")
}
