#!/usr/bin/env Rscript

# ============================================================================
# Self-contained joint refinement algorithm for the binary probit
# independent-mixture factor model
#
# This file contains the actual refinement code.  It assumes the pretraining
# file has been sourced, because the refinement uses the same univariate mixture
# fitting utilities:
#
#   source("binary_probit_pretraining_algorithm_commented.R")
#   source("binary_probit_refinement_algorithm_commented.R")
#
# Refined model:
#
#   X_ij | f_i, lambda_j ~ Bernoulli(Phi(lambda_j' f_i))
#   f_ih independently follows a Gaussian mixture.
#
# Probit scale normalization:
#
#   Psi = I.
#
# Refinement objective:
#
#   L(F, Lambda, mixture params)
#     = sum_ij log Bernoulli{X_ij; Phi(lambda_j' f_i)}
#       + sum_ih log p_h(f_ih).
#
# Refinement algorithm:
#
#   1. Start from pretrained F and mixtures.
#   2. Pre-estimate Lambda by itemwise probit regressions X_j ~ 0 + F.
#   3. Repeat:
#        a. Update each f_i by L-BFGS-B under the binary probit likelihood
#           plus independent mixture log prior.
#        b. Update Lambda by itemwise probit regressions given the new F.
#        c. Refit each marginal mixture.  Either:
#             - keep G_h fixed from the current fit, or
#             - choose G_h by BIC over 1,...,G_max.
#        d. Optionally prune small or weakly separated mixture components,
#           then refit the pruned mixtures with fixed component counts.
#   4. Return refined F, Lambda, mixtures, labels, and the objective history.
# ============================================================================

script_dir <- local({
  frames <- sys.frames()
  files <- vapply(
    frames,
    function(frame) {
      if (!is.null(frame$ofile)) frame$ofile else NA_character_
    },
    character(1)
  )
  files <- files[!is.na(files)]
  if (length(files) > 0L) dirname(normalizePath(tail(files, 1L))) else getwd()
})

if (!exists("fit_binary_probit_pretraining", mode = "function") ||
    !exists("parallel_lapply", mode = "function")) {
  source(file.path(script_dir, "binary_probit_pretraining_algorithm_commented.R"))
}

# ----------------------------------------------------------------------------
# Objective functions
# ----------------------------------------------------------------------------

binary_probit_loglik <- function(X, F_hat, Lambda, alpha = NULL) {
  # Conditional binary likelihood under
  # P(X_ij = 1 | f_i) = Phi(alpha_j + lambda_j' f_i).
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  if (is.null(alpha)) alpha <- rep(0, ncol(X))
  alpha <- as.numeric(alpha)

  # Clamp probabilities away from zero to avoid log(0) in separated cases.
  eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  sum(X * log(p1) + (1 - X) * log(p0))
}

mixture_prior_loglik <- function(F_hat, mixture_fits) {
  # Independent factor prior: add the marginal mixture log density for each
  # coordinate h and subject i.
  sum(vapply(seq_len(ncol(F_hat)), function(h) {
    sum(log_dmix_1d(F_hat[, h], mixture_fits[[h]]))
  }, numeric(1)))
}

joint_binary_probit_objective <- function(
    X,
    F_hat,
    Lambda,
    alpha = NULL,
    mixture_fits,
    mixture_prior_weight = 1) {
  # Complete-data/MAP objective for the current latent scores and parameters.
  binary_probit_loglik(X, F_hat, Lambda, alpha = alpha) +
    mixture_prior_weight * mixture_prior_loglik(F_hat, mixture_fits)
}

coerce_lambda_l1_penalty <- function(lambda_l1_penalty, H) {
  penalty <- as.numeric(lambda_l1_penalty)
  if (length(penalty) == 0L || any(!is.finite(penalty))) {
    stop("lambda_l1_penalty must contain finite numeric values.")
  }
  if (length(penalty) == 1L) {
    rep(penalty, H)
  } else if (length(penalty) == H) {
    penalty
  } else {
    stop("lambda_l1_penalty must be a scalar or have length H.")
  }
}

lambda_laplace_log_prior <- function(Lambda, lambda_l1_penalty = 0) {
  # Kernel of independent Laplace priors on loading entries.  Constants are
  # omitted because they do not affect MAP optimization for fixed penalty.
  Lambda <- as.matrix(Lambda)
  penalty <- coerce_lambda_l1_penalty(lambda_l1_penalty, ncol(Lambda))
  -sum(sweep(abs(Lambda), 2L, penalty, "*"))
}

lambda_column_spike_slab_penalties <- function(
    Lambda,
    spike_penalty = 1,
    slab_penalty = 0.05,
    slab_prior = 0.5,
    effective_p = 1) {
  # Empirical-Bayes column spike/slab lasso.  Each loading column is assigned a
  # posterior slab probability under a two-Laplace mixture, then its next
  # surrogate-regression penalty is E[lambda_h | current Lambda_h].
  Lambda <- as.matrix(Lambda)
  if (!is.finite(spike_penalty) || !is.finite(slab_penalty) ||
      spike_penalty <= 0 || slab_penalty <= 0) {
    stop("Spike and slab penalties must be positive finite numbers.")
  }
  if (spike_penalty < slab_penalty) {
    stop("spike_penalty should be greater than or equal to slab_penalty.")
  }
  if (!is.finite(slab_prior) || slab_prior <= 0 || slab_prior >= 1) {
    stop("slab_prior must lie strictly between 0 and 1.")
  }
  if (!is.finite(effective_p) || effective_p < 0) {
    stop("effective_p must be a nonnegative finite number.")
  }

  col_l1 <- colSums(abs(Lambda))
  col_l2 <- sqrt(colSums(Lambda^2))
  logit_slab <- qlogis(slab_prior) +
    effective_p * log(slab_penalty / spike_penalty) +
    (spike_penalty - slab_penalty) * col_l1
  logit_slab <- pmax(pmin(logit_slab, 35), -35)
  slab_prob <- plogis(logit_slab)
  penalty <- slab_prob * slab_penalty + (1 - slab_prob) * spike_penalty

  list(
    penalty = penalty,
    slab_prob = slab_prob,
    col_l1 = col_l1,
    col_l2 = col_l2
  )
}

mixture_parameter_log_prior <- function(
    mixture_fits,
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1) {
  # Kernel of the same MAP priors used by the mixture M-step:
  #   pi_h ~ Dirichlet(alpha)
  #   sigma_hg^2 ~ Inv-Gamma(shape, scale)
  #   mu_hg | sigma_hg^2 ~ N(mu0, sigma_hg^2 / kappa0)
  # Constants independent of the fitted parameters are omitted.
  sum(vapply(mixture_fits, function(fit) {
    G <- length(fit$pi)
    alpha <- rep(weight_prior_alpha, G)
    if (length(weight_prior_alpha) == G) alpha <- weight_prior_alpha
    pi_log_prior <- sum((alpha - 1) * log(pmax(fit$pi, 1e-12)))

    var <- pmax(fit$var, 1e-12)
    var_log_prior <- sum(-(var_prior_shape + 1) * log(var) - var_prior_scale / var)
    mu_log_prior <- if (mu_prior_kappa > 0) {
      sum(-0.5 * mu_prior_kappa * (fit$mu - mu_prior_mean)^2 / var)
    } else {
      0
    }

    pi_log_prior + var_log_prior + mu_log_prior
  }, numeric(1)))
}

posterior_binary_probit_objective <- function(
    X,
    F_hat,
    Lambda,
    alpha = NULL,
    mixture_fits,
    mixture_prior_weight = 1,
    lambda_l1_penalty = 0,
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1) {
  # MAP posterior kernel:
  #   log p(X | F, Lambda)
  #   + w log p(F | mixture parameters)
  #   + log p(mixture parameters)
  #   + log p(Lambda).
  joint_binary_probit_objective(
    X,
    F_hat,
    Lambda,
    alpha = NULL,
    mixture_fits,
    alpha = alpha,
    mixture_prior_weight = mixture_prior_weight
  ) +
    mixture_parameter_log_prior(
      mixture_fits,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha
    ) +
    lambda_laplace_log_prior(Lambda, lambda_l1_penalty = lambda_l1_penalty)
}

# ----------------------------------------------------------------------------
# Loading update: itemwise probit regressions
# ----------------------------------------------------------------------------

probit_item_negloglik <- function(y, F_hat, beta, estimate_intercept = TRUE) {
  if (isTRUE(estimate_intercept)) {
    alpha <- beta[1L]
    loading <- beta[-1L]
  } else {
    alpha <- 0
    loading <- beta
  }
  eta <- as.numeric(alpha + F_hat %*% loading)
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  -sum(y * log(p1) + (1 - y) * log(p0))
}

probit_item_negloglik_grad <- function(y, F_hat, beta, estimate_intercept = TRUE) {
  if (isTRUE(estimate_intercept)) {
    alpha <- beta[1L]
    loading <- beta[-1L]
  } else {
    alpha <- 0
    loading <- beta
  }
  eta <- as.numeric(alpha + F_hat %*% loading)
  phi <- dnorm(eta)
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)

  # This is the negative of the binary probit score with respect to beta.
  score_eta <- y * phi / p1 - (1 - y) * phi / p0
  grad_loading <- as.numeric(crossprod(F_hat, score_eta))
  if (isTRUE(estimate_intercept)) {
    -c(sum(score_eta), grad_loading)
  } else {
    -grad_loading
  }
}

fit_probit_lasso_item <- function(
    y,
    F_hat,
    alpha_init = NULL,
    beta_init = NULL,
    lambda_l1_penalty = 0,
    estimate_intercept = TRUE,
    maxit = 200L,
    tol = 1e-6) {
  H <- ncol(F_hat)
  beta_loading <- if (is.null(beta_init)) rep(0, H) else as.numeric(beta_init)
  if (length(beta_loading) != H || any(!is.finite(beta_loading))) beta_loading <- rep(0, H)
  alpha <- if (isTRUE(estimate_intercept) && !is.null(alpha_init) && is.finite(alpha_init)) {
    as.numeric(alpha_init)
  } else {
    0
  }
  beta <- if (isTRUE(estimate_intercept)) c(alpha, beta_loading) else beta_loading
  lambda_l1_penalty <- coerce_lambda_l1_penalty(lambda_l1_penalty, H)
  penalty <- if (isTRUE(estimate_intercept)) c(0, lambda_l1_penalty) else lambda_l1_penalty

  objective <- function(b) {
    probit_item_negloglik(y, F_hat, b, estimate_intercept = estimate_intercept) +
      sum(penalty * abs(b))
  }

  old_obj <- objective(beta)

  for (iter in seq_len(maxit)) {
    grad <- probit_item_negloglik_grad(y, F_hat, beta, estimate_intercept = estimate_intercept)
    step <- 1
    accepted <- FALSE

    # Backtracking line search around the proximal soft-thresholding step.
    for (bt in seq_len(30L)) {
      beta_trial <- soft_threshold(beta - step * grad, step * penalty)
      trial_obj <- objective(beta_trial)

      if (is.finite(trial_obj) && trial_obj <= old_obj + 1e-10) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }

    if (!accepted) break
    max_change <- max(abs(beta_trial - beta))
    beta <- beta_trial

    if (max_change <= tol * (1 + max(abs(beta)))) break
    old_obj <- trial_obj
  }

  if (isTRUE(estimate_intercept)) {
    list(alpha = beta[1L], lambda = beta[-1L])
  } else {
    list(alpha = 0, lambda = beta)
  }
}

update_binary_probit_loadings_glm <- function(
    X,
    F_hat,
    alpha_init = NULL,
    Lambda_init = NULL,
    lambda_l1_penalty = 0,
    estimate_intercept = TRUE,
    lasso_maxit = 200L,
    lasso_tol = 1e-6,
    parallel = FALSE,
    workers = NULL) {
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  p <- ncol(X)
  H <- ncol(F_hat)
  lambda_l1_penalty <- coerce_lambda_l1_penalty(lambda_l1_penalty, H)
  Lambda <- matrix(0, p, H)
  alpha <- rep(0, p)
  colnames(Lambda) <- paste0("factor_", seq_len(H))
  rownames(Lambda) <- colnames(X)
  names(alpha) <- colnames(X)

  # Reuse the same design matrix for every item; only the binary response
  # vector changes from column to column.
  dat <- as.data.frame(F_hat)
  names(dat) <- paste0("factor_", seq_len(H))
  form <- if (isTRUE(estimate_intercept)) {
    as.formula(paste("y ~", paste(names(dat), collapse = " + ")))
  } else {
    as.formula(paste("y ~ 0 +", paste(names(dat), collapse = " + ")))
  }

  fit_one_item <- function(j) {
    beta_init <- if (!is.null(Lambda_init)) Lambda_init[j, ] else NULL
    alpha_j_init <- if (!is.null(alpha_init)) alpha_init[j] else NULL

    if (any(lambda_l1_penalty > 0)) {
      return(fit_probit_lasso_item(
        y = X[, j],
        F_hat = F_hat,
        alpha_init = alpha_j_init,
        beta_init = beta_init,
        lambda_l1_penalty = lambda_l1_penalty,
        estimate_intercept = estimate_intercept,
        maxit = lasso_maxit,
        tol = lasso_tol
      ))
    }

    dat_j <- dat
    dat_j$y <- X[, j]
    fit_j <- tryCatch(
      suppressWarnings(glm(form, data = dat_j, family = binomial(link = "probit"))),
      error = function(e) NULL
    )

    if (!is.null(fit_j) && all(is.finite(coef(fit_j)))) {
      coef_j <- coef(fit_j)
      if (isTRUE(estimate_intercept)) {
        list(alpha = unname(coef_j["(Intercept)"]), lambda = unname(coef_j[names(dat)]))
      } else {
        list(alpha = 0, lambda = unname(coef_j[names(dat)]))
      }
    } else if (!is.null(beta_init)) {
      # Separation or convergence failures can happen for binary data.  In that
      # case, keep the previous loading row rather than injecting NA values.
      list(alpha = if (!is.null(alpha_j_init)) alpha_j_init else 0, lambda = beta_init)
    } else {
      list(alpha = 0, lambda = rep(0, H))
    }
  }

  rows <- parallel_lapply(
    seq_len(p),
    fit_one_item,
    parallel = parallel,
    workers = workers
  )
  alpha[] <- vapply(rows, `[[`, numeric(1), "alpha")
  Lambda[,] <- do.call(rbind, lapply(rows, `[[`, "lambda"))

  if (isTRUE(estimate_intercept)) {
    list(alpha = alpha, Lambda = Lambda)
  } else {
    Lambda
  }
}

# ----------------------------------------------------------------------------
# Factor-score update: optimize each f_i given Lambda and mixture priors
# ----------------------------------------------------------------------------

single_factor_log_prior_grad <- function(x, fit) {
  # Compute both log p_h(x) and d/dx log p_h(x).  The gradient is a
  # responsibility-weighted average of component Gaussian score functions.
  log_components <- log(pmax(fit$pi, 1e-12)) +
    dnorm(x, mean = fit$mu, sd = sqrt(fit$var), log = TRUE)
  m <- max(log_components)
  w <- exp(log_components - m)
  resp <- w / sum(w)
  log_density <- m + log(sum(w))

  # Gradient of log sum_g pi_g phi(x; mu_g, var_g) with respect to x.
  grad <- sum(resp * (-(x - fit$mu) / fit$var))

  list(log_density = log_density, grad = grad)
}

single_factor_conditional_prior_grad <- function(x, fit, resp) {
  # Conditional/EM-style prior contribution for one factor coordinate.  The
  # responsibilities are held fixed during the current factor-score update, so
  # the mixture prior becomes a weighted Gaussian quadratic rather than a
  # log-sum of Gaussian components.
  resp <- as.numeric(resp)
  resp <- resp / sum(resp)
  var <- pmax(fit$var, 1e-8)

  log_density <- sum(resp * (
    log(pmax(fit$pi, 1e-12)) +
      dnorm(x, mean = fit$mu, sd = sqrt(var), log = TRUE)
  ))
  grad <- sum(resp * (-(x - fit$mu) / var))

  list(log_density = log_density, grad = grad)
}

mixture_responsibility_list <- function(F_hat, mixture_fits, hard = FALSE) {
  # Compute r_ihg for every factor h.  In hard mode, each row is replaced by
  # the MAP component indicator; in soft mode, the full posterior weights are
  # retained for a smoother EM-like update.
  lapply(seq_len(ncol(F_hat)), function(h) {
    resp_h <- mixture_responsibilities(F_hat[, h], mixture_fits[[h]])
    if (isTRUE(hard)) {
      z <- max.col(resp_h, ties.method = "first")
      hard_resp <- matrix(0, nrow = nrow(resp_h), ncol = ncol(resp_h))
      hard_resp[cbind(seq_len(nrow(resp_h)), z)] <- 1
      resp_h <- hard_resp
    }
    resp_h
  })
}

update_one_factor_score <- function(
    x_i,
    f_init,
    Lambda,
    alpha = NULL,
    mixture_fits,
    responsibility_i = NULL,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    mixture_prior_weight = 1,
    maxit = 50L) {
  # Optimize one subject's H-dimensional factor vector while holding Lambda and
  # the mixture prior fixed.
  factor_update <- match.arg(factor_update)
  x_i <- as.numeric(x_i)
  f_init <- as.numeric(f_init)
  Lambda <- as.matrix(Lambda)
  if (is.null(alpha)) alpha <- rep(0, length(x_i))
  alpha <- as.numeric(alpha)
  H <- length(f_init)

  objective <- function(f) {
    # Binary probit contribution for subject i.
    eta <- as.numeric(alpha + Lambda %*% f)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    ll <- sum(x_i * log(p1) + (1 - x_i) * log(p0))

    # Independent mixture prior contribution across factor coordinates.
    lp <- 0
    for (h in seq_len(H)) {
      if (factor_update == "marginal") {
        lp <- lp + single_factor_log_prior_grad(f[h], mixture_fits[[h]])$log_density
      } else {
        lp <- lp + single_factor_conditional_prior_grad(
          f[h],
          mixture_fits[[h]],
          responsibility_i[[h]]
        )$log_density
      }
    }

    -(ll + mixture_prior_weight * lp)
  }

  gradient <- function(f) {
    # Analytic gradient keeps the per-subject optimization reasonably fast.
    eta <- as.numeric(alpha + Lambda %*% f)
    phi <- dnorm(eta)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)

    # d/d eta of the Bernoulli-probit log likelihood.
    score_eta <- x_i * phi / p1 - (1 - x_i) * phi / p0
    grad <- as.numeric(crossprod(Lambda, score_eta))

    for (h in seq_len(H)) {
      # Add the hth marginal mixture-prior gradient.
      if (factor_update == "marginal") {
        grad[h] <- grad[h] + mixture_prior_weight *
          single_factor_log_prior_grad(f[h], mixture_fits[[h]])$grad
      } else {
        grad[h] <- grad[h] + mixture_prior_weight *
          single_factor_conditional_prior_grad(
            f[h],
            mixture_fits[[h]],
            responsibility_i[[h]]
          )$grad
      }
    }

    -grad
  }

  # L-BFGS-B is unconstrained here, but robust for small H and allows bounds to
  # be added later if desired.
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

update_factor_scores_joint_map <- function(
    X,
    F_hat,
    Lambda,
    alpha = NULL,
    mixture_fits,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    mixture_prior_weight = 1,
    maxit_per_subject = 50L,
    parallel = FALSE,
    workers = NULL,
    verbose = FALSE) {
  factor_update <- match.arg(factor_update)
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  if (is.null(alpha)) alpha <- rep(0, ncol(X))
  alpha <- as.numeric(alpha)
  F_new <- F_hat
  responsibility_list <- NULL
  if (factor_update != "marginal") {
    responsibility_list <- mixture_responsibility_list(
      F_hat,
      mixture_fits,
      hard = factor_update == "conditional_hard"
    )
  }

  rows <- parallel_lapply(seq_len(nrow(X)), function(i) {
    # The factor-score update separates by subject conditional on Lambda and
    # the current mixture parameters.
    responsibility_i <- NULL
    if (!is.null(responsibility_list)) {
      responsibility_i <- lapply(responsibility_list, function(resp_h) resp_h[i, ])
    }
    out_i <- update_one_factor_score(
      x_i = X[i, ],
      f_init = F_hat[i, ],
      Lambda = Lambda,
      alpha = alpha,
      mixture_fits = mixture_fits,
      responsibility_i = responsibility_i,
      factor_update = factor_update,
      mixture_prior_weight = mixture_prior_weight,
      maxit = maxit_per_subject
    )
    if (verbose && i %% 100 == 0L) {
      message("  optimized factor scores for ", i, " subjects.")
    }
    out_i
  }, parallel = parallel, workers = workers)

  F_new[,] <- do.call(rbind, rows)

  F_new
}

# ----------------------------------------------------------------------------
# Factor-scale normalization
# ----------------------------------------------------------------------------

normalize_refinement_factor_scale <- function(
    F_hat,
    Lambda,
    mixture_fits = NULL,
    target_scale = 1,
    scale_method = c("sd", "rms"),
    min_scale = 1e-6) {
  # The binary probit likelihood is unchanged by the column-wise transformation
  # F_h <- a_h F_h and Lambda_h <- Lambda_h / a_h.  We use that invariance to
  # keep refined factor coordinates on a stable scale, and transform the current
  # mixture means/variances to the same new factor scale.
  scale_method <- match.arg(scale_method)
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  H <- ncol(F_hat)

  current_scale <- vapply(seq_len(H), function(h) {
    if (scale_method == "sd") {
      sd(F_hat[, h])
    } else {
      sqrt(mean(F_hat[, h]^2))
    }
  }, numeric(1))
  current_scale[!is.finite(current_scale) | current_scale < min_scale] <- 1

  multiplier <- target_scale / current_scale
  F_new <- sweep(F_hat, 2L, multiplier, "*")
  Lambda_new <- sweep(Lambda, 2L, multiplier, "/")

  mixture_new <- mixture_fits
  if (!is.null(mixture_new)) {
    for (h in seq_len(H)) {
      mixture_new[[h]]$mu <- mixture_new[[h]]$mu * multiplier[h]
      mixture_new[[h]]$var <- pmax(mixture_new[[h]]$var * multiplier[h]^2, 1e-8)
    }
  }

  list(
    F_hat = F_new,
    Lambda = Lambda_new,
    mixture_fits = mixture_new,
    scale_before = current_scale,
    multiplier = multiplier
  )
}

normalize_refinement_factor_location <- function(
    F_hat,
    Lambda,
    alpha = NULL,
    mixture_fits = NULL) {
  # The probit index is invariant to F_h <- F_h - m_h and
  # alpha <- alpha + Lambda m.  Enforce centered refined factor coordinates so
  # item intercepts carry item difficulty rather than factor-location drift.
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  if (is.null(alpha)) alpha <- rep(0, nrow(Lambda))
  alpha <- as.numeric(alpha)

  location <- colMeans(F_hat)
  F_new <- sweep(F_hat, 2L, location, "-")
  alpha_new <- alpha + as.numeric(Lambda %*% location)

  mixture_new <- mixture_fits
  if (!is.null(mixture_new)) {
    for (h in seq_len(ncol(F_hat))) {
      mixture_new[[h]]$mu <- mixture_new[[h]]$mu - location[h]
    }
  }

  list(
    F_hat = F_new,
    alpha = alpha_new,
    mixture_fits = mixture_new,
    location = location
  )
}

# ----------------------------------------------------------------------------
# Mixture update during refinement
# ----------------------------------------------------------------------------

update_mixture_fits_fixed_G <- function(
    F_hat,
    mixture_fits,
    n_starts = 3L,
    max_iter = 20L,
    min_var = 1e-3,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    parallel = FALSE,
    workers = NULL) {
  mixture_update <- match.arg(mixture_update)
  parallel_lapply(seq_len(ncol(F_hat)), function(h) {
    # Keep the current number of components for factor h and warm-start EM from
    # the current mixture parameters.
    fit_gmm_1d(
      F_hat[, h],
      G = length(mixture_fits[[h]]$pi),
      n_starts = n_starts,
      max_iter = max_iter,
      min_var = min_var,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      init = mixture_fits[[h]]
    )
  }, parallel = parallel, workers = workers)
}

update_mixture_fits_refinement <- function(
    F_hat,
    mixture_fits,
    G_selection = c("fixed", "bic"),
    G_max = 5L,
    n_starts = 3L,
    max_iter = 20L,
    mixture_penalty_multiplier = 1,
    min_var = 1e-3,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    parallel = FALSE,
    workers = NULL) {
  G_selection <- match.arg(G_selection)
  mixture_update <- match.arg(mixture_update)

  if (G_selection == "fixed") {
    # Useful for controlled comparisons where pretraining's G_h should not
    # change during refinement.
    return(update_mixture_fits_fixed_G(
      F_hat = F_hat,
      mixture_fits = mixture_fits,
      n_starts = n_starts,
      max_iter = max_iter,
      min_var = min_var,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      parallel = parallel,
      workers = workers
    ))
  }

  parallel_lapply(seq_len(ncol(F_hat)), function(h) {
    # Let the refined factor scores correct earlier mis-clustering by
    # reselecting G_h with BIC at the current orientation.
    select_gmm_bic(
      F_hat[, h],
      G_max = G_max,
      n_starts = n_starts,
      max_iter = max_iter,
      mixture_penalty_multiplier = mixture_penalty_multiplier,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      previous_fit = mixture_fits[[h]]
    )
  }, parallel = parallel, workers = workers)
}

merge_two_refinement_components <- function(fit, a, b) {
  # Moment-match two univariate Gaussian mixture components into one component.
  # The merged fit is used as an EM initialization, not as the final estimate.
  keep <- setdiff(seq_along(fit$pi), c(a, b))
  weight <- fit$pi[a] + fit$pi[b]
  mean_ab <- (fit$pi[a] * fit$mu[a] + fit$pi[b] * fit$mu[b]) / weight
  second_ab <- (
    fit$pi[a] * (fit$var[a] + fit$mu[a]^2) +
      fit$pi[b] * (fit$var[b] + fit$mu[b]^2)
  ) / weight
  var_ab <- max(second_ab - mean_ab^2, 1e-3)

  pi_new <- c(fit$pi[keep], weight)
  mu_new <- c(fit$mu[keep], mean_ab)
  var_new <- c(fit$var[keep], var_ab)
  ord <- order(mu_new)

  list(
    pi = pi_new[ord] / sum(pi_new),
    mu = mu_new[ord],
    var = var_new[ord],
    loglik = NA_real_
  )
}

nearest_refinement_component_pair <- function(fit, g) {
  # Attach a tiny component to the nearest component by center location.
  others <- setdiff(seq_along(fit$pi), g)
  nearest <- others[which.min(abs(fit$mu[others] - fit$mu[g]))]
  sort(c(g, nearest))
}

closest_refinement_adjacent_pair <- function(fit) {
  # Separation is measured in pooled within-component standard-deviation units.
  ord <- order(fit$mu)
  fit <- list(pi = fit$pi[ord], mu = fit$mu[ord], var = fit$var[ord])
  pooled_var <- pmax((fit$var[-length(fit$var)] + fit$var[-1L]) / 2, 1e-8)
  sep <- diff(fit$mu) / sqrt(pooled_var)
  which_min <- which.min(sep)

  list(pair = c(which_min, which_min + 1L), separation = sep[which_min])
}

prune_refinement_gmm_1d <- function(
    fit,
    min_weight = 0,
    min_separation = 0,
    max_iter = 20L) {
  # Merge components that are either too small or too close to be credible
  # separate latent profiles.  Setting both thresholds to zero disables pruning.
  fit <- list(pi = fit$pi, mu = fit$mu, var = fit$var, loglik = fit$loglik)
  ord <- order(fit$mu)
  fit$pi <- fit$pi[ord] / sum(fit$pi)
  fit$mu <- fit$mu[ord]
  fit$var <- fit$var[ord]

  for (iter in seq_len(max_iter)) {
    if (length(fit$pi) <= 1L) break

    small <- which(fit$pi < min_weight)
    if (length(small) > 0L) {
      pair <- nearest_refinement_component_pair(fit, small[which.min(fit$pi[small])])
      fit <- merge_two_refinement_components(fit, pair[1], pair[2])
      next
    }

    closest <- closest_refinement_adjacent_pair(fit)
    if (is.finite(closest$separation) && closest$separation < min_separation) {
      fit <- merge_two_refinement_components(fit, closest$pair[1], closest$pair[2])
      next
    }

    break
  }

  fit
}

prune_mixture_fits_refinement <- function(
    F_hat,
    mixture_fits,
    min_weight = 0,
    min_separation = 0,
    max_iter = 20L,
    n_starts = 3L,
    mixture_max_iter = 20L,
    min_var = 1e-3,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    parallel = FALSE,
    workers = NULL) {
  mixture_update <- match.arg(mixture_update)
  # Apply the pruning rule factor by factor.  If any component count changes,
  # refit the pruned mixture with fixed G so the retained parameters are EM
  # estimates at the current factor scores rather than raw moment matches.
  G_before <- vapply(mixture_fits, function(z) length(z$pi), integer(1))
  pruned <- parallel_lapply(mixture_fits, prune_refinement_gmm_1d,
                            min_weight = min_weight,
                            min_separation = min_separation,
                            max_iter = max_iter,
                            parallel = parallel,
                            workers = workers)
  G_after_moment <- vapply(pruned, function(z) length(z$pi), integer(1))
  changed <- G_after_moment != G_before

  if (any(changed)) {
    pruned[changed] <- parallel_lapply(which(changed), function(h) {
      fit_gmm_1d(
        F_hat[, h],
        G = length(pruned[[h]]$pi),
        n_starts = n_starts,
        max_iter = mixture_max_iter,
        min_var = min_var,
        mixture_update = mixture_update,
        mu_prior_mean = mu_prior_mean,
        mu_prior_kappa = mu_prior_kappa,
        var_prior_shape = var_prior_shape,
        var_prior_scale = var_prior_scale,
        weight_prior_alpha = weight_prior_alpha,
        init = pruned[[h]]
      )
    }, parallel = parallel, workers = workers)
  }

  list(
    mixture_fits = pruned,
    G_before = G_before,
    G_after = vapply(pruned, function(z) length(z$pi), integer(1)),
    n_pruned = sum(G_before) - sum(G_after_moment)
  )
}

profile_select_one_factor_G <- function(
    h,
    X,
    F_hat,
    Lambda,
    mixture_fits,
    G_max = 5L,
    n_starts = 3L,
    mixture_max_iter = 20L,
    mixture_penalty_multiplier = 1,
    profile_mixture_penalty_multiplier = mixture_penalty_multiplier,
    lambda_l1_penalty = 0,
    lasso_maxit = 200L,
    lasso_tol = 1e-6,
    maxit_per_subject = 20L,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    mixture_prior_weight = 1,
    min_mixture_var = 1e-3,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    normalize_factor_scale = TRUE,
    normalize_factor_location = TRUE,
    factor_scale_target = 1,
    factor_scale_method = c("sd", "rms"),
    parallel_inner = FALSE,
    workers = NULL) {
  # Profile candidate G_h values through a short local update of the full data
  # model.  This is deliberately local: it starts from the current refinement
  # state and asks which candidate G_h gives the best nearby joint objective.
  factor_update <- match.arg(factor_update)
  mixture_update <- match.arg(mixture_update)
  factor_scale_method <- match.arg(factor_scale_method)
  n <- nrow(F_hat)
  candidates <- seq_len(G_max)
  best <- NULL

  for (g in candidates) {
    init_h <- if (length(mixture_fits[[h]]$pi) == g) mixture_fits[[h]] else NULL
    candidate_mixtures <- mixture_fits
    candidate_mixtures[[h]] <- fit_gmm_1d(
      F_hat[, h],
      G = g,
      n_starts = n_starts,
      max_iter = mixture_max_iter,
      min_var = min_mixture_var,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      init = init_h
    )

    candidate_F <- update_factor_scores_joint_map(
      X = X,
      F_hat = F_hat,
      Lambda = Lambda,
      alpha = alpha,
      mixture_fits = candidate_mixtures,
      factor_update = factor_update,
      mixture_prior_weight = mixture_prior_weight,
      maxit_per_subject = maxit_per_subject,
      parallel = parallel_inner,
      workers = workers,
      verbose = FALSE
    )

    if (isTRUE(normalize_factor_scale)) {
      scaled <- normalize_refinement_factor_scale(
        F_hat = candidate_F,
        Lambda = Lambda,
        mixture_fits = candidate_mixtures,
        target_scale = factor_scale_target,
        scale_method = factor_scale_method
      )
      candidate_F <- scaled$F_hat
      candidate_Lambda_init <- scaled$Lambda
      candidate_mixtures <- scaled$mixture_fits
    } else {
      candidate_Lambda_init <- Lambda
    }

    candidate_loading <- update_binary_probit_loadings_glm(
      X = X,
      F_hat = candidate_F,
      alpha_init = alpha,
      Lambda_init = candidate_Lambda_init,
      lambda_l1_penalty = lambda_l1_penalty,
      estimate_intercept = estimate_intercept,
      lasso_maxit = lasso_maxit,
      lasso_tol = lasso_tol,
      parallel = parallel_inner,
      workers = workers
    )
    candidate_alpha <- if (isTRUE(estimate_intercept)) candidate_loading$alpha else rep(0, ncol(X))
    candidate_Lambda <- if (isTRUE(estimate_intercept)) candidate_loading$Lambda else candidate_loading

    # Refit mixtures at the profiled factor scores while keeping all candidate
    # component counts fixed.  This lets the prior term respond to the data-fit
    # update without doing a full restart of the global algorithm.
    candidate_mixtures <- update_mixture_fits_fixed_G(
      F_hat = candidate_F,
      mixture_fits = candidate_mixtures,
      n_starts = n_starts,
      max_iter = mixture_max_iter,
      min_var = min_mixture_var,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      parallel = parallel_inner,
      workers = workers
    )

    df_h <- 3 * g - 1
    penalty_h <- profile_mixture_penalty_multiplier * 0.5 * df_h * log(n)
    joint <- joint_binary_probit_objective(
      X = X,
      F_hat = candidate_F,
      Lambda = candidate_Lambda,
      alpha = candidate_alpha,
      mixture_fits = candidate_mixtures,
      mixture_prior_weight = mixture_prior_weight
    )
    score <- joint - penalty_h

    if (is.null(best) || score > best$profile_score) {
      best <- list(
        G = g,
        profile_score = score,
        joint_objective = joint,
        penalty = penalty_h,
        F_hat = candidate_F,
        alpha = candidate_alpha,
        Lambda = candidate_Lambda,
        mixture_fits = candidate_mixtures
      )
    }
  }

  best
}

profile_select_mixture_fits_refinement <- function(
    X,
    F_hat,
    Lambda,
    mixture_fits,
    G_max = 5L,
    n_starts = 3L,
    mixture_max_iter = 20L,
    mixture_penalty_multiplier = 1,
    profile_mixture_penalty_multiplier = mixture_penalty_multiplier,
    lambda_l1_penalty = 0,
    lasso_maxit = 200L,
    lasso_tol = 1e-6,
    maxit_per_subject = 20L,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    mixture_prior_weight = 1,
    min_mixture_var = 1e-3,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    normalize_factor_scale = TRUE,
    factor_scale_target = 1,
    factor_scale_method = c("sd", "rms"),
    parallel = FALSE,
    workers = NULL) {
  factor_update <- match.arg(factor_update)
  mixture_update <- match.arg(mixture_update)
  factor_scale_method <- match.arg(factor_scale_method)
  profile_history <- vector("list", ncol(F_hat))

  for (h in seq_len(ncol(F_hat))) {
    selected <- profile_select_one_factor_G(
      h = h,
      X = X,
      F_hat = F_hat,
      Lambda = Lambda,
      alpha = alpha,
      mixture_fits = mixture_fits,
      G_max = G_max,
      n_starts = n_starts,
      mixture_max_iter = mixture_max_iter,
      mixture_penalty_multiplier = mixture_penalty_multiplier,
      profile_mixture_penalty_multiplier = profile_mixture_penalty_multiplier,
      lambda_l1_penalty = lambda_l1_penalty,
      estimate_intercept = estimate_intercept,
      lasso_maxit = lasso_maxit,
      lasso_tol = lasso_tol,
      maxit_per_subject = maxit_per_subject,
      factor_update = factor_update,
      mixture_prior_weight = mixture_prior_weight,
      min_mixture_var = min_mixture_var,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      normalize_factor_scale = normalize_factor_scale,
      factor_scale_target = factor_scale_target,
      factor_scale_method = factor_scale_method,
      parallel_inner = parallel,
      workers = workers
    )

    # Commit the locally profiled state before moving to the next factor, giving
    # later factors the benefit of earlier component-count decisions.
    F_hat <- selected$F_hat
    alpha <- selected$alpha
    Lambda <- selected$Lambda
    mixture_fits <- selected$mixture_fits
    profile_history[[h]] <- data.frame(
      factor = h,
      selected_G = selected$G,
      profile_score = selected$profile_score,
      joint_objective = selected$joint_objective,
      penalty = selected$penalty
    )
  }

  list(
    F_hat = F_hat,
    alpha = alpha,
    Lambda = Lambda,
    mixture_fits = mixture_fits,
    profile_history = do.call(rbind, profile_history)
  )
}

# ----------------------------------------------------------------------------
# Main joint refinement algorithm
# ----------------------------------------------------------------------------

fit_binary_probit_refinement <- function(
    X,
    pretrain_fit,
    n_refine_iter = 5L,
    maxit_per_subject = 50L,
    n_mix_starts = 3L,
    mixture_max_iter = 20L,
    G_selection = c("fixed", "bic", "profile_bic"),
    G_max = NULL,
    mixture_penalty_multiplier = 1,
    profile_mixture_penalty_multiplier = NULL,
    profile_G_maxit_per_subject = 10L,
    profile_G_every = 1L,
    prune_mixtures = FALSE,
    pruning_min_weight = 0,
    pruning_min_separation = 0,
    pruning_max_iter = 20L,
    pruning_start_iter = 1L,
    pruning_every = 1L,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    min_mixture_var = 1e-3,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    mixture_prior_weight = 1,
    estimate_intercept = TRUE,
    preestimate_loadings = TRUE,
    lambda_l1_penalty = 0,
    lambda_column_spike_slab = FALSE,
    lambda_spike_penalty = 1,
    lambda_slab_penalty = 0.05,
    lambda_slab_prior = 0.5,
    lambda_column_effective_p = 1,
    lasso_maxit = 200L,
    lasso_tol = 1e-6,
    normalize_factor_scale = TRUE,
    normalize_factor_location = TRUE,
    factor_scale_target = 1,
    factor_scale_method = c("sd", "rms"),
    objective_tolerance = 1e-5,
    min_refine_iter = 1L,
    stopping_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    return_best_refinement_iteration = FALSE,
    refinement_selection_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    store_refinement_trace = FALSE,
    parallel = FALSE,
    workers = NULL,
    verbose = TRUE,
    ...) {
  G_selection <- match.arg(G_selection)
  factor_update <- match.arg(factor_update)
  mixture_update <- match.arg(mixture_update)
  factor_scale_method <- match.arg(factor_scale_method)
  stopping_objective <- match.arg(stopping_objective)
  refinement_selection_objective <- match.arg(refinement_selection_objective)
  if (is.null(profile_mixture_penalty_multiplier)) {
    profile_mixture_penalty_multiplier <- mixture_penalty_multiplier
  }
  X <- as.matrix(X)
  workers <- resolve_workers(workers)
  min_refine_iter <- as.integer(min_refine_iter)
  prune_mixtures <- isTRUE(prune_mixtures)
  pruning_start_iter <- as.integer(pruning_start_iter)
  pruning_every <- as.integer(pruning_every)
  mixture_prior_weight <- as.numeric(mixture_prior_weight)
  estimate_intercept <- isTRUE(estimate_intercept)
  if (!is.finite(mixture_prior_weight) || mixture_prior_weight < 0) {
    stop("mixture_prior_weight must be a nonnegative finite number.")
  }
  lambda_column_spike_slab <- isTRUE(lambda_column_spike_slab)

  # Initialize refinement from the output of the pretraining algorithm.
  F_hat <- pretrain_fit$F_hat
  Lambda <- pretrain_fit$Lambda_hat
  alpha <- if (estimate_intercept && !is.null(pretrain_fit$alpha_hat)) {
    as.numeric(pretrain_fit$alpha_hat)
  } else {
    rep(0, ncol(X))
  }
  mixture_fits <- pretrain_fit$mixture_fits

  lambda_penalty_info <- if (lambda_column_spike_slab) {
    lambda_column_spike_slab_penalties(
      Lambda = Lambda,
      spike_penalty = lambda_spike_penalty,
      slab_penalty = lambda_slab_penalty,
      slab_prior = lambda_slab_prior,
      effective_p = lambda_column_effective_p
    )
  } else {
    list(
      penalty = coerce_lambda_l1_penalty(lambda_l1_penalty, ncol(Lambda)),
      slab_prob = rep(NA_real_, ncol(Lambda)),
      col_l1 = colSums(abs(Lambda)),
      col_l2 = sqrt(colSums(Lambda^2))
    )
  }
  lambda_penalty_history <- list()

  if (is.null(G_max)) {
    # If no search limit is supplied, use the largest pretrained component
    # count.  Passing a larger G_max allows new splits during BIC refinement.
    G_max <- max(vapply(mixture_fits, function(z) length(z$pi), integer(1)))
  }

  if (isTRUE(preestimate_loadings)) {
    # This is the requested ordering: finish pretraining, then estimate Lambda
    # given the final pretrained factors before joint refinement begins.
    if (verbose) message("Pre-estimating probit loadings before joint refinement.")
    loading_fit <- update_binary_probit_loadings_glm(
      X = X,
      F_hat = F_hat,
      alpha_init = alpha,
      Lambda_init = Lambda,
      lambda_l1_penalty = lambda_penalty_info$penalty,
      estimate_intercept = estimate_intercept,
      lasso_maxit = lasso_maxit,
      lasso_tol = lasso_tol,
      parallel = parallel,
      workers = workers
    )
    if (estimate_intercept) {
      alpha <- loading_fit$alpha
      Lambda <- loading_fit$Lambda
    } else {
      Lambda <- loading_fit
      alpha <- rep(0, ncol(X))
    }
    if (lambda_column_spike_slab) {
      lambda_penalty_info <- lambda_column_spike_slab_penalties(
        Lambda = Lambda,
        spike_penalty = lambda_spike_penalty,
        slab_penalty = lambda_slab_penalty,
        slab_prior = lambda_slab_prior,
        effective_p = lambda_column_effective_p
      )
    }
  }

  if (isTRUE(normalize_factor_location)) {
    located <- normalize_refinement_factor_location(
      F_hat = F_hat,
      Lambda = Lambda,
      alpha = alpha,
      mixture_fits = mixture_fits
    )
    F_hat <- located$F_hat
    alpha <- located$alpha
    mixture_fits <- located$mixture_fits
  }

  # Iteration 0 is the objective after loading pre-estimation but before any
  # joint updates of F, Lambda, or mixture parameters.
  history <- vector("list", n_refine_iter + 1L)
  initial_binary_loglik <- binary_probit_loglik(X, F_hat, Lambda, alpha = alpha)
  initial_mixture_loglik <- mixture_prior_loglik(F_hat, mixture_fits)
  initial_mixture_parameter_logprior <- mixture_parameter_log_prior(
    mixture_fits,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha
  )
  initial_lambda_logprior <- lambda_laplace_log_prior(
    Lambda,
    lambda_l1_penalty = lambda_penalty_info$penalty
  )
  lambda_penalty_history[[1L]] <- data.frame(
    iteration = 0L,
    factor = seq_len(ncol(Lambda)),
    column_l1 = lambda_penalty_info$col_l1,
    column_l2 = lambda_penalty_info$col_l2,
    slab_prob = lambda_penalty_info$slab_prob,
    l1_penalty = lambda_penalty_info$penalty
  )
  initial_joint_objective <- initial_binary_loglik +
    mixture_prior_weight * initial_mixture_loglik
  initial_posterior_objective <- initial_joint_objective +
    initial_mixture_parameter_logprior +
    initial_lambda_logprior
  history[[1L]] <- data.frame(
    iteration = 0L,
    binary_loglik = initial_binary_loglik,
    mixture_loglik = initial_mixture_loglik,
    mixture_parameter_logprior = initial_mixture_parameter_logprior,
    lambda_logprior = initial_lambda_logprior,
    joint_objective = initial_joint_objective,
    posterior_objective = initial_posterior_objective,
    objective_improvement = NA_real_,
    relative_objective_improvement = NA_real_,
    relative_objective_change = NA_real_,
    G_hat_before_prune = paste(vapply(mixture_fits, function(z) length(z$pi), integer(1)),
                               collapse = ","),
    n_pruned_components = 0L,
    max_factor_scale_before_normalize = max(apply(F_hat, 2L, sd)),
    max_factor_scale_after_normalize = max(apply(F_hat, 2L, sd)),
    G_hat = paste(vapply(mixture_fits, function(z) length(z$pi), integer(1)),
                  collapse = ","),
    factor_update_seconds = NA_real_,
    location_normalize_seconds = NA_real_,
    normalize_seconds = NA_real_,
    lambda_update_seconds = NA_real_,
    mixture_update_seconds = NA_real_,
    pruning_seconds = NA_real_,
    objective_seconds = NA_real_,
    iteration_seconds = NA_real_
  )

  converged <- FALSE
  n_completed <- 0L
  profile_G_history <- list()
  refinement_trace <- NULL
  keep_refinement_trace <- isTRUE(store_refinement_trace) || isTRUE(return_best_refinement_iteration)
  if (isTRUE(keep_refinement_trace)) {
    refinement_trace <- vector("list", n_refine_iter + 1L)
    refinement_trace[[1L]] <- list(
      iteration = 0L,
      F_hat = F_hat,
      alpha_hat = alpha,
      Lambda_hat = Lambda,
      mixture_fits = mixture_fits,
      G_hat = vapply(mixture_fits, function(z) length(z$pi), integer(1))
    )
  }

  for (iter in seq_len(n_refine_iter)) {
    iter_start <- Sys.time()
    if (verbose) message("Joint refinement iteration ", iter, " of ", n_refine_iter, ".")

    # Step 1: update latent factor scores subject by subject.
    factor_update_start <- Sys.time()
    F_hat <- update_factor_scores_joint_map(
      X = X,
      F_hat = F_hat,
      Lambda = Lambda,
      alpha = alpha,
      mixture_fits = mixture_fits,
      factor_update = factor_update,
      mixture_prior_weight = mixture_prior_weight,
      maxit_per_subject = maxit_per_subject,
      parallel = parallel,
      workers = workers,
      verbose = FALSE
    )
    factor_update_seconds <- as.numeric(difftime(Sys.time(), factor_update_start, units = "secs"))

    location_normalize_start <- Sys.time()
    if (isTRUE(normalize_factor_location)) {
      located <- normalize_refinement_factor_location(
        F_hat = F_hat,
        Lambda = Lambda,
        alpha = alpha,
        mixture_fits = mixture_fits
      )
      F_hat <- located$F_hat
      alpha <- located$alpha
      mixture_fits <- located$mixture_fits
    }
    location_normalize_seconds <- as.numeric(difftime(Sys.time(), location_normalize_start, units = "secs"))

    normalize_start <- Sys.time()
    scale_before_normalize <- apply(F_hat, 2L, sd)
    scale_after_normalize <- scale_before_normalize
    if (isTRUE(normalize_factor_scale)) {
      scaled <- normalize_refinement_factor_scale(
        F_hat = F_hat,
        Lambda = Lambda,
        mixture_fits = mixture_fits,
        target_scale = factor_scale_target,
        scale_method = factor_scale_method
      )
      F_hat <- scaled$F_hat
      Lambda <- scaled$Lambda
      mixture_fits <- scaled$mixture_fits
      scale_before_normalize <- scaled$scale_before
      scale_after_normalize <- apply(F_hat, 2L, sd)
    }
    normalize_seconds <- as.numeric(difftime(Sys.time(), normalize_start, units = "secs"))

    # Step 2: update loadings given the refined factors.
    lambda_update_start <- Sys.time()
    if (lambda_column_spike_slab) {
      lambda_penalty_info <- lambda_column_spike_slab_penalties(
        Lambda = Lambda,
        spike_penalty = lambda_spike_penalty,
        slab_penalty = lambda_slab_penalty,
        slab_prior = lambda_slab_prior,
        effective_p = lambda_column_effective_p
      )
    }
    loading_fit <- update_binary_probit_loadings_glm(
      X = X,
      F_hat = F_hat,
      alpha_init = alpha,
      Lambda_init = Lambda,
      lambda_l1_penalty = lambda_penalty_info$penalty,
      estimate_intercept = estimate_intercept,
      lasso_maxit = lasso_maxit,
      lasso_tol = lasso_tol,
      parallel = parallel,
      workers = workers
    )
    if (estimate_intercept) {
      alpha <- loading_fit$alpha
      Lambda <- loading_fit$Lambda
    } else {
      Lambda <- loading_fit
      alpha <- rep(0, ncol(X))
    }
    if (lambda_column_spike_slab) {
      lambda_penalty_info <- lambda_column_spike_slab_penalties(
        Lambda = Lambda,
        spike_penalty = lambda_spike_penalty,
        slab_penalty = lambda_slab_penalty,
        slab_prior = lambda_slab_prior,
        effective_p = lambda_column_effective_p
      )
    }
    lambda_update_seconds <- as.numeric(difftime(Sys.time(), lambda_update_start, units = "secs"))

    # Step 3: update marginal mixture profiles.  The usual BIC mode selects G_h
    # from the marginal factor distribution only.  The profiled mode lets the
    # full binary-probit data model respond locally before scoring candidates.
    mixture_update_start <- Sys.time()
    if (G_selection == "profile_bic" && iter %% as.integer(profile_G_every) == 0L) {
      profiled <- profile_select_mixture_fits_refinement(
        X = X,
        F_hat = F_hat,
        Lambda = Lambda,
        alpha = alpha,
        mixture_fits = mixture_fits,
        G_max = G_max,
        n_starts = n_mix_starts,
        mixture_max_iter = mixture_max_iter,
        mixture_penalty_multiplier = mixture_penalty_multiplier,
        profile_mixture_penalty_multiplier = profile_mixture_penalty_multiplier,
        lambda_l1_penalty = lambda_penalty_info$penalty,
        lasso_maxit = lasso_maxit,
        lasso_tol = lasso_tol,
        maxit_per_subject = profile_G_maxit_per_subject,
        factor_update = factor_update,
        mixture_prior_weight = mixture_prior_weight,
        min_mixture_var = min_mixture_var,
        mixture_update = mixture_update,
        mu_prior_mean = mu_prior_mean,
        mu_prior_kappa = mu_prior_kappa,
        var_prior_shape = var_prior_shape,
        var_prior_scale = var_prior_scale,
        weight_prior_alpha = weight_prior_alpha,
        normalize_factor_scale = normalize_factor_scale,
        factor_scale_target = factor_scale_target,
        factor_scale_method = factor_scale_method,
        parallel = parallel,
        workers = workers
      )
      F_hat <- profiled$F_hat
      alpha <- profiled$alpha
      Lambda <- profiled$Lambda
      mixture_fits <- profiled$mixture_fits
      profiled$profile_history$iteration <- iter
      profile_G_history[[length(profile_G_history) + 1L]] <- profiled$profile_history
    } else {
      mixture_fits <- update_mixture_fits_refinement(
        F_hat = F_hat,
        mixture_fits = mixture_fits,
        G_selection = if (G_selection == "profile_bic") "fixed" else G_selection,
        G_max = G_max,
        n_starts = n_mix_starts,
        max_iter = mixture_max_iter,
        mixture_penalty_multiplier = mixture_penalty_multiplier,
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
    }
    mixture_update_seconds <- as.numeric(difftime(Sys.time(), mixture_update_start, units = "secs"))

    G_before_prune <- vapply(mixture_fits, function(z) length(z$pi), integer(1))
    n_pruned_components <- 0L
    pruning_seconds <- 0
    if (prune_mixtures &&
        iter >= pruning_start_iter &&
        pruning_every > 0L &&
        iter %% pruning_every == 0L) {
      pruning_start <- Sys.time()
      pruned <- prune_mixture_fits_refinement(
        F_hat = F_hat,
        mixture_fits = mixture_fits,
        min_weight = pruning_min_weight,
        min_separation = pruning_min_separation,
        max_iter = pruning_max_iter,
        n_starts = n_mix_starts,
        mixture_max_iter = mixture_max_iter,
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
      mixture_fits <- pruned$mixture_fits
      n_pruned_components <- pruned$n_pruned
      pruning_seconds <- as.numeric(difftime(Sys.time(), pruning_start, units = "secs"))
    }

    # Track both pieces of the objective so we can see whether improvements
    # come from binary prediction, better mixture fit, or both.
    objective_start <- Sys.time()
    current_binary_loglik <- binary_probit_loglik(X, F_hat, Lambda, alpha = alpha)
    current_mixture_loglik <- mixture_prior_loglik(F_hat, mixture_fits)
    current_mixture_parameter_logprior <- mixture_parameter_log_prior(
      mixture_fits,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha
    )
    current_lambda_logprior <- lambda_laplace_log_prior(
      Lambda,
      lambda_l1_penalty = lambda_penalty_info$penalty
    )
    lambda_penalty_history[[iter + 1L]] <- data.frame(
      iteration = iter,
      factor = seq_len(ncol(Lambda)),
      column_l1 = lambda_penalty_info$col_l1,
      column_l2 = lambda_penalty_info$col_l2,
      slab_prob = lambda_penalty_info$slab_prob,
      l1_penalty = lambda_penalty_info$penalty
    )
    current_objective <- current_binary_loglik +
      mixture_prior_weight * current_mixture_loglik
    current_posterior_objective <- current_objective +
      current_mixture_parameter_logprior +
      current_lambda_logprior

    previous_stopping_value <- history[[iter]][[stopping_objective]]
    current_stopping_value <- switch(
      stopping_objective,
      binary_loglik = current_binary_loglik,
      mixture_loglik = current_mixture_loglik,
      joint_objective = current_objective,
      posterior_objective = current_posterior_objective
    )
    objective_improvement <- current_stopping_value - previous_stopping_value
    relative_improvement <- objective_improvement / (1 + abs(previous_stopping_value))
    objective_seconds <- as.numeric(difftime(Sys.time(), objective_start, units = "secs"))
    iteration_seconds <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

    history[[iter + 1L]] <- data.frame(
      iteration = iter,
      binary_loglik = current_binary_loglik,
      mixture_loglik = current_mixture_loglik,
      mixture_parameter_logprior = current_mixture_parameter_logprior,
      lambda_logprior = current_lambda_logprior,
      joint_objective = current_objective,
      posterior_objective = current_posterior_objective,
      objective_improvement = objective_improvement,
      relative_objective_improvement = relative_improvement,
      # Backward-compatible column name; this is now signed improvement, not
      # absolute change.
      relative_objective_change = relative_improvement,
      G_hat_before_prune = paste(G_before_prune, collapse = ","),
      n_pruned_components = n_pruned_components,
      max_factor_scale_before_normalize = max(scale_before_normalize),
      max_factor_scale_after_normalize = max(scale_after_normalize),
      G_hat = paste(vapply(mixture_fits, function(z) length(z$pi), integer(1)),
                    collapse = ","),
      factor_update_seconds = factor_update_seconds,
      location_normalize_seconds = location_normalize_seconds,
      normalize_seconds = normalize_seconds,
      lambda_update_seconds = lambda_update_seconds,
      mixture_update_seconds = mixture_update_seconds,
      pruning_seconds = pruning_seconds,
      objective_seconds = objective_seconds,
      iteration_seconds = iteration_seconds
    )

    if (isTRUE(keep_refinement_trace)) {
      refinement_trace[[iter + 1L]] <- list(
        iteration = iter,
        F_hat = F_hat,
        alpha_hat = alpha,
        Lambda_hat = Lambda,
        mixture_fits = mixture_fits,
        G_hat = vapply(mixture_fits, function(z) length(z$pi), integer(1))
      )
    }

    n_completed <- iter
    if (!is.null(objective_tolerance) &&
        is.finite(objective_tolerance) &&
        iter >= min_refine_iter &&
        is.finite(relative_improvement) &&
        relative_improvement <= objective_tolerance) {
      converged <- TRUE
      if (verbose) {
        message(
          "Stopping refinement: relative ",
          stopping_objective,
          " improvement ",
          signif(relative_improvement, 3),
          " <= ",
          objective_tolerance,
          "."
        )
      }
      break
    }
  }

  history_out <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  trace_out <- if (isTRUE(keep_refinement_trace)) {
    refinement_trace[!vapply(refinement_trace, is.null, logical(1))]
  } else {
    NULL
  }

  selected_refinement_iteration <- tail(history_out$iteration, 1L)
  if (isTRUE(return_best_refinement_iteration)) {
    best_row <- which.max(history_out[[refinement_selection_objective]])
    selected_refinement_iteration <- history_out$iteration[best_row]
    selected_snapshot <- trace_out[[which(vapply(
      trace_out,
      function(z) identical(as.integer(z$iteration), as.integer(selected_refinement_iteration)),
      logical(1)
    ))[1L]]]

    F_hat <- selected_snapshot$F_hat
    alpha <- selected_snapshot$alpha_hat
    Lambda <- selected_snapshot$Lambda_hat
    mixture_fits <- selected_snapshot$mixture_fits
    lambda_penalty_info <- if (lambda_column_spike_slab) {
      lambda_column_spike_slab_penalties(
        Lambda = Lambda,
        spike_penalty = lambda_spike_penalty,
        slab_penalty = lambda_slab_penalty,
        slab_prior = lambda_slab_prior,
        effective_p = lambda_column_effective_p
      )
    } else {
      list(
        penalty = coerce_lambda_l1_penalty(lambda_l1_penalty, ncol(Lambda)),
        slab_prob = rep(NA_real_, ncol(Lambda)),
        col_l1 = colSums(abs(Lambda)),
        col_l2 = sqrt(colSums(Lambda^2))
      )
    }
  }

  # Convert the returned marginal mixture responsibilities into hard labels.
  responsibilities <- lapply(seq_len(ncol(F_hat)), function(h) {
    mixture_responsibilities(F_hat[, h], mixture_fits[[h]])
  })
  class_map <- sapply(responsibilities, max.col, ties.method = "first")
  if (ncol(F_hat) == 1L) class_map <- matrix(class_map, ncol = 1L)
  colnames(class_map) <- paste0("factor_", seq_len(ncol(F_hat)))

  # Return the same object shape as pretraining, but with refined quantities.
  out <- pretrain_fit
  out$model <- "joint_refined_binary_probit_independent_mixture_factor"
  out$F_hat <- F_hat
  out$alpha_hat <- alpha
  out$Lambda_hat <- Lambda
  out$Psi_hat <- diag(1, ncol(X))
  out$mixture_fits <- mixture_fits
  out$G_hat <- vapply(mixture_fits, function(z) length(z$pi), integer(1))
  out$responsibilities <- responsibilities
  out$class_map <- class_map
  out$profile_id <- binary_profile_id(class_map)
  out$joint_refinement <- list(
    n_refine_iter = n_refine_iter,
    n_completed = n_completed,
    converged = converged,
    objective_tolerance = objective_tolerance,
    min_refine_iter = min_refine_iter,
    stopping_objective = stopping_objective,
    return_best_refinement_iteration = isTRUE(return_best_refinement_iteration),
    refinement_selection_objective = refinement_selection_objective,
    selected_refinement_iteration = selected_refinement_iteration,
    maxit_per_subject = maxit_per_subject,
    mixture_max_iter = mixture_max_iter,
    G_selection = G_selection,
    G_max = G_max,
    mixture_penalty_multiplier = mixture_penalty_multiplier,
    profile_mixture_penalty_multiplier = profile_mixture_penalty_multiplier,
    profile_G_maxit_per_subject = profile_G_maxit_per_subject,
    profile_G_every = profile_G_every,
    prune_mixtures = prune_mixtures,
    pruning_min_weight = pruning_min_weight,
    pruning_min_separation = pruning_min_separation,
    pruning_max_iter = pruning_max_iter,
    pruning_start_iter = pruning_start_iter,
    pruning_every = pruning_every,
    factor_update = factor_update,
    min_mixture_var = min_mixture_var,
    mixture_update = mixture_update,
    mixture_prior = list(
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha
    ),
    mixture_prior_weight = mixture_prior_weight,
    estimate_intercept = estimate_intercept,
    profile_G_history = if (length(profile_G_history)) do.call(rbind, profile_G_history) else NULL,
    lambda_l1_penalty = lambda_l1_penalty,
    lambda_column_spike_slab = list(
      enabled = lambda_column_spike_slab,
      spike_penalty = lambda_spike_penalty,
      slab_penalty = lambda_slab_penalty,
      slab_prior = lambda_slab_prior,
      effective_p = lambda_column_effective_p,
      final_column_l1 = colSums(abs(Lambda)),
      final_column_l2 = sqrt(colSums(Lambda^2)),
      final_slab_prob = lambda_penalty_info$slab_prob,
      final_l1_penalty = lambda_penalty_info$penalty,
      history = if (length(lambda_penalty_history)) do.call(rbind, lambda_penalty_history) else NULL
    ),
    lasso_maxit = lasso_maxit,
    lasso_tol = lasso_tol,
    normalize_factor_scale = normalize_factor_scale,
    normalize_factor_location = isTRUE(normalize_factor_location),
    factor_scale_target = factor_scale_target,
    factor_scale_method = factor_scale_method,
    preestimate_loadings = isTRUE(preestimate_loadings),
    store_refinement_trace = isTRUE(store_refinement_trace),
    refinement_trace = if (isTRUE(store_refinement_trace)) {
      trace_out
    } else {
      NULL
    },
    history = history_out,
    psi = "identity",
    parallel = list(enabled = isTRUE(parallel), workers = workers)
  )
  out
}

# ----------------------------------------------------------------------------
# Convenience wrapper: pretrain, then refine
# ----------------------------------------------------------------------------

fit_binary_probit_pretrain_then_refine <- function(
    X,
    H = NULL,
    H_max = min(10L, nrow(as.matrix(X)) - 1L, ncol(as.matrix(X))),
    G_max = 5L,
    pretrain_G_selection = c("bic", "fixed"),
    G_fixed = NULL,
    n_aug_iter = 4L,
    z_update = c("sample", "expectation"),
    n_refine_iter = 5L,
    refine_G_selection = c("fixed", "bic"),
    refine_G_max = G_max,
    mixture_penalty_multiplier = 1,
    profile_mixture_penalty_multiplier = NULL,
    prune_mixtures = FALSE,
    pruning_min_weight = 0,
    pruning_min_separation = 0,
    pruning_max_iter = 20L,
    pruning_start_iter = 1L,
    pruning_every = 1L,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    min_mixture_var = 1e-3,
    mixture_prior_weight = 1,
    lambda_l1_penalty = 0,
    lambda_column_spike_slab = FALSE,
    lambda_spike_penalty = 1,
    lambda_slab_penalty = 0.05,
    lambda_slab_prior = 0.5,
    lambda_column_effective_p = 1,
    normalize_factor_scale = TRUE,
    normalize_factor_location = TRUE,
    factor_scale_target = 1,
    factor_scale_method = c("sd", "rms"),
    objective_tolerance = 1e-5,
    min_refine_iter = 1L,
    stopping_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    return_best_refinement_iteration = FALSE,
    refinement_selection_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    parallel = FALSE,
    workers = NULL,
    seed = 20260715L,
    verbose = TRUE,
    ...) {
  z_update <- match.arg(z_update)
  pretrain_G_selection <- match.arg(pretrain_G_selection)
  refine_G_selection <- match.arg(refine_G_selection)
  factor_update <- match.arg(factor_update)
  factor_scale_method <- match.arg(factor_scale_method)
  stopping_objective <- match.arg(stopping_objective)
  refinement_selection_objective <- match.arg(refinement_selection_objective)

  pretrain_fit <- fit_binary_probit_pretraining(
    X = X,
    H = H,
    H_max = H_max,
    pretrain_G_selection = pretrain_G_selection,
    G_max = G_max,
    G_fixed = G_fixed,
    n_aug_iter = n_aug_iter,
    z_update = z_update,
    mixture_penalty_multiplier = mixture_penalty_multiplier,
    parallel = parallel,
    workers = workers,
    seed = seed,
    verbose = verbose,
    ...
  )

  refine_fit <- fit_binary_probit_refinement(
    X = X,
    pretrain_fit = pretrain_fit,
    n_refine_iter = n_refine_iter,
    G_selection = refine_G_selection,
    G_max = refine_G_max,
    mixture_penalty_multiplier = mixture_penalty_multiplier,
    profile_mixture_penalty_multiplier = profile_mixture_penalty_multiplier,
    prune_mixtures = prune_mixtures,
    pruning_min_weight = pruning_min_weight,
    pruning_min_separation = pruning_min_separation,
    pruning_max_iter = pruning_max_iter,
    pruning_start_iter = pruning_start_iter,
    pruning_every = pruning_every,
    factor_update = factor_update,
    min_mixture_var = min_mixture_var,
    mixture_prior_weight = mixture_prior_weight,
    lambda_l1_penalty = lambda_l1_penalty,
    lambda_column_spike_slab = lambda_column_spike_slab,
    lambda_spike_penalty = lambda_spike_penalty,
    lambda_slab_penalty = lambda_slab_penalty,
    lambda_slab_prior = lambda_slab_prior,
    lambda_column_effective_p = lambda_column_effective_p,
    normalize_factor_scale = normalize_factor_scale,
    normalize_factor_location = normalize_factor_location,
    factor_scale_target = factor_scale_target,
    factor_scale_method = factor_scale_method,
    objective_tolerance = objective_tolerance,
    min_refine_iter = min_refine_iter,
    stopping_objective = stopping_objective,
    return_best_refinement_iteration = return_best_refinement_iteration,
    refinement_selection_objective = refinement_selection_objective,
    parallel = parallel,
    workers = workers,
    verbose = verbose,
    ...
  )

  list(pretrain_fit = pretrain_fit, refine_fit = refine_fit)
}

summarize_binary_probit_refinement <- function(fit) {
  hist <- fit$joint_refinement$history

  data.frame(
    stage = "joint_refinement",
    H_hat = fit$H,
    G_hat = paste(fit$G_hat, collapse = ","),
    G_selection = fit$joint_refinement$G_selection,
    psi = fit$joint_refinement$psi,
    lambda_l1_penalty = fit$joint_refinement$lambda_l1_penalty,
    converged = fit$joint_refinement$converged,
    n_completed = fit$joint_refinement$n_completed,
    objective_tolerance = fit$joint_refinement$objective_tolerance,
    prune_mixtures = fit$joint_refinement$prune_mixtures,
    pruning_min_weight = fit$joint_refinement$pruning_min_weight,
    pruning_min_separation = fit$joint_refinement$pruning_min_separation,
    pruning_start_iter = fit$joint_refinement$pruning_start_iter,
    pruning_every = fit$joint_refinement$pruning_every,
    normalize_factor_scale = fit$joint_refinement$normalize_factor_scale,
    factor_scale_target = fit$joint_refinement$factor_scale_target,
    factor_scale_method = fit$joint_refinement$factor_scale_method,
    initial_joint_objective = hist$joint_objective[1L],
    final_joint_objective = tail(hist$joint_objective, 1L),
    objective_gain = tail(hist$joint_objective, 1L) - hist$joint_objective[1L],
    stringsAsFactors = FALSE
  )
}
