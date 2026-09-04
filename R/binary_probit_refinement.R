#!/usr/bin/env Rscript

# ============================================================================
# Self-contained joint refinement algorithm for the binary probit
# independent-mixture factor model
#
# This file contains the actual refinement code.  It assumes the pretraining
# file has been sourced, because the refinement uses the same univariate mixture
# fitting utilities:
#
#   source("R/binary_probit_pretraining.R")
#   source("R/binary_probit_refinement.R")
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
#        c. Refit the fixed-G marginal mixture distribution for each factor.
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
  source(file.path(script_dir, "binary_probit_pretraining.R"))
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
  # Accept either one global loading penalty or one penalty per factor column.
  # The simulation and IFEval scripts use this for the sparse itemwise probit
  # loading regressions in refinement.
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
    X = X,
    F_hat = F_hat,
    Lambda = Lambda,
    alpha = alpha,
    mixture_fits = mixture_fits,
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
  # Smooth negative log likelihood for one item.  beta contains
  # (alpha_j, lambda_j) when intercepts are estimated and only lambda_j
  # otherwise.  The L1 penalty is handled by proximal steps in
  # fit_probit_lasso_item(), not inside this smooth objective.
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
  # Pipeline role: update one item's alpha_j and loading row lambda_j given
  # current factor scores.  This is proximal-gradient probit regression:
  # gradient step on the binary likelihood, soft-threshold loadings, no
  # shrinkage on the intercept.
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

fit_probit_glmnet_lasso_item <- function(
    y,
    F_hat,
    alpha_init = NULL,
    beta_init = NULL,
    lambda_l1_penalty = 0,
    estimate_intercept = TRUE,
    maxit = 200L,
    tol = 1e-6,
    standardize = FALSE) {
  # Pipeline role: faster exact probit-lasso item update using glmnet >= 4,
  # which accepts GLM family objects such as binomial(link = "probit").  glmnet
  # minimizes average negative log likelihood, so the simulation's summed-loss
  # penalty lambda * |beta| is passed as lambda / n.
  H <- ncol(F_hat)
  lambda_l1_penalty <- coerce_lambda_l1_penalty(lambda_l1_penalty, H)

  fallback <- function() {
    fit_probit_lasso_item(
      y = y,
      F_hat = F_hat,
      alpha_init = alpha_init,
      beta_init = beta_init,
      lambda_l1_penalty = lambda_l1_penalty,
      estimate_intercept = estimate_intercept,
      maxit = maxit,
      tol = tol
    )
  }

  if (!requireNamespace("glmnet", quietly = TRUE)) return(fallback())

  penalty_mean <- mean(lambda_l1_penalty)
  if (!is.finite(penalty_mean) || penalty_mean <= 0) return(fallback())
  penalty_factor <- lambda_l1_penalty / penalty_mean
  lambda_glmnet <- penalty_mean / length(y)

  fit <- tryCatch(
    glmnet::glmnet(
      x = as.matrix(F_hat),
      y = as.numeric(y),
      family = binomial(link = "probit"),
      alpha = 1,
      lambda = lambda_glmnet,
      standardize = isTRUE(standardize),
      intercept = isTRUE(estimate_intercept),
      penalty.factor = penalty_factor,
      thresh = tol,
      maxit = max(100000L, as.integer(maxit))
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(fallback())

  coef_mat <- tryCatch(
    as.matrix(stats::coef(fit, s = lambda_glmnet)),
    error = function(e) NULL
  )
  if (is.null(coef_mat)) return(fallback())

  alpha <- 0
  if (isTRUE(estimate_intercept)) {
    if (!"(Intercept)" %in% rownames(coef_mat)) return(fallback())
    alpha <- unname(coef_mat["(Intercept)", 1L])
  }

  loading_rows <- setdiff(rownames(coef_mat), "(Intercept)")
  if (length(loading_rows) < H) return(fallback())
  loading <- as.numeric(coef_mat[loading_rows[seq_len(H)], 1L])

  if (!is.finite(alpha) || any(!is.finite(loading))) return(fallback())
  list(alpha = alpha, lambda = loading)
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
    lasso_backend = c("proximal", "glmnet"),
    glmnet_standardize = FALSE,
    parallel = FALSE,
    workers = NULL) {
  # Pipeline role: itemwise loading/intercept update in refinement.  Conditional
  # on F, every item is an independent probit regression, so this step can be
  # parallelized over item columns.
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  p <- ncol(X)
  H <- ncol(F_hat)
  lasso_backend <- match.arg(lasso_backend)
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
      if (identical(lasso_backend, "glmnet")) {
        return(fit_probit_glmnet_lasso_item(
          y = X[, j],
          F_hat = F_hat,
          alpha_init = alpha_j_init,
          beta_init = beta_init,
          lambda_l1_penalty = lambda_l1_penalty,
          estimate_intercept = estimate_intercept,
          maxit = lasso_maxit,
          tol = lasso_tol,
          standardize = glmnet_standardize
        ))
      }
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

moment_corrected_mixture_responsibility_list <- function(
    F_hat,
    factor_var_diag,
    mixture_fits,
    hard = FALSE,
    min_var = 1e-8) {
  # Variational/Laplace soft classification.  If f_ih | X is approximated by
  # N(m_ih, v_ih), classify using E_q[log N(f_ih; mu_hg, sigma_hg^2)] rather
  # than pretending the point estimate m_ih is observed.
  F_hat <- as.matrix(F_hat)
  factor_var_diag <- as.matrix(factor_var_diag)
  lapply(seq_len(ncol(F_hat)), function(h) {
    fit_h <- mixture_fits[[h]]
    var_h <- pmax(fit_h$var, min_var)
    log_resp <- vapply(seq_along(fit_h$pi), function(g) {
      log(pmax(fit_h$pi[g], 1e-12)) -
        0.5 * log(2 * pi * var_h[g]) -
        0.5 * ((F_hat[, h] - fit_h$mu[g])^2 + factor_var_diag[, h]) / var_h[g]
    }, numeric(nrow(F_hat)))
    log_den <- row_log_sum_exp(log_resp)
    resp_h <- exp(log_resp - log_den)
    if (isTRUE(hard)) {
      z <- max.col(resp_h, ties.method = "first")
      hard_resp <- matrix(0, nrow = nrow(resp_h), ncol = ncol(resp_h))
      hard_resp[cbind(seq_len(nrow(resp_h)), z)] <- 1
      resp_h <- hard_resp
    }
    resp_h
  })
}

factor_posterior_diag_variance <- function(
    X,
    F_hat,
    Lambda,
    alpha = NULL,
    mixture_fits,
    responsibility_list = NULL,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    mixture_prior_weight = 1,
    min_var = 1e-6,
    max_var = 1e6,
    parallel = FALSE,
    workers = NULL) {
  # Diagonal Laplace approximation for Var(f_i | X_i).  The binary probit
  # likelihood contributes item information Lambda' W_i Lambda.  The mixture
  # prior contributes the local expected Gaussian precision sum_g gamma_ihg /
  # sigma_hg^2.  We keep only the diagonal for speed and stability.
  factor_update <- match.arg(factor_update)
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  if (is.null(alpha)) alpha <- rep(0, ncol(X))
  alpha <- as.numeric(alpha)
  H <- ncol(F_hat)

  rows <- parallel_lapply(seq_len(nrow(X)), function(i) {
    eta <- as.numeric(alpha + Lambda %*% F_hat[i, ])
    phi <- dnorm(eta)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    y <- X[i, ]

    # Observed information for log Phi(eta) or log Phi(-eta).  Clamp to avoid
    # numerical trouble in nearly separated probit items.
    w1 <- phi / p1
    w0 <- phi / p0
    info_eta <- y * w1 * (eta + w1) + (1 - y) * w0 * (w0 - eta)
    info_eta[!is.finite(info_eta) | info_eta < 0] <- 0
    precision <- as.numeric(crossprod(Lambda^2, info_eta))

    for (h in seq_len(H)) {
      if (!is.null(responsibility_list)) {
        resp_h <- responsibility_list[[h]][i, ]
      } else {
        resp_h <- mixture_responsibilities(F_hat[i, h], mixture_fits[[h]])
        resp_h <- as.numeric(resp_h[1L, ])
      }
      resp_h <- resp_h / sum(resp_h)
      precision[h] <- precision[h] + mixture_prior_weight *
        sum(resp_h / pmax(mixture_fits[[h]]$var, min_var))
    }

    1 / pmin(pmax(precision, 1 / max_var), 1 / min_var)
  }, parallel = parallel, workers = workers)

  out <- do.call(rbind, rows)
  colnames(out) <- colnames(F_hat)
  out
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
  # Pipeline role: update one subject's H-dimensional factor score by MAP while
  # holding Lambda, alpha, and the mixture prior fixed.
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
  # Pipeline role: update all subject factor scores.  This is usually the
  # slowest refinement block in large n,p simulations because it solves n small
  # H-dimensional optimizations, each evaluating all p item probabilities.
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

update_gmm_1d_from_responsibilities <- function(
    x,
    resp,
    min_var = 1e-3,
    min_weight = 1e-4,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1) {
  # Weighted M-step for one marginal Gaussian mixture when component
  # responsibilities have already been estimated.  This shares the same MLE/MAP
  # formulas as fit_gmm_1d(), but skips the inner EM loop.
  mixture_update <- match.arg(mixture_update)
  x <- as.numeric(x)
  resp <- as.matrix(resp)
  n <- length(x)
  G <- ncol(resp)
  if (nrow(resp) != n || G < 1L) {
    stop("responsibility matrix must have one row per factor score.")
  }

  resp[!is.finite(resp) | resp < 0] <- 0
  row_total <- rowSums(resp)
  bad_rows <- !is.finite(row_total) | row_total <= 0
  if (any(bad_rows)) {
    resp[bad_rows, ] <- 1 / G
    row_total[bad_rows] <- 1
  }
  resp <- resp / row_total

  nk <- colSums(resp) + 1e-12
  xbar <- colSums(resp * x) / nk
  centered <- sweep(matrix(x, nrow = n, ncol = G), 2L, xbar, "-")
  ss <- colSums(resp * centered^2)

  if (mixture_update == "map") {
    alpha <- rep(weight_prior_alpha, G)
    if (length(weight_prior_alpha) == G) alpha <- weight_prior_alpha
    if (any(alpha < 1)) stop("weight_prior_alpha must be >= 1 for MAP weights.")
    pi_g <- pmax(nk + alpha - 1, min_weight)
    pi_g <- pi_g / sum(pi_g)

    kappa0 <- pmax(mu_prior_kappa, 0)
    kappa_n <- kappa0 + nk
    mu_g <- (kappa0 * mu_prior_mean + nk * xbar) / kappa_n
    shape_n <- var_prior_shape + nk / 2
    scale_n <- var_prior_scale + 0.5 * ss +
      (kappa0 * nk * (xbar - mu_prior_mean)^2) / (2 * kappa_n)
    var_g <- pmax(scale_n / (shape_n + 1), min_var)
  } else {
    pi_g <- pmax(nk / n, min_weight)
    pi_g <- pi_g / sum(pi_g)
    mu_g <- xbar
    var_g <- pmax(ss / nk, min_var)
  }

  fit <- list(
    pi = pi_g,
    mu = mu_g,
    var = var_g,
    G = G,
    converged = TRUE,
    mixture_update = mixture_update,
    fixed_responsibility_mstep = TRUE,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha
  )
  fit$loglik <- sum(log_dmix_1d(x, fit))

  ord <- order(fit$mu)
  fit$pi <- fit$pi[ord]
  fit$mu <- fit$mu[ord]
  fit$var <- fit$var[ord]
  fit
}

update_gmm_1d_from_posterior_moments <- function(
    m,
    v,
    resp,
    min_var = 1e-3,
    min_weight = 1e-4,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1) {
  # Moment-corrected M-step for q(f_ih | X_i) ~= N(m_ih, v_ih).  This is the
  # deterministic analogue of Gibbs using sampled latent factors: mixture
  # updates receive posterior first and second moments rather than point MAP
  # scores alone.
  mixture_update <- match.arg(mixture_update)
  m <- as.numeric(m)
  v <- pmax(as.numeric(v), 0)
  resp <- as.matrix(resp)
  n <- length(m)
  G <- ncol(resp)
  if (length(v) != n || nrow(resp) != n || G < 1L) {
    stop("posterior moments and responsibilities must have compatible sizes.")
  }

  resp[!is.finite(resp) | resp < 0] <- 0
  row_total <- rowSums(resp)
  bad_rows <- !is.finite(row_total) | row_total <= 0
  if (any(bad_rows)) {
    resp[bad_rows, ] <- 1 / G
    row_total[bad_rows] <- 1
  }
  resp <- resp / row_total

  nk <- colSums(resp) + 1e-12
  xbar <- colSums(resp * m) / nk
  centered <- sweep(matrix(m, nrow = n, ncol = G), 2L, xbar, "-")
  ss <- colSums(resp * (centered^2 + v))

  if (mixture_update == "map") {
    alpha <- rep(weight_prior_alpha, G)
    if (length(weight_prior_alpha) == G) alpha <- weight_prior_alpha
    if (any(alpha < 1)) stop("weight_prior_alpha must be >= 1 for MAP weights.")
    pi_g <- pmax(nk + alpha - 1, min_weight)
    pi_g <- pi_g / sum(pi_g)

    kappa0 <- pmax(mu_prior_kappa, 0)
    kappa_n <- kappa0 + nk
    mu_g <- (kappa0 * mu_prior_mean + nk * xbar) / kappa_n
    shape_n <- var_prior_shape + nk / 2
    scale_n <- var_prior_scale + 0.5 * ss +
      (kappa0 * nk * (xbar - mu_prior_mean)^2) / (2 * kappa_n)
    var_g <- pmax(scale_n / (shape_n + 1), min_var)
  } else {
    pi_g <- pmax(nk / n, min_weight)
    pi_g <- pi_g / sum(pi_g)
    mu_g <- xbar
    var_g <- pmax(ss / nk, min_var)
  }

  fit <- list(
    pi = pi_g,
    mu = mu_g,
    var = var_g,
    G = G,
    converged = TRUE,
    mixture_update = mixture_update,
    posterior_moment_mstep = TRUE,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha
  )
  fit$loglik <- sum(log_dmix_1d(m, fit))

  ord <- order(fit$mu)
  fit$pi <- fit$pi[ord]
  fit$mu <- fit$mu[ord]
  fit$var <- fit$var[ord]
  fit
}

update_mixture_fits_fixed_G <- function(
    F_hat,
    mixture_fits,
    responsibility_list = NULL,
    factor_var_diag = NULL,
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
  if (!is.null(factor_var_diag)) factor_var_diag <- as.matrix(factor_var_diag)
  parallel_lapply(seq_len(ncol(F_hat)), function(h) {
    if (!is.null(responsibility_list)) {
      if (!is.null(factor_var_diag)) {
        return(update_gmm_1d_from_posterior_moments(
          m = F_hat[, h],
          v = factor_var_diag[, h],
          resp = responsibility_list[[h]],
          min_var = min_var,
          mixture_update = mixture_update,
          mu_prior_mean = mu_prior_mean,
          mu_prior_kappa = mu_prior_kappa,
          var_prior_shape = var_prior_shape,
          var_prior_scale = var_prior_scale,
          weight_prior_alpha = weight_prior_alpha
        ))
      }
      # One weighted M-step using responsibilities fixed before the current
      # factor-score update.  This is a fast classification/soft-EM variant:
      # clusters are estimated first, then the mixture parameters are refreshed
      # from the updated factor scores without a free inner EM relabeling step.
      return(update_gmm_1d_from_responsibilities(
        F_hat[, h],
        resp = responsibility_list[[h]],
        min_var = min_var,
        mixture_update = mixture_update,
        mu_prior_mean = mu_prior_mean,
        mu_prior_kappa = mu_prior_kappa,
        var_prior_shape = var_prior_shape,
        var_prior_scale = var_prior_scale,
        weight_prior_alpha = weight_prior_alpha
      ))
    }

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
    responsibility_list = NULL,
    factor_var_diag = NULL,
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
  # Pipeline role: refresh the independent marginal mixture prior after factor
  # scores move.  Fixed component counts are preserved factor by factor.
  update_mixture_fits_fixed_G(
    F_hat = F_hat,
    mixture_fits = mixture_fits,
    responsibility_list = responsibility_list,
    factor_var_diag = factor_var_diag,
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
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    min_mixture_var = 1e-3,
    mixture_update = c("map", "mle"),
    mixture_refit = c("em", "fixed_responsibility_mstep", "posterior_moment_mstep"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    mixture_prior_weight = 1,
    estimate_intercept = TRUE,
    preestimate_loadings = TRUE,
    lambda_l1_penalty = 0,
    lasso_maxit = 200L,
    lasso_tol = 1e-6,
    lasso_backend = c("proximal", "glmnet"),
    glmnet_standardize = FALSE,
    normalize_factor_scale = TRUE,
    normalize_factor_location = TRUE,
    factor_scale_target = 1,
    factor_scale_method = c("sd", "rms"),
    objective_tolerance = 1e-5,
    min_refine_iter = 1L,
    stopping_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    enforce_monotone_refinement = TRUE,
    monotone_tolerance = 1e-8,
    return_best_refinement_iteration = FALSE,
    refinement_selection_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    require_mixture_convergence_for_stop = FALSE,
    store_refinement_trace = FALSE,
    parallel = FALSE,
    workers = NULL,
    verbose = TRUE,
    ...) {
  # Pipeline role: Stage 3 joint MAP refinement.  Starting from pretraining,
  # alternate factor-score MAP updates, itemwise probit loading/intercept
  # updates, factor scale/location normalization, and fixed-G mixture updates.
  factor_update <- match.arg(factor_update)
  mixture_update <- match.arg(mixture_update)
  mixture_refit <- match.arg(mixture_refit)
  factor_scale_method <- match.arg(factor_scale_method)
  stopping_objective <- match.arg(stopping_objective)
  refinement_selection_objective <- match.arg(refinement_selection_objective)
  lasso_backend <- match.arg(lasso_backend)
  X <- as.matrix(X)
  workers <- resolve_workers(workers)
  min_refine_iter <- as.integer(min_refine_iter)
  mixture_prior_weight <- as.numeric(mixture_prior_weight)
  estimate_intercept <- isTRUE(estimate_intercept)
  enforce_monotone_refinement <- isTRUE(enforce_monotone_refinement)
  monotone_tolerance <- as.numeric(monotone_tolerance)
  if (!is.finite(monotone_tolerance) || monotone_tolerance < 0) monotone_tolerance <- 0
  if (!is.finite(mixture_prior_weight) || mixture_prior_weight < 0) {
    stop("mixture_prior_weight must be a nonnegative finite number.")
  }

  # Initialize refinement from the output of the pretraining algorithm.
  F_hat <- pretrain_fit$F_hat
  Lambda <- pretrain_fit$Lambda_hat
  alpha <- if (estimate_intercept && !is.null(pretrain_fit$alpha_hat)) {
    as.numeric(pretrain_fit$alpha_hat)
  } else {
    rep(0, ncol(X))
  }
  mixture_fits <- pretrain_fit$mixture_fits

  lambda_l1_penalty_vec <- coerce_lambda_l1_penalty(lambda_l1_penalty, ncol(Lambda))
  lambda_penalty_history <- list()

  if (isTRUE(preestimate_loadings)) {
    # This is the requested ordering: finish pretraining, then estimate Lambda
    # given the final pretrained factors before joint refinement begins.
    if (verbose) message("Pre-estimating probit loadings before joint refinement.")
    loading_fit <- update_binary_probit_loadings_glm(
      X = X,
      F_hat = F_hat,
      alpha_init = alpha,
      Lambda_init = Lambda,
      lambda_l1_penalty = lambda_l1_penalty_vec,
      estimate_intercept = estimate_intercept,
      lasso_maxit = lasso_maxit,
      lasso_tol = lasso_tol,
      lasso_backend = lasso_backend,
      glmnet_standardize = glmnet_standardize,
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
    lambda_l1_penalty = lambda_l1_penalty_vec
  )
  lambda_penalty_history[[1L]] <- data.frame(
    iteration = 0L,
    factor = seq_len(ncol(Lambda)),
    column_l1 = colSums(abs(Lambda)),
    column_l2 = sqrt(colSums(Lambda^2)),
    l1_penalty = lambda_l1_penalty_vec
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
    max_factor_scale_before_normalize = max(apply(F_hat, 2L, sd)),
    max_factor_scale_after_normalize = max(apply(F_hat, 2L, sd)),
    G_hat = paste(vapply(mixture_fits, function(z) length(z$pi), integer(1)),
                  collapse = ","),
    factor_update_seconds = NA_real_,
    location_normalize_seconds = NA_real_,
    normalize_seconds = NA_real_,
    lambda_update_seconds = NA_real_,
    mixture_update_seconds = NA_real_,
    factor_posterior_var_median = NA_real_,
    all_mixtures_converged = all(vapply(mixture_fits, function(z) isTRUE(z$converged), logical(1))),
    iteration_rejected_by_monotone_guard = FALSE,
    objective_seconds = NA_real_,
    iteration_seconds = NA_real_
  )

  converged <- FALSE
  n_completed <- 0L
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
    previous_F_hat <- F_hat
    previous_alpha <- alpha
    previous_Lambda <- Lambda
    previous_mixture_fits <- mixture_fits
    fixed_mixture_responsibilities <- NULL
    if (mixture_refit == "fixed_responsibility_mstep") {
      fixed_mixture_responsibilities <- mixture_responsibility_list(
        F_hat,
        mixture_fits,
        hard = factor_update == "conditional_hard"
      )
    }

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
    loading_fit <- update_binary_probit_loadings_glm(
      X = X,
      F_hat = F_hat,
      alpha_init = alpha,
      Lambda_init = Lambda,
      lambda_l1_penalty = lambda_l1_penalty_vec,
      estimate_intercept = estimate_intercept,
      lasso_maxit = lasso_maxit,
      lasso_tol = lasso_tol,
      lasso_backend = lasso_backend,
      glmnet_standardize = glmnet_standardize,
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
    lambda_update_seconds <- as.numeric(difftime(Sys.time(), lambda_update_start, units = "secs"))

    factor_var_diag_for_mixture <- NULL
    mixture_refit_responsibilities <- fixed_mixture_responsibilities
    if (mixture_refit == "posterior_moment_mstep") {
      point_responsibilities <- mixture_responsibility_list(
        F_hat,
        mixture_fits,
        hard = factor_update == "conditional_hard"
      )
      factor_var_diag_for_mixture <- factor_posterior_diag_variance(
        X = X,
        F_hat = F_hat,
        Lambda = Lambda,
        alpha = alpha,
        mixture_fits = mixture_fits,
        responsibility_list = point_responsibilities,
        factor_update = factor_update,
        mixture_prior_weight = mixture_prior_weight,
        min_var = min_mixture_var,
        parallel = parallel,
        workers = workers
      )
      mixture_refit_responsibilities <- moment_corrected_mixture_responsibility_list(
        F_hat = F_hat,
        factor_var_diag = factor_var_diag_for_mixture,
        mixture_fits = mixture_fits,
        hard = FALSE,
        min_var = min_mixture_var
      )
    }

    # Step 3: update marginal mixture profiles, keeping the pretrained component
    # count fixed for each factor.
    mixture_update_start <- Sys.time()
    mixture_fits <- update_mixture_fits_refinement(
      F_hat = F_hat,
      mixture_fits = mixture_fits,
      responsibility_list = mixture_refit_responsibilities,
      factor_var_diag = factor_var_diag_for_mixture,
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
    mixture_update_seconds <- as.numeric(difftime(Sys.time(), mixture_update_start, units = "secs"))
    all_mixtures_converged <- all(vapply(mixture_fits, function(z) isTRUE(z$converged), logical(1)))

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
      lambda_l1_penalty = lambda_l1_penalty_vec
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

    iteration_rejected_by_monotone_guard <- FALSE
    if (enforce_monotone_refinement &&
        is.finite(objective_improvement) &&
        objective_improvement < -monotone_tolerance) {
      iteration_rejected_by_monotone_guard <- TRUE
      rejected_drop <- -objective_improvement
      F_hat <- previous_F_hat
      alpha <- previous_alpha
      Lambda <- previous_Lambda
      mixture_fits <- previous_mixture_fits
      all_mixtures_converged <- all(vapply(mixture_fits, function(z) isTRUE(z$converged), logical(1)))

      previous_row <- history[[iter]]
      current_binary_loglik <- previous_row$binary_loglik
      current_mixture_loglik <- previous_row$mixture_loglik
      current_mixture_parameter_logprior <- previous_row$mixture_parameter_logprior
      current_lambda_logprior <- previous_row$lambda_logprior
      current_objective <- previous_row$joint_objective
      current_posterior_objective <- previous_row$posterior_objective
      current_stopping_value <- previous_stopping_value
      objective_improvement <- 0
      relative_improvement <- 0

      if (verbose) {
        message(
          "Rejecting refinement iteration ", iter,
          ": ", stopping_objective, " decreased by ",
          signif(rejected_drop, 3),
          "."
        )
      }
    }
    lambda_penalty_history[[iter + 1L]] <- data.frame(
      iteration = iter,
      factor = seq_len(ncol(Lambda)),
      column_l1 = colSums(abs(Lambda)),
      column_l2 = sqrt(colSums(Lambda^2)),
      l1_penalty = lambda_l1_penalty_vec
    )

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
      max_factor_scale_before_normalize = max(scale_before_normalize),
      max_factor_scale_after_normalize = max(scale_after_normalize),
      G_hat = paste(vapply(mixture_fits, function(z) length(z$pi), integer(1)),
                    collapse = ","),
      factor_update_seconds = factor_update_seconds,
      location_normalize_seconds = location_normalize_seconds,
      normalize_seconds = normalize_seconds,
      lambda_update_seconds = lambda_update_seconds,
      mixture_update_seconds = mixture_update_seconds,
      factor_posterior_var_median = if (!is.null(factor_var_diag_for_mixture)) {
        median(factor_var_diag_for_mixture, na.rm = TRUE)
      } else {
        NA_real_
      },
      all_mixtures_converged = all_mixtures_converged,
      iteration_rejected_by_monotone_guard = iteration_rejected_by_monotone_guard,
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
    if (isTRUE(iteration_rejected_by_monotone_guard)) {
      if (verbose) message("Stopping refinement after monotone-guard rejection.")
      break
    }
    if (!is.null(objective_tolerance) &&
        is.finite(objective_tolerance) &&
        iter >= min_refine_iter &&
        is.finite(relative_improvement) &&
        relative_improvement <= objective_tolerance &&
        (!isTRUE(require_mixture_convergence_for_stop) || all_mixtures_converged)) {
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
    enforce_monotone_refinement = enforce_monotone_refinement,
    monotone_tolerance = monotone_tolerance,
    return_best_refinement_iteration = isTRUE(return_best_refinement_iteration),
    refinement_selection_objective = refinement_selection_objective,
    selected_refinement_iteration = selected_refinement_iteration,
    maxit_per_subject = maxit_per_subject,
    mixture_max_iter = mixture_max_iter,
    G_selection = "fixed",
    factor_update = factor_update,
    min_mixture_var = min_mixture_var,
    mixture_update = mixture_update,
    mixture_refit = mixture_refit,
    mixture_prior = list(
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha
    ),
    mixture_prior_weight = mixture_prior_weight,
    estimate_intercept = estimate_intercept,
    lambda_l1_penalty = lambda_l1_penalty,
    lambda_penalty_history = if (length(lambda_penalty_history)) do.call(rbind, lambda_penalty_history) else NULL,
    lasso_maxit = lasso_maxit,
    lasso_tol = lasso_tol,
    lasso_backend = lasso_backend,
    glmnet_standardize = isTRUE(glmnet_standardize),
    normalize_factor_scale = normalize_factor_scale,
    normalize_factor_location = isTRUE(normalize_factor_location),
    factor_scale_target = factor_scale_target,
    factor_scale_method = factor_scale_method,
    preestimate_loadings = isTRUE(preestimate_loadings),
    store_refinement_trace = isTRUE(store_refinement_trace),
    require_mixture_convergence_for_stop = isTRUE(require_mixture_convergence_for_stop),
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
    G_fixed = NULL,
    n_aug_iter = 4L,
    z_update = c("sample", "expectation"),
    n_refine_iter = 5L,
    n_mix_starts = 3L,
    mixture_max_iter = 20L,
    mixture_update = c("map", "mle"),
    mixture_refit = c("em", "fixed_responsibility_mstep", "posterior_moment_mstep"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    refine_mu_prior_mean = NULL,
    refine_mu_prior_kappa = NULL,
    refine_var_prior_shape = NULL,
    refine_var_prior_scale = NULL,
    refine_weight_prior_alpha = NULL,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    min_mixture_var = 1e-3,
    mixture_prior_weight = 1,
    lambda_l1_penalty = 0,
    lasso_backend = c("proximal", "glmnet"),
    glmnet_standardize = FALSE,
    normalize_factor_scale = TRUE,
    normalize_factor_location = TRUE,
    factor_scale_target = 1,
    factor_scale_method = c("sd", "rms"),
    objective_tolerance = 1e-5,
    min_refine_iter = 1L,
    stopping_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    enforce_monotone_refinement = TRUE,
    monotone_tolerance = 1e-8,
    return_best_refinement_iteration = FALSE,
    refinement_selection_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    require_mixture_convergence_for_stop = FALSE,
    parallel = FALSE,
    workers = NULL,
    seed = 20260715L,
    verbose = TRUE,
    ...) {
  z_update <- match.arg(z_update)
  factor_update <- match.arg(factor_update)
  mixture_update <- match.arg(mixture_update)
  mixture_refit <- match.arg(mixture_refit)
  factor_scale_method <- match.arg(factor_scale_method)
  stopping_objective <- match.arg(stopping_objective)
  refinement_selection_objective <- match.arg(refinement_selection_objective)
  lasso_backend <- match.arg(lasso_backend)
  if (is.null(refine_mu_prior_mean)) refine_mu_prior_mean <- mu_prior_mean
  if (is.null(refine_mu_prior_kappa)) refine_mu_prior_kappa <- mu_prior_kappa
  if (is.null(refine_var_prior_shape)) refine_var_prior_shape <- var_prior_shape
  if (is.null(refine_var_prior_scale)) refine_var_prior_scale <- var_prior_scale
  if (is.null(refine_weight_prior_alpha)) refine_weight_prior_alpha <- weight_prior_alpha

  pretrain_fit <- fit_binary_probit_pretraining(
    X = X,
    H = H,
    H_max = H_max,
    G_fixed = G_fixed,
    n_aug_iter = n_aug_iter,
    z_update = z_update,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha,
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
    mixture_prior_weight = mixture_prior_weight,
    lambda_l1_penalty = lambda_l1_penalty,
    lasso_backend = lasso_backend,
    glmnet_standardize = glmnet_standardize,
    normalize_factor_scale = normalize_factor_scale,
    normalize_factor_location = normalize_factor_location,
    factor_scale_target = factor_scale_target,
    factor_scale_method = factor_scale_method,
    objective_tolerance = objective_tolerance,
    min_refine_iter = min_refine_iter,
    stopping_objective = stopping_objective,
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
    normalize_factor_scale = fit$joint_refinement$normalize_factor_scale,
    factor_scale_target = fit$joint_refinement$factor_scale_target,
    factor_scale_method = fit$joint_refinement$factor_scale_method,
    initial_joint_objective = hist$joint_objective[1L],
    final_joint_objective = tail(hist$joint_objective, 1L),
    objective_gain = tail(hist$joint_objective, 1L) - hist$joint_objective[1L],
    stringsAsFactors = FALSE
  )
}
