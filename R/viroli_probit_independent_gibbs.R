#!/usr/bin/env Rscript

# Probit-augmented Viroli-style independent factor analysis Gibbs sampler.
#
# Model:
#   X_ij = 1{Z_ij > 0}
#   Z_i | f_i ~ N(alpha + Lambda f_i, I_p)
#   f_ih | c_ih = g ~ N(mu_hg, sig2_hg), independently over h.
#
# This differs from the joint-mixture Gibbs comparator because it samples one
# mixture allocation per factor coordinate.  The allocation cost is O(n sum_h
# G_h), not O(n prod_h G_h).  Each Gibbs draw is normalized so that the current
# marginal mixture for every factor has mean zero and variance one; alpha and
# Lambda are transformed simultaneously, so the probit linear predictor is
# unchanged.

viroli_sample_dirichlet <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

viroli_sample_inv_gamma <- function(shape, rate) {
  1 / rgamma(1L, shape = shape, rate = rate)
}

viroli_rmvnorm_chol <- function(mean, cov) {
  mean + as.numeric(t(chol(cov)) %*% rnorm(length(mean)))
}

viroli_rinvgauss <- function(mu, lambda) {
  # Michael-Schucany-Haas sampler for IG(mu, lambda).
  mu <- as.numeric(mu)
  lambda <- as.numeric(lambda)
  if (!is.finite(mu) || !is.finite(lambda) || mu <= 0 || lambda <= 0) return(NA_real_)
  y <- rnorm(1L)^2
  x <- mu + (mu^2 * y) / (2 * lambda) -
    (mu / (2 * lambda)) * sqrt(4 * mu * lambda * y + mu^2 * y^2)
  if (runif(1L) <= mu / (mu + x)) x else mu^2 / x
}

viroli_sample_binary_Z <- function(X, F, Lambda, alpha) {
  mean_mat <- sweep(F %*% t(Lambda), 2L, alpha, "+")
  Z <- matrix(NA_real_, nrow(X), ncol(X))
  one <- X == 1
  Z[one] <- rtruncnorm_binary_vec(mean_mat[one], 1, 0, Inf)
  Z[!one] <- rtruncnorm_binary_vec(mean_mat[!one], 1, -Inf, 0)
  Z
}

viroli_sample_alpha_lambda <- function(
    Z,
    F,
    tau_intercept = 5,
    tau_lambda = 1.5,
    parallel = FALSE,
    workers = NULL) {
  n <- nrow(Z)
  p <- ncol(Z)
  H <- ncol(F)
  W <- cbind(1, F)
  prior_prec <- diag(c(1 / tau_intercept^2, rep(1 / tau_lambda^2, H)), H + 1L)
  V <- solve(crossprod(W) + prior_prec)
  chol_V <- chol(V)
  WtZ <- crossprod(W, Z)

  rows <- parallel_lapply(seq_len(p), function(j) {
    m <- V %*% WtZ[, j]
    as.numeric(m + t(chol_V) %*% rnorm(H + 1L))
  }, parallel = parallel, workers = workers)
  draw_mat <- do.call(rbind, rows)

  list(alpha = draw_mat[, 1L], Lambda = draw_mat[, -1L, drop = FALSE])
}

viroli_coerce_loading_penalty <- function(lambda_l1_penalty, H) {
  penalty <- as.numeric(lambda_l1_penalty)
  if (!length(penalty) || any(!is.finite(penalty)) || any(penalty < 0)) {
    stop("lambda_l1_penalty must contain nonnegative finite values.")
  }
  if (length(penalty) == 1L) penalty <- rep(penalty, H)
  if (length(penalty) != H) {
    stop("lambda_l1_penalty must be a scalar or have length H.")
  }
  penalty
}

viroli_sample_laplace_local_variances <- function(
    Lambda,
    lambda_l1_penalty,
    min_local_var = 1e-8,
    max_local_var = 1e8,
    beta_floor = 1e-8) {
  # Bayesian lasso hierarchy:
  #   lambda_jh | tau_jh ~ N(0, tau_jh),
  #   tau_jh ~ Exp(lambda_l1_h^2 / 2).
  # Therefore 1/tau_jh | lambda_jh is inverse-Gaussian with
  # mean lambda_l1_h / |lambda_jh| and shape lambda_l1_h^2.
  Lambda <- as.matrix(Lambda)
  H <- ncol(Lambda)
  penalty <- viroli_coerce_loading_penalty(lambda_l1_penalty, H)
  local_var <- matrix(Inf, nrow(Lambda), H)

  for (h in seq_len(H)) {
    lambda_h <- penalty[h]
    if (lambda_h <= 0) next
    shape_h <- lambda_h^2
    for (j in seq_len(nrow(Lambda))) {
      mu_inv <- lambda_h / max(abs(Lambda[j, h]), beta_floor)
      inv_tau <- viroli_rinvgauss(mu_inv, shape_h)
      if (!is.finite(inv_tau) || inv_tau <= 0) inv_tau <- 1 / min_local_var
      local_var[j, h] <- pmin(pmax(1 / inv_tau, min_local_var), max_local_var)
    }
  }

  local_var
}

viroli_sample_alpha_lambda_laplace <- function(
    Z,
    F,
    local_var,
    tau_intercept = 5,
    parallel = FALSE,
    workers = NULL) {
  Z <- as.matrix(Z)
  F <- as.matrix(F)
  p <- ncol(Z)
  H <- ncol(F)
  W <- cbind(1, F)
  XtX <- crossprod(W)
  WtZ <- crossprod(W, Z)

  rows <- parallel_lapply(seq_len(p), function(j) {
    prior_prec <- diag(c(1 / tau_intercept^2, 1 / pmax(local_var[j, ], 1e-12)), H + 1L)
    V <- solve(XtX + prior_prec)
    m <- V %*% WtZ[, j]
    as.numeric(m + t(chol(V)) %*% rnorm(H + 1L))
  }, parallel = parallel, workers = workers)
  draw_mat <- do.call(rbind, rows)

  list(alpha = draw_mat[, 1L], Lambda = draw_mat[, -1L, drop = FALSE])
}

viroli_soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

viroli_fit_gaussian_lasso_item <- function(
    z,
    F,
    alpha_init = NULL,
    beta_init = NULL,
    lambda_l1_penalty = 0,
    maxit = 200L,
    tol = 1e-6) {
  # Sparse MAP version of the Viroli loading update for one item:
  #   min_{alpha_j, lambda_j} 0.5 ||z_j - alpha_j - F lambda_j||^2
  #                            + sum_h lambda_l1_h |lambda_jh|.
  # The intercept is never penalized.
  F <- as.matrix(F)
  H <- ncol(F)
  penalty <- viroli_coerce_loading_penalty(lambda_l1_penalty, H)
  beta <- if (is.null(beta_init)) rep(0, H) else as.numeric(beta_init)
  if (length(beta) != H || any(!is.finite(beta))) beta <- rep(0, H)
  alpha <- if (!is.null(alpha_init) && is.finite(alpha_init)) as.numeric(alpha_init) else mean(z)
  x2 <- pmax(colSums(F^2), 1e-12)

  for (iter in seq_len(maxit)) {
    old <- c(alpha, beta)

    # Given beta, the unpenalized intercept has a closed-form update.
    alpha <- mean(z - as.numeric(F %*% beta))
    residual <- z - alpha - as.numeric(F %*% beta)

    for (h in seq_len(H)) {
      partial <- residual + F[, h] * beta[h]
      raw <- sum(F[, h] * partial) / x2[h]
      beta_new <- viroli_soft_threshold(raw, penalty[h] / x2[h])
      residual <- partial - F[, h] * beta_new
      beta[h] <- beta_new
    }

    if (max(abs(c(alpha, beta) - old)) <= tol * (1 + max(abs(old)))) break
  }

  list(alpha = alpha, lambda = beta)
}

viroli_map_alpha_lambda_laplace <- function(
    Z,
    F,
    alpha_init = NULL,
    Lambda_init = NULL,
    lambda_l1_penalty = 0,
    lasso_maxit = 200L,
    lasso_tol = 1e-6,
    parallel = FALSE,
    workers = NULL) {
  Z <- as.matrix(Z)
  F <- as.matrix(F)
  p <- ncol(Z)
  H <- ncol(F)
  penalty <- viroli_coerce_loading_penalty(lambda_l1_penalty, H)

  rows <- parallel_lapply(seq_len(p), function(j) {
    viroli_fit_gaussian_lasso_item(
      z = Z[, j],
      F = F,
      alpha_init = if (!is.null(alpha_init)) alpha_init[j] else NULL,
      beta_init = if (!is.null(Lambda_init)) Lambda_init[j, ] else NULL,
      lambda_l1_penalty = penalty,
      maxit = lasso_maxit,
      tol = lasso_tol
    )
  }, parallel = parallel, workers = workers)

  Lambda <- do.call(rbind, lapply(rows, `[[`, "lambda"))
  alpha <- vapply(rows, `[[`, numeric(1L), "alpha")
  list(alpha = alpha, Lambda = Lambda)
}

sample_viroli_factor_allocations <- function(F, pi_mat, mu_mat, sig2_mat, G) {
  F <- as.matrix(F)
  n <- nrow(F)
  H <- ncol(F)
  C <- matrix(NA_integer_, n, H)

  for (h in seq_len(H)) {
    Gh <- G[h]
    log_prob <- matrix(NA_real_, n, Gh)
    for (g in seq_len(Gh)) {
      log_prob[, g] <- log(pi_mat[h, g] + 1e-300) +
        dnorm(F[, h], mean = mu_mat[h, g], sd = sqrt(pmax(sig2_mat[h, g], 1e-8)), log = TRUE)
    }
    log_prob <- log_prob - apply(log_prob, 1L, max)
    prob <- exp(log_prob)
    prob <- prob / rowSums(prob)
    cumprob <- t(apply(prob, 1L, cumsum))
    C[, h] <- 1L + rowSums(runif(n) > cumprob)
    C[, h] <- pmin(pmax(C[, h], 1L), Gh)
  }

  C
}

sample_viroli_factor_scores <- function(
    Z,
    alpha,
    Lambda,
    C,
    mu_mat,
    sig2_mat,
    min_var = 0.05,
    parallel = FALSE,
    workers = NULL) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  H <- ncol(Lambda)
  LtL <- crossprod(Lambda)
  Ltz <- sweep(Z, 2L, alpha, "-") %*% Lambda
  rows <- parallel_lapply(seq_len(n), function(i) {
    prior_mean <- mu_mat[cbind(seq_len(H), C[i, ])]
    prior_var <- pmax(sig2_mat[cbind(seq_len(H), C[i, ])], min_var)
    prior_prec <- diag(1 / prior_var, H)
    V <- solve(LtL + prior_prec)
    m <- V %*% (Ltz[i, ] + prior_prec %*% prior_mean)
    viroli_rmvnorm_chol(as.numeric(m), V)
  }, parallel = parallel, workers = workers)
  F <- do.call(rbind, rows)

  F
}

viroli_scalar_effective_sample_size <- function(x, max_lag = NULL) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L || !is.finite(var(x)) || var(x) <= 0) return(NA_real_)
  if (is.null(max_lag)) max_lag <- min(n - 1L, 250L)
  ac <- as.numeric(acf(x, lag.max = max_lag, plot = FALSE, demean = TRUE)$acf[-1L])
  ac <- ac[is.finite(ac)]
  if (!length(ac)) return(n)
  first_nonpositive <- which(ac <= 0)[1L]
  use_lags <- if (is.na(first_nonpositive)) seq_along(ac) else seq_len(first_nonpositive - 1L)
  if (!length(use_lags)) return(n)
  tau <- 1 + 2 * sum(ac[use_lags])
  if (!is.finite(tau) || tau <= 0) return(NA_real_)
  min(n, n / tau)
}

viroli_compute_parameter_ess_table <- function(trace_list) {
  make_table <- function(group, mat) {
    if (is.null(mat) || !length(mat)) return(NULL)
    mat <- as.matrix(mat)
    data.frame(
      parameter_group = group,
      parameter = colnames(mat),
      ess = apply(mat, 2L, viroli_scalar_effective_sample_size),
      n_draws = nrow(mat),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, Filter(Negate(is.null), list(
    make_table("alpha", trace_list$alpha),
    make_table("lambda", trace_list$lambda),
    make_table("pi", trace_list$pi),
    make_table("mu", trace_list$mu),
    make_table("sig2", trace_list$sig2)
  )))
  rownames(out) <- NULL
  out
}

viroli_summarize_parameter_ess <- function(ess_table) {
  if (is.null(ess_table) || !nrow(ess_table)) {
    return(data.frame(
      gibbs_min_parameter_ess = NA_real_,
      gibbs_median_parameter_ess = NA_real_,
      gibbs_mean_parameter_ess = NA_real_,
      gibbs_min_alpha_ess = NA_real_,
      gibbs_median_alpha_ess = NA_real_,
      gibbs_min_lambda_ess = NA_real_,
      gibbs_median_lambda_ess = NA_real_,
      gibbs_min_pi_ess = NA_real_,
      gibbs_median_pi_ess = NA_real_,
      gibbs_min_mu_ess = NA_real_,
      gibbs_median_mu_ess = NA_real_,
      gibbs_min_sig2_ess = NA_real_,
      gibbs_median_sig2_ess = NA_real_
    ))
  }
  group_value <- function(group, FUN) {
    z <- ess_table$ess[ess_table$parameter_group == group]
    if (!length(z) || all(!is.finite(z))) return(NA_real_)
    FUN(z, na.rm = TRUE)
  }
  data.frame(
    gibbs_min_parameter_ess = if (all(!is.finite(ess_table$ess))) NA_real_ else min(ess_table$ess, na.rm = TRUE),
    gibbs_median_parameter_ess = median(ess_table$ess, na.rm = TRUE),
    gibbs_mean_parameter_ess = mean(ess_table$ess, na.rm = TRUE),
    gibbs_min_alpha_ess = group_value("alpha", min),
    gibbs_median_alpha_ess = group_value("alpha", median),
    gibbs_min_lambda_ess = group_value("lambda", min),
    gibbs_median_lambda_ess = group_value("lambda", median),
    gibbs_min_pi_ess = group_value("pi", min),
    gibbs_median_pi_ess = group_value("pi", median),
    gibbs_min_mu_ess = group_value("mu", min),
    gibbs_median_mu_ess = group_value("mu", median),
    gibbs_min_sig2_ess = group_value("sig2", min),
    gibbs_median_sig2_ess = group_value("sig2", median)
  )
}

sample_viroli_mixture_parameters <- function(
    F,
    C,
    G,
    alpha_dirichlet = 1,
    mu0 = 0,
    kappa0 = 0.05,
    a0 = 3,
    b0 = 2,
    min_var = 0.05) {
  F <- as.matrix(F)
  n <- nrow(F)
  H <- ncol(F)
  G_max <- max(G)
  pi_mat <- matrix(NA_real_, H, G_max)
  mu_mat <- matrix(NA_real_, H, G_max)
  sig2_mat <- matrix(NA_real_, H, G_max)

  for (h in seq_len(H)) {
    Gh <- G[h]
    counts <- tabulate(C[, h], nbins = Gh)
    pi_h <- viroli_sample_dirichlet(counts + alpha_dirichlet)
    pi_mat[h, seq_len(Gh)] <- pi_h

    for (g in seq_len(Gh)) {
      idx <- which(C[, h] == g)
      n_g <- length(idx)
      if (n_g > 0L) {
        x <- F[idx, h]
        xbar <- mean(x)
        ss <- sum((x - xbar)^2)
      } else {
        xbar <- 0
        ss <- 0
      }

      kappa_n <- kappa0 + n_g
      a_n <- a0 + n_g / 2
      b_n <- b0 + 0.5 * ss + (kappa0 * n_g * (xbar - mu0)^2) / (2 * kappa_n)
      sig2_g <- max(viroli_sample_inv_gamma(a_n, b_n), min_var)
      mu_n <- (kappa0 * mu0 + n_g * xbar) / kappa_n
      mu_g <- rnorm(1L, mu_n, sqrt(sig2_g / kappa_n))

      mu_mat[h, g] <- mu_g
      sig2_mat[h, g] <- sig2_g
    }
  }

  list(pi = pi_mat, mu = mu_mat, sig2 = sig2_mat)
}

sort_viroli_components <- function(C, pi_mat, mu_mat, sig2_mat, G) {
  H <- length(G)
  for (h in seq_len(H)) {
    Gh <- G[h]
    ord <- order(mu_mat[h, seq_len(Gh)])
    inv_ord <- integer(Gh)
    inv_ord[ord] <- seq_len(Gh)
    C[, h] <- inv_ord[C[, h]]
    pi_mat[h, seq_len(Gh)] <- pi_mat[h, ord]
    mu_mat[h, seq_len(Gh)] <- mu_mat[h, ord]
    sig2_mat[h, seq_len(Gh)] <- sig2_mat[h, ord]
  }
  list(C = C, pi = pi_mat, mu = mu_mat, sig2 = sig2_mat)
}

normalize_viroli_draw <- function(F, alpha, Lambda, pi_mat, mu_mat, sig2_mat, G, min_scale = 1e-4) {
  H <- ncol(F)
  marginal_mean <- numeric(H)
  marginal_scale <- numeric(H)

  for (h in seq_len(H)) {
    Gh <- G[h]
    pi_h <- pi_mat[h, seq_len(Gh)]
    pi_h <- pi_h / sum(pi_h)
    mu_h <- mu_mat[h, seq_len(Gh)]
    sig2_h <- sig2_mat[h, seq_len(Gh)]
    m_h <- sum(pi_h * mu_h)
    second_h <- sum(pi_h * (sig2_h + mu_h^2))
    v_h <- max(second_h - m_h^2, min_scale^2)
    s_h <- sqrt(v_h)
    marginal_mean[h] <- m_h
    marginal_scale[h] <- s_h

    F[, h] <- (F[, h] - m_h) / s_h
    mu_mat[h, seq_len(Gh)] <- (mu_h - m_h) / s_h
    sig2_mat[h, seq_len(Gh)] <- sig2_h / s_h^2
    Lambda[, h] <- Lambda[, h] * s_h
  }

  alpha <- alpha + as.numeric(Lambda %*% (marginal_mean / marginal_scale))

  list(
    F = F,
    alpha = alpha,
    Lambda = Lambda,
    pi = pi_mat,
    mu = mu_mat,
    sig2 = sig2_mat,
    marginal_mean_before = marginal_mean,
    marginal_scale_before = marginal_scale
  )
}

initialize_viroli_probit_state <- function(X, H, G, seed = 1L) {
  set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X)
  G <- normalize_G_fixed(G, H)

  alpha0 <- qnorm((colSums(X) + 0.5) / (n + 1))
  Z <- initialize_binary_Z(X, seed = seed, alpha = alpha0)
  svd_out <- svd_scores_from_Z(Z, H = H, center_Z = TRUE)
  F <- scale(svd_out$S)
  if (any(!is.finite(F))) F <- matrix(rnorm(n * H), n, H)

  loading_fit <- update_working_loadings_with_intercept(Z, F, loading_penalty = 0)
  alpha <- loading_fit$alpha
  Lambda <- loading_fit$Lambda

  G_max <- max(G)
  C <- matrix(NA_integer_, n, H)
  pi_mat <- matrix(NA_real_, H, G_max)
  mu_mat <- matrix(NA_real_, H, G_max)
  sig2_mat <- matrix(NA_real_, H, G_max)

  for (h in seq_len(H)) {
    Gh <- G[h]
    if (n >= Gh) {
      km <- try(kmeans(F[, h], centers = Gh, nstart = 5, iter.max = 30), silent = TRUE)
      C[, h] <- if (inherits(km, "try-error")) sample.int(Gh, n, replace = TRUE) else km$cluster
    } else {
      C[, h] <- sample.int(Gh, n, replace = TRUE)
    }
    for (g in seq_len(Gh)) {
      idx <- which(C[, h] == g)
      pi_mat[h, g] <- max(length(idx) / n, 1e-4)
      mu_mat[h, g] <- if (length(idx)) mean(F[idx, h]) else 0
      sig2_mat[h, g] <- if (length(idx) > 1L) max(var(F[idx, h]), 0.05) else 1
    }
    pi_mat[h, seq_len(Gh)] <- pi_mat[h, seq_len(Gh)] / sum(pi_mat[h, seq_len(Gh)])
  }

  sorted <- sort_viroli_components(C, pi_mat, mu_mat, sig2_mat, G)
  list(Z = Z, F = as.matrix(F), C = sorted$C, alpha = alpha, Lambda = Lambda,
       pi = sorted$pi, mu = sorted$mu, sig2 = sorted$sig2)
}

fit_viroli_probit_independent_gibbs <- function(
    X,
    H,
    G,
    n_iter = 2000L,
    burn = 1000L,
    thin = 1L,
    tau_lambda = 1.5,
    tau_intercept = 5,
    alpha_dirichlet = 1,
    mu0 = 0,
    kappa0 = 0.05,
    a0 = 3,
    b0 = 2,
    min_var = 0.05,
    normalize_each_draw = TRUE,
    min_scale = 1e-4,
    lambda_l1_penalty = 0,
    laplace_min_local_var = 1e-8,
    laplace_max_local_var = 1e8,
    parallel = FALSE,
    workers = NULL,
    compute_parameter_ess = TRUE,
    seed = 1L,
    verbose = TRUE) {
  set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  G <- normalize_G_fixed(G, H)
  G_max <- max(G)

  state <- initialize_viroli_probit_state(X, H = H, G = G, seed = seed)
  Z <- state$Z
  F <- state$F
  C <- state$C
  alpha <- state$alpha
  Lambda <- state$Lambda
  pi_mat <- state$pi
  mu_mat <- state$mu
  sig2_mat <- state$sig2
  lambda_l1_penalty <- viroli_coerce_loading_penalty(lambda_l1_penalty, H)
  use_laplace_loading_prior <- any(lambda_l1_penalty > 0)
  local_loading_var <- matrix(tau_lambda^2, p, H)

  keep_target <- max(0L, floor((n_iter - burn) / max(thin, 1L)))
  keep_F <- matrix(0, n, H)
  keep_alpha <- numeric(p)
  keep_Lambda <- matrix(0, p, H)
  keep_pi <- matrix(0, H, G_max)
  keep_mu <- matrix(0, H, G_max)
  keep_sig2 <- matrix(0, H, G_max)
  n_keep <- 0L
  if (isTRUE(compute_parameter_ess) && keep_target > 0L) {
    alpha_trace <- matrix(NA_real_, keep_target, p)
    colnames(alpha_trace) <- paste0("alpha_", seq_len(p))
    lambda_trace <- matrix(NA_real_, keep_target, p * H)
    colnames(lambda_trace) <- as.vector(outer(paste0("item", seq_len(p)), paste0("factor", seq_len(H)), paste, sep = "_"))
    pi_trace <- matrix(NA_real_, keep_target, H * G_max)
    colnames(pi_trace) <- as.vector(outer(paste0("factor", seq_len(H)), paste0("component", seq_len(G_max)), paste, sep = "_"))
    mu_trace <- matrix(NA_real_, keep_target, H * G_max)
    colnames(mu_trace) <- colnames(pi_trace)
    sig2_trace <- matrix(NA_real_, keep_target, H * G_max)
    colnames(sig2_trace) <- colnames(pi_trace)
  } else {
    alpha_trace <- lambda_trace <- pi_trace <- mu_trace <- sig2_trace <- NULL
  }

  history <- data.frame(
    iteration = seq_len(n_iter),
    occupied_factor_classes = NA_integer_,
    max_abs_mixture_mean_before_normalize = NA_real_,
    max_abs_log_mixture_scale_before_normalize = NA_real_,
    z_sample_seconds = NA_real_,
    factor_sample_seconds = NA_real_,
    class_sample_seconds = NA_real_,
    mixture_weight_seconds = NA_real_,
    mixture_parameter_seconds = NA_real_,
    regression_seconds = NA_real_,
    normalization_seconds = NA_real_,
    keep_draw_seconds = NA_real_,
    iteration_seconds = NA_real_
  )

  t0 <- proc.time()[["elapsed"]]
  for (iter in seq_len(n_iter)) {
    iter_start <- Sys.time()

    z_start <- Sys.time()
    Z <- viroli_sample_binary_Z(X, F, Lambda, alpha)
    history$z_sample_seconds[iter] <- as.numeric(difftime(Sys.time(), z_start, units = "secs"))

    factor_start <- Sys.time()
    F <- sample_viroli_factor_scores(
      Z,
      alpha,
      Lambda,
      C,
      mu_mat,
      sig2_mat,
      min_var = min_var,
      parallel = parallel,
      workers = workers
    )
    history$factor_sample_seconds[iter] <- as.numeric(difftime(Sys.time(), factor_start, units = "secs"))

    class_start <- Sys.time()
    C <- sample_viroli_factor_allocations(F, pi_mat, mu_mat, sig2_mat, G)
    history$class_sample_seconds[iter] <- as.numeric(difftime(Sys.time(), class_start, units = "secs"))

    mixture_start <- Sys.time()
    mix <- sample_viroli_mixture_parameters(
      F = F,
      C = C,
      G = G,
      alpha_dirichlet = alpha_dirichlet,
      mu0 = mu0,
      kappa0 = kappa0,
      a0 = a0,
      b0 = b0,
      min_var = min_var
    )
    pi_mat <- mix$pi
    mu_mat <- mix$mu
    sig2_mat <- mix$sig2
    sorted <- sort_viroli_components(C, pi_mat, mu_mat, sig2_mat, G)
    C <- sorted$C
    pi_mat <- sorted$pi
    mu_mat <- sorted$mu
    sig2_mat <- sorted$sig2
    history$mixture_parameter_seconds[iter] <- as.numeric(difftime(Sys.time(), mixture_start, units = "secs"))

    regression_start <- Sys.time()
    if (isTRUE(use_laplace_loading_prior)) {
      local_loading_var <- viroli_sample_laplace_local_variances(
        Lambda = Lambda,
        lambda_l1_penalty = lambda_l1_penalty,
        min_local_var = laplace_min_local_var,
        max_local_var = laplace_max_local_var
      )
      regression_draw <- viroli_sample_alpha_lambda_laplace(
        Z = Z,
        F = F,
        local_var = local_loading_var,
        tau_intercept = tau_intercept,
        parallel = parallel,
        workers = workers
      )
    } else {
      regression_draw <- viroli_sample_alpha_lambda(
        Z = Z,
        F = F,
        tau_intercept = tau_intercept,
        tau_lambda = tau_lambda,
        parallel = parallel,
        workers = workers
      )
    }
    alpha <- regression_draw$alpha
    Lambda <- regression_draw$Lambda
    history$regression_seconds[iter] <- as.numeric(difftime(Sys.time(), regression_start, units = "secs"))

    normalization_start <- Sys.time()
    if (isTRUE(normalize_each_draw)) {
      normalized <- normalize_viroli_draw(
        F = F,
        alpha = alpha,
        Lambda = Lambda,
        pi_mat = pi_mat,
        mu_mat = mu_mat,
        sig2_mat = sig2_mat,
        G = G,
        min_scale = min_scale
      )
      F <- normalized$F
      alpha <- normalized$alpha
      Lambda <- normalized$Lambda
      pi_mat <- normalized$pi
      mu_mat <- normalized$mu
      sig2_mat <- normalized$sig2
      history$max_abs_mixture_mean_before_normalize[iter] <-
        max(abs(normalized$marginal_mean_before), na.rm = TRUE)
      history$max_abs_log_mixture_scale_before_normalize[iter] <-
        max(abs(log(normalized$marginal_scale_before)), na.rm = TRUE)
    }
    history$normalization_seconds[iter] <- as.numeric(difftime(Sys.time(), normalization_start, units = "secs"))

    keep_start <- Sys.time()
    if (iter > burn && ((iter - burn) %% thin == 0L)) {
      keep_F <- keep_F + F
      keep_alpha <- keep_alpha + alpha
      keep_Lambda <- keep_Lambda + Lambda
      keep_pi <- keep_pi + replace(pi_mat, is.na(pi_mat), 0)
      keep_mu <- keep_mu + replace(mu_mat, is.na(mu_mat), 0)
      keep_sig2 <- keep_sig2 + replace(sig2_mat, is.na(sig2_mat), 0)
      n_keep <- n_keep + 1L
      if (isTRUE(compute_parameter_ess) && n_keep <= keep_target) {
        alpha_trace[n_keep, ] <- alpha
        lambda_trace[n_keep, ] <- as.numeric(Lambda)
        pi_trace[n_keep, ] <- as.numeric(pi_mat)
        mu_trace[n_keep, ] <- as.numeric(mu_mat)
        sig2_trace[n_keep, ] <- as.numeric(sig2_mat)
      }
    }
    history$keep_draw_seconds[iter] <- as.numeric(difftime(Sys.time(), keep_start, units = "secs"))

    history$occupied_factor_classes[iter] <- sum(vapply(seq_len(H), function(h) {
      length(unique(C[, h]))
    }, integer(1L)))
    history$iteration_seconds[iter] <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

    if (isTRUE(verbose) && (iter == 1L || iter %% max(1L, floor(n_iter / 5)) == 0L)) {
      cat(sprintf(
        "Viroli probit Gibbs iter %d/%d | occupied factor classes=%d/%d\n",
        iter, n_iter, history$occupied_factor_classes[iter], sum(G)
      ))
    }
  }

  if (n_keep > 0L) {
    F_hat <- keep_F / n_keep
    alpha_hat <- keep_alpha / n_keep
    Lambda_hat <- keep_Lambda / n_keep
    pi_hat <- keep_pi / n_keep
    mu_hat <- keep_mu / n_keep
    sig2_hat <- keep_sig2 / n_keep
  } else {
    F_hat <- F
    alpha_hat <- alpha
    Lambda_hat <- Lambda
    pi_hat <- pi_mat
    mu_hat <- mu_mat
    sig2_hat <- sig2_mat
  }

  if (isTRUE(compute_parameter_ess) && n_keep > 0L) {
    trace_list <- list(
      alpha = alpha_trace[seq_len(n_keep), , drop = FALSE],
      lambda = lambda_trace[seq_len(n_keep), , drop = FALSE],
      pi = pi_trace[seq_len(n_keep), , drop = FALSE],
      mu = mu_trace[seq_len(n_keep), , drop = FALSE],
      sig2 = sig2_trace[seq_len(n_keep), , drop = FALSE]
    )
    ess_table <- viroli_compute_parameter_ess_table(trace_list)
    ess_summary <- viroli_summarize_parameter_ess(ess_table)
  } else {
    ess_table <- NULL
    ess_summary <- viroli_summarize_parameter_ess(NULL)
  }

  if (isTRUE(normalize_each_draw)) {
    posterior_mean_normalized <- normalize_viroli_draw(
      F = F_hat,
      alpha = alpha_hat,
      Lambda = Lambda_hat,
      pi_mat = pi_hat,
      mu_mat = mu_hat,
      sig2_mat = sig2_hat,
      G = G,
      min_scale = min_scale
    )
    F_hat <- posterior_mean_normalized$F
    alpha_hat <- posterior_mean_normalized$alpha
    Lambda_hat <- posterior_mean_normalized$Lambda
    pi_hat <- posterior_mean_normalized$pi
    mu_hat <- posterior_mean_normalized$mu
    sig2_hat <- posterior_mean_normalized$sig2
  }

  mixture_fits <- lapply(seq_len(H), function(h) {
    Gh <- G[h]
    list(
      pi = pi_hat[h, seq_len(Gh)],
      mu = mu_hat[h, seq_len(Gh)],
      var = sig2_hat[h, seq_len(Gh)],
      G = Gh,
      mixture_update = "gibbs"
    )
  })

  list(
    F_hat = F_hat,
    alpha_hat = alpha_hat,
    Lambda_hat = Lambda_hat,
    C = C,
    pi = pi_hat,
    mu = mu_hat,
    sig2 = sig2_hat,
    mixture_fits = mixture_fits,
    history = history,
    ess_table = ess_table,
    ess_summary = ess_summary,
    n_keep = n_keep,
    G = G,
    normalize_each_draw = normalize_each_draw,
    lambda_l1_penalty = lambda_l1_penalty,
    loading_prior = if (isTRUE(use_laplace_loading_prior)) "bayesian_lasso_scale_mixture" else "normal",
    seconds = proc.time()[["elapsed"]] - t0,
    seconds_per_iter = (proc.time()[["elapsed"]] - t0) / n_iter
  )
}
