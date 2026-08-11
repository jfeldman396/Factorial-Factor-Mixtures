#!/usr/bin/env Rscript

# ============================================================================
# Self-contained pretraining algorithm for the binary probit
# independent-mixture factor model
#
# This file contains the actual fitting code for the pretraining stage.  It is
# intentionally written in straightforward sequential R so that the algorithm is
# easy to read and check.
#
# Observed binary model:
#
#   X_ij = 1{Z_ij > 0}
#
# Working latent-Gaussian factor model:
#
#   Z_ij = alpha_j + lambda_j' f_i + epsilon_ij
#
# Independent mixture factor model:
#
#   f_ih independently follows sum_g pi_hg N(mu_hg, sigma_hg^2).
#
# Pretraining algorithm:
#
#   repeat augmentation:
#       1. Center the current latent Z by the current alpha estimate and
#          compute SVD scores S.
#       2. Rotate S by an orthogonal matrix R.
#       3. At the current orientation, fit a univariate Gaussian mixture to
#          each factor coordinate.  Either choose G_h by BIC or keep a supplied
#          fixed G_h.
#       4. Improve R using pairwise Givens rotations.  Each pairwise angle is
#          selected by grid search followed by local one-dimensional
#          optimization.
#       5. Estimate alpha and working loadings Lambda from Z ~ alpha + F Lambda'.
#       6. Fix Psi = I, then sample or average Z | X, F, alpha, Lambda, Psi.
#   end
#
# Main exported functions in this file:
#
#   simulate_binary_probit_independent_mixture_factor_model()
#   fit_binary_probit_pretraining()
#   summarize_binary_probit_pretraining()
#   summarize_binary_mixture_profiles()
# ============================================================================

# ----------------------------------------------------------------------------
# Basic numerical utilities
# ----------------------------------------------------------------------------

row_log_sum_exp <- function(A) {
  m <- apply(A, 1L, max)
  m + log(rowSums(exp(A - m)))
}

soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

random_orthogonal <- function(H) {
  qr_out <- qr(matrix(rnorm(H * H), H, H))
  Q <- qr.Q(qr_out)
  if (det(Q) < 0) Q[, 1L] <- -Q[, 1L]
  Q
}

project_to_orthogonal <- function(R) {
  dec <- svd(R)
  dec$u %*% t(dec$v)
}

resolve_workers <- function(workers = NULL) {
  # Pick a conservative default: leave one core free when possible.
  if (!is.null(workers)) return(max(1L, as.integer(workers)))
  n_cores <- parallel::detectCores(logical = TRUE)
  if (!is.finite(n_cores) || is.na(n_cores)) return(1L)
  max(1L, n_cores - 1L)
}

parallel_lapply <- function(X, FUN, ..., parallel = FALSE, workers = NULL) {
  # Small dependency-free parallel map.  On Unix-like systems, mclapply avoids
  # explicit cluster export.  On Windows, use a temporary PSOCK cluster.
  workers <- resolve_workers(workers)
  if (!isTRUE(parallel) || workers <= 1L || length(X) <= 1L) {
    return(lapply(X, FUN, ...))
  }

  if (.Platform$OS.type == "unix") {
    return(parallel::mclapply(
      X,
      FUN,
      ...,
      mc.cores = workers,
      mc.preschedule = TRUE
    ))
  }

  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::parLapply(cl, X, FUN, ...)
}

# ----------------------------------------------------------------------------
# Simulation utilities
# ----------------------------------------------------------------------------

standardize_mixture_parameters <- function(par) {
  # Normalize input weights and check that every component has a valid scale.
  stopifnot(length(par$pi) == length(par$mu), length(par$mu) == length(par$sd))

  pi <- par$pi / sum(par$pi)
  mu <- par$mu
  sd <- par$sd
  if (any(sd <= 0)) stop("Mixture standard deviations must be positive.")

  # Standardize the marginal mixture distribution so every factor coordinate
  # has mean zero and variance one before it enters the factor model.
  marginal_mean <- sum(pi * mu)
  marginal_var <- sum(pi * (sd^2 + mu^2)) - marginal_mean^2
  if (!is.finite(marginal_var) || marginal_var <= 0) {
    stop("The mixture must have positive marginal variance.")
  }

  # Order by the standardized means to make component labels reproducible.
  mu_std <- (mu - marginal_mean) / sqrt(marginal_var)
  sd_std <- sd / sqrt(marginal_var)
  ord <- order(mu_std)

  list(pi = pi[ord], mu = mu_std[ord], sd = sd_std[ord])
}

sample_standardized_mixture <- function(n, par) {
  # Draw mixture labels first, then draw factor values conditional on labels.
  par <- standardize_mixture_parameters(par)
  G <- length(par$pi)
  component <- sample.int(G, n, replace = TRUE, prob = par$pi)
  x <- rnorm(n, mean = par$mu[component], sd = par$sd[component])

  list(x = x, component = component, parameters = par)
}

make_block_sparse_loadings <- function(
    p,
    H,
    block_sizes = NULL,
    primary_range = c(0.8, 1.4),
    factor_strength = 1,
    random_signs = TRUE,
    cross_loading_prob = 0.1,
    cross_loading_range = c(0.04, 0.18)) {
  # If no block sizes are supplied, split variables as evenly as possible
  # across the H latent factors.
  if (is.null(block_sizes)) {
    block_sizes <- rep(floor(p / H), H)
    block_sizes[seq_len(p - sum(block_sizes))] <-
      block_sizes[seq_len(p - sum(block_sizes))] + 1L
  }
  if (sum(block_sizes) != p) stop("block_sizes must sum to p.")

  Lambda <- matrix(0, p, H)
  block_id <- rep(seq_len(H), times = block_sizes)

  for (j in seq_len(p)) {
    # Every variable has one primary loading in its assigned block.
    h <- block_id[j]
    sign_j <- if (random_signs) sample(c(-1, 1), 1L) else 1
    Lambda[j, h] <- sign_j * runif(1, primary_range[1], primary_range[2]) *
      factor_strength

    for (k in setdiff(seq_len(H), h)) {
      # Weak cross-loadings make the design less artificially block diagonal.
      if (runif(1) < cross_loading_prob) {
        Lambda[j, k] <- sample(c(-1, 1), 1L) *
          runif(1, cross_loading_range[1], cross_loading_range[2])
      }
    }
  }

  list(Lambda = Lambda, block_id = block_id, block_sizes = block_sizes)
}

make_very_sparse_loadings <- function(
    p,
    H,
    block_sizes = NULL,
    active_primary_fraction = 0.35,
    min_primary_loadings_per_factor = 4L,
    primary_range = c(0.9, 1.5),
    factor_strength = 1,
    random_signs = TRUE,
    cross_loading_prob = 0.01,
    cross_loading_range = c(0.03, 0.10)) {
  # A deliberately sparse loading design.  Only a subset of items in each block
  # has a primary loading; the remaining rows are noise variables with all-zero
  # loadings unless they receive a rare weak cross-loading.
  if (is.null(block_sizes)) {
    block_sizes <- rep(floor(p / H), H)
    block_sizes[seq_len(p - sum(block_sizes))] <-
      block_sizes[seq_len(p - sum(block_sizes))] + 1L
  }
  if (sum(block_sizes) != p) stop("block_sizes must sum to p.")

  Lambda <- matrix(0, p, H)
  block_id <- rep(seq_len(H), times = block_sizes)

  for (h in seq_len(H)) {
    block_rows <- which(block_id == h)
    n_active <- ceiling(length(block_rows) * active_primary_fraction)
    n_active <- max(min_primary_loadings_per_factor, n_active)
    n_active <- min(length(block_rows), n_active)
    active_rows <- sample(block_rows, n_active)

    for (j in active_rows) {
      sign_j <- if (random_signs) sample(c(-1, 1), 1L) else 1
      Lambda[j, h] <- sign_j * runif(1, primary_range[1], primary_range[2]) *
        factor_strength
    }
  }

  for (j in seq_len(p)) {
    primary_h <- block_id[j]
    for (k in setdiff(seq_len(H), primary_h)) {
      # Sparse weak cross-loadings keep the design from being perfectly simple
      # while preserving a genuinely sparse true Lambda.
      if (runif(1) < cross_loading_prob) {
        Lambda[j, k] <- sample(c(-1, 1), 1L) *
          runif(1, cross_loading_range[1], cross_loading_range[2])
      }
    }
  }

  list(Lambda = Lambda, block_id = block_id, block_sizes = block_sizes)
}

make_block_sparse_multisigned_loadings <- function(
    p,
    H,
    block_sizes = NULL,
    active_primary_fraction = 0.40,
    min_primary_loadings_per_factor = 4L,
    primary_range = c(1.0, 1.7),
    factor_strength = 1,
    random_signs = TRUE,
    cross_factors_per_item = 2L,
    cross_loading_range = c(0.25, 0.65),
    cross_factor_mode = c("neighbor", "random")) {
  # A sparse-but-interpretable design.  Each factor owns a block of items, but
  # only a subset of rows in each block are active.  Active rows have a strong
  # primary loading plus signed secondary loadings on a small number of other
  # factors.  Inactive rows remain exactly zero, so the design keeps genuine
  # row-level sparsity while allowing positive and negative multi-factor items.
  cross_factor_mode <- match.arg(cross_factor_mode)
  if (is.null(block_sizes)) {
    block_sizes <- rep(floor(p / H), H)
    block_sizes[seq_len(p - sum(block_sizes))] <-
      block_sizes[seq_len(p - sum(block_sizes))] + 1L
  }
  if (sum(block_sizes) != p) stop("block_sizes must sum to p.")

  Lambda <- matrix(0, p, H)
  block_id <- rep(seq_len(H), times = block_sizes)
  n_cross <- min(as.integer(cross_factors_per_item), max(0L, H - 1L))

  neighboring_factors <- function(h) {
    if (H <= 1L || n_cross == 0L) return(integer(0))
    offsets <- as.vector(rbind(seq_len(H - 1L), -seq_len(H - 1L)))
    out <- integer(0)
    for (offset in offsets) {
      candidate <- ((h - 1L + offset) %% H) + 1L
      if (!(candidate %in% out) && candidate != h) out <- c(out, candidate)
      if (length(out) == n_cross) break
    }
    out
  }

  for (h in seq_len(H)) {
    block_rows <- which(block_id == h)
    n_active <- ceiling(length(block_rows) * active_primary_fraction)
    n_active <- max(min_primary_loadings_per_factor, n_active)
    n_active <- min(length(block_rows), n_active)
    active_rows <- sample(block_rows, n_active)

    for (j in active_rows) {
      primary_sign <- if (random_signs) sample(c(-1, 1), 1L) else 1
      Lambda[j, h] <- primary_sign *
        runif(1, primary_range[1], primary_range[2]) * factor_strength

      cross_factors <- if (cross_factor_mode == "neighbor") {
        neighboring_factors(h)
      } else {
        sample(setdiff(seq_len(H), h), n_cross)
      }

      for (k in cross_factors) {
        Lambda[j, k] <- sample(c(-1, 1), 1L) *
          runif(1, cross_loading_range[1], cross_loading_range[2]) *
          factor_strength
      }
    }
  }

  list(Lambda = Lambda, block_id = block_id, block_sizes = block_sizes)
}

make_block_cross_loadings <- function(
    p,
    H,
    block_sizes = NULL,
    primary_range = c(0.8, 1.3),
    factor_strength = 1,
    random_signs = TRUE,
    cross_factors_per_item = 2L,
    cross_loading_range = c(0.25, 0.55),
    cross_factor_mode = c("neighbor", "random")) {
  # A structured cross-loading design.  Every item has one primary loading in
  # its assigned block plus a fixed number of cross-factor loadings.  With the
  # default cross_factors_per_item = 2, each row loads on three factors total
  # whenever H >= 3.
  cross_factor_mode <- match.arg(cross_factor_mode)
  if (is.null(block_sizes)) {
    block_sizes <- rep(floor(p / H), H)
    block_sizes[seq_len(p - sum(block_sizes))] <-
      block_sizes[seq_len(p - sum(block_sizes))] + 1L
  }
  if (sum(block_sizes) != p) stop("block_sizes must sum to p.")

  Lambda <- matrix(0, p, H)
  block_id <- rep(seq_len(H), times = block_sizes)
  n_cross <- min(as.integer(cross_factors_per_item), max(0L, H - 1L))

  neighboring_factors <- function(h) {
    if (H <= 1L || n_cross == 0L) return(integer(0))
    offsets <- as.vector(rbind(seq_len(H - 1L), -seq_len(H - 1L)))
    out <- integer(0)
    for (offset in offsets) {
      candidate <- ((h - 1L + offset) %% H) + 1L
      if (!(candidate %in% out) && candidate != h) out <- c(out, candidate)
      if (length(out) == n_cross) break
    }
    out
  }

  for (j in seq_len(p)) {
    h <- block_id[j]
    sign_j <- if (random_signs) sample(c(-1, 1), 1L) else 1
    Lambda[j, h] <- sign_j * runif(1, primary_range[1], primary_range[2]) *
      factor_strength

    cross_factors <- if (cross_factor_mode == "neighbor") {
      neighboring_factors(h)
    } else {
      sample(setdiff(seq_len(H), h), n_cross)
    }

    for (k in cross_factors) {
      Lambda[j, k] <- sample(c(-1, 1), 1L) *
        runif(1, cross_loading_range[1], cross_loading_range[2])
    }
  }

  list(Lambda = Lambda, block_id = block_id, block_sizes = block_sizes)
}

simulate_binary_probit_independent_mixture_factor_model <- function(
    n = 500L,
    p = 80L,
    mixture_params,
    block_sizes = NULL,
    loading_design = c("block_sparse", "very_sparse", "block_sparse_multisigned", "block_cross"),
    active_primary_fraction = 0.35,
    min_primary_loadings_per_factor = 4L,
    cross_factors_per_item = 2L,
    cross_factor_mode = c("neighbor", "random"),
    primary_range = c(0.8, 1.4),
    factor_strength = 1,
    random_loading_signs = TRUE,
    cross_loading_prob = 0.1,
    cross_loading_range = c(0.04, 0.18),
    noise_sd = 1,
    rotate_loadings = TRUE,
    seed = 1L) {
  set.seed(seed)
  loading_design <- match.arg(loading_design)
  cross_factor_mode <- match.arg(cross_factor_mode)
  H <- length(mixture_params)

  # Generate a loading matrix with either dense primary blocks or a deliberately
  # sparse subset of active indicators per factor.
  loading_out <- switch(
    loading_design,
    block_sparse = make_block_sparse_loadings(
      p = p,
      H = H,
      block_sizes = block_sizes,
      primary_range = primary_range,
      factor_strength = factor_strength,
      random_signs = random_loading_signs,
      cross_loading_prob = cross_loading_prob,
      cross_loading_range = cross_loading_range
    ),
    very_sparse = make_very_sparse_loadings(
      p = p,
      H = H,
      block_sizes = block_sizes,
      active_primary_fraction = active_primary_fraction,
      min_primary_loadings_per_factor = min_primary_loadings_per_factor,
      primary_range = primary_range,
      factor_strength = factor_strength,
      random_signs = random_loading_signs,
      cross_loading_prob = cross_loading_prob,
      cross_loading_range = cross_loading_range
    ),
    block_sparse_multisigned = make_block_sparse_multisigned_loadings(
      p = p,
      H = H,
      block_sizes = block_sizes,
      active_primary_fraction = active_primary_fraction,
      min_primary_loadings_per_factor = min_primary_loadings_per_factor,
      primary_range = primary_range,
      factor_strength = factor_strength,
      random_signs = random_loading_signs,
      cross_factors_per_item = cross_factors_per_item,
      cross_loading_range = cross_loading_range,
      cross_factor_mode = cross_factor_mode
    ),
    block_cross = make_block_cross_loadings(
      p = p,
      H = H,
      block_sizes = block_sizes,
      primary_range = primary_range,
      factor_strength = factor_strength,
      random_signs = random_loading_signs,
      cross_factors_per_item = cross_factors_per_item,
      cross_loading_range = cross_loading_range,
      cross_factor_mode = cross_factor_mode
    )
  )

  Lambda <- loading_out$Lambda
  F <- matrix(NA_real_, n, H)
  component <- matrix(NA_integer_, n, H)
  standardized_params <- vector("list", H)

  # Draw independent non-Gaussian factor coordinates and their true labels.
  for (h in seq_len(H)) {
    draw_h <- sample_standardized_mixture(n, mixture_params[[h]])
    F[, h] <- draw_h$x
    component[, h] <- draw_h$component
    standardized_params[[h]] <- draw_h$parameters
  }

  if (isTRUE(rotate_loadings)) {
    # Rotating the loadings hides the block structure from the SVD basis while
    # preserving the same signal subspace.  The estimation algorithm must then
    # recover a rotation toward the non-Gaussian factor axes.
    Q <- random_orthogonal(H)
    Lambda <- Lambda %*% Q
  } else {
    Q <- diag(H)
  }

  if (length(noise_sd) == 1L) noise_sd <- rep(noise_sd, p)
  E <- matrix(rnorm(n * p), n, p)
  E <- sweep(E, 2L, noise_sd, "*")

  # Create latent probit variables and threshold them at zero.
  Z_latent <- F %*% t(Lambda) + E
  X_binary <- 1L * (Z_latent > 0)

  list(
    X_binary = X_binary,
    Z_latent = Z_latent,
    F = F,
    Lambda = Lambda,
    loading_rotation = Q,
    component = component,
    mixture_params = standardized_params,
    Psi = diag(noise_sd^2, p),
    block_id = loading_out$block_id,
    block_sizes = loading_out$block_sizes,
    prevalence = colMeans(X_binary)
  )
}

# ----------------------------------------------------------------------------
# Univariate Gaussian mixture fitting and BIC selection
# ----------------------------------------------------------------------------

log_dmix_1d <- function(x, fit) {
  x <- as.numeric(x)
  G <- length(fit$pi)

  # Evaluate every component density on the log scale for numerical stability.
  log_components <- vapply(seq_len(G), function(g) {
    log(pmax(fit$pi[g], 1e-12)) +
      dnorm(x, mean = fit$mu[g], sd = sqrt(fit$var[g]), log = TRUE)
  }, numeric(length(x)))

  # The one-component case is already a vector; otherwise use log-sum-exp.
  if (G == 1L) return(log_components)
  row_log_sum_exp(log_components)
}

mixture_responsibilities <- function(x, fit) {
  x <- as.numeric(x)
  G <- length(fit$pi)

  # Posterior probabilities are proportional to pi_g phi(x; mu_g, var_g).
  log_resp <- vapply(seq_len(G), function(g) {
    log(pmax(fit$pi[g], 1e-12)) +
      dnorm(x, mean = fit$mu[g], sd = sqrt(fit$var[g]), log = TRUE)
  }, numeric(length(x)))

  if (G == 1L) {
    return(matrix(1, nrow = length(x), ncol = 1L))
  }

  # Normalize rows on the log scale to avoid underflow in separated mixtures.
  log_den <- row_log_sum_exp(log_resp)
  exp(log_resp - log_den)
}

fit_gmm_1d <- function(
    x,
    G = 2L,
    n_starts = 8L,
    max_iter = 20L,
    tol = 1e-8,
    min_var = 1e-3,
    min_weight = 1e-4,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    init = NULL) {
  mixture_update <- match.arg(mixture_update)
  x <- as.numeric(x)
  n <- length(x)
  G <- as.integer(G)
  if (G < 1L || n <= G) stop("Require 1 <= G < length(x).")

  sample_var <- var(x)
  if (!is.finite(sample_var) || sample_var <= 0) {
    stop("The input to fit_gmm_1d must have positive variance.")
  }

  update_component_parameters <- function(resp, pi_old = NULL) {
    # M-step for a univariate Gaussian mixture.  In MLE mode this is the usual
    # weighted mean/variance update.  In MAP mode it is the conjugate
    # Normal-Inverse-Gamma update:
    #
    #   mu_g | sigma_g^2 ~ N(m0, sigma_g^2 / kappa0)
    #   sigma_g^2        ~ Inv-Gamma(a0, b0)
    #
    # The b0 / sigma_g^2 term makes tiny variances unattractive, which is the
    # regularization we want for the G = 3 collapse pathology.
    nk <- colSums(resp) + 1e-12

    if (mixture_update == "map") {
      alpha <- rep(weight_prior_alpha, G)
      if (length(weight_prior_alpha) == G) alpha <- weight_prior_alpha
      if (any(alpha < 1)) stop("weight_prior_alpha must be >= 1 for MAP weights.")
      pi_g <- pmax(nk + alpha - 1, min_weight)
      pi_g <- pi_g / sum(pi_g)
    } else {
      pi_g <- pmax(nk / n, min_weight)
      pi_g <- pi_g / sum(pi_g)
    }

    xbar <- colSums(resp * x) / nk
    centered <- sweep(matrix(x, nrow = n, ncol = G), 2L, xbar, "-")
    ss <- colSums(resp * centered^2)

    if (mixture_update == "map") {
      kappa0 <- pmax(mu_prior_kappa, 0)
      kappa_n <- kappa0 + nk
      mu_g <- (kappa0 * mu_prior_mean + nk * xbar) / kappa_n
      shape_n <- var_prior_shape + nk / 2
      scale_n <- var_prior_scale + 0.5 * ss +
        (kappa0 * nk * (xbar - mu_prior_mean)^2) / (2 * kappa_n)
      var_g <- pmax(scale_n / (shape_n + 1), min_var)
    } else {
      mu_g <- xbar
      var_g <- pmax(ss / nk, min_var)
    }

    list(pi = pi_g, mu = mu_g, var = var_g)
  }

  if (G == 1L) {
    # With one Gaussian component, use the closed-form M-step.
    updated <- update_component_parameters(matrix(1, nrow = n, ncol = 1L))
    fit <- list(
      pi = 1,
      mu = updated$mu,
      var = updated$var,
      G = 1L,
      converged = TRUE
    )
    fit$loglik <- sum(log_dmix_1d(x, fit))
    fit$mixture_update <- mixture_update
    return(fit)
  }

  best <- NULL

  for (s in seq_len(n_starts)) {
    # Start 1 can warm-start from the previous iteration if the component count
    # matches.  Otherwise we use deterministic quantiles plus random restarts.
    use_previous <- s == 1L && !is.null(init) && length(init$pi) == G

    if (use_previous) {
      pi_g <- init$pi
      mu_g <- init$mu
      var_g <- init$var
    } else if (s == 1L) {
      mu_g <- as.numeric(quantile(x, probs = (seq_len(G) - 0.5) / G,
                                  names = FALSE))
      pi_g <- rep(1 / G, G)
      var_g <- rep(sample_var, G)
    } else {
      mu_g <- sample(x, size = G, replace = FALSE)
      pi_g <- rep(1 / G, G)
      var_g <- rep(sample_var, G)
    }

    pi_g <- pmax(pi_g, min_weight)
    pi_g <- pi_g / sum(pi_g)
    var_g <- pmax(var_g, min_var)
    old_loglik <- -Inf
    converged <- FALSE

    for (iter in seq_len(max_iter)) {
      # E-step: compute posterior component probabilities for every x_i.
      log_resp <- vapply(seq_len(G), function(g) {
        log(pi_g[g]) + dnorm(x, mean = mu_g[g], sd = sqrt(var_g[g]), log = TRUE)
      }, numeric(n))

      log_den <- row_log_sum_exp(log_resp)
      loglik <- sum(log_den)
      resp <- exp(log_resp - log_den)

      # M-step: update weights, means, and variances using either MLE or MAP.
      updated <- update_component_parameters(resp)
      pi_g <- updated$pi
      mu_g <- updated$mu
      var_g <- updated$var

      # Stop when the log likelihood is no longer changing materially.
      if (is.finite(old_loglik) &&
          abs(loglik - old_loglik) <= tol * (1 + abs(old_loglik))) {
        converged <- TRUE
        break
      }
      old_loglik <- loglik
    }

    fit <- list(
      pi = pi_g,
      mu = mu_g,
      var = var_g,
      G = G,
      converged = converged,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha
    )
    fit$loglik <- sum(log_dmix_1d(x, fit))

    # Sort labels by mean so component 1 is the leftmost cluster, etc.
    ord <- order(fit$mu)
    fit$pi <- fit$pi[ord]
    fit$mu <- fit$mu[ord]
    fit$var <- fit$var[ord]

    if (is.null(best) || fit$loglik > best$loglik) best <- fit
  }

  best
}

select_gmm_bic <- function(
    x,
    G_max = 5L,
    n_starts = 8L,
    max_iter = 20L,
    mixture_penalty_multiplier = 1,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    previous_fit = NULL) {
  mixture_update <- match.arg(mixture_update)
  x <- as.numeric(x)
  n <- length(x)
  G_max <- min(as.integer(G_max), max(1L, floor(n / 10)))

  fits <- vector("list", G_max)
  bic <- numeric(G_max)

  for (G in seq_len(G_max)) {
    # Fit each candidate G.  If the previous fit had this same G, use it as a
    # warm start so refinement does not ignore useful pretraining information.
    init_G <- NULL
    if (!is.null(previous_fit) && length(previous_fit$pi) == G) {
      init_G <- previous_fit
    }

    # Univariate Gaussian mixture degrees of freedom:
    # (G - 1) weights + G means + G variances = 3G - 1.
    fits[[G]] <- fit_gmm_1d(
      x,
      G = G,
      n_starts = n_starts,
      max_iter = max_iter,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      init = init_G
    )
    n_parameters <- 3 * G - 1
    bic[G] <- -2 * fits[[G]]$loglik +
      mixture_penalty_multiplier * n_parameters * log(n)
  }

  # Smaller BIC is better; return both the winner and the full BIC path.
  best_G <- which.min(bic)
  best <- fits[[best_G]]
  best$bic <- bic[best_G]
  best$all_bic <- bic
  best$all_fits <- fits
  best
}

fit_column_mixtures_bic <- function(
    F,
    G_max = 5L,
    n_starts = 8L,
    max_iter = 20L,
    mixture_penalty_multiplier = 1,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    previous_fits = NULL,
    parallel = FALSE,
    workers = NULL) {
  mixture_update <- match.arg(mixture_update)
  H <- ncol(F)

  # Because the model assumes independent factor coordinates, each marginal
  # mixture can be fit separately.
  parallel_lapply(seq_len(H), function(h) {
    prev <- if (!is.null(previous_fits)) previous_fits[[h]] else NULL
    select_gmm_bic(
      F[, h],
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
      previous_fit = prev
    )
  }, parallel = parallel, workers = workers)
}

normalize_G_fixed <- function(G_fixed, H) {
  # Accept either one common G for all factors or one G_h per factor.
  if (is.null(G_fixed)) {
    stop("G_fixed must be supplied when pretrain_G_selection = 'fixed'.")
  }

  if (length(G_fixed) == 1L) G_fixed <- rep(G_fixed, H)
  if (length(G_fixed) != H) {
    stop("G_fixed must have length 1 or length H.")
  }

  G_fixed <- as.integer(G_fixed)
  if (any(G_fixed < 1L)) stop("All fixed component counts must be positive.")
  G_fixed
}

fit_column_mixtures_fixed_G <- function(
    F,
    G_fixed,
    n_starts = 8L,
    max_iter = 20L,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    previous_fits = NULL,
    parallel = FALSE,
    workers = NULL) {
  mixture_update <- match.arg(mixture_update)
  H <- ncol(F)
  G_fixed <- normalize_G_fixed(G_fixed, H)

  parallel_lapply(seq_len(H), function(h) {
    # Fit the hth marginal mixture using the user-specified component count.
    prev <- if (!is.null(previous_fits) &&
                length(previous_fits[[h]]$pi) == G_fixed[h]) {
      previous_fits[[h]]
    } else {
      NULL
    }

    fit_gmm_1d(
      F[, h],
      G = G_fixed[h],
      n_starts = n_starts,
      max_iter = max_iter,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      init = prev
    )
  }, parallel = parallel, workers = workers)
}

fit_column_mixtures <- function(
    F,
    G_selection = c("bic", "fixed"),
    G_max = 5L,
    G_fixed = NULL,
    n_starts = 8L,
    max_iter = 20L,
    mixture_penalty_multiplier = 1,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    previous_fits = NULL,
    parallel = FALSE,
    workers = NULL) {
  G_selection <- match.arg(G_selection)
  mixture_update <- match.arg(mixture_update)

  if (G_selection == "fixed") {
    return(fit_column_mixtures_fixed_G(
      F = F,
      G_fixed = G_fixed,
      n_starts = n_starts,
      max_iter = max_iter,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      previous_fits = previous_fits,
      parallel = parallel,
      workers = workers
    ))
  }

  fit_column_mixtures_bic(
    F = F,
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
    previous_fits = previous_fits,
    parallel = parallel,
    workers = workers
  )
}

mixture_loglik_total <- function(F, fits) {
  sum(vapply(seq_len(ncol(F)), function(h) {
    sum(log_dmix_1d(F[, h], fits[[h]]))
  }, numeric(1)))
}

mixture_bic_score <- function(F, fits, mixture_penalty_multiplier = 1) {
  n <- nrow(F)

  # This is the rotation objective: marginal mixture log likelihood minus the
  # BIC-style complexity penalty for all selected univariate mixtures.
  penalty <- sum(vapply(fits, function(fit) {
    G <- length(fit$pi)
    0.5 * mixture_penalty_multiplier * (3 * G - 1) * log(n)
  }, numeric(1)))
  mixture_loglik_total(F, fits) - penalty
}

# ----------------------------------------------------------------------------
# Givens-rotation mixture-ICA step
# ----------------------------------------------------------------------------

apply_pair_rotation <- function(F, a, b, theta) {
  # Right-multiply columns a and b of F by a two-dimensional Givens rotation.
  ca <- cos(theta)
  sa <- sin(theta)
  Fa <- F[, a]
  Fb <- F[, b]
  F[, a] <- ca * Fa + sa * Fb
  F[, b] <- -sa * Fa + ca * Fb
  F
}

pair_rotation_matrix <- function(H, a, b, theta) {
  # Build the full H x H matrix corresponding to the same pairwise rotation.
  G <- diag(H)
  ca <- cos(theta)
  sa <- sin(theta)
  G[a, a] <- ca
  G[b, b] <- ca
  G[a, b] <- -sa
  G[b, a] <- sa
  G
}

pair_negative_loglik <- function(theta, Fa, Fb, fit_a, fit_b) {
  # Rotate only the selected pair and score the two affected marginal mixtures.
  # All other coordinates are unchanged, so they need not be recomputed.
  ca <- cos(theta)
  sa <- sin(theta)
  Fa_new <- ca * Fa + sa * Fb
  Fb_new <- -sa * Fa + ca * Fb

  -sum(log_dmix_1d(Fa_new, fit_a)) -
    sum(log_dmix_1d(Fb_new, fit_b))
}

optimize_pair_angle <- function(
    Fa,
    Fb,
    fit_a,
    fit_b,
    interval = c(-pi / 2, pi / 2),
    grid_size = 21L) {
  # The pairwise mixture likelihood is not guaranteed to be globally unimodal.
  # We first search a coarse grid, then locally optimize around the best grid
  # point.  This is cheap because it is only one-dimensional.
  grid <- seq(interval[1], interval[2], length.out = grid_size)
  values <- vapply(grid, pair_negative_loglik, numeric(1),
                   Fa = Fa, Fb = Fb, fit_a = fit_a, fit_b = fit_b)

  # Restrict the local optimizer to the grid neighborhood of the best angle.
  k <- which.min(values)
  local_interval <- if (k == 1L) {
    grid[1:2]
  } else if (k == length(grid)) {
    grid[(length(grid) - 1L):length(grid)]
  } else {
    grid[c(k - 1L, k + 1L)]
  }

  local <- optimize(
    f = pair_negative_loglik,
    interval = local_interval,
    Fa = Fa,
    Fb = Fb,
    fit_a = fit_a,
    fit_b = fit_b,
    tol = 1e-8
  )

  # Keep the no-rotation angle if the local optimizer does not improve it.
  zero_value <- pair_negative_loglik(0, Fa, Fb, fit_a, fit_b)
  if (zero_value <= local$objective) {
    list(theta = 0, objective = zero_value)
  } else {
    list(theta = local$minimum, objective = local$objective)
  }
}

make_disjoint_pairs <- function(H, pairing = c("random", "adjacent"), seed = NULL) {
  # Construct one matching of disjoint factor pairs.  Disjoint pair rotations
  # commute, so their one-dimensional angle searches can be run in parallel.
  pairing <- match.arg(pairing)
  if (H < 2L) return(matrix(integer(0), ncol = 2L))

  index <- seq_len(H)
  if (pairing == "random") {
    if (!is.null(seed)) set.seed(seed)
    index <- sample(index)
  }

  n_pairs <- floor(H / 2)
  matrix(index[seq_len(2L * n_pairs)], ncol = 2L, byrow = TRUE)
}

disjoint_mixture_sweep <- function(
    S,
    R,
    fits,
    grid_size = 21L,
    pairing = c("random", "adjacent"),
    seed = NULL,
    parallel = FALSE,
    workers = NULL) {
  pairing <- match.arg(pairing)
  workers <- resolve_workers(workers)
  H <- ncol(S)
  F <- S %*% R
  if (H < 2L) return(list(F = F, R = R, angles = numeric(0), pairs = matrix(integer(0), ncol = 2L)))

  pairs <- make_disjoint_pairs(H, pairing = pairing, seed = seed)
  angles <- parallel_lapply(
    seq_len(nrow(pairs)),
    function(r) {
      a <- pairs[r, 1L]
      b <- pairs[r, 2L]
      optimize_pair_angle(
        Fa = F[, a],
        Fb = F[, b],
        fit_a = fits[[a]],
        fit_b = fits[[b]],
        grid_size = grid_size
      )$theta
    },
    parallel = parallel,
    workers = workers
  )
  angles <- unlist(angles, use.names = FALSE)

  # Apply the accepted disjoint rotations as one block update.  Since no pair
  # shares a coordinate, these rotations commute and do not need sequential
  # re-scoring within this sweep.
  G_all <- diag(H)
  for (r in seq_len(nrow(pairs))) {
    theta <- angles[r]
    if (abs(theta) > 0) {
      a <- pairs[r, 1L]
      b <- pairs[r, 2L]
      ca <- cos(theta)
      sa <- sin(theta)
      G_all[a, a] <- ca
      G_all[b, b] <- ca
      G_all[a, b] <- -sa
      G_all[b, a] <- sa
    }
  }

  R <- project_to_orthogonal(R %*% G_all)
  list(F = S %*% R, R = R, angles = angles, pairs = pairs)
}

mixture_profile_scores <- function(F, fits) {
  # Summarize each subject's posterior mixture profile on each axis.  Because
  # F remains orthogonal under rotations, value correlations are not useful for
  # finding coupled axes; dependence in the inferred mixture profiles is.
  F <- as.matrix(F)
  H <- ncol(F)
  scores <- matrix(0, nrow(F), H)
  for (h in seq_len(H)) {
    resp <- mixture_responsibilities(F[, h], fits[[h]])
    score_h <- as.numeric(resp %*% fits[[h]]$mu)
    sd_h <- sd(score_h)
    if (is.finite(sd_h) && sd_h > 1e-10) {
      scores[, h] <- (score_h - mean(score_h)) / sd_h
    }
  }
  scores
}

select_promising_pairs <- function(
    F,
    fits,
    max_pairs = NULL,
    fraction = 0.15,
    min_score = 0) {
  H <- ncol(F)
  if (H < 2L) return(matrix(integer(0), ncol = 2L))

  all_pairs <- t(combn(H, 2L))
  profile_scores <- mixture_profile_scores(F, fits)
  pair_scores <- apply(all_pairs, 1L, function(pair) {
    z <- suppressWarnings(cor(profile_scores[, pair[1L]], profile_scores[, pair[2L]]))
    if (is.finite(z)) abs(z) else 0
  })

  if (is.null(max_pairs)) {
    max_pairs <- max(H, ceiling(fraction * nrow(all_pairs)))
  }
  max_pairs <- max(1L, min(as.integer(max_pairs), nrow(all_pairs)))
  keep_order <- order(pair_scores, decreasing = TRUE)
  keep_order <- keep_order[pair_scores[keep_order] >= min_score]
  if (length(keep_order) == 0L) return(matrix(integer(0), ncol = 2L))
  all_pairs[head(keep_order, max_pairs), , drop = FALSE]
}

selected_pair_mixture_sweep <- function(S, R, fits, pairs, grid_size = 21L) {
  H <- ncol(S)
  F <- S %*% R
  if (H < 2L || is.null(pairs) || nrow(pairs) == 0L) {
    return(list(F = F, R = R, angles = numeric(0), pairs = matrix(integer(0), ncol = 2L)))
  }

  angles <- numeric(nrow(pairs))
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]
    b <- pairs[r, 2L]
    opt <- optimize_pair_angle(
      Fa = F[, a],
      Fb = F[, b],
      fit_a = fits[[a]],
      fit_b = fits[[b]],
      grid_size = grid_size
    )
    theta <- opt$theta
    if (abs(theta) > 0) {
      F <- apply_pair_rotation(F, a, b, theta)
      R <- R %*% pair_rotation_matrix(H, a, b, theta)
    }
    angles[r] <- theta
  }

  R <- project_to_orthogonal(R)
  list(F = S %*% R, R = R, angles = angles, pairs = pairs)
}

jacobi_mixture_sweep <- function(S, R, fits, grid_size = 21L) {
  H <- ncol(S)
  F <- S %*% R
  if (H < 2L) return(list(F = F, R = R, angles = numeric(0)))

  # A sweep updates each coordinate pair once while holding mixture parameters
  # fixed.  The mixtures are re-estimated after the sweep.
  pairs <- t(combn(H, 2L))
  angles <- numeric(nrow(pairs))

  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1L]
    b <- pairs[r, 2L]

    # Choose the best angle for this pair under the current marginal mixtures.
    opt <- optimize_pair_angle(
      Fa = F[, a],
      Fb = F[, b],
      fit_a = fits[[a]],
      fit_b = fits[[b]],
      grid_size = grid_size
    )

    theta <- opt$theta
    if (abs(theta) > 0) {
      F <- apply_pair_rotation(F, a, b, theta)
      R <- R %*% pair_rotation_matrix(H, a, b, theta)
    }
    angles[r] <- theta
  }

  # Numerical projection removes floating-point drift from repeated rotations.
  R <- project_to_orthogonal(R)
  list(F = S %*% R, R = R, angles = angles, pairs = pairs)
}

estimate_mixture_ica_unknown_G <- function(
    S,
    G_selection = c("bic", "fixed"),
    G_max = 5L,
    G_fixed = NULL,
    n_random_starts = 1L,
    max_outer = 4L,
    n_mix_starts = 3L,
    mixture_max_iter = 20L,
    mixture_penalty_multiplier = 1,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    grid_size = 21L,
    rotation_sweep = c("full", "disjoint", "multi_disjoint", "promising", "hybrid"),
    disjoint_pairing = c("random", "adjacent"),
    n_disjoint_rounds = 4L,
    promising_max_pairs = NULL,
    promising_fraction = 0.15,
    promising_min_score = 0,
    seed = 1L,
    parallel = FALSE,
    workers = NULL,
    verbose = TRUE) {
  G_selection <- match.arg(G_selection)
  mixture_update <- match.arg(mixture_update)
  rotation_sweep <- match.arg(rotation_sweep)
  disjoint_pairing <- match.arg(disjoint_pairing)
  n_disjoint_rounds <- max(1L, as.integer(n_disjoint_rounds))
  workers <- resolve_workers(workers)
  set.seed(seed)
  S <- as.matrix(S)
  H <- ncol(S)

  if (G_selection == "fixed") {
    G_fixed <- normalize_G_fixed(G_fixed, H)
  }

  # Use identity plus random orthogonal starts to reduce local optimum risk.
  starts <- list(identity = diag(H))
  for (s in seq_len(n_random_starts)) {
    starts[[paste0("random_", s)]] <- random_orthogonal(H)
  }

  fit_one_start <- function(start_index, use_parallel_mixtures = FALSE) {
    set.seed(seed + 100L * start_index)
    start_name <- names(starts)[start_index]
    if (verbose) message("  rotation start: ", start_name)
    R <- starts[[start_name]]
    F <- S %*% R

    # Initial mixture fits at the starting orientation.
    fits <- fit_column_mixtures(
      F,
      G_selection = G_selection,
      G_max = G_max,
      G_fixed = G_fixed,
      n_starts = n_mix_starts,
      max_iter = mixture_max_iter,
      mixture_penalty_multiplier = mixture_penalty_multiplier,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      parallel = use_parallel_mixtures,
      workers = workers
    )

    refit_current <- function(F_current, previous_fits) {
      fit_column_mixtures(
        F_current,
        G_selection = G_selection,
        G_max = G_max,
        G_fixed = G_fixed,
        n_starts = n_mix_starts,
        max_iter = mixture_max_iter,
        mixture_penalty_multiplier = mixture_penalty_multiplier,
        mixture_update = mixture_update,
        mu_prior_mean = mu_prior_mean,
        mu_prior_kappa = mu_prior_kappa,
        var_prior_shape = var_prior_shape,
        var_prior_scale = var_prior_scale,
        weight_prior_alpha = weight_prior_alpha,
        previous_fits = previous_fits,
        parallel = use_parallel_mixtures,
        workers = workers
      )
    }

    for (outer in seq_len(max_outer)) {
      # Alternate rotation updates with re-estimating the one-dimensional
      # mixtures.  Component counts are either fixed by G_fixed or reselected
      # by BIC, depending on G_selection.
      if (rotation_sweep == "full") {
        sweep_out <- jacobi_mixture_sweep(
          S = S,
          R = R,
          fits = fits,
          grid_size = grid_size
        )
        R <- sweep_out$R
        F <- sweep_out$F
        fits <- refit_current(F, fits)
      } else if (rotation_sweep == "disjoint") {
        sweep_out <- disjoint_mixture_sweep(
          S = S,
          R = R,
          fits = fits,
          grid_size = grid_size,
          pairing = disjoint_pairing,
          seed = seed + 10000L * start_index + outer,
          parallel = !use_parallel_mixtures && isTRUE(parallel),
          workers = workers
        )
        R <- sweep_out$R
        F <- sweep_out$F
        fits <- refit_current(F, fits)
      } else if (rotation_sweep %in% c("multi_disjoint", "hybrid")) {
        for (round_id in seq_len(n_disjoint_rounds)) {
          sweep_out <- disjoint_mixture_sweep(
            S = S,
            R = R,
            fits = fits,
            grid_size = grid_size,
            pairing = disjoint_pairing,
            seed = seed + 10000L * start_index + 100L * outer + round_id,
            parallel = !use_parallel_mixtures && isTRUE(parallel),
            workers = workers
          )
          R <- sweep_out$R
          F <- sweep_out$F
          fits <- refit_current(F, fits)
        }
        if (rotation_sweep == "hybrid") {
          promising_pairs <- select_promising_pairs(
            F,
            fits,
            max_pairs = promising_max_pairs,
            fraction = promising_fraction,
            min_score = promising_min_score
          )
          sweep_out <- selected_pair_mixture_sweep(
            S = S,
            R = R,
            fits = fits,
            pairs = promising_pairs,
            grid_size = grid_size
          )
          R <- sweep_out$R
          F <- sweep_out$F
          fits <- refit_current(F, fits)
        }
      } else if (rotation_sweep == "promising") {
        promising_pairs <- select_promising_pairs(
          F,
          fits,
          max_pairs = promising_max_pairs,
          fraction = promising_fraction,
          min_score = promising_min_score
        )
        sweep_out <- selected_pair_mixture_sweep(
          S = S,
          R = R,
          fits = fits,
          pairs = promising_pairs,
          grid_size = grid_size
        )
        R <- sweep_out$R
        F <- sweep_out$F
        fits <- refit_current(F, fits)
      } else {
        sweep_out <- disjoint_mixture_sweep(
          S = S,
          R = R,
          fits = fits,
          grid_size = grid_size,
          pairing = disjoint_pairing,
          seed = seed + 10000L * start_index + outer,
          parallel = !use_parallel_mixtures && isTRUE(parallel),
          workers = workers
        )
        R <- sweep_out$R
        F <- sweep_out$F
        fits <- refit_current(F, fits)
      }
    }

    # Rotation starts are compared by the profiled marginal mixture likelihood.
    # When G is fixed, there is no model-selection penalty to apply here.
    loglik <- mixture_loglik_total(F, fits)
    penalized_score <- mixture_bic_score(F, fits, mixture_penalty_multiplier)
    list(
      F_hat = F,
      R = R,
      fits = fits,
      G_hat = vapply(fits, function(z) length(z$pi), integer(1)),
      loglik = loglik,
      penalized_score = penalized_score,
      G_selection = G_selection,
      mixture_update = mixture_update,
      rotation_sweep = rotation_sweep,
      disjoint_pairing = disjoint_pairing,
      n_disjoint_rounds = n_disjoint_rounds
    )
  }

  # Parallelize over rotation starts when there is more than one start.  If
  # there is only one start, use workers for the independent marginal mixtures.
  use_start_parallel <- isTRUE(parallel) && workers > 1L && length(starts) > 1L
  results <- parallel_lapply(
    seq_along(starts),
    function(s) fit_one_start(s, use_parallel_mixtures = !use_start_parallel && isTRUE(parallel)),
    parallel = use_start_parallel,
    workers = workers
  )
  names(results) <- names(starts)

  # Return the best rotation start and retain scores for diagnostics.  The
  # rotation initialization itself is likelihood based; BIC-style penalties are
  # only used inside the optional unknown-G mixture fitting path.
  scores <- vapply(results, function(z) z$loglik, numeric(1))
  best <- which.max(scores)
  out <- results[[best]]
  out$start_name <- names(results)[best]
  out$all_start_scores <- scores
  out$all_start_loglik <- scores
  out$all_start_penalized_scores <- vapply(results, function(z) z$penalized_score, numeric(1))
  out
}

# ----------------------------------------------------------------------------
# Binary probit Z augmentation and working loading updates
# ----------------------------------------------------------------------------

rtruncnorm_binary_vec <- function(mean, sd, lower, upper) {
  # Inverse-CDF sampling from N(mean, sd^2) truncated to [lower, upper].
  mean <- as.numeric(mean)
  sd <- as.numeric(sd)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)

  if (length(sd) == 1L) sd <- rep(sd, length(mean))
  if (length(lower) == 1L) lower <- rep(lower, length(mean))
  if (length(upper) == 1L) upper <- rep(upper, length(mean))

  # Transform a uniform draw on the truncated CDF interval back to the normal.
  a <- pnorm((lower - mean) / sd)
  b <- pnorm((upper - mean) / sd)
  width <- pmax(b - a, .Machine$double.eps)
  u <- a + runif(length(mean)) * width
  mean + sd * qnorm(pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps))
}

initialize_binary_intercepts <- function(X, offset = 0.5) {
  X <- as.matrix(X)
  n_j <- colSums(!is.na(X))
  p_hat <- (colSums(X, na.rm = TRUE) + offset) / (n_j + 2 * offset)
  as.numeric(qnorm(pmin(pmax(p_hat, 1e-6), 1 - 1e-6)))
}

initialize_binary_Z <- function(X, seed = NULL, alpha = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  Z <- matrix(NA_real_, nrow(X), ncol(X))
  if (is.null(alpha)) alpha <- rep(0, ncol(X))
  alpha <- as.numeric(alpha)
  if (length(alpha) != ncol(X) || any(!is.finite(alpha))) alpha <- rep(0, ncol(X))

  # X = 1 means Z is positive; X = 0 means Z is negative.
  mean_mat <- matrix(rep(alpha, each = nrow(X)), nrow(X), ncol(X))
  Z[X == 1] <- rtruncnorm_binary_vec(mean_mat[X == 1], 1, 0, Inf)
  Z[X == 0] <- rtruncnorm_binary_vec(mean_mat[X == 0], 1, -Inf, 0)
  dimnames(Z) <- dimnames(X)
  Z
}

truncnorm_binary_moments_vec <- function(mean, sd, lower, upper) {
  # Closed-form first two moments of a truncated Gaussian.  The pretraining
  # algorithm can use E[Z | X, F, Lambda] instead of a random draw.
  a <- (lower - mean) / sd
  b <- (upper - mean) / sd
  Z <- pmax(pnorm(b) - pnorm(a), .Machine$double.eps)
  pa <- dnorm(a)
  pb <- dnorm(b)
  Ez_std <- (pa - pb) / Z
  Vz_std <- 1 + (a * pa - b * pb) / Z - Ez_std^2
  list(mean = mean + sd * Ez_std, var = pmax(sd^2 * Vz_std, 1e-8))
}

sample_binary_Z_given_model <- function(X, F_hat, Lambda, Psi, alpha = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  if (is.null(alpha)) alpha <- rep(0, ncol(X))
  alpha <- as.numeric(alpha)
  psi <- if (is.matrix(Psi)) diag(Psi) else as.numeric(Psi)
  psi <- pmax(psi, 1e-8)

  # Conditional mean of every latent Z_ij under the working factor layer.
  mean_mat <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  Z <- matrix(NA_real_, nrow(X), ncol(X))

  for (j in seq_len(ncol(X))) {
    # Each item column has its own probit residual variance psi_j.
    lower <- ifelse(X[, j] == 1, 0, -Inf)
    upper <- ifelse(X[, j] == 1, Inf, 0)
    Z[, j] <- rtruncnorm_binary_vec(mean_mat[, j], sqrt(psi[j]), lower, upper)
  }

  dimnames(Z) <- dimnames(X)
  Z
}

expected_binary_Z_given_model <- function(X, F_hat, Lambda, Psi, alpha = NULL) {
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  if (is.null(alpha)) alpha <- rep(0, ncol(X))
  alpha <- as.numeric(alpha)
  psi <- if (is.matrix(Psi)) diag(Psi) else as.numeric(Psi)
  psi <- pmax(psi, 1e-8)

  # Same conditional distribution as the sampling function, but replace each
  # latent Z_ij by its conditional expectation.
  mean_mat <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  Z_mean <- matrix(NA_real_, nrow(X), ncol(X))

  for (j in seq_len(ncol(X))) {
    lower <- ifelse(X[, j] == 1, 0, -Inf)
    upper <- ifelse(X[, j] == 1, Inf, 0)
    Z_mean[, j] <- truncnorm_binary_moments_vec(
      mean = mean_mat[, j],
      sd = sqrt(psi[j]),
      lower = lower,
      upper = upper
    )$mean
  }

  dimnames(Z_mean) <- dimnames(X)
  Z_mean
}

update_working_loadings_no_intercept <- function(
    Z,
    F_hat,
    loading_penalty = 0,
    min_residual_var = 1e-4) {
  Z <- as.matrix(Z)
  F_hat <- as.matrix(F_hat)
  n <- nrow(Z)
  p <- ncol(Z)

  # F_hat is scaled like sqrt(n) times an orthonormal SVD basis, and the
  # rotation is orthogonal.  Therefore crossprod(F_hat) / n is approximately I.
  # With that normalization, soft-thresholding the least-squares coefficients is
  # the coordinate-wise lasso solution for each item regression Z_j ~ F.
  Lambda_ls <- crossprod(Z, F_hat) / n
  Lambda <- soft_threshold(Lambda_ls, loading_penalty)
  fitted <- F_hat %*% t(Lambda)
  residual <- Z - fitted

  # The residual variances are diagnostic only in this script.  The probit
  # augmentation below always uses Psi = I for scale identification.
  psi <- pmax(colMeans(residual^2), min_residual_var)

  rownames(Lambda) <- colnames(Z)
  colnames(Lambda) <- paste0("factor_", seq_len(ncol(F_hat)))

  list(
    Lambda = Lambda,
    Lambda_ls = Lambda_ls,
    loading_update = if (loading_penalty > 0) "lasso_soft_threshold" else "least_squares",
    loading_penalty = loading_penalty,
    Psi = diag(psi, p),
    fitted = fitted,
    residual = residual
  )
}

update_working_loadings_with_intercept <- function(
    Z,
    F_hat,
    loading_penalty = 0,
    min_residual_var = 1e-4) {
  Z <- as.matrix(Z)
  F_hat <- as.matrix(F_hat)
  p <- ncol(Z)
  H <- ncol(F_hat)

  alpha <- colMeans(Z)
  F_center <- sweep(F_hat, 2L, colMeans(F_hat), "-")
  Z_center <- sweep(Z, 2L, alpha, "-")
  XtX <- crossprod(F_center) + diag(1e-8, H)
  Lambda_ls <- t(solve(XtX, crossprod(F_center, Z_center)))
  Lambda <- soft_threshold(Lambda_ls, loading_penalty)
  alpha <- as.numeric(colMeans(Z - F_hat %*% t(Lambda)))
  fitted <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  residual <- Z - fitted
  psi <- pmax(colMeans(residual^2), min_residual_var)

  rownames(Lambda) <- colnames(Z)
  colnames(Lambda) <- paste0("factor_", seq_len(H))

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

binary_probit_observed_loglik <- function(X, F_hat, Lambda, alpha = NULL) {
  # Observed binary probit log likelihood under
  # P(X_ij = 1 | f_i) = Phi(lambda_j' f_i).  This is used as a coherent
  # pretraining diagnostic/stopping target because it scores the actual binary
  # data rather than the augmented Z reconstruction.
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  if (is.null(alpha)) alpha <- rep(0, ncol(X))
  alpha <- as.numeric(alpha)

  eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)

  sum(X * log(p1) + (1 - X) * log(p0))
}

prune_factors_by_probit_drop <- function(
    X,
    F_hat,
    Lambda,
    mixture_fits = NULL,
    class_map = NULL,
    responsibilities = NULL,
    loading_threshold = 0.10,
    min_active_loadings = 4L,
    min_loading_energy_share = 0,
    likelihood_penalty_multiplier = 1,
    min_keep = 1L) {
  # Factor pruning after a few Z-augmentation updates.  The goal is to avoid
  # choosing H from the crude initial Z.  We first run an intentionally generous
  # factor model, then keep factors whose estimated loadings are active and
  # whose removal causes a meaningful drop in the observed binary probit
  # likelihood.
  X <- as.matrix(X)
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  H <- ncol(F_hat)
  n <- nrow(X)

  full_ll <- binary_probit_observed_loglik(X, F_hat, Lambda)
  active_count <- colSums(abs(Lambda) > loading_threshold)
  loading_energy <- colSums(Lambda^2)
  loading_energy_share <- loading_energy / pmax(sum(loading_energy), 1e-12)
  ll_drop <- numeric(H)

  for (h in seq_len(H)) {
    keep_h <- setdiff(seq_len(H), h)
    minus_ll <- if (length(keep_h) == 0L) {
      binary_probit_observed_loglik(
        X,
        matrix(0, nrow = nrow(F_hat), ncol = 1L),
        matrix(0, nrow = ncol(X), ncol = 1L)
      )
    } else {
      binary_probit_observed_loglik(
        X,
        F_hat[, keep_h, drop = FALSE],
        Lambda[, keep_h, drop = FALSE]
      )
    }
    ll_drop[h] <- full_ll - minus_ll
  }

  # A factor should survive only if it has enough active loadings and improves
  # the data likelihood by more than a BIC-like active-support penalty.
  penalty <- likelihood_penalty_multiplier * active_count * log(n)
  keep <- active_count >= min_active_loadings &
    loading_energy_share >= min_loading_energy_share &
    (2 * ll_drop) > penalty

  # If the rule is too aggressive, retain the strongest factors by likelihood
  # evidence so the algorithm can continue.
  min_keep <- max(1L, min(as.integer(min_keep), H))
  if (sum(keep) < min_keep) {
    rank_score <- 2 * ll_drop - penalty
    keep[order(rank_score, loading_energy, decreasing = TRUE)[seq_len(min_keep)]] <- TRUE
  }

  kept <- which(keep)
  dropped <- which(!keep)

  out <- list(
    F_hat = F_hat[, kept, drop = FALSE],
    Lambda = Lambda[, kept, drop = FALSE],
    mixture_fits = if (!is.null(mixture_fits)) mixture_fits[kept] else NULL,
    class_map = if (!is.null(class_map)) class_map[, kept, drop = FALSE] else NULL,
    responsibilities = if (!is.null(responsibilities)) responsibilities[kept] else NULL,
    kept = kept,
    dropped = dropped,
    diagnostics = data.frame(
      factor = seq_len(H),
      active_loadings = active_count,
      loading_energy = loading_energy,
      loading_energy_share = loading_energy_share,
      probit_loglik_drop = ll_drop,
      penalty = penalty,
      keep = keep
    )
  )

  if (!is.null(out$class_map)) {
    colnames(out$class_map) <- paste0("factor_", seq_len(ncol(out$F_hat)))
  }
  colnames(out$Lambda) <- paste0("factor_", seq_len(ncol(out$F_hat)))
  out
}

svd_scores_from_Z <- function(Z, H, center_Z = FALSE) {
  Z <- as.matrix(Z)
  # The left singular vectors estimate the signal column space of latent Z.
  Z_work <- if (isTRUE(center_Z)) sweep(Z, 2L, colMeans(Z), "-") else Z
  dec <- svd(Z_work, nu = H, nv = H)
  S <- sqrt(nrow(Z_work)) * dec$u[, seq_len(H), drop = FALSE]

  list(
    S = S,
    singular_values = dec$d,
    U = dec$u,
    V = dec$v,
    center_Z = center_Z,
    column_center = if (isTRUE(center_Z)) colMeans(Z) else rep(0, ncol(Z))
  )
}

select_num_factors_bai_ng <- function(
    X,
    H_max = min(20L, nrow(X) - 1L, ncol(X)),
    H_min = 1L,
    criterion = c("ICp2", "ICp1", "ICp3")) {
  criterion <- match.arg(criterion)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  H_max <- min(as.integer(H_max), n - 1L, p)

  dec <- svd(X, nu = H_max, nv = H_max)
  values <- numeric(H_max)

  for (H in seq_len(H_max)) {
    # Reconstruct the rank-H approximation and score residual variance plus
    # the Bai-Ng complexity penalty.
    X_hat <- dec$u[, seq_len(H), drop = FALSE] %*%
      diag(dec$d[seq_len(H)], H, H) %*%
      t(dec$v[, seq_len(H), drop = FALSE])
    residual_var <- mean((X - X_hat)^2)
    penalty <- switch(
      criterion,
      ICp1 = (n + p) / (n * p) * log(n * p / (n + p)),
      ICp2 = (n + p) / (n * p) * log(min(n, p)),
      ICp3 = log(min(n, p)) / min(n, p)
    )
    values[H] <- log(residual_var) + H * penalty
  }

  H_grid <- seq_len(H_max)
  keep <- H_grid >= H_min
  # Select the rank with the smallest information criterion value.
  H_hat <- H_grid[keep][which.min(values[keep])]

  list(H_hat = H_hat, criterion = criterion, values = values)
}

binary_profile_id <- function(class_map) {
  apply(class_map, 1L, paste, collapse = "-")
}

summarize_binary_mixture_profiles <- function(fit) {
  do.call(rbind, lapply(seq_len(fit$H), function(h) {
    # Summarize the hard labels and fitted marginal mixture parameters for
    # each factor axis.
    cls <- fit$class_map[, h]
    tab <- tabulate(cls, nbins = fit$G_hat[h])
    data.frame(
      factor = h,
      group = seq_len(fit$G_hat[h]),
      n = tab,
      prop = tab / length(cls),
      mean = fit$mixture_fits[[h]]$mu,
      sd = sqrt(fit$mixture_fits[[h]]$var),
      weight = fit$mixture_fits[[h]]$pi
    )
  }))
}

# ----------------------------------------------------------------------------
# Main pretraining algorithm
# ----------------------------------------------------------------------------

fit_binary_probit_pretraining <- function(
    X,
    H = NULL,
    H_max = min(10L, nrow(as.matrix(X)) - 1L, ncol(as.matrix(X))),
    factor_criterion = "ICp2",
    H_selection_strategy = c("initial_bai_ng", "overfit_prune"),
    H_prune_after_iter = 4L,
    H_prune_loading_threshold = 0.10,
    H_prune_min_active_loadings = 4L,
    H_prune_min_loading_energy_share = 0,
    H_prune_likelihood_penalty_multiplier = 1,
    H_prune_min_keep = 1L,
    pretrain_G_selection = c("bic", "fixed"),
    G_max = 5L,
    G_fixed = NULL,
    n_aug_iter = 4L,
    z_update = c("sample", "expectation"),
    center_Z_for_svd = FALSE,
    n_random_starts = 1L,
    max_outer = 4L,
    n_mix_starts = 3L,
    mixture_max_iter = 20L,
    mixture_penalty_multiplier = 1,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 0.3,
    weight_prior_alpha = 1,
    estimate_intercept = TRUE,
    loading_penalty = 0,
    pretrain_objective = c(
      "full_data_loglik",
      "lpretrain",
      "probit_loglik",
      "mixture_loglik",
      "combined",
      "penalized_score",
      "negative_reconstruction_mse"
    ),
    pretrain_objective_tolerance = NULL,
    pretrain_objective_patience = 0L,
    min_aug_iter = 2L,
    grid_size = 21L,
    rotation_sweep = c("full", "disjoint", "multi_disjoint", "promising", "hybrid"),
    disjoint_pairing = c("random", "adjacent"),
    n_disjoint_rounds = 4L,
    promising_max_pairs = NULL,
    promising_fraction = 0.15,
    promising_min_score = 0,
    store_iteration_trace = FALSE,
    return_best_iteration = FALSE,
    parallel = FALSE,
    workers = NULL,
    seed = 20260715L,
    verbose = TRUE,
    ...) {
  pretrain_G_selection <- match.arg(pretrain_G_selection)
  H_selection_strategy <- match.arg(H_selection_strategy)
  z_update <- match.arg(z_update)
  pretrain_objective <- match.arg(pretrain_objective)
  mixture_update <- match.arg(mixture_update)
  rotation_sweep <- match.arg(rotation_sweep)
  disjoint_pairing <- match.arg(disjoint_pairing)
  n_disjoint_rounds <- max(1L, as.integer(n_disjoint_rounds))
  workers <- resolve_workers(workers)
  min_aug_iter <- max(1L, as.integer(min_aug_iter))
  pretrain_objective_patience <- max(0L, as.integer(pretrain_objective_patience))
  estimate_intercept <- isTRUE(estimate_intercept)

  # Validate the observed data and add default item names for readable output.
  X <- as.matrix(X)
  if (!all(X %in% c(0, 1))) stop("X must be a binary matrix with entries 0/1.")
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))

  alpha_hat <- if (estimate_intercept) initialize_binary_intercepts(X) else rep(0, ncol(X))

  if (z_update == "expectation") {
    # Deterministic initialization at the mean of a standard half-normal.
    Z <- initialize_binary_Z(X, seed = seed, alpha = alpha_hat)
    Z[X == 1] <- pmax(alpha_hat[col(Z)][X == 1], 0) + sqrt(2 / pi)
    Z[X == 0] <- pmin(alpha_hat[col(Z)][X == 0], 0) - sqrt(2 / pi)
  } else {
    # Stochastic initialization from the standard probit truncation.
    Z <- initialize_binary_Z(X, seed = seed, alpha = alpha_hat)
  }

  H_selection <- NULL
  H_pruning <- NULL
  H_pruned <- FALSE
  if (is.null(H)) {
    Z_for_selection <- if (isTRUE(center_Z_for_svd)) {
      sweep(Z, 2L, colMeans(Z), "-")
    } else {
      Z
    }
    H_selection_initial <- select_num_factors_bai_ng(
      Z_for_selection,
      H_max = H_max,
      H_min = 1L,
      criterion = factor_criterion
    )

    if (H_selection_strategy == "overfit_prune") {
      # Start deliberately large, then prune after a few Z updates.  This avoids
      # making an irreversible low-rank choice from the crude initial Z.
      H <- min(as.integer(H_max), nrow(X) - 1L, ncol(X))
      H_selection <- H_selection_initial
      H_selection$H_initial_bai_ng <- H_selection_initial$H_hat
      H_selection$H_start <- H
      H_selection$strategy <- H_selection_strategy
      if (verbose) {
        message(
          "Starting with H = ", H,
          " for overfit-prune pretraining; initial ",
          factor_criterion,
          " would choose H = ", H_selection_initial$H_hat,
          "."
        )
      }
    } else {
      H_selection <- H_selection_initial
      H_selection$strategy <- H_selection_strategy
      H <- H_selection$H_hat
      if (verbose) message("Selected H = ", H, " by ", factor_criterion, ".")
    }
  }

  H <- as.integer(H)
  if (pretrain_G_selection == "fixed") {
    G_fixed <- normalize_G_fixed(G_fixed, H)
  }
  history <- vector("list", n_aug_iter)
  iteration_trace <- if (isTRUE(store_iteration_trace)) vector("list", n_aug_iter) else NULL
  candidate_snapshots <- if (isTRUE(return_best_iteration)) vector("list", n_aug_iter) else NULL
  current <- NULL
  converged <- FALSE
  completed_iter <- 0L
  best_pretrain_objective <- -Inf
  best_pretrain_iteration <- NA_integer_
  no_improvement_count <- 0L

  for (iter in seq_len(n_aug_iter)) {
    iter_start <- Sys.time()
    if (verbose) message("Pretraining augmentation iteration ", iter, " of ", n_aug_iter, ".")

    # Step 1: estimate the factor subspace from the current latent Z.
    svd_start <- Sys.time()
    svd_out <- svd_scores_from_Z(Z, H = H, center_Z = center_Z_for_svd)
    svd_seconds <- as.numeric(difftime(Sys.time(), svd_start, units = "secs"))
    S <- svd_out$S

    # Step 2: rotate the SVD basis toward independent non-Gaussian mixture
    # coordinates.  The marginal mixture counts are either selected by BIC or
    # fixed at the supplied G_fixed values.
    rotation_start <- Sys.time()
    rotation_out <- estimate_mixture_ica_unknown_G(
      S = S,
      G_selection = pretrain_G_selection,
      G_max = G_max,
      G_fixed = G_fixed,
      n_random_starts = n_random_starts,
      max_outer = max_outer,
      n_mix_starts = n_mix_starts,
      mixture_max_iter = mixture_max_iter,
      mixture_penalty_multiplier = mixture_penalty_multiplier,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      grid_size = grid_size,
      rotation_sweep = rotation_sweep,
      disjoint_pairing = disjoint_pairing,
      n_disjoint_rounds = n_disjoint_rounds,
      promising_max_pairs = promising_max_pairs,
      promising_fraction = promising_fraction,
      promising_min_score = promising_min_score,
      seed = seed + 1000L + iter,
      parallel = parallel,
      workers = workers,
      verbose = verbose
    )
    rotation_seconds <- as.numeric(difftime(Sys.time(), rotation_start, units = "secs"))

    F_hat <- rotation_out$F_hat

    # Convert marginal mixture responsibilities into hard component labels.
    responsibilities <- lapply(seq_len(H), function(h) {
      mixture_responsibilities(F_hat[, h], rotation_out$fits[[h]])
    })
    class_map <- sapply(responsibilities, max.col, ties.method = "first")
    if (H == 1L) class_map <- matrix(class_map, ncol = 1L)
    colnames(class_map) <- paste0("factor_", seq_len(H))

    # Step 3: estimate working loadings from the current latent Z and factors.
    loading_start <- Sys.time()
    loading_out <- if (estimate_intercept) {
      update_working_loadings_with_intercept(
        Z = Z,
        F_hat = F_hat,
        loading_penalty = loading_penalty
      )
    } else {
      update_working_loadings_no_intercept(
        Z = Z,
        F_hat = F_hat,
        loading_penalty = loading_penalty
      )
    }
    alpha_hat <- if (!is.null(loading_out$alpha)) loading_out$alpha else rep(0, ncol(X))
    loading_seconds <- as.numeric(difftime(Sys.time(), loading_start, units = "secs"))

    H_before_prune <- H
    pruned_factors <- integer(0)
    prune_diagnostics <- NULL
    if (H_selection_strategy == "overfit_prune" &&
        !isTRUE(H_pruned) &&
        iter >= as.integer(H_prune_after_iter)) {
      prune_out <- prune_factors_by_probit_drop(
        X = X,
        F_hat = F_hat,
        Lambda = loading_out$Lambda,
        mixture_fits = rotation_out$fits,
        class_map = class_map,
        responsibilities = responsibilities,
        loading_threshold = H_prune_loading_threshold,
        min_active_loadings = H_prune_min_active_loadings,
        min_loading_energy_share = H_prune_min_loading_energy_share,
        likelihood_penalty_multiplier = H_prune_likelihood_penalty_multiplier,
        min_keep = H_prune_min_keep
      )

      F_hat <- prune_out$F_hat
      rotation_out$F_hat <- F_hat
      rotation_out$fits <- prune_out$mixture_fits
      rotation_out$G_hat <- vapply(rotation_out$fits, function(z) length(z$pi), integer(1))
      rotation_out$loglik <- mixture_loglik_total(F_hat, rotation_out$fits)
      rotation_out$penalized_score <- mixture_bic_score(
        F_hat,
        rotation_out$fits,
        mixture_penalty_multiplier
      )
      loading_out$Lambda <- prune_out$Lambda
      loading_out$Lambda_ls <- loading_out$Lambda_ls[, prune_out$kept, drop = FALSE]
      loading_out$fitted <- sweep(F_hat %*% t(loading_out$Lambda), 2L, alpha_hat, "+")
      loading_out$residual <- Z - loading_out$fitted
      responsibilities <- prune_out$responsibilities
      class_map <- prune_out$class_map
      H <- ncol(F_hat)
      if (pretrain_G_selection == "fixed") {
        G_fixed <- normalize_G_fixed(G_fixed[prune_out$kept], H)
      }
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
      if (verbose) {
        message(
          "Pruned H from ", H_before_prune,
          " to ", H,
          " at augmentation iteration ", iter,
          "; dropped axes: ",
          if (length(pruned_factors)) paste(pruned_factors, collapse = ",") else "none",
          "."
        )
      }
    }

    # Step 4: fix the probit residual scale to the identity.  The unconstrained
    # residual variance is saved only as a diagnostic and is never used to
    # sample or average the next Z.
    objective_start <- Sys.time()
    Psi_for_update <- diag(1, ncol(X))
    reconstruction_mse <- mean(loading_out$residual^2)
    probit_loglik <- binary_probit_observed_loglik(
      X = X,
      F_hat = F_hat,
      Lambda = loading_out$Lambda,
      alpha = alpha_hat
    )
    full_data_loglik <- probit_loglik + rotation_out$loglik
    pretrain_objective_value <- switch(
      pretrain_objective,
      full_data_loglik = full_data_loglik,
      lpretrain = full_data_loglik,
      probit_loglik = probit_loglik,
      combined = rotation_out$penalized_score - 0.5 * length(X) * log(pmax(reconstruction_mse, 1e-12)),
      penalized_score = rotation_out$penalized_score,
      mixture_loglik = rotation_out$loglik,
      negative_reconstruction_mse = -reconstruction_mse
    )
    previous_objective <- if (iter > 1L) history[[iter - 1L]]$pretrain_objective else NA_real_
    relative_objective_change <- if (iter > 1L && is.finite(previous_objective)) {
      abs(pretrain_objective_value - previous_objective) / (1 + abs(previous_objective))
    } else {
      NA_real_
    }
    previous_best_objective <- best_pretrain_objective
    objective_improvement_from_best <- if (is.finite(previous_best_objective)) {
      pretrain_objective_value - previous_best_objective
    } else {
      Inf
    }
    relative_best_objective_improvement <- if (is.finite(previous_best_objective)) {
      objective_improvement_from_best / (1 + abs(previous_best_objective))
    } else {
      Inf
    }
    meaningful_improvement <- is.finite(pretrain_objective_value) &&
      (
        !is.finite(previous_best_objective) ||
          (
            objective_improvement_from_best > 0 &&
              (
                is.null(pretrain_objective_tolerance) ||
                  !is.finite(pretrain_objective_tolerance) ||
                  relative_best_objective_improvement > pretrain_objective_tolerance
              )
          )
      )
    if (meaningful_improvement) {
      best_pretrain_objective <- pretrain_objective_value
      best_pretrain_iteration <- iter
      no_improvement_count <- 0L
    } else {
      no_improvement_count <- no_improvement_count + 1L
    }
    objective_seconds <- as.numeric(difftime(Sys.time(), objective_start, units = "secs"))

    history[[iter]] <- data.frame(
      iteration = iter,
      H = H,
      G_hat = paste(rotation_out$G_hat, collapse = ","),
      probit_loglik = probit_loglik,
      full_data_loglik = full_data_loglik,
      mixture_loglik = rotation_out$loglik,
      penalized_score = rotation_out$penalized_score,
      reconstruction_mse = reconstruction_mse,
      pretrain_objective_name = pretrain_objective,
      pretrain_objective = pretrain_objective_value,
      relative_objective_change = relative_objective_change,
      best_pretrain_objective = best_pretrain_objective,
      best_pretrain_iteration = best_pretrain_iteration,
      objective_improvement_from_best = objective_improvement_from_best,
      relative_best_objective_improvement = relative_best_objective_improvement,
      no_improvement_count = no_improvement_count,
      converged = FALSE,
      mean_psi = mean(diag(Psi_for_update)),
      H_before_prune = H_before_prune,
      H_after_prune = H,
      n_pruned_factors = length(pruned_factors),
      pretrain_G_selection = pretrain_G_selection,
      mixture_update = mixture_update,
      svd_seconds = svd_seconds,
      rotation_seconds = rotation_seconds,
      loading_seconds = loading_seconds,
      objective_seconds = objective_seconds,
      z_update_seconds = NA_real_,
      iteration_seconds = NA_real_
    )
    completed_iter <- iter

    # Store a full snapshot from this augmentation iteration.  The final
    # snapshot is returned after the loop.
    current <- list(
      model = "binary_probit_pretraining_independent_mixture_factor",
      X = X,
      S = S,
      svd_fit = svd_out,
      R = rotation_out$R,
      F_hat = F_hat,
      mixture_fits = rotation_out$fits,
      G_hat = rotation_out$G_hat,
      rotation_fit = rotation_out,
      responsibilities = responsibilities,
      class_map = class_map,
      profile_id = binary_profile_id(class_map),
      Lambda_hat = loading_out$Lambda,
      alpha_hat = alpha_hat,
      Lambda_ls = loading_out$Lambda_ls,
      loading_update = loading_out$loading_update,
      loading_penalty = loading_out$loading_penalty,
      Psi_hat = Psi_for_update,
      Psi_hat_unconstrained = loading_out$Psi,
      fitted = loading_out$fitted,
      residual = loading_out$residual
    )
    if (isTRUE(return_best_iteration)) {
      # Keep the full returned-object state for this iteration, so a sampled-Z
      # run can return the best-scoring iterate instead of whichever iterate
      # happened to occur last.
      candidate_snapshots[[iter]] <- current
    }

    if (isTRUE(store_iteration_trace)) {
      # Keep the iteration-level quantities needed for diagnostics and plots.
      # This is opt-in because storing every F_hat/class_map/Lambda snapshot can
      # be nontrivial for large n and p.
      iteration_trace[[iter]] <- list(
        iteration = iter,
        S = S,
        R = rotation_out$R,
        F_hat = F_hat,
        mixture_fits = rotation_out$fits,
        G_hat = rotation_out$G_hat,
        responsibilities = responsibilities,
        class_map = class_map,
        Lambda_hat = loading_out$Lambda,
        alpha_hat = alpha_hat,
        reconstruction_mse = reconstruction_mse,
        probit_loglik = probit_loglik,
        full_data_loglik = full_data_loglik,
        mixture_loglik = rotation_out$loglik,
        penalized_score = rotation_out$penalized_score,
        pretrain_objective = pretrain_objective_value,
        relative_objective_change = relative_objective_change,
        H_before_prune = H_before_prune,
        H_after_prune = H,
        pruned_factors = pruned_factors,
        prune_diagnostics = prune_diagnostics
      )
    }

    if (verbose) {
      message(
        "  G = [", paste(rotation_out$G_hat, collapse = ", "),
        "]; reconstruction MSE = ", round(reconstruction_mse, 4),
        "; pretraining objective = ", round(pretrain_objective_value, 4)
      )
    }

    previous_change_converged <- !is.null(pretrain_objective_tolerance) &&
      is.finite(pretrain_objective_tolerance) &&
      pretrain_objective_patience == 0L &&
      iter >= min_aug_iter &&
      is.finite(relative_objective_change) &&
      relative_objective_change <= pretrain_objective_tolerance
    patience_converged <- !is.null(pretrain_objective_tolerance) &&
      is.finite(pretrain_objective_tolerance) &&
      pretrain_objective_patience > 0L &&
      iter >= min_aug_iter &&
      no_improvement_count >= pretrain_objective_patience

    if (previous_change_converged || patience_converged) {
      converged <- TRUE
      history[[iter]]$converged <- TRUE
      history[[iter]]$iteration_seconds <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
      if (verbose) {
        if (patience_converged) {
          message(
            "Stopping pretraining: no best-score improvement larger than ",
            signif(pretrain_objective_tolerance, 4),
            " for ", no_improvement_count,
            " consecutive iterations."
          )
        } else {
          message(
            "Stopping pretraining: relative objective change ",
            signif(relative_objective_change, 4),
            " <= tolerance ",
            signif(pretrain_objective_tolerance, 4),
            "."
          )
        }
      }
      break
    }

    if (z_update == "sample") {
      # Step 5a: stochastic DA update, Z ~ p(Z | X, F, Lambda, Psi).
      z_update_start <- Sys.time()
      Z <- sample_binary_Z_given_model(
        X = X,
        F_hat = F_hat,
        Lambda = loading_out$Lambda,
        alpha = alpha_hat,
        Psi = Psi_for_update,
        seed = seed + 2000L + iter
      )
      history[[iter]]$z_update_seconds <- as.numeric(difftime(Sys.time(), z_update_start, units = "secs"))
    } else {
      # Step 5b: deterministic EM-style update, Z <- E[Z | X, F, Lambda, Psi].
      z_update_start <- Sys.time()
      Z <- expected_binary_Z_given_model(
        X = X,
        F_hat = F_hat,
        Lambda = loading_out$Lambda,
        alpha = alpha_hat,
        Psi = Psi_for_update
      )
      history[[iter]]$z_update_seconds <- as.numeric(difftime(Sys.time(), z_update_start, units = "secs"))
    }
    history[[iter]]$iteration_seconds <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
  }

  history_out <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  selected_pretraining_iteration <- completed_iter
  if (isTRUE(return_best_iteration)) {
    selectable <- rep(TRUE, nrow(history_out))
    if (H_selection_strategy == "overfit_prune" && !is.null(H_pruning)) {
      # Once H has been pruned, do not return a pre-pruning snapshot with the
      # intentionally over-complete factor dimension.
      selectable <- history_out$iteration >= H_pruning$iteration
    }
    best_row <- which(selectable)[which.max(history_out$pretrain_objective[selectable])]
    selected_pretraining_iteration <- history_out$iteration[best_row]
    current <- candidate_snapshots[[selected_pretraining_iteration]]
  }

  # Attach run-level metadata and the final latent Z used to stop the loop.
  current$Z_last_sampled <- Z
  current$H <- H
  current$H_selection <- H_selection
  current$H_selection_strategy <- H_selection_strategy
  current$H_pruning <- H_pruning
  current$history <- history_out
  current$iteration_trace <- iteration_trace[!vapply(iteration_trace, is.null, logical(1))]
  current$pretraining_converged <- converged
  current$pretraining_completed_iter <- completed_iter
  current$selected_pretraining_iteration <- selected_pretraining_iteration
  current$return_best_iteration <- isTRUE(return_best_iteration)
  current$pretrain_objective <- pretrain_objective
  current$pretrain_objective_tolerance <- pretrain_objective_tolerance
  current$pretrain_objective_patience <- pretrain_objective_patience
  current$min_aug_iter <- min_aug_iter
  current$z_update <- z_update
  current$mixture_max_iter <- mixture_max_iter
  current$fix_psi_identity <- TRUE
  current$estimate_intercept <- estimate_intercept
  current$pretrain_G_selection <- pretrain_G_selection
  current$G_fixed <- G_fixed
  current$mixture_update <- mixture_update
  current$rotation_sweep <- rotation_sweep
  current$disjoint_pairing <- disjoint_pairing
  current$n_disjoint_rounds <- n_disjoint_rounds
  current$promising_max_pairs <- promising_max_pairs
  current$promising_fraction <- promising_fraction
  current$promising_min_score <- promising_min_score
  current$mixture_prior <- list(
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha
  )
  current$parallel <- list(enabled = isTRUE(parallel), workers = workers)
  current$call <- match.call()
  current
}

summarize_binary_probit_pretraining <- function(fit) {
  data.frame(
    stage = "pretraining",
    H_hat = fit$H,
    G_hat = paste(fit$G_hat, collapse = ","),
    pretrain_G_selection = fit$pretrain_G_selection,
    z_update = fit$z_update,
    psi = "identity",
    final_reconstruction_mse = tail(fit$history$reconstruction_mse, 1L),
    final_probit_loglik = tail(fit$history$probit_loglik, 1L),
    final_mixture_loglik = tail(fit$history$mixture_loglik, 1L),
    final_full_data_loglik = tail(fit$history$full_data_loglik, 1L),
    stringsAsFactors = FALSE
  )
}
