#!/usr/bin/env Rscript

# Hybrid probit EM-SVD signal estimator.
#
# This keeps the fast EM structure used by fit_probit_signal_em_svd(), but
# replaces the hard rank-H M-step with singular-value soft-thresholding:
#
#   W_c = U diag(d) V'
#   L   = U diag((d - tau)_+) V'.
#
# The shrinkage level tau can be selected by held-out probit log likelihood.
# A rank cap such as H+5 keeps the SVD inexpensive while allowing a few extra
# directions to enter before the final rank-H projection used downstream.

probit_ifa_truncated_mean_missing <- function(X, alpha, L) {
  X <- as.matrix(X)
  eta <- sweep(as.matrix(L), 2L, as.numeric(alpha), "+")
  observed <- !is.na(X)
  W <- eta
  if (any(observed)) {
    s <- 2 * X - 1
    W[observed] <- eta[observed] +
      s[observed] * probit_ifa_inverse_mills(s[observed] * eta[observed])
  }
  W
}

fill_binary_matrix_for_svd_start <- function(X, epsilon = 1e-4) {
  # The Zhang-Chen-Li SVD initializer is defined for complete binary matrices.
  # For held-out diagnostics, fill masked entries with their training item mean
  # so the initializer does not peek at validation responses.
  X <- as.matrix(X)
  observed <- !is.na(X)
  col_mean <- colMeans(X, na.rm = TRUE)
  global_mean <- mean(X[observed])
  if (!is.finite(global_mean)) global_mean <- 0.5
  col_mean[!is.finite(col_mean)] <- global_mean
  col_mean <- pmin(pmax(col_mean, epsilon), 1 - epsilon)

  X_fill <- X
  miss <- which(!observed, arr.ind = TRUE)
  if (nrow(miss)) X_fill[miss] <- col_mean[miss[, 2L]]
  X_fill
}

zcl_svd_probit_signal_start <- function(
    X,
    H,
    epsilon = 1e-4,
    threshold_multiplier = 1.01,
    include_intercept_rank = TRUE) {
  # SVD initializer inspired by Zhang, Chen, and Li's note on exploratory item
  # factor analysis by SVD:
  #   1. low-rank SVD approximation of the response/probability matrix,
  #   2. clip fitted probabilities,
  #   3. map through the probit inverse link,
  #   4. center by item means and retain rank H.
  X_fill <- fill_binary_matrix_for_svd_start(X, epsilon = epsilon)
  n <- nrow(X_fill)
  p <- ncol(X_fill)
  H <- min(as.integer(H), n, p)
  if (H < 1L) stop("H must be positive.")

  dec_prob <- svd(X_fill, nu = min(n, p), nv = min(n, p))
  d <- dec_prob$d
  spectral_rank <- sum(d >= threshold_multiplier * sqrt(n))
  base_rank <- if (isTRUE(include_intercept_rank)) H + 1L else H
  r_tilde <- min(max(base_rank, spectral_rank, 1L), length(d), n, p)

  U <- dec_prob$u[, seq_len(r_tilde), drop = FALSE]
  V <- dec_prob$v[, seq_len(r_tilde), drop = FALSE]
  P_hat <- U %*% diag(d[seq_len(r_tilde)], r_tilde, r_tilde) %*% t(V)
  P_hat <- pmin(pmax(P_hat, epsilon), 1 - epsilon)

  M <- qnorm(P_hat)
  alpha <- colMeans(M)
  Mc <- sweep(M, 2L, alpha, "-")
  dec_signal <- svd(Mc, nu = H, nv = H)
  L <- dec_signal$u[, seq_len(H), drop = FALSE] %*%
    diag(dec_signal$d[seq_len(H)], H, H) %*%
    t(dec_signal$v[, seq_len(H), drop = FALSE])
  L <- sweep(L, 2L, colMeans(L), "-")

  list(
    alpha = as.numeric(alpha),
    L = L,
    probability_rank = r_tilde,
    spectral_rank = spectral_rank,
    probability_svd_values = d,
    signal_svd = dec_signal,
    epsilon = epsilon,
    threshold_multiplier = threshold_multiplier
  )
}

soft_svd_centered_projection <- function(W, shrinkage = 0, rank_cap = NULL) {
  W <- as.matrix(W)
  alpha <- colMeans(W)
  Wc <- sweep(W, 2L, alpha, "-")
  max_rank <- min(dim(Wc))
  if (!is.null(rank_cap)) max_rank <- min(max_rank, as.integer(rank_cap))
  if (max_rank < 1L) {
    return(list(alpha = as.numeric(alpha), L = 0 * Wc, svd = NULL, rank = 0L))
  }

  dec <- svd(Wc, nu = max_rank, nv = max_rank)
  d <- dec$d[seq_len(max_rank)]
  d_soft <- pmax(d - shrinkage, 0)
  keep <- which(d_soft > 0)
  if (!length(keep)) {
    L <- matrix(0, nrow(Wc), ncol(Wc), dimnames = dimnames(Wc))
  } else {
    L <- dec$u[, keep, drop = FALSE] %*%
      diag(d_soft[keep], length(keep), length(keep)) %*%
      t(dec$v[, keep, drop = FALSE])
    dimnames(L) <- dimnames(Wc)
  }
  L <- sweep(L, 2L, colMeans(L), "-")
  list(alpha = as.numeric(alpha), L = L, svd = dec, rank = length(keep))
}

rank_project_signal <- function(L, H) {
  L <- sweep(as.matrix(L), 2L, colMeans(L), "-")
  H_eff <- min(as.integer(H), min(dim(L)))
  dec <- svd(L, nu = H_eff, nv = H_eff)
  if (H_eff < 1L) return(0 * L)
  out <- dec$u[, seq_len(H_eff), drop = FALSE] %*%
    diag(dec$d[seq_len(H_eff)], H_eff, H_eff) %*%
    t(dec$v[, seq_len(H_eff), drop = FALSE])
  sweep(out, 2L, colMeans(out), "-")
}

fit_probit_signal_em_svd_soft_one <- function(
    X,
    H,
    shrinkage,
    rank_cap = NULL,
    max_iter = 50L,
    tol_loglik = 1e-5,
    tol_L = 1e-4,
    alpha_init = NULL,
    L_init = NULL,
    verbose = FALSE) {
  if (!exists("probit_ifa_inverse_mills", mode = "function", inherits = TRUE)) {
    stop("probit_ifa_inverse_mills() must be sourced before fit_probit_signal_em_svd_soft_one().")
  }
  if (!exists("initialize_binary_intercepts", mode = "function", inherits = TRUE)) {
    stop("initialize_binary_intercepts() must be sourced before fit_probit_signal_em_svd_soft_one().")
  }

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  H <- as.integer(H)
  if (H < 1L || H > min(n, p)) stop("H must be between 1 and min(n, p).")

  alpha <- if (is.null(alpha_init)) initialize_binary_intercepts(X) else as.numeric(alpha_init)
  if (length(alpha) != p || any(!is.finite(alpha))) alpha <- initialize_binary_intercepts(X)
  L <- if (is.null(L_init)) matrix(0, n, p) else as.matrix(L_init)
  if (!identical(dim(L), dim(X))) stop("L_init must have the same dimensions as X.")
  L <- sweep(L, 2L, colMeans(L), "-")

  old_loglik <- probit_signal_loglik(X, sweep(L, 2L, alpha, "+"))
  if (!is.finite(old_loglik)) old_loglik <- -Inf
  history <- vector("list", max_iter)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    iter_start <- Sys.time()
    e_step_start <- Sys.time()
    W <- probit_ifa_truncated_mean_missing(X, alpha, L)
    e_step_seconds <- as.numeric(difftime(Sys.time(), e_step_start, units = "secs"))
    projection_start <- Sys.time()
    projection <- soft_svd_centered_projection(W, shrinkage = shrinkage, rank_cap = rank_cap)
    alpha_new <- projection$alpha
    L_new <- projection$L
    projection_seconds <- as.numeric(difftime(Sys.time(), projection_start, units = "secs"))
    objective_start <- Sys.time()
    eta_new <- sweep(L_new, 2L, alpha_new, "+")
    loglik <- probit_signal_loglik(X, eta_new)
    objective_seconds <- as.numeric(difftime(Sys.time(), objective_start, units = "secs"))

    rel_loglik <- abs(loglik - old_loglik) / (1 + abs(old_loglik))
    if (!is.finite(rel_loglik)) rel_loglik <- Inf
    rel_L <- sqrt(sum((L_new - L)^2)) / (1 + sqrt(sum(L^2)))
    history[[iter]] <- data.frame(
      iteration = iter,
      probit_loglik = loglik,
      relative_loglik_change = rel_loglik,
      relative_L_change = rel_L,
      active_rank = projection$rank,
      shrinkage = shrinkage,
      e_step_seconds = e_step_seconds,
      projection_seconds = projection_seconds,
      objective_seconds = objective_seconds,
      iteration_seconds = as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
    )

    if (isTRUE(verbose)) {
      message(
        "  EM-SVD soft iter ", iter,
        ": probit loglik=", round(loglik, 3),
        ", rank=", projection$rank,
        ", tau=", signif(shrinkage, 4),
        ", rel ll=", signif(rel_loglik, 3),
        ", rel L=", signif(rel_L, 3)
      )
    }

    alpha <- alpha_new
    L <- L_new
    old_loglik <- loglik

    if (rel_loglik <= tol_loglik && rel_L <= tol_L) {
      converged <- TRUE
      break
    }
  }

  history <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  eta_hat <- sweep(L, 2L, alpha, "+")
  can <- canonicalize_probit_signal(eta_hat)
  list(
    method = "em_svd_soft",
    eta_hat = can$eta_hat,
    alpha_hat = can$alpha_hat,
    L_hat = can$L_hat,
    loglik = probit_signal_loglik(X, can$eta_hat),
    converged = converged,
    iterations = nrow(history),
    elapsed_sec = NA_real_,
    rank = matrix_rank_estimate(can$L_hat),
    tuning = list(
      H = H,
      shrinkage = shrinkage,
      rank_cap = rank_cap,
      max_iter = max_iter,
      tol_loglik = tol_loglik,
      tol_L = tol_L
    ),
    history = history
  )
}

initial_soft_svd_shrinkage_reference <- function(
    X,
    H,
    reference = c("dHplus1", "dH", "sqrt_np")) {
  reference <- match.arg(reference)
  X <- as.matrix(X)
  alpha <- initialize_binary_intercepts(X)
  W0 <- probit_ifa_truncated_mean_missing(
    X = X,
    alpha = alpha,
    L = matrix(0, nrow(X), ncol(X))
  )
  Wc <- sweep(W0, 2L, colMeans(W0), "-")
  d <- svd(Wc, nu = 0, nv = 0)$d
  if (reference == "sqrt_np") return(sqrt(nrow(X) * ncol(X)))
  idx <- if (reference == "dH") as.integer(H) else as.integer(H) + 1L
  if (idx <= length(d)) return(d[idx])
  tail(d, 1L)
}

fit_probit_signal_em_svd_soft_path <- function(
    X,
    H,
    shrinkage_grid = NULL,
    shrinkage_ratio_grid = c(0, 0.1, 0.2, 0.35, 0.5, 0.75, 1),
    shrinkage_reference = c("dHplus1", "dH", "sqrt_np"),
    validation_mask = NULL,
    rank_cap = NULL,
    alpha_init = NULL,
    L_init = NULL,
    max_iter = 50L,
    tol_loglik = 1e-5,
    tol_L = 1e-4,
    return_all_fits = TRUE,
    verbose = FALSE) {
  shrinkage_reference <- match.arg(shrinkage_reference)
  X <- as.matrix(X)
  X_train <- X
  if (!is.null(validation_mask)) X_train[as.matrix(validation_mask)] <- NA

  ref <- initial_soft_svd_shrinkage_reference(X_train, H = H, reference = shrinkage_reference)
  if (is.null(shrinkage_grid)) shrinkage_grid <- ref * shrinkage_ratio_grid
  shrinkage_grid <- sort(unique(as.numeric(shrinkage_grid)))
  shrinkage_grid <- shrinkage_grid[is.finite(shrinkage_grid) & shrinkage_grid >= 0]
  if (!length(shrinkage_grid)) stop("shrinkage_grid is empty after filtering.")

  path <- vector("list", length(shrinkage_grid))
  fits <- vector("list", length(shrinkage_grid))
  total_start <- Sys.time()

  for (ell in seq_along(shrinkage_grid)) {
    tic <- Sys.time()
    fit <- fit_probit_signal_em_svd_soft_one(
      X = X_train,
      H = H,
      shrinkage = shrinkage_grid[ell],
      rank_cap = rank_cap,
      max_iter = max_iter,
      tol_loglik = tol_loglik,
      tol_L = tol_L,
      alpha_init = alpha_init,
      L_init = L_init,
      verbose = verbose
    )
    fit$elapsed_sec <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
    fit$method <- "em_svd_soft"
    fit$tuning$shrinkage_reference <- shrinkage_reference
    fit$tuning$shrinkage_reference_value <- ref
    fit$tuning$shrinkage_ratio <- shrinkage_grid[ell] / max(ref, .Machine$double.eps)

    valid_ll <- probit_signal_loglik(X, fit$eta_hat, validation_mask)
    valid_n <- if (is.null(validation_mask)) nrow(X) * ncol(X) else sum(validation_mask & !is.na(X))
    path[[ell]] <- data.frame(
      shrinkage_index = ell,
      shrinkage = shrinkage_grid[ell],
      shrinkage_reference = ref,
      shrinkage_ratio = shrinkage_grid[ell] / max(ref, .Machine$double.eps),
      selected_by_validation = FALSE,
      validation_loglik = valid_ll,
      validation_loglik_per_response = valid_ll / max(valid_n, 1L),
      full_loglik = probit_signal_loglik(X, fit$eta_hat),
      estimated_rank = fit$rank,
      converged = fit$converged,
      iterations = fit$iterations,
      elapsed_sec = fit$elapsed_sec
    )
    fits[[ell]] <- if (isTRUE(return_all_fits)) fit else NULL
  }

  path <- do.call(rbind, path)
  best_idx <- which.max(path$validation_loglik_per_response)
  if (!length(best_idx) || !is.finite(path$validation_loglik_per_response[best_idx])) {
    best_idx <- which.max(path$full_loglik)
  }
  path$selected_by_validation[best_idx] <- TRUE
  best <- if (isTRUE(return_all_fits)) fits[[best_idx]] else NULL
  if (is.null(best)) stop("return_all_fits must be TRUE to return the selected fitted signal.")
  best$method <- "em_svd_soft_cv"
  best$elapsed_sec <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  best$tuning$selected_shrinkage_index <- best_idx
  best$tuning$selection_rule <- "heldout_probit_loglik_per_response"

  L_rankH <- rank_project_signal(best$L_hat, H = H)
  eta_rankH <- sweep(L_rankH, 2L, best$alpha_hat, "+")
  rankH <- list(
    method = "em_svd_soft_cv_rankH",
    eta_hat = eta_rankH,
    alpha_hat = best$alpha_hat,
    L_hat = L_rankH,
    loglik = probit_signal_loglik(X, eta_rankH),
    converged = best$converged,
    iterations = best$iterations,
    elapsed_sec = best$elapsed_sec,
    tuning = c(best$tuning, list(debias = "rank_H_svd_projection")),
    history = best$history
  )

  list(best = best, rankH = rankH, path = path, fits = fits)
}
