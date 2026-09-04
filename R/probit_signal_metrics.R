#!/usr/bin/env Rscript

# Utilities for comparing first-stage probit signal estimators.
#
# The target is the centered signal
#
#   L0,c = M_n eta0,
#
# where eta0 = 1 alpha0' + F0 Lambda0'.  This removes the arbitrary constant
# split between item intercepts and low-rank signal before comparing methods.

canonicalize_probit_signal <- function(eta_hat) {
  eta_hat <- as.matrix(eta_hat)
  alpha_hat <- colMeans(eta_hat)
  L_hat <- sweep(eta_hat, 2L, alpha_hat, "-")
  list(alpha_hat = as.numeric(alpha_hat), L_hat = L_hat, eta_hat = eta_hat)
}

make_centered_signal_truth <- function(eta0) {
  eta0 <- as.matrix(eta0)
  alpha0_c <- colMeans(eta0)
  L0_c <- sweep(eta0, 2L, alpha0_c, "-")
  list(alpha0_c = as.numeric(alpha0_c), L0_c = L0_c, eta0 = eta0)
}

probit_signal_loglik <- function(X, eta, mask = NULL) {
  X <- as.matrix(X)
  eta <- as.matrix(eta)
  if (!identical(dim(X), dim(eta))) stop("X and eta must have the same dimensions.")
  observed <- !is.na(X)
  if (!is.null(mask)) observed <- observed & as.matrix(mask)
  if (!any(observed)) return(NA_real_)
  s <- 2 * X - 1
  sum(pnorm((s * eta)[observed], log.p = TRUE))
}

matrix_rank_estimate <- function(M, tol = NULL) {
  M <- as.matrix(M)
  d <- svd(M, nu = 0, nv = 0)$d
  if (!length(d)) return(0L)
  if (is.null(tol)) tol <- max(d) * 1e-8
  sum(d > tol)
}

principal_angle_errors <- function(L_hat, L0_c, H) {
  L_hat <- as.matrix(L_hat)
  L0_c <- as.matrix(L0_c)
  if (!identical(dim(L_hat), dim(L0_c))) stop("L_hat and L0_c must have the same dimensions.")
  H <- min(as.integer(H), min(dim(L_hat)), min(dim(L0_c)))
  if (H < 1L) {
    return(list(sinTheta_op = NA_real_, sinTheta_fro = NA_real_, cosines = numeric(0)))
  }

  U_hat <- svd(L_hat, nu = H, nv = 0)$u[, seq_len(H), drop = FALSE]
  U0 <- svd(L0_c, nu = H, nv = 0)$u[, seq_len(H), drop = FALSE]
  cosines <- svd(crossprod(U_hat, U0), nu = 0, nv = 0)$d
  cosines <- pmin(pmax(cosines, 0), 1)
  sins <- sqrt(pmax(1 - cosines^2, 0))
  list(
    sinTheta_op = max(sins),
    sinTheta_fro = sqrt(sum(sins^2)),
    cosines = cosines
  )
}

probit_signal_metrics <- function(fit, truth, X, H, validation_mask = NULL) {
  method <- if (!is.null(fit$method)) fit$method else NA_character_
  n <- nrow(X)
  p <- ncol(X)

  unavailable <- is.null(fit$eta_hat) ||
    length(fit$eta_hat) == 1L && is.na(fit$eta_hat[1L]) ||
    !is.matrix(fit$eta_hat)

  if (isTRUE(unavailable)) {
    return(data.frame(
      method = method,
      relative_frobenius_error = NA_real_,
      per_entry_mse = NA_real_,
      relative_operator_error = NA_real_,
      signal_correlation = NA_real_,
      intercept_rmse = NA_real_,
      sinTheta_op = NA_real_,
      sinTheta_fro = NA_real_,
      estimated_rank = NA_integer_,
      full_loglik = NA_real_,
      heldout_loglik = NA_real_,
      heldout_loglik_per_response = NA_real_,
      converged = isTRUE(fit$converged),
      iterations = if (!is.null(fit$iterations)) fit$iterations else NA_integer_,
      elapsed_sec = if (!is.null(fit$elapsed_sec)) fit$elapsed_sec else NA_real_
    ))
  }

  eta_hat <- as.matrix(fit$eta_hat)
  can <- canonicalize_probit_signal(eta_hat)
  L_hat <- can$L_hat
  alpha_hat <- can$alpha_hat
  L0_c <- truth$L0_c
  alpha0_c <- truth$alpha0_c
  err <- L_hat - L0_c
  sv0 <- svd(L0_c, nu = 0, nv = 0)$d
  sv_err <- svd(err, nu = 0, nv = 0)$d
  denom_fro <- sqrt(sum(L0_c^2))
  denom_op <- if (length(sv0) >= H) sv0[H] else tail(sv0, 1L)
  angles <- principal_angle_errors(L_hat, L0_c, H)
  heldout_ll <- probit_signal_loglik(X, eta_hat, validation_mask)
  heldout_n <- if (is.null(validation_mask)) n * p else sum(as.matrix(validation_mask) & !is.na(X))

  data.frame(
    method = method,
    relative_frobenius_error = sqrt(sum(err^2)) / max(denom_fro, .Machine$double.eps),
    per_entry_mse = sum(err^2) / (n * p),
    relative_operator_error = max(sv_err) / max(denom_op, .Machine$double.eps),
    signal_correlation = suppressWarnings(cor(as.vector(L_hat), as.vector(L0_c))),
    intercept_rmse = sqrt(mean((alpha_hat - alpha0_c)^2)),
    sinTheta_op = angles$sinTheta_op,
    sinTheta_fro = angles$sinTheta_fro,
    estimated_rank = matrix_rank_estimate(L_hat),
    full_loglik = probit_signal_loglik(X, eta_hat),
    heldout_loglik = heldout_ll,
    heldout_loglik_per_response = heldout_ll / max(heldout_n, 1L),
    converged = isTRUE(fit$converged),
    iterations = if (!is.null(fit$iterations)) fit$iterations else NA_integer_,
    elapsed_sec = if (!is.null(fit$elapsed_sec)) fit$elapsed_sec else NA_real_
  )
}

probit_signal_singular_values <- function(fit, truth, H, n_extra = 5L) {
  r <- min(min(dim(truth$L0_c)), H + n_extra)
  truth_sv <- svd(truth$L0_c, nu = 0, nv = 0)$d
  if (is.null(fit$L_hat) || !is.matrix(fit$L_hat)) {
    fit_sv <- rep(NA_real_, r)
  } else {
    fit_sv <- svd(fit$L_hat, nu = 0, nv = 0)$d
  }
  data.frame(
    method = if (!is.null(fit$method)) fit$method else NA_character_,
    singular_index = seq_len(r),
    truth_singular_value = truth_sv[seq_len(r)],
    fitted_singular_value = fit_sv[seq_len(r)]
  )
}

probit_signal_sanity_checks <- function(fit, X) {
  if (is.null(fit$eta_hat) || !is.matrix(fit$eta_hat)) {
    return(data.frame(
      method = if (!is.null(fit$method)) fit$method else NA_character_,
      centered_L_max_abs_colmean = NA_real_,
      eta_reconstruction_max_abs_error = NA_real_,
      finite_eta = FALSE,
      correct_dimensions = FALSE
    ))
  }
  eta <- as.matrix(fit$eta_hat)
  can <- canonicalize_probit_signal(eta)
  eta_rec <- sweep(can$L_hat, 2L, can$alpha_hat, "+")
  data.frame(
    method = if (!is.null(fit$method)) fit$method else NA_character_,
    centered_L_max_abs_colmean = max(abs(colMeans(can$L_hat))),
    eta_reconstruction_max_abs_error = max(abs(eta - eta_rec)),
    finite_eta = all(is.finite(eta)),
    correct_dimensions = identical(dim(eta), dim(X))
  )
}
