#!/usr/bin/env Rscript

# Final canonical normalization for fitted probit factor-mixture models.
#
# These transformations are for reporting and cross-run comparison after
# convergence.  They leave the probit linear predictor
#   alpha_j + f_i' lambda_j
# unchanged up to floating-point error.

mixture_marginal_location_scale <- function(mixture_fits, min_scale = 1e-8) {
  H <- length(mixture_fits)
  mean <- numeric(H)
  scale <- numeric(H)

  for (h in seq_len(H)) {
    fit <- mixture_fits[[h]]
    pi_h <- fit$pi / sum(fit$pi)
    mu_h <- fit$mu
    var_h <- fit$var
    mean[h] <- sum(pi_h * mu_h)
    second <- sum(pi_h * (var_h + mu_h^2))
    var <- second - mean[h]^2
    scale[h] <- sqrt(max(var, min_scale^2))
  }

  list(mean = mean, scale = scale)
}

canonical_normalize_factor_parameters <- function(
    F_hat,
    Lambda,
    alpha,
    mixture_fits,
    min_scale = 1e-8,
    sign_rule = c("largest_loading_positive", "none"),
    order_rule = c("loading_l2", "none")) {
  sign_rule <- match.arg(sign_rule)
  order_rule <- match.arg(order_rule)
  F_hat <- as.matrix(F_hat)
  Lambda <- as.matrix(Lambda)
  alpha <- as.numeric(alpha)
  H <- ncol(F_hat)
  if (ncol(Lambda) != H) stop("F_hat and Lambda must have the same number of factors.")
  if (length(alpha) != nrow(Lambda)) stop("alpha must have length nrow(Lambda).")
  if (length(mixture_fits) != H) stop("mixture_fits must have length H.")

  eta_before <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")

  loc_scale <- mixture_marginal_location_scale(mixture_fits, min_scale = min_scale)
  location <- loc_scale$mean
  scale <- loc_scale$scale

  alpha_new <- alpha + as.numeric(Lambda %*% location)
  F_new <- sweep(sweep(F_hat, 2L, location, "-"), 2L, scale, "/")
  Lambda_new <- sweep(Lambda, 2L, scale, "*")
  mixture_new <- mixture_fits
  for (h in seq_len(H)) {
    mixture_new[[h]]$mu <- (mixture_new[[h]]$mu - location[h]) / scale[h]
    mixture_new[[h]]$var <- pmax(mixture_new[[h]]$var / scale[h]^2, min_scale^2)
    ord <- order(mixture_new[[h]]$mu)
    mixture_new[[h]]$pi <- mixture_new[[h]]$pi[ord]
    mixture_new[[h]]$mu <- mixture_new[[h]]$mu[ord]
    mixture_new[[h]]$var <- mixture_new[[h]]$var[ord]
  }

  signs <- rep(1, H)
  if (sign_rule == "largest_loading_positive") {
    for (h in seq_len(H)) {
      j <- which.max(abs(Lambda_new[, h]))
      if (length(j) == 1L && is.finite(Lambda_new[j, h]) && Lambda_new[j, h] < 0) {
        signs[h] <- -1
      }
    }
    F_new <- sweep(F_new, 2L, signs, "*")
    Lambda_new <- sweep(Lambda_new, 2L, signs, "*")
    for (h in seq_len(H)) {
      if (signs[h] < 0) {
        mixture_new[[h]]$mu <- -mixture_new[[h]]$mu
        ord <- order(mixture_new[[h]]$mu)
        mixture_new[[h]]$pi <- mixture_new[[h]]$pi[ord]
        mixture_new[[h]]$mu <- mixture_new[[h]]$mu[ord]
        mixture_new[[h]]$var <- mixture_new[[h]]$var[ord]
      }
    }
  }

  permutation <- seq_len(H)
  if (order_rule == "loading_l2") {
    permutation <- order(sqrt(colSums(Lambda_new^2)), decreasing = TRUE)
    F_new <- F_new[, permutation, drop = FALSE]
    Lambda_new <- Lambda_new[, permutation, drop = FALSE]
    mixture_new <- mixture_new[permutation]
    signs <- signs[permutation]
    location <- location[permutation]
    scale <- scale[permutation]
  }

  colnames(F_new) <- paste0("factor_", seq_len(H))
  colnames(Lambda_new) <- paste0("factor_", seq_len(H))

  eta_after <- sweep(F_new %*% t(Lambda_new), 2L, alpha_new, "+")
  list(
    F_hat = F_new,
    Lambda = Lambda_new,
    alpha = alpha_new,
    mixture_fits = mixture_new,
    location_before = location,
    scale_before = scale,
    signs = signs,
    permutation = permutation,
    max_abs_eta_difference = max(abs(eta_before - eta_after)),
    rms_eta_difference = sqrt(mean((eta_before - eta_after)^2))
  )
}

canonical_normalize_refined_fit <- function(
    fit,
    min_scale = 1e-8,
    sign_rule = c("largest_loading_positive", "none"),
    order_rule = c("loading_l2", "none")) {
  norm <- canonical_normalize_factor_parameters(
    F_hat = fit$F_hat,
    Lambda = fit$Lambda_hat,
    alpha = fit$alpha_hat,
    mixture_fits = fit$mixture_fits,
    min_scale = min_scale,
    sign_rule = sign_rule,
    order_rule = order_rule
  )

  fit$F_hat <- norm$F_hat
  fit$Lambda_hat <- norm$Lambda
  fit$alpha_hat <- norm$alpha
  fit$mixture_fits <- norm$mixture_fits
  fit$canonical_normalization <- norm[
    c(
      "location_before",
      "scale_before",
      "signs",
      "permutation",
      "max_abs_eta_difference",
      "rms_eta_difference"
    )
  ]
  fit
}
