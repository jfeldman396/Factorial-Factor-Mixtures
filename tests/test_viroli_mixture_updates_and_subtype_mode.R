#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this test with Rscript.")
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."))

source(file.path(repo_root, "R", "binary_probit_pretraining.R"))
source(file.path(repo_root, "R", "canonical_factor_normalization.R"))
source(file.path(repo_root, "R", "viroli_probit_independent_gibbs.R"))
source(file.path(repo_root, "R", "component_profile_recovery.R"))

# Check the normal-inverse-gamma and Dirichlet full conditionals against their
# analytic posterior means, including an empty component.
set.seed(20260924)
F <- matrix(c(-1.4, -0.8, -0.2, 0.5, 1.1, 1.7), ncol = 1L)
C <- matrix(c(1L, 1L, 1L, 2L, 2L, 2L), ncol = 1L)
G <- 3L
alpha_dirichlet <- 1.7
mu0 <- 0.15
kappa0 <- 0.6
a0 <- 4.2
b0 <- 1.3
n_draw <- 30000L

posterior_moments <- lapply(seq_len(G), function(g) {
  x <- F[C[, 1L] == g, 1L]
  n_g <- length(x)
  xbar <- if (n_g) mean(x) else 0
  ss <- if (n_g) sum((x - xbar)^2) else 0
  kappa_n <- kappa0 + n_g
  a_n <- a0 + n_g / 2
  b_n <- b0 + 0.5 * ss +
    (kappa0 * n_g * (xbar - mu0)^2) / (2 * kappa_n)
  c(
    mu = (kappa0 * mu0 + n_g * xbar) / kappa_n,
    sig2 = b_n / (a_n - 1)
  )
})
posterior_moments <- do.call(rbind, posterior_moments)
expected_pi <- (tabulate(C[, 1L], nbins = G) + alpha_dirichlet) /
  (nrow(F) + G * alpha_dirichlet)

draws <- replicate(n_draw, {
  draw <- sample_viroli_mixture_parameters(
    F = F,
    C = C,
    G = G,
    alpha_dirichlet = alpha_dirichlet,
    mu0 = mu0,
    kappa0 = kappa0,
    a0 = a0,
    b0 = b0,
    min_var = 1e-12
  )
  c(draw$pi[1L, ], draw$mu[1L, ], draw$sig2[1L, ])
})
empirical_pi <- rowMeans(draws[seq_len(G), , drop = FALSE])
empirical_mu <- rowMeans(draws[G + seq_len(G), , drop = FALSE])
empirical_sig2 <- rowMeans(draws[2L * G + seq_len(G), , drop = FALSE])

stopifnot(
  max(abs(empirical_pi - expected_pi)) < 0.006,
  max(abs(empirical_mu - posterior_moments[, "mu"])) < 0.025,
  max(abs(empirical_sig2 - posterior_moments[, "sig2"])) < 0.015
)

# The joint posterior profile mode need not equal the vector of marginal modes.
# Here (1, 2) occurs most often, while the two marginal modes combine to (2, 2).
profile_trace <- matrix(c(3L, 3L, 3L, 2L, 2L, 4L, 4L), nrow = 1L)
profile_draws <- viroli_component_profile_from_id(as.integer(profile_trace), c(2L, 2L))
component_counts <- lapply(seq_len(2L), function(h) {
  matrix(tabulate(profile_draws[, h], nbins = 2L), nrow = 1L)
})
posterior <- viroli_finalize_component_posterior(
  component_counts = component_counts,
  profile_trace = profile_trace,
  G = c(2L, 2L),
  n_keep = ncol(profile_trace)
)
stopifnot(
  identical(unname(posterior$marginal_mode), matrix(c(2L, 2L), nrow = 1L)),
  identical(unname(posterior$profile_mode), matrix(c(1L, 2L), nrow = 1L)),
  posterior$profile_mode_id == 3L,
  isTRUE(all.equal(posterior$profile_mode_probability, 3 / 7))
)

# Hard subtype accuracy must use the posterior modal profile, while probability
# scores continue to use the marginal posterior allocation probabilities.
truth <- rbind(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L))
soft <- lapply(seq_len(2L), function(h) {
  out <- matrix(0.4, nrow(truth), 2L)
  out[cbind(seq_len(nrow(truth)), truth[, h])] <- 0.6
  out
})
recovery <- component_profile_recovery(
  true_component = truth,
  estimated_responsibilities = soft,
  estimated_component = truth
)
stopifnot(
  recovery$component_profile_hamming_accuracy == 1,
  recovery$component_profile_exact_accuracy == 1,
  isTRUE(all.equal(recovery$mean_true_component_probability, 0.6)),
  recovery$component_brier_score > 0
)

# End-to-end smoke test: retained aligned Gibbs allocations must produce valid
# marginal probabilities and one posterior modal joint profile per observation.
set.seed(19)
n <- 24L
p <- 18L
H <- 2L
G_smoke <- c(2L, 2L)
F_smoke <- cbind(rnorm(n, -0.7, 0.5), rnorm(n, 0.7, 0.5))
Lambda_smoke <- matrix(rnorm(p * H, 0, 0.7), p, H)
X_smoke <- 1L * (
  matrix(rnorm(n * p), n, p) + F_smoke %*% t(Lambda_smoke) > 0
)
fit <- fit_viroli_probit_independent_gibbs(
  X = X_smoke,
  H = H,
  G = G_smoke,
  n_iter = 30L,
  burn = 10L,
  thin = 2L,
  compute_parameter_ess = FALSE,
  parallel = FALSE,
  seed = 23L,
  verbose = FALSE
)
stopifnot(
  fit$n_keep == 10L,
  identical(dim(fit$component_profile_mode), c(n, H)),
  length(fit$component_probabilities) == H,
  all(vapply(fit$component_probabilities, function(x) {
    max(abs(rowSums(x) - 1)) < 1e-12
  }, logical(1L))),
  all(fit$component_profile_mode_probability > 0 &
        fit$component_profile_mode_probability <= 1)
)

cat("Viroli mixture-update and posterior modal-profile tests passed.\n")
