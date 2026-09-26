#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this test with Rscript.")
repo_root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."))

source(file.path(repo_root, "R", "binary_probit_pretraining.R"))
source(file.path(repo_root, "R", "canonical_factor_normalization.R"))
source(file.path(repo_root, "R", "viroli_probit_independent_gibbs.R"))

if (!requireNamespace("clue", quietly = TRUE)) {
  stop("Package 'clue' is required for this test.")
}

set.seed(42)
n <- 20L
p <- 30L
H <- 3L
G <- c(2L, 3L, 2L)

F_reference <- matrix(rnorm(n * H), n, H)
Lambda_reference <- matrix(rnorm(p * H), p, H)
C_reference <- cbind(
  sample.int(2L, n, replace = TRUE),
  sample.int(3L, n, replace = TRUE),
  sample.int(2L, n, replace = TRUE)
)
pi_reference <- rbind(c(0.4, 0.6, NA), c(0.2, 0.5, 0.3), c(0.55, 0.45, NA))
mu_reference <- rbind(c(-1.0, 1.0, NA), c(-1.4, 0.0, 1.4), c(-0.8, 0.8, NA))
sig2_reference <- rbind(c(0.3, 0.4, NA), c(0.2, 0.5, 0.3), c(0.6, 0.4, NA))

# Swap the two G=2 factors and independently reverse two orientations.
current_order <- c(3L, 2L, 1L)
current_sign <- c(-1, -1, 1)
F_current <- sweep(F_reference[, current_order, drop = FALSE], 2L, current_sign, "*")
Lambda_current <- sweep(Lambda_reference[, current_order, drop = FALSE], 2L, current_sign, "*")
C_current <- C_reference[, current_order, drop = FALSE]
pi_current <- pi_reference[current_order, , drop = FALSE]
mu_current <- sweep(mu_reference[current_order, , drop = FALSE], 1L, current_sign, "*")
sig2_current <- sig2_reference[current_order, , drop = FALSE]

# Match the sampler's increasing-component-mean convention before retention.
current_sorted <- sort_viroli_components(
  C_current,
  pi_current,
  mu_current,
  sig2_current,
  G
)

aligned <- align_viroli_draw_to_reference(
  reference_Lambda = Lambda_reference,
  F = F_current,
  Lambda = Lambda_current,
  C = current_sorted$C,
  pi_mat = current_sorted$pi,
  mu_mat = current_sorted$mu,
  sig2_mat = current_sorted$sig2,
  G = G
)

stopifnot(
  isTRUE(all.equal(aligned$F, F_reference, tolerance = 1e-10)),
  isTRUE(all.equal(aligned$Lambda, Lambda_reference, tolerance = 1e-10)),
  identical(unname(aligned$C), unname(C_reference)),
  isTRUE(all.equal(aligned$pi, pi_reference, tolerance = 1e-10)),
  isTRUE(all.equal(aligned$mu, mu_reference, tolerance = 1e-10)),
  isTRUE(all.equal(aligned$sig2, sig2_reference, tolerance = 1e-10)),
  aligned$alignment$n_permuted_factors == 2L,
  aligned$alignment$n_sign_flips == 2L,
  max(abs(F_current %*% t(Lambda_current) - aligned$F %*% t(aligned$Lambda))) < 1e-10
)

cat("Viroli retained-draw alignment test passed.\n")
