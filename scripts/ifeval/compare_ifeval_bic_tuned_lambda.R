#!/usr/bin/env Rscript

# Compare BIC-style criteria for the selected IFEval independent-mixture probit
# fit and an ordinary binary probit factor baseline fit with the same tuned
# entrywise loading penalty.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
source(file.path(script_dir, "math500_intercept_imfm_fit.R"))

matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  file.path(bundle_root, "data", "openeval_ifeval_only_binary_matrix.csv")
)
mixture_fit_path <- Sys.getenv(
  "MIXTURE_FIT_PATH",
  file.path(
    bundle_root,
    "results",
    "reproduced_openeval_ifeval_H3_G3_interpretation",
    "openeval_H3_G3_fit.rds"
  )
)
ordinary_fit_path <- Sys.getenv(
  "ORDINARY_FIT_PATH",
  file.path(
    bundle_root,
    "results",
    "reproduced_openeval_ifeval_ordinary_probit_H3_lambda10_visualization",
    "ordinary_probit_full_fit.rds"
  )
)
previous_ordinary_fit_path <- Sys.getenv(
  "PREVIOUS_ORDINARY_FIT_PATH",
  file.path(
    bundle_root,
    "results",
    "reproduced_openeval_ifeval_ordinary_probit_H3_visualization",
    "ordinary_probit_full_fit.rds"
  )
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(bundle_root, "results", "bic_comparisons")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_binary_matrix <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- raw[[1L]]
  X <- X[, colSums(is.na(X)) == 0L, drop = FALSE]
  X <- X[, colSums(X == 1) > 0L & colSums(X == 0) > 0L, drop = FALSE]
  X
}

ordinary_binary_loglik <- function(X, fit) {
  eta <- sweep(fit$F_hat %*% t(fit$Lambda), 2L, fit$alpha, "+")
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  sum(X * log(p1) + (1 - X) * log(p0))
}

ordinary_binary_loglik_parts <- function(X, F_hat, Lambda, alpha) {
  eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  sum(X * log(p1) + (1 - X) * log(p0))
}

make_row <- function(model, version, lambda_l1, H, G, binary_loglik, Lambda,
                     n, p, N, mixture_df = 0L, zero_threshold = 1e-6) {
  nnz_lambda <- sum(abs(Lambda) > zero_threshold)
  effective_df <- p + n * H + nnz_lambda + mixture_df
  data.frame(
    model = model,
    version = version,
    lambda_l1 = lambda_l1,
    H = H,
    G = G,
    binary_loglik = binary_loglik,
    loglik_per_response = binary_loglik / N,
    nnz_lambda = nnz_lambda,
    mixture_df = mixture_df,
    effective_df = effective_df,
    bic = -2 * binary_loglik + log(N) * effective_df,
    stringsAsFactors = FALSE
  )
}

X <- read_binary_matrix(matrix_path)
n <- nrow(X)
p <- ncol(X)
N <- length(X)

mixture_fit <- readRDS(mixture_fit_path)
H <- ncol(mixture_fit$F_hat)
G <- mixture_fit$G_hat[1L]
mixture_ll <- binary_probit_loglik_alpha(
  X,
  mixture_fit$F_hat,
  mixture_fit$Lambda_hat,
  mixture_fit$alpha_hat
)
mixture_lambda_thresholded <- mixture_fit$Lambda_hat
mixture_lambda_thresholded[abs(mixture_lambda_thresholded) < 0.5] <- 0
mixture_thresholded_ll <- binary_probit_loglik_alpha(
  X,
  mixture_fit$F_hat,
  mixture_lambda_thresholded,
  mixture_fit$alpha_hat
)
mixture_df <- sum(3L * mixture_fit$G_hat - 1L)

ordinary_fit <- readRDS(ordinary_fit_path)
ordinary_ll <- ordinary_binary_loglik(X, ordinary_fit)
ordinary_lambda_thresholded <- ordinary_fit$Lambda
ordinary_lambda_thresholded[abs(ordinary_lambda_thresholded) < 0.5] <- 0
ordinary_thresholded_ll <- ordinary_binary_loglik_parts(
  X,
  ordinary_fit$F_hat,
  ordinary_lambda_thresholded,
  ordinary_fit$alpha
)

rows <- list(
  make_row(
    model = "independent-mixture probit",
    version = "full penalized",
    lambda_l1 = 10,
    H = H,
    G = G,
    binary_loglik = mixture_ll,
    Lambda = mixture_fit$Lambda_hat,
    n = n,
    p = p,
    N = N,
    mixture_df = mixture_df
  ),
  make_row(
    model = "independent-mixture probit",
    version = "hard-thresholded |lambda| >= 0.5",
    lambda_l1 = 10,
    H = H,
    G = G,
    binary_loglik = mixture_thresholded_ll,
    Lambda = mixture_lambda_thresholded,
    n = n,
    p = p,
    N = N,
    mixture_df = mixture_df
  ),
  make_row(
    model = "ordinary binary probit",
    version = "full penalized",
    lambda_l1 = 10,
    H = ncol(ordinary_fit$F_hat),
    G = NA_integer_,
    binary_loglik = ordinary_ll,
    Lambda = ordinary_fit$Lambda,
    n = n,
    p = p,
    N = N,
    mixture_df = 0L
  ),
  make_row(
    model = "ordinary binary probit",
    version = "hard-thresholded |lambda| >= 0.5",
    lambda_l1 = 10,
    H = ncol(ordinary_fit$F_hat),
    G = NA_integer_,
    binary_loglik = ordinary_thresholded_ll,
    Lambda = ordinary_lambda_thresholded,
    n = n,
    p = p,
    N = N,
    mixture_df = 0L
  )
)

if (file.exists(previous_ordinary_fit_path)) {
  previous_ordinary_fit <- readRDS(previous_ordinary_fit_path)
  previous_ordinary_ll <- ordinary_binary_loglik(X, previous_ordinary_fit)
  rows[[length(rows) + 1L]] <- make_row(
    model = "ordinary binary probit previous",
    version = "full penalized",
    lambda_l1 = 2,
    H = ncol(previous_ordinary_fit$F_hat),
    G = NA_integer_,
    binary_loglik = previous_ordinary_ll,
    Lambda = previous_ordinary_fit$Lambda,
    n = n,
    p = p,
    N = N,
    mixture_df = 0L
  )
}

comparison <- do.call(rbind, rows)
comparison$delta_bic_vs_best <- comparison$bic - min(comparison$bic)

out_path <- file.path(out_dir, "ifeval_H3_tuned_lambda_bic_comparison.csv")
write.csv(comparison, out_path, row.names = FALSE)

cat("\nBIC-style comparison for IFEval H=3 fits:\n")
print(comparison)
cat("\nOutput saved in: ", normalizePath(out_path), "\n", sep = "")
