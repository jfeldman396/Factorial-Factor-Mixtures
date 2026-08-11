#!/usr/bin/env Rscript

# Ordinary binary probit factor-model rank selection on an LLM-by-item matrix.
# This is the Gaussian-factor baseline:
#
#   X_ij | f_i ~ Bernoulli(Phi(alpha_j + lambda_j' f_i)),
#   f_i      ~ N(0, I_H).
#
# Rank H is selected by held-out cell predictive log likelihood.

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
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(bundle_root, "results", "openeval_ifeval_cv_H1_10_ordinary")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parse_H_grid <- function(x, default) {
  if (is.na(x) || !nzchar(x)) return(default)
  x <- gsub("[[:space:]]+", "", x)
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    bounds <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq(bounds[1], bounds[2]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

clip01_obs <- function(p, n) {
  eps <- 0.5 / max(1L, n)
  pmin(pmax(p, eps), 1 - eps)
}

read_binary_matrix <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- raw[[1L]]

  complete <- colSums(is.na(X)) == 0L
  nonconstant <- apply(X[, complete, drop = FALSE], 2L, function(z) {
    length(unique(z[is.finite(z)])) > 1L
  })
  X <- X[, complete, drop = FALSE][, nonconstant, drop = FALSE]
  colnames(X) <- paste0("item_", seq_len(ncol(X)))
  X
}

binary_loglik_masked <- function(X, W, F_hat, Lambda, alpha) {
  eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  ll <- X * log(p1) + (1 - X) * log(p0)
  sum(ll[W])
}

score_heldout <- function(X, heldout, fit) {
  eta <- sweep(fit$F_hat %*% t(fit$Lambda), 2L, fit$alpha, "+")
  prob <- pmin(pmax(pnorm(eta), 1e-12), 1 - 1e-12)
  y <- X[heldout]
  phat <- prob[heldout]
  pred <- as.numeric(phat >= 0.5)
  data.frame(
    heldout_loglik_per_response = mean(y * log(phat) + (1 - y) * log(1 - phat)),
    heldout_accuracy = mean(pred == y),
    heldout_brier = mean((phat - y)^2),
    n_heldout = length(y)
  )
}

sample_Z_missing <- function(X, W, F_hat = NULL, Lambda = NULL, alpha = NULL) {
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(F_hat)) {
    n_obs <- colSums(W)
    p_obs <- colSums(X * W) / pmax(n_obs, 1L)
    alpha <- qnorm(clip01_obs(p_obs, n_obs))
    eta <- matrix(rep(alpha, each = n), n, p, dimnames = dimnames(X))
  } else {
    eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  }

  Z <- eta + matrix(rnorm(n * p), n, p)
  dimnames(Z) <- dimnames(X)
  for (j in seq_len(p)) {
    obs <- W[, j]
    if (!any(obs)) next
    mu <- eta[obs, j]
    y <- X[obs, j]
    lower <- ifelse(y == 1, 0, -Inf)
    upper <- ifelse(y == 1, Inf, 0)
    Z[obs, j] <- rtruncnorm_binary_vec(mu, 1, lower, upper)
  }
  Z
}

initialize_ordinary_from_Z <- function(Z, W, H, loading_penalty = 0.05) {
  svd_out <- svd_scores_from_Z(Z, H = H, center_Z = TRUE)
  F_hat <- svd_out$S
  load <- update_working_loadings_missing(Z, W, F_hat, loading_penalty)
  list(F_hat = F_hat, Lambda = load$Lambda, alpha = load$alpha)
}

update_working_loadings_missing <- function(Z, W, F_hat, loading_penalty = 0.05) {
  p <- ncol(Z)
  H <- ncol(F_hat)
  alpha <- numeric(p)
  Lambda <- matrix(0, p, H, dimnames = list(colnames(Z), paste0("factor_", seq_len(H))))

  for (j in seq_len(p)) {
    obs <- W[, j]
    if (sum(obs) <= H + 1L) {
      alpha[j] <- mean(Z[obs, j])
      next
    }
    D <- cbind(1, F_hat[obs, , drop = FALSE])
    y <- Z[obs, j]
    XtX <- crossprod(D) + diag(c(0, rep(1e-5, H)), H + 1L)
    coef <- tryCatch(solve(XtX, crossprod(D, y)), error = function(e) rep(0, H + 1L))
    alpha[j] <- coef[1L]
    Lambda[j, ] <- soft_threshold(coef[-1L], loading_penalty)
  }

  list(alpha = alpha, Lambda = Lambda)
}

update_loadings_probit_missing <- function(X, W, F_hat, Lambda_init, alpha_init,
                                           lambda_l1_penalty, workers) {
  p <- ncol(X)
  H <- ncol(F_hat)
  rows <- parallel_lapply(seq_len(p), function(j) {
    obs <- W[, j]
    y <- X[obs, j]
    if (length(unique(y)) < 2L) {
      a <- qnorm(clip01_obs(mean(y), length(y)))
      return(c(a, rep(0, H)))
    }
    fit_probit_lasso_item_alpha(
      y = y,
      F_hat = F_hat[obs, , drop = FALSE],
      theta_init = c(alpha_init[j], Lambda_init[j, ]),
      lambda_l1_penalty = lambda_l1_penalty,
      maxit = 120L,
      tol = 1e-5
    )
  }, parallel = workers > 1L, workers = workers)

  mat <- do.call(rbind, rows)
  list(
    alpha = as.numeric(mat[, 1L]),
    Lambda = matrix(
      mat[, -1L, drop = FALSE],
      nrow = p,
      ncol = H,
      dimnames = list(colnames(X), paste0("factor_", seq_len(H)))
    )
  )
}

update_one_factor_score_gaussian_prior <- function(x_i, obs_i, f_init, Lambda,
                                                   alpha, prior_weight = 1,
                                                   maxit = 50L) {
  y <- as.numeric(x_i[obs_i])
  L <- Lambda[obs_i, , drop = FALSE]
  a <- alpha[obs_i]

  objective <- function(f) {
    eta <- a + as.numeric(L %*% f)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    ll <- sum(y * log(p1) + (1 - y) * log(p0))
    -ll + 0.5 * prior_weight * sum(f^2)
  }

  gradient <- function(f) {
    eta <- a + as.numeric(L %*% f)
    phi <- dnorm(eta)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    score_eta <- y * phi / p1 - (1 - y) * phi / p0
    -as.numeric(crossprod(L, score_eta)) + prior_weight * f
  }

  opt <- tryCatch(
    optim(
      par = as.numeric(f_init),
      fn = objective,
      gr = gradient,
      method = "L-BFGS-B",
      control = list(maxit = maxit, factr = 1e7)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || !all(is.finite(opt$par))) f_init else opt$par
}

update_factor_scores_gaussian_missing <- function(X, W, F_hat, Lambda, alpha,
                                                  prior_weight, workers) {
  H <- ncol(F_hat)
  rows <- parallel_lapply(seq_len(nrow(X)), function(i) {
    update_one_factor_score_gaussian_prior(
      x_i = X[i, ],
      obs_i = W[i, ],
      f_init = F_hat[i, ],
      Lambda = Lambda,
      alpha = alpha,
      prior_weight = prior_weight,
      maxit = 50L
    )
  }, parallel = workers > 1L, workers = workers)

  out <- matrix(unlist(rows), ncol = H, byrow = TRUE)
  rownames(out) <- rownames(F_hat)
  colnames(out) <- colnames(F_hat)
  out
}

fit_ordinary_missing_aware_probit_factor <- function(
    X,
    W,
    H,
    fold,
    workers,
    n_aug_iter = 4L,
    n_refine_iter = 5L,
    pretrain_loading_penalty = 0.05,
    lambda_l1_penalty = 2,
    prior_weight = 1,
    seed = 20260724L) {
  set.seed(seed + 10000L * fold + 100L * H)
  Z <- sample_Z_missing(X, W)
  history <- data.frame()
  current <- NULL

  for (iter in seq_len(n_aug_iter)) {
    current <- initialize_ordinary_from_Z(Z, W, H, pretrain_loading_penalty)
    Z <- sample_Z_missing(X, W, current$F_hat, current$Lambda, current$alpha)
    history <- rbind(history, data.frame(
      stage = "pretrain",
      iter = iter,
      train_loglik_per_response = binary_loglik_masked(
        X, W, current$F_hat, current$Lambda, current$alpha
      ) / sum(W)
    ))
  }

  F_hat <- current$F_hat
  Lambda <- current$Lambda
  alpha <- current$alpha

  load <- update_loadings_probit_missing(
    X, W, F_hat, Lambda, alpha, lambda_l1_penalty, workers
  )
  Lambda <- load$Lambda
  alpha <- load$alpha

  prev_score <- -Inf
  for (iter in seq_len(n_refine_iter)) {
    F_hat <- update_factor_scores_gaussian_missing(
      X, W, F_hat, Lambda, alpha, prior_weight, workers
    )
    load <- update_loadings_probit_missing(
      X, W, F_hat, Lambda, alpha, lambda_l1_penalty, workers
    )
    Lambda <- load$Lambda
    alpha <- load$alpha

    train_ll <- binary_loglik_masked(X, W, F_hat, Lambda, alpha) / sum(W)
    history <- rbind(history, data.frame(
      stage = "refine",
      iter = iter,
      train_loglik_per_response = train_ll
    ))
    if (iter >= 3L && is.finite(prev_score) && abs(train_ll - prev_score) < 1e-4) break
    prev_score <- train_ll
  }

  list(F_hat = F_hat, Lambda = Lambda, alpha = alpha, history = history)
}

X_full <- read_binary_matrix(matrix_path)
max_feasible_H <- min(nrow(X_full) - 1L, ncol(X_full))
H_grid <- parse_H_grid(Sys.getenv("H_GRID"), default = 1:min(10L, max_feasible_H))
H_grid <- H_grid[H_grid <= max_feasible_H]
K_folds <- as.integer(Sys.getenv("K_FOLDS", "3"))
workers <- as.integer(Sys.getenv("WORKERS", "8"))
n_aug_iter <- as.integer(Sys.getenv("AUG_ITER", "4"))
n_refine_iter <- as.integer(Sys.getenv("REFINE_ITER", "5"))
lambda_l1_penalty <- as.numeric(Sys.getenv("LAMBDA_L1", "2"))

set.seed(20260724L)
fold_id <- matrix(
  sample.int(K_folds, length(X_full), replace = TRUE),
  nrow = nrow(X_full),
  ncol = ncol(X_full),
  dimnames = dimnames(X_full)
)

score_rows <- list()
history_rows <- list()

for (H in H_grid) {
  for (fold in seq_len(K_folds)) {
    message("Ordinary probit CV: H=", H, ", fold=", fold, "/", K_folds)
    heldout <- fold_id == fold
    W <- !heldout
    start_time <- proc.time()[["elapsed"]]
    fit <- fit_ordinary_missing_aware_probit_factor(
      X = X_full,
      W = W,
      H = H,
      fold = fold,
      workers = workers,
      n_aug_iter = n_aug_iter,
      n_refine_iter = n_refine_iter,
      lambda_l1_penalty = lambda_l1_penalty
    )
    elapsed <- proc.time()[["elapsed"]] - start_time
    sc <- score_heldout(X_full, heldout, fit)
    train_ll <- binary_loglik_masked(X_full, W, fit$F_hat, fit$Lambda, fit$alpha) / sum(W)

    score_rows[[length(score_rows) + 1L]] <- cbind(
      method = "ordinary_probit_factor",
      H = H,
      fold = fold,
      sc,
      train_loglik_per_response = train_ll,
      fit_seconds = elapsed
    )
    history_rows[[length(history_rows) + 1L]] <- cbind(
      method = "ordinary_probit_factor",
      H = H,
      fold = fold,
      fit$history
    )

    write.csv(
      do.call(rbind, score_rows),
      file.path(out_dir, "ordinary_probit_factor_fold_scores_partial.csv"),
      row.names = FALSE
    )
    write.csv(
      do.call(rbind, history_rows),
      file.path(out_dir, "ordinary_probit_factor_histories_partial.csv"),
      row.names = FALSE
    )
    saveRDS(fit, file.path(out_dir, sprintf("ordinary_fit_H%02d_fold%02d.rds", H, fold)))
  }
}

fold_scores <- do.call(rbind, score_rows)
write.csv(fold_scores, file.path(out_dir, "ordinary_probit_factor_fold_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, history_rows), file.path(out_dir, "ordinary_probit_factor_histories.csv"), row.names = FALSE)

summary_scores <- aggregate(
  heldout_loglik_per_response ~ H,
  data = fold_scores,
  FUN = function(z) c(mean = mean(z), sd = sd(z), se = sd(z) / sqrt(length(z)))
)
summary_scores <- data.frame(
  H = summary_scores$H,
  mean_heldout_loglik_per_response = summary_scores$heldout_loglik_per_response[, "mean"],
  sd_heldout_loglik_per_response = summary_scores$heldout_loglik_per_response[, "sd"],
  se_heldout_loglik_per_response = summary_scores$heldout_loglik_per_response[, "se"]
)
acc_summary <- aggregate(
  heldout_accuracy ~ H,
  data = fold_scores,
  FUN = function(z) c(mean = mean(z), sd = sd(z), se = sd(z) / sqrt(length(z)))
)
summary_scores$mean_heldout_accuracy <- acc_summary$heldout_accuracy[, "mean"]
summary_scores$se_heldout_accuracy <- acc_summary$heldout_accuracy[, "se"]
summary_scores$mean_fit_seconds <- aggregate(fit_seconds ~ H, fold_scores, mean)$fit_seconds

best <- summary_scores[which.max(summary_scores$mean_heldout_loglik_per_response), ]
threshold <- best$mean_heldout_loglik_per_response - best$se_heldout_loglik_per_response
summary_scores$best_H <- summary_scores$H == best$H
summary_scores$selected_by_one_se <- summary_scores$H == min(summary_scores$H[
  summary_scores$mean_heldout_loglik_per_response >= threshold
])
summary_scores$one_se_threshold <- threshold
write.csv(summary_scores, file.path(out_dir, "ordinary_probit_H_summary.csv"), row.names = FALSE)

png(file.path(out_dir, "ordinary_probit_H_lineplot.png"), width = 1500, height = 950, res = 160)
plot(
  summary_scores$H,
  summary_scores$mean_heldout_loglik_per_response,
  type = "b",
  pch = 19,
  col = "#4A5568",
  ylim = range(
    summary_scores$mean_heldout_loglik_per_response -
      summary_scores$se_heldout_loglik_per_response,
    summary_scores$mean_heldout_loglik_per_response +
      summary_scores$se_heldout_loglik_per_response
  ),
  xlab = "number of factors H",
  ylab = "held-out log likelihood per response",
  main = "Ordinary probit factor H selection, OpenEval"
)
arrows(
  summary_scores$H,
  summary_scores$mean_heldout_loglik_per_response - summary_scores$se_heldout_loglik_per_response,
  summary_scores$H,
  summary_scores$mean_heldout_loglik_per_response + summary_scores$se_heldout_loglik_per_response,
  angle = 90,
  code = 3,
  length = 0.035,
  col = "#4A5568"
)
abline(h = threshold, col = "#BC4749", lty = 2)
dev.off()

cat("\nData dimensions after complete/nonconstant filtering:\n")
cat("models n=", nrow(X_full), ", items p=", ncol(X_full), "\n", sep = "")
cat("\nOrdinary probit CV summary:\n")
print(summary_scores)
cat("\nOutputs saved in: ", normalizePath(out_dir), "\n", sep = "")
