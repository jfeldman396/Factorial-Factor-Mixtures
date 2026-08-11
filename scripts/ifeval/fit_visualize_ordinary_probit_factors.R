#!/usr/bin/env Rscript

# Fit and visualize the ordinary binary probit factor model on the IFEval
# complete/nonconstant item matrix. This reproducible copy defaults to H=3.

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
  "MMLU_MATRIX",
  file.path(bundle_root, "data", "openeval_ifeval_only_binary_matrix.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(bundle_root, "results", "reproduced_openeval_ifeval_ordinary_probit_H3_visualization")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

H <- as.integer(Sys.getenv("H", "3"))
workers <- as.integer(Sys.getenv("WORKERS", "8"))
n_aug_iter <- as.integer(Sys.getenv("AUG_ITER", "4"))
n_refine_iter <- as.integer(Sys.getenv("REFINE_ITER", "10"))
lambda_l1_penalty <- as.numeric(Sys.getenv("LAMBDA_L1", "2"))

read_mmlu_binary_matrix <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  X0 <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(X0) <- "numeric"
  rownames(X0) <- raw[[1L]]

  complete <- colSums(is.na(X0)) == 0L
  nonconstant <- apply(X0[, complete, drop = FALSE], 2L, function(z) {
    length(unique(z[is.finite(z)])) > 1L
  })
  X <- X0[, complete, drop = FALSE][, nonconstant, drop = FALSE]
  colnames(X) <- paste0("item_", seq_len(ncol(X)))
  storage.mode(X) <- "numeric"
  X
}

binary_loglik_full <- function(X, F_hat, Lambda, alpha) {
  eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  sum(X * log(p1) + (1 - X) * log(p0))
}

sample_Z_full <- function(X, F_hat = NULL, Lambda = NULL, alpha = NULL) {
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(F_hat)) {
    alpha <- qnorm(clip01(colMeans(X), n))
    eta <- matrix(rep(alpha, each = n), n, p, dimnames = dimnames(X))
  } else {
    eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  }

  Z <- matrix(NA_real_, n, p, dimnames = dimnames(X))
  for (j in seq_len(p)) {
    lower <- ifelse(X[, j] == 1, 0, -Inf)
    upper <- ifelse(X[, j] == 1, Inf, 0)
    Z[, j] <- rtruncnorm_binary_vec(eta[, j], 1, lower, upper)
  }
  Z
}

update_working_loadings_full <- function(Z, F_hat, loading_penalty = 0.05) {
  alpha <- colMeans(Z)
  Zc <- sweep(Z, 2L, alpha, "-")
  Lambda_ls <- crossprod(Zc, F_hat) / nrow(Z)
  Lambda <- soft_threshold(Lambda_ls, loading_penalty)
  rownames(Lambda) <- colnames(Z)
  colnames(Lambda) <- paste0("factor_", seq_len(ncol(F_hat)))
  list(alpha = alpha, Lambda = Lambda)
}

update_one_factor_score_gaussian_prior <- function(x_i, f_init, Lambda, alpha,
                                                   prior_weight = 1,
                                                   maxit = 50L) {
  y <- as.numeric(x_i)

  objective <- function(f) {
    eta <- alpha + as.numeric(Lambda %*% f)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    ll <- sum(y * log(p1) + (1 - y) * log(p0))
    -ll + 0.5 * prior_weight * sum(f^2)
  }

  gradient <- function(f) {
    eta <- alpha + as.numeric(Lambda %*% f)
    phi <- dnorm(eta)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    score_eta <- y * phi / p1 - (1 - y) * phi / p0
    -as.numeric(crossprod(Lambda, score_eta)) + prior_weight * f
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

update_factor_scores_gaussian <- function(X, F_hat, Lambda, alpha, workers) {
  H <- ncol(F_hat)
  rows <- parallel_lapply(seq_len(nrow(X)), function(i) {
    update_one_factor_score_gaussian_prior(
      x_i = X[i, ],
      f_init = F_hat[i, ],
      Lambda = Lambda,
      alpha = alpha
    )
  }, parallel = workers > 1L, workers = workers)
  out <- matrix(unlist(rows), ncol = H, byrow = TRUE)
  rownames(out) <- rownames(F_hat)
  colnames(out) <- colnames(F_hat)
  out
}

fit_ordinary_probit_factor_full <- function(X, H, workers, n_aug_iter,
                                            n_refine_iter, lambda_l1_penalty) {
  set.seed(20260727L)
  Z <- sample_Z_full(X)

  current <- NULL
  history <- data.frame()
  for (iter in seq_len(n_aug_iter)) {
    svd_out <- svd_scores_from_Z(Z, H = H, center_Z = TRUE)
    F_hat <- svd_out$S
    load <- update_working_loadings_full(Z, F_hat, loading_penalty = 0.05)
    Z <- sample_Z_full(X, F_hat, load$Lambda, load$alpha)
    current <- list(F_hat = F_hat, Lambda = load$Lambda, alpha = load$alpha)
    history <- rbind(history, data.frame(
      stage = "pretrain",
      iter = iter,
      loglik_per_response = binary_loglik_full(X, F_hat, load$Lambda, load$alpha) / length(X)
    ))
  }

  F_hat <- current$F_hat
  Lambda <- current$Lambda
  alpha <- current$alpha

  load <- update_binary_probit_loadings_glm_alpha(
    X = X,
    F_hat = F_hat,
    Lambda_init = Lambda,
    alpha_init = alpha,
    lambda_l1_penalty = lambda_l1_penalty,
    lasso_maxit = 120L,
    lasso_tol = 1e-5,
    parallel = workers > 1L,
    workers = workers
  )
  Lambda <- load$Lambda
  alpha <- load$alpha

  prev <- -Inf
  for (iter in seq_len(n_refine_iter)) {
    F_hat <- update_factor_scores_gaussian(X, F_hat, Lambda, alpha, workers)
    load <- update_binary_probit_loadings_glm_alpha(
      X = X,
      F_hat = F_hat,
      Lambda_init = Lambda,
      alpha_init = alpha,
      lambda_l1_penalty = lambda_l1_penalty,
      lasso_maxit = 120L,
      lasso_tol = 1e-5,
      parallel = workers > 1L,
      workers = workers
    )
    Lambda <- load$Lambda
    alpha <- load$alpha
    ll <- binary_loglik_full(X, F_hat, Lambda, alpha) / length(X)
    history <- rbind(history, data.frame(
      stage = "refine",
      iter = iter,
      loglik_per_response = ll
    ))
    if (iter >= 3L && abs(ll - prev) < 1e-4) break
    prev <- ll
  }

  list(F_hat = F_hat, Lambda = Lambda, alpha = alpha, history = history)
}

label_points <- function(x, y, labels, cex = 0.55) {
  text(x, y, labels = labels, pos = 3, cex = cex, col = "#2D3748")
}

X <- read_mmlu_binary_matrix(matrix_path)
fit <- fit_ordinary_probit_factor_full(
  X = X,
  H = H,
  workers = workers,
  n_aug_iter = n_aug_iter,
  n_refine_iter = n_refine_iter,
  lambda_l1_penalty = lambda_l1_penalty
)

acc <- rowMeans(X)

# Orient factor 1 to overall accuracy; orient factor 2 so its largest loading
# is positive.  Signs remain arbitrary, but this makes the plot reproducible.
if (cor(fit$F_hat[, 1L], acc) < 0) {
  fit$F_hat[, 1L] <- -fit$F_hat[, 1L]
  fit$Lambda[, 1L] <- -fit$Lambda[, 1L]
}
if (H >= 2L) {
  idx <- which.max(abs(fit$Lambda[, 2L]))
  if (fit$Lambda[idx, 2L] < 0) {
    fit$F_hat[, 2L] <- -fit$F_hat[, 2L]
    fit$Lambda[, 2L] <- -fit$Lambda[, 2L]
  }
}

scores <- data.frame(
  model = rownames(X),
  accuracy = acc,
  fit$F_hat,
  check.names = FALSE
)
write.csv(scores, file.path(out_dir, "ordinary_probit_factor_scores.csv"), row.names = FALSE)
write.csv(fit$Lambda, file.path(out_dir, "ordinary_probit_lambda.csv"), row.names = TRUE)
write.csv(fit$history, file.path(out_dir, "ordinary_probit_fit_history.csv"), row.names = FALSE)
saveRDS(fit, file.path(out_dir, "ordinary_probit_full_fit.rds"))

png(file.path(out_dir, "ordinary_probit_factor_scatter_labeled.png"), width = 1500, height = 1100, res = 160)
par(mar = c(5, 5, 4, 2))
pal <- colorRampPalette(c("#B91C1C", "#F59E0B", "#2563EB"))(100)
col <- pal[pmax(1L, pmin(100L, as.integer(cut(acc, breaks = 100, labels = FALSE))))]
plot(
  fit$F_hat[, 1L],
  fit$F_hat[, 2L],
  pch = 19,
  col = col,
  xlab = "ordinary probit factor 1",
  ylab = "ordinary probit factor 2",
  main = "Ordinary probit factor scores, IFEval"
)
abline(h = 0, v = 0, lty = 3, col = "gray60")
label_points(fit$F_hat[, 1L], fit$F_hat[, 2L], rownames(X), cex = 0.45)
legend("bottomright", legend = c("lower accuracy", "higher accuracy"),
       col = c("#B91C1C", "#2563EB"), pch = 19, bty = "n")
dev.off()

png(file.path(out_dir, "ordinary_probit_factor_scores_by_llm.png"), width = 1700, height = 1050, res = 160)
ord <- order(acc)
ylim <- range(fit$F_hat)
plot(
  seq_along(ord),
  fit$F_hat[ord, 1L],
  type = "b",
  pch = 19,
  col = "#2B6CB0",
  ylim = ylim,
  xaxt = "n",
  xlab = "LLM, ordered by empirical accuracy",
  ylab = "factor score",
  main = "Ordinary probit factor scores by LLM"
)
if (H >= 2L) lines(seq_along(ord), fit$F_hat[ord, 2L], type = "b", pch = 17, col = "#C05621")
axis(1, at = seq_along(ord), labels = rownames(X)[ord], las = 2, cex.axis = 0.45)
abline(h = 0, lty = 3, col = "gray60")
legend("topleft", legend = paste0("factor ", seq_len(H)), col = c("#2B6CB0", "#C05621")[seq_len(H)],
       pch = c(19, 17)[seq_len(H)], lty = 1, bty = "n")
dev.off()

png(file.path(out_dir, "ordinary_probit_lambda_heatmap.png"), width = 1250, height = 1000, res = 160)
op <- par(no.readonly = TRUE)
par(mar = c(5, 5, 4, 6))
L <- fit$Lambda
zlim <- max(abs(L))
image(
  x = seq_len(ncol(L)),
  y = seq_len(nrow(L)),
  z = t(L[nrow(L):1L, , drop = FALSE]),
  col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  zlim = c(-zlim, zlim),
  axes = FALSE,
  xlab = "factor",
  ylab = "item",
  main = "Ordinary probit loadings"
)
axis(1, at = seq_len(ncol(L)), labels = seq_len(ncol(L)))
axis(2, at = pretty(seq_len(nrow(L))), labels = rev(pretty(seq_len(nrow(L)))))
box()
par(op)
dev.off()

cat("\nFitted ordinary probit factor model on full IFEval matrix.\n")
cat("n=", nrow(X), ", p=", ncol(X), ", H=", H, "\n", sep = "")
cat("Final loglik per response=", tail(fit$history$loglik_per_response, 1L), "\n", sep = "")
cat("Outputs saved in: ", normalizePath(out_dir), "\n", sep = "")
