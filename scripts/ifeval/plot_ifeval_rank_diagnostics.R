#!/usr/bin/env Rscript

# IFEval rank diagnostic:
#   singular-value shelf from an intercept-only probit augmentation.
#
# This is a descriptive complement to random-cell held-out likelihood tuning.
# It is intentionally lightweight: no model fitting over H is required.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
source(file.path(repo_root, "R", "binary_probit_pretraining.R"))

matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  file.path(repo_root, "data", "ifeval", "openeval_ifeval_only_binary_matrix.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(repo_root, "results", "reproduced_openeval_ifeval_rank_diagnostics")
)
H_max <- as.integer(Sys.getenv("H_MAX", "10"))
center_augmented_Z <- tolower(Sys.getenv("CENTER_AUGMENTED_Z", "TRUE")) %in% c("true", "t", "1", "yes")

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

make_intercept_augmented_Z <- function(X) {
  # Use the same probit intercept initialization as pretraining, then replace
  # each latent Z_ij by E[Z_ij | X_ij, alpha_j] under the intercept-only model.
  alpha <- initialize_binary_intercepts(X)
  F0 <- matrix(0, nrow(X), 1L)
  Lambda0 <- matrix(0, ncol(X), 1L)
  Z <- expected_binary_Z_given_model(
    X = X,
    F_hat = F0,
    Lambda = Lambda0,
    Psi = diag(1, ncol(X)),
    alpha = alpha
  )
  list(Z = Z, alpha = alpha)
}

spectral_rank_table <- function(Z, H_max, center = TRUE) {
  # Singular values of the augmented latent matrix. This is a quick spectral
  # diagnostic used to see whether there is a visible singular-value shelf.
  Z_work <- if (isTRUE(center)) sweep(Z, 2L, colMeans(Z), "-") else Z
  n <- nrow(Z_work)
  p <- ncol(Z_work)
  H_max <- min(as.integer(H_max), n - 1L, p)
  dec <- svd(Z_work, nu = 0, nv = 0)
  d <- dec$d
  H <- seq_len(H_max)

  data.frame(
    H = H,
    singular_value = d[H],
    singular_value_ratio_to_next = d[H] / c(d[H[-1L]], if (length(d) > H_max) d[H_max + 1L] else NA_real_)
  )
}

plot_rank_diagnostics <- function(rank_df, out_dir) {
  png(file.path(out_dir, "ifeval_singular_value_shelf.png"), width = 1500, height = 850, res = 160)
  op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

  plot(rank_df$H, rank_df$singular_value, type = "b", pch = 19, lwd = 2,
       col = "#2B6CB0", xlab = "candidate rank H", ylab = "singular value",
       main = "IFEval Singular-Value Shelf")
  grid(col = "gray85")

  plot(rank_df$H, rank_df$singular_value_ratio_to_next, type = "b", pch = 19, lwd = 2,
       col = "#C05621", xlab = "candidate rank H", ylab = "d_H / d_{H+1}",
       main = "Adjacent Singular-Value Ratio")
  grid(col = "gray85")
  par(op)
  dev.off()
}

X <- read_binary_matrix(matrix_path)
aug <- make_intercept_augmented_Z(X)
rank_df <- spectral_rank_table(aug$Z, H_max = H_max, center = center_augmented_Z)

singular_out <- rank_df[, c("H", "singular_value", "singular_value_ratio_to_next"), drop = FALSE]
write.csv(singular_out, file.path(out_dir, "ifeval_singular_value_shelf.csv"), row.names = FALSE)
write.csv(
  data.frame(item = colnames(X), intercept_only_alpha = aug$alpha),
  file.path(out_dir, "ifeval_intercept_only_alpha_for_rank_diagnostics.csv"),
  row.names = FALSE
)

plot_rank_diagnostics(rank_df, out_dir)

cat("\nIFEval rank diagnostics written to:\n", out_dir, "\n", sep = "")
