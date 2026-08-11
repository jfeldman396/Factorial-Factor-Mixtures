#!/usr/bin/env Rscript

# Generate binary probit data from an H=3, G=3 DGP calibrated to the selected
# IFEval independent-mixture probit fit.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

get_env <- function(name, default, FUN = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  FUN(value)
}

fit_path <- get_env(
  "FIT_PATH",
  file.path(
    bundle_root,
    "results",
    "reproduced_openeval_ifeval_lambda_sparsity_tuning",
    "ifeval_H3_G3_lambda_l1_10_fit.rds"
  ),
  as.character
)
real_matrix_path <- get_env(
  "REAL_MATRIX_PATH",
  file.path(bundle_root, "data", "openeval_ifeval_only_binary_matrix.csv"),
  as.character
)
out_dir <- get_env(
  "OUT_DIR",
  file.path(bundle_root, "results", "ifeval_like_H3_G3_simulated_dgp"),
  as.character
)
sim_n <- get_env("SIM_N", NA_integer_, as.integer)
seed <- get_env("SEED", 20260809L, as.integer)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fit <- readRDS(fit_path)
if (!identical(as.integer(fit$H), 3L) || any(as.integer(fit$G_hat) != 3L)) {
  stop("Expected the selected IFEval fit to have H=3 and G=3 for every factor.")
}

X_real <- fit$X
if (file.exists(real_matrix_path)) {
  raw_real <- read.csv(real_matrix_path, check.names = FALSE)
  X_from_csv <- as.matrix(raw_real[, -1L, drop = FALSE])
  storage.mode(X_from_csv) <- "numeric"
  rownames(X_from_csv) <- raw_real[[1L]]
  common_items <- intersect(colnames(X_real), colnames(X_from_csv))
  X_from_csv <- X_from_csv[, common_items, drop = FALSE]
  X_from_csv <- X_from_csv[, colSums(is.na(X_from_csv)) == 0L, drop = FALSE]
  X_from_csv <- X_from_csv[, colSums(X_from_csv == 1) > 0L & colSums(X_from_csv == 0) > 0L, drop = FALSE]
  X_real <- X_from_csv[, colnames(fit$X), drop = FALSE]
}

n_real <- nrow(X_real)
p <- ncol(X_real)
if (is.na(sim_n)) sim_n <- n_real

alpha <- as.numeric(fit$alpha_hat[colnames(X_real)])
names(alpha) <- colnames(X_real)
Lambda <- fit$Lambda_hat[colnames(X_real), , drop = FALSE]
H <- ncol(Lambda)

set.seed(seed)
F_true <- matrix(NA_real_, sim_n, H)
component <- matrix(NA_integer_, sim_n, H)
colnames(F_true) <- paste0("F", seq_len(H))
colnames(component) <- paste0("F", seq_len(H))

mixture_rows <- list()
for (h in seq_len(H)) {
  mix <- fit$mixture_fits[[h]]
  pi_h <- mix$pi / sum(mix$pi)
  component[, h] <- sample.int(length(pi_h), size = sim_n, replace = TRUE, prob = pi_h)
  F_true[, h] <- rnorm(sim_n, mean = mix$mu[component[, h]], sd = sqrt(mix$var[component[, h]]))
  mixture_rows[[h]] <- data.frame(
    factor = paste0("F", h),
    component = seq_along(pi_h),
    weight = pi_h,
    mean = mix$mu,
    variance = mix$var,
    sd = sqrt(mix$var)
  )
}
mixture_params <- do.call(rbind, mixture_rows)

eta <- sweep(F_true %*% t(Lambda), 2L, alpha, "+")
prob <- pnorm(eta)
Z <- eta + matrix(rnorm(sim_n * p), sim_n, p)
X_sim <- 1L * (Z > 0)
colnames(X_sim) <- colnames(X_real)
rownames(X_sim) <- sprintf("ifeval_like_sim_%03d", seq_len(sim_n))
colnames(prob) <- colnames(X_real)
rownames(prob) <- rownames(X_sim)
colnames(Z) <- colnames(X_real)
rownames(Z) <- rownames(X_sim)
rownames(F_true) <- rownames(X_sim)
rownames(component) <- rownames(X_sim)

write_matrix_with_id <- function(x, path, id_name) {
  out <- data.frame(id = rownames(x), x, check.names = FALSE)
  names(out)[1L] <- id_name
  write.csv(out, path, row.names = FALSE)
}

write_matrix_with_id(X_sim, file.path(out_dir, "ifeval_like_H3_G3_binary_matrix.csv"), "sim_subject_id")
write_matrix_with_id(round(prob, 8), file.path(out_dir, "ifeval_like_H3_G3_probabilities.csv"), "sim_subject_id")
write_matrix_with_id(round(Z, 8), file.path(out_dir, "ifeval_like_H3_G3_latent_Z.csv"), "sim_subject_id")
write_matrix_with_id(round(F_true, 8), file.path(out_dir, "ifeval_like_H3_G3_true_factors.csv"), "sim_subject_id")
write_matrix_with_id(component, file.path(out_dir, "ifeval_like_H3_G3_true_components.csv"), "sim_subject_id")

item_params <- data.frame(
  item_id = colnames(X_real),
  alpha = alpha,
  F1 = Lambda[, 1L],
  F2 = Lambda[, 2L],
  F3 = Lambda[, 3L],
  max_abs_loading = apply(abs(Lambda), 1L, max),
  n_abs_gt_0_5 = rowSums(abs(Lambda) > 0.5),
  empirical_accuracy_real = colMeans(X_real),
  empirical_accuracy_sim = colMeans(X_sim),
  mean_prob_sim = colMeans(prob),
  check.names = FALSE
)
write.csv(item_params, file.path(out_dir, "ifeval_like_H3_G3_dgp_item_parameters.csv"), row.names = FALSE)
write.csv(mixture_params, file.path(out_dir, "ifeval_like_H3_G3_dgp_mixture_parameters.csv"), row.names = FALSE)

profile_id <- apply(component, 1L, paste, collapse = "-")
profile_summary <- as.data.frame(table(profile_id), stringsAsFactors = FALSE)
names(profile_summary) <- c("profile", "n_subjects")
profile_summary$prop_subjects <- profile_summary$n_subjects / sim_n
profile_summary <- profile_summary[order(-profile_summary$n_subjects), ]
write.csv(profile_summary, file.path(out_dir, "ifeval_like_H3_G3_profile_summary.csv"), row.names = FALSE)

component_summary <- do.call(rbind, lapply(seq_len(H), function(h) {
  tab <- table(factor(component[, h], levels = seq_len(3L)))
  data.frame(
    factor = paste0("F", h),
    component = seq_len(3L),
    n_subjects = as.integer(tab),
    prop_subjects = as.numeric(tab) / sim_n
  )
}))
write.csv(component_summary, file.path(out_dir, "ifeval_like_H3_G3_component_summary.csv"), row.names = FALSE)

summary <- data.frame(
  fit_path = fit_path,
  seed = seed,
  n_real = n_real,
  n_sim = sim_n,
  p = p,
  H = H,
  G = 3L,
  real_mean_accuracy = mean(X_real),
  sim_mean_accuracy = mean(X_sim),
  real_item_accuracy_sd = sd(colMeans(X_real)),
  sim_item_accuracy_sd = sd(colMeans(X_sim)),
  real_subject_accuracy_sd = sd(rowMeans(X_real)),
  sim_subject_accuracy_sd = sd(rowMeans(X_sim)),
  item_accuracy_correlation = suppressWarnings(cor(colMeans(X_real), colMeans(X_sim))),
  lambda_nnz = sum(abs(Lambda) > 1e-12),
  lambda_n_abs_gt_0_5 = sum(abs(Lambda) > 0.5),
  lambda_density = mean(abs(Lambda) > 1e-12),
  mean_abs_lambda = mean(abs(Lambda)),
  median_abs_lambda = median(abs(Lambda)),
  n_profiles_occupied = nrow(profile_summary)
)
write.csv(summary, file.path(out_dir, "ifeval_like_H3_G3_simulation_summary.csv"), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  acc_df <- rbind(
    data.frame(type = "real IFEval", accuracy = colMeans(X_real), unit = "item"),
    data.frame(type = "simulated", accuracy = colMeans(X_sim), unit = "item"),
    data.frame(type = "real IFEval", accuracy = rowMeans(X_real), unit = "subject/model"),
    data.frame(type = "simulated", accuracy = rowMeans(X_sim), unit = "subject/model")
  )
  p_acc <- ggplot(acc_df, aes(x = accuracy, fill = type)) +
    geom_histogram(position = "identity", alpha = 0.55, bins = 30) +
    facet_wrap(~unit, scales = "free_y") +
    labs(
      title = "IFEval-like H=3, G=3 simulated data",
      subtitle = "Generated from fitted IFEval alpha, Lambda, and marginal mixture parameters",
      x = "accuracy",
      y = "count",
      fill = NULL
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, "ifeval_like_H3_G3_accuracy_histograms.png"), p_acc, width = 8, height = 4.8, dpi = 200)

  lambda_df <- data.frame(
    item_index = rep(seq_len(p), times = H),
    factor = rep(paste0("F", seq_len(H)), each = p),
    lambda = as.numeric(Lambda)
  )
  p_lam <- ggplot(lambda_df, aes(x = factor, y = item_index, fill = lambda)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
    scale_y_reverse() +
    labs(
      title = "DGP loading matrix from selected IFEval H=3, G=3 fit",
      x = "factor",
      y = "item",
      fill = "lambda"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, "ifeval_like_H3_G3_dgp_lambda_heatmap.png"), p_lam, width = 5.5, height = 7, dpi = 200)

  factor_df <- data.frame(F_true, component)
  names(factor_df) <- c(paste0("F", seq_len(H)), paste0("C", seq_len(H)))
  dens_rows <- do.call(rbind, lapply(seq_len(H), function(h) {
    data.frame(
      factor = paste0("F", h),
      value = F_true[, h],
      component = factor(component[, h])
    )
  }))
  p_fac <- ggplot(dens_rows, aes(x = value, fill = component)) +
    geom_density(alpha = 0.45) +
    facet_wrap(~factor, scales = "free") +
    labs(
      title = "Simulated factor distributions",
      x = "factor value",
      y = "density",
      fill = "component"
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, "ifeval_like_H3_G3_factor_densities.png"), p_fac, width = 8, height = 4.8, dpi = 200)
}

cat("Wrote IFEval-like H=3, G=3 simulated DGP data to:", out_dir, "\n")
print(summary, row.names = FALSE)
