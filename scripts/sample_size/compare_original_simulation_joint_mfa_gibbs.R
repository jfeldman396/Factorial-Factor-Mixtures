#!/usr/bin/env Rscript

# Compare the independent marginal-mixture probit factor algorithm against two
# Gibbs samplers under the canonical sample-size simulation DGP:
#   1. an unrestricted joint-mixture factor prior with K = prod_h G_h classes;
#   2. a Viroli-style independent-mixture factor prior with one univariate
#      mixture per factor coordinate.
#
# Baseline model:
#
#   eta_i | c_i = k ~ N(mu_k, diag(sigma_k^2)),  k = 1,...,K
#   Z_i | eta_i   ~ N(alpha + Lambda eta_i, I_p)
#   X_ij          = 1{Z_ij > 0}
#
# with K = G^H for the unrestricted joint baseline.  The simulation DGP and
# loading design aliases are centralized in R/sample_size_dgp.R so that "Sparse"
# and "Cross" mean the same thing across drivers, diagnostics, and plots.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
file_arg <- gsub("~\\+~", " ", file_arg)
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "../.."))
source(file.path(repo_root, "R", "binary_probit_pretraining.R"))
source(file.path(repo_root, "R", "binary_probit_refinement.R"))
source(file.path(repo_root, "R", "probit_ifa_em_svd_pretraining.R"))
source(file.path(repo_root, "R", "viroli_probit_independent_gibbs.R"))
source(file.path(repo_root, "R", "sample_size_dgp.R"))

get_env <- function(name, default, FUN = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  FUN(value)
}

split_csv <- function(x) {
  x <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  x[nzchar(x)]
}

parse_int_csv <- function(x) as.integer(split_csv(x))
parse_num_csv <- function(x) as.numeric(split_csv(x))
parse_g_count_string <- function(x) parse_int_csv(gsub("-", ",", x, fixed = TRUE))

parse_g_configs <- function(x) {
  configs <- trimws(strsplit(x, ";", fixed = TRUE)[[1L]])
  configs <- configs[nzchar(configs)]
  lapply(configs, parse_int_csv)
}

normalize_G_counts <- function(G, H) {
  G <- as.integer(G)
  if (length(G) == 1L) G <- rep(G, H)
  if (length(G) != H) stop("Each G configuration must have length 1 or length H.")
  if (any(!is.finite(G)) || any(G < 1L)) stop("All component counts must be positive integers.")
  G
}

format_G_config <- function(G) paste(as.integer(G), collapse = "-")

parse_np_grid <- function(x) {
  parts <- split_csv(x)
  out <- do.call(rbind, lapply(parts, function(part) {
    z <- strsplit(part, ":", fixed = TRUE)[[1L]]
    if (length(z) != 3L) stop("NP_GRID entries must look like label:n:p")
    data.frame(np_regime = z[1L], n = as.integer(z[2L]), p = as.integer(z[3L]))
  }))
  rownames(out) <- NULL
  out
}

out_dir <- get_env(
  "OUT_DIR",
  file.path(dirname(script_dir), "results", "original_simulation_joint_mfa_comparison")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seed_base <- get_env("SEED", 20260731L, as.integer)
rep_values <- get_env("REP_VALUES", 1L, parse_int_csv)
H_values <- get_env("H_VALUES", get_env("H_TRUE", 4L, as.integer), parse_int_csv)
G_values <- get_env("G_VALUES", get_env("G_TRUE", 2L, as.integer), parse_int_csv)
G_configs_env <- Sys.getenv("G_CONFIGS", unset = "")
G_configs_input <- if (nzchar(G_configs_env)) parse_g_configs(G_configs_env) else NULL
separations <- get_env("SEPARATIONS", 1.0, parse_num_csv)
mixture_param_mode <- get_env("MIXTURE_PARAM_MODE", "equal", as.character)
mixture_variance_mode <- get_env("MIXTURE_VARIANCE_MODE", "unequal", as.character)
intercept_mode <- get_env("INTERCEPT_MODE", "none", as.character)
intercept_sd <- get_env("INTERCEPT_SD", 0.45, as.numeric)
intercept_block_span <- get_env("INTERCEPT_BLOCK_SPAN", 1.6, as.numeric)
intercept_clip <- get_env("INTERCEPT_CLIP", 1.75, as.numeric)
loading_sign_mode <- get_env("LOADING_SIGN_MODE", "block", as.character)
if (!loading_sign_mode %in% c("block", "smoke")) {
  stop("LOADING_SIGN_MODE must be either 'block' or 'smoke'.")
}
alignment_mode <- get_env("ALIGNMENT_MODE", "factors", as.character)
if (!alignment_mode %in% c("factors", "loadings")) {
  stop("ALIGNMENT_MODE must be either 'factors' or 'loadings'.")
}
dgp_order_mode <- get_env("DGP_ORDER_MODE", "driver", as.character)
if (!dgp_order_mode %in% c("driver", "viroli_smoke_exact")) {
  stop("DGP_ORDER_MODE must be either 'driver' or 'viroli_smoke_exact'.")
}
dgp_seed_mode <- get_env("DGP_SEED_MODE", "driver", as.character)
if (!dgp_seed_mode %in% c("driver", "viroli_smoke")) {
  stop("DGP_SEED_MODE must be either 'driver' or 'viroli_smoke'.")
}
g3_raw_sd <- get_env("G3_RAW_SD", c(0.25, 0.50, 0.75), parse_num_csv)
if (length(g3_raw_sd) != 3L || any(!is.finite(g3_raw_sd)) || any(g3_raw_sd <= 0)) {
  stop("G3_RAW_SD must contain three positive numbers, e.g. G3_RAW_SD='0.25,0.50,0.30'.")
}
loading_designs <- get_env(
  "LOADING_DESIGNS",
  c("balanced_moderate_few_positive_cross", "balanced_moderate_dense_signed_cross"),
  split_csv
)
loading_designs <- normalize_sample_size_loading_design(loading_designs)
block_size_mode <- get_env("BLOCK_SIZE_MODE", "balanced", as.character)
np_settings <- get_env(
  "NP_GRID",
  data.frame(np_regime = "smoke_wide", n = 80L, p = 120L),
  parse_np_grid
)

run_ours <- get_env("RUN_OURS", TRUE, as.logical)
run_joint_mfa <- get_env("RUN_JOINT_MFA", TRUE, as.logical)
run_viroli <- get_env("RUN_VIROLI", run_joint_mfa, as.logical)
resume_existing <- get_env("RESUME_EXISTING", TRUE, as.logical)

ours_convergence_summary <- function(fit) {
  pre <- fit$pretrain_fit
  ref <- fit$refine_fit$joint_refinement
  rotation_converged <- if (!is.null(pre$rotation_converged)) isTRUE(pre$rotation_converged) else NA
  rotation_completed_outer <- if (!is.null(pre$rotation_completed_outer)) pre$rotation_completed_outer else NA_integer_
  refinement_all_mixtures_converged <- if (!is.null(ref$history) &&
                                           "all_mixtures_converged" %in% names(ref$history)) {
    tail(ref$history$all_mixtures_converged, 1L)
  } else {
    NA
  }
  data.frame(
    pretraining_converged = isTRUE(pre$pretraining_converged),
    pretraining_completed_iter = if (!is.null(pre$pretraining_completed_iter)) pre$pretraining_completed_iter else NA_integer_,
    selected_pretraining_iteration = if (!is.null(pre$selected_pretraining_iteration)) pre$selected_pretraining_iteration else NA_integer_,
    pretraining_return_best_iteration = isTRUE(pre$return_best_iteration),
    rotation_converged = rotation_converged,
    rotation_completed_outer = rotation_completed_outer,
    refinement_converged = isTRUE(ref$converged),
    refinement_completed_iter = if (!is.null(ref$n_completed)) ref$n_completed else NA_integer_,
    selected_refinement_iteration = if (!is.null(ref$selected_refinement_iteration)) ref$selected_refinement_iteration else NA_integer_,
    refinement_return_best_iteration = isTRUE(ref$return_best_refinement_iteration),
    refinement_stopping_objective = if (!is.null(ref$stopping_objective)) ref$stopping_objective else NA_character_,
    refinement_all_mixtures_converged = refinement_all_mixtures_converged,
    refinement_monotone_guard_triggered = if (!is.null(ref$history) &&
                                               "iteration_rejected_by_monotone_guard" %in% names(ref$history)) {
      any(ref$history$iteration_rejected_by_monotone_guard, na.rm = TRUE)
    } else {
      NA
    },
    stringsAsFactors = FALSE
  )
}

pretrain_aug_iter <- get_env("PRETRAIN_AUG_ITER", 4L, as.integer)
# The core estimator uses deterministic low-rank probit EM-SVD pretraining:
# estimate the rank-H signal L, take its SVD, then rotate the left singular
# vectors toward independent marginal mixtures.  sampled_z is retained only as
# a legacy ablation.
ours_pretraining_method <- get_env("OURS_PRETRAINING_METHOD", "em_svd", as.character)
if (!ours_pretraining_method %in% c("sampled_z", "em_svd")) {
  stop("OURS_PRETRAINING_METHOD must be either 'sampled_z' or 'em_svd'.")
}
pretrain_z_update <- get_env("PRETRAIN_Z_UPDATE", "sample", as.character)
if (!pretrain_z_update %in% c("sample", "expectation")) {
  stop("PRETRAIN_Z_UPDATE must be either 'sample' or 'expectation'.")
}
pretrain_min_aug_iter <- get_env("PRETRAIN_MIN_AUG_ITER", 2L, as.integer)
pretrain_objective <- get_env("PRETRAIN_OBJECTIVE", "full_data_loglik", as.character)
pretrain_objective_tolerance <- get_env("PRETRAIN_OBJECTIVE_TOLERANCE", NA_real_, as.numeric)
if (is.na(pretrain_objective_tolerance)) pretrain_objective_tolerance <- NULL
pretrain_objective_patience <- get_env("PRETRAIN_OBJECTIVE_PATIENCE", 0L, as.integer)
pretrain_return_best_iteration <- get_env("PRETRAIN_RETURN_BEST_ITERATION", FALSE, as.logical)
center_Z_for_svd <- get_env("CENTER_Z_FOR_SVD", intercept_mode != "none", as.logical)
refine_iter <- get_env("REFINE_ITER", 4L, as.integer)
refine_min_iter <- get_env("REFINE_MIN_ITER", 1L, as.integer)
refine_objective_tolerance <- get_env("REFINE_OBJECTIVE_TOLERANCE", 1e-5, as.numeric)
refine_stopping_objective <- get_env("REFINE_STOPPING_OBJECTIVE", "posterior_objective", as.character)
refine_return_best_iteration <- get_env("REFINE_RETURN_BEST_ITERATION", FALSE, as.logical)
refine_selection_objective <- get_env("REFINE_SELECTION_OBJECTIVE", refine_stopping_objective, as.character)
factor_update <- get_env("FACTOR_UPDATE", "marginal", as.character)
lambda_l1_penalty <- get_env("LAMBDA_L1_PENALTY", 5, as.numeric)
pretrain_loading_penalty <- get_env("PRETRAIN_LOADING_PENALTY", 0, as.numeric)
lasso_backend <- get_env("LASSO_BACKEND", "proximal", as.character)
if (!lasso_backend %in% c("glmnet", "proximal")) {
  stop("LASSO_BACKEND must be either 'glmnet' or 'proximal'.")
}
glmnet_standardize <- get_env("GLMNET_STANDARDIZE", FALSE, as.logical)
mixture_max_iter <- get_env("MIXTURE_MAX_ITER", 20L, as.integer)
mixture_update <- get_env("MIXTURE_UPDATE", "map", as.character)
mixture_refit <- get_env("MIXTURE_REFIT", "em", as.character)
mu_prior_mean <- get_env("MU_PRIOR_MEAN", 0, as.numeric)
mu_prior_kappa <- get_env("MU_PRIOR_KAPPA", 0.01, as.numeric)
var_prior_shape <- get_env("VAR_PRIOR_SHAPE", 2, as.numeric)
var_prior_scale <- get_env("VAR_PRIOR_SCALE", 1.5, as.numeric)
weight_prior_alpha <- get_env("WEIGHT_PRIOR_ALPHA", 1, as.numeric)
refine_mu_prior_mean <- get_env("REFINE_MU_PRIOR_MEAN", mu_prior_mean, as.numeric)
refine_mu_prior_kappa <- get_env("REFINE_MU_PRIOR_KAPPA", mu_prior_kappa, as.numeric)
refine_var_prior_shape <- get_env("REFINE_VAR_PRIOR_SHAPE", var_prior_shape, as.numeric)
refine_var_prior_scale <- get_env("REFINE_VAR_PRIOR_SCALE", var_prior_scale, as.numeric)
refine_weight_prior_alpha <- get_env("REFINE_WEIGHT_PRIOR_ALPHA", weight_prior_alpha, as.numeric)
min_mixture_var <- get_env("MIN_MIXTURE_VAR", 1e-3, as.numeric)
em_svd_iter <- get_env("EM_SVD_ITER", pretrain_aug_iter, as.integer)
em_svd_tol_loglik <- get_env("EM_SVD_TOL_LOGLIK", 0, as.numeric)
em_svd_tol_L <- get_env("EM_SVD_TOL_L", NA_real_, as.numeric)
if (is.na(em_svd_tol_L)) em_svd_tol_L <- NULL
em_svd_init_method <- get_env("EM_SVD_INIT", "intercept_only", as.character)
if (!em_svd_init_method %in% c("intercept_only", "viroli_svd", "both")) {
  stop("EM_SVD_INIT must be one of 'intercept_only', 'viroli_svd', or 'both'.")
}
em_svd_init_z <- get_env("EM_SVD_INIT_Z", "sample", as.character)
if (!em_svd_init_z %in% c("sample", "expectation")) {
  stop("EM_SVD_INIT_Z must be either 'sample' or 'expectation'.")
}
em_svd_random_starts <- get_env("EM_SVD_RANDOM_STARTS", 0L, as.integer)
em_svd_random_start_scale <- get_env("EM_SVD_RANDOM_START_SCALE", 0.05, as.numeric)
rotation_random_starts <- get_env("ROTATION_RANDOM_STARTS", 1L, as.integer)
rotation_max_outer <- get_env("ROTATION_ITER", pretrain_aug_iter, as.integer)
rotation_n_mix_starts <- get_env("ROTATION_N_MIX_STARTS", 3L, as.integer)
rotation_grid_size <- get_env("ROTATION_GRID_SIZE", 21L, as.integer)
rotation_sweep <- get_env("ROTATION_SWEEP", "full", as.character)
rotation_loading_l1_penalty <- get_env("ROTATION_LOADING_L1_PENALTY", 0, as.numeric)
rotation_objective_tolerance <- get_env("ROTATION_OBJECTIVE_TOLERANCE", 1e-4, as.numeric)
rotation_min_outer <- get_env("ROTATION_MIN_OUTER", 2L, as.integer)
rotation_require_mixture_convergence <- get_env("ROTATION_REQUIRE_MIXTURE_CONVERGENCE", TRUE, as.logical)
refine_require_mixture_convergence <- get_env("REFINE_REQUIRE_MIXTURE_CONVERGENCE", FALSE, as.logical)
refine_enforce_monotone <- get_env("REFINE_ENFORCE_MONOTONE", TRUE, as.logical)
refine_monotone_tolerance <- get_env("REFINE_MONOTONE_TOLERANCE", 1e-8, as.numeric)

mfa_iter <- get_env("MFA_ITER", 2000L, as.integer)
mfa_burn <- get_env("MFA_BURN", 1000L, as.integer)
mfa_thin <- get_env("MFA_THIN", 1L, as.integer)
mfa_tau_lambda <- get_env("MFA_TAU_LAMBDA", 1.5, as.numeric)
mfa_tau_intercept <- get_env("MFA_TAU_INTERCEPT", 5, as.numeric)
mfa_alpha_dirichlet_override <- get_env("MFA_ALPHA_DIRICHLET", NA_real_, as.numeric)
mfa_stop_on_stability <- get_env("MFA_STOP_ON_STABILITY", FALSE, as.logical)
mfa_min_iter <- get_env("MFA_MIN_ITER", max(50L, floor(mfa_iter / 3)), as.integer)
mfa_check_window <- get_env("MFA_CHECK_WINDOW", max(25L, floor(mfa_iter / 10)), as.integer)
mfa_lambda_stability_tol <- get_env("MFA_LAMBDA_STABILITY_TOL", 0.01, as.numeric)
mfa_intercept_stability_tol <- get_env("MFA_INTERCEPT_STABILITY_TOL", 0.01, as.numeric)
mfa_occupied_stability_tol <- get_env("MFA_OCCUPIED_STABILITY_TOL", 1L, as.integer)
mfa_normalize_scale <- get_env("MFA_NORMALIZE_SCALE", TRUE, as.logical)
mfa_min_scale <- get_env("MFA_MIN_SCALE", 1e-4, as.numeric)
mfa_verbose <- get_env("MFA_VERBOSE", TRUE, as.logical)
mfa_compute_parameter_ess <- get_env("MFA_COMPUTE_PARAMETER_ESS", TRUE, as.logical)
viroli_iter <- get_env("VIROLI_ITER", mfa_iter, as.integer)
viroli_burn <- get_env("VIROLI_BURN", mfa_burn, as.integer)
viroli_thin <- get_env("VIROLI_THIN", mfa_thin, as.integer)
viroli_tau_lambda <- get_env("VIROLI_TAU_LAMBDA", mfa_tau_lambda, as.numeric)
viroli_tau_intercept <- get_env("VIROLI_TAU_INTERCEPT", mfa_tau_intercept, as.numeric)
viroli_lambda_l1_penalty <- get_env("VIROLI_LAMBDA_L1_PENALTY", 0, as.numeric)
viroli_method_name <- get_env("VIROLI_METHOD_NAME", "viroli_independent_factor_gibbs", as.character)
viroli_alpha_dirichlet <- get_env("VIROLI_ALPHA_DIRICHLET", 1, as.numeric)
viroli_normalize_each_draw <- get_env("VIROLI_NORMALIZE_EACH_DRAW", TRUE, as.logical)
viroli_min_scale <- get_env("VIROLI_MIN_SCALE", mfa_min_scale, as.numeric)
viroli_compute_parameter_ess <- get_env("VIROLI_COMPUTE_PARAMETER_ESS", mfa_compute_parameter_ess, as.logical)
viroli_verbose <- get_env("VIROLI_VERBOSE", mfa_verbose, as.logical)
viroli_seed_override <- get_env("VIROLI_SEED", NA_real_, as.numeric)
write_iteration_histories <- get_env("WRITE_ITERATION_HISTORIES", FALSE, as.logical)
write_parameter_tables <- get_env("WRITE_PARAMETER_TABLES", TRUE, as.logical)
write_parameter_table_max_K <- get_env("WRITE_PARAMETER_TABLE_MAX_K", 5000L, as.integer)
parallel_ours <- get_env("PARALLEL_OURS", FALSE, as.logical)
parallel_gibbs <- get_env("PARALLEL_GIBBS", FALSE, as.logical)
parallel_workers <- get_env("PARALLEL_WORKERS", NULL, as.integer)
max_joint_alignment_K <- get_env("MAX_JOINT_ALIGNMENT_K", 5000L, as.integer)
max_joint_parameter_K <- get_env("MAX_JOINT_PARAMETER_K", 5000L, as.integer)
max_joint_profile_ari_K <- get_env("MAX_JOINT_PROFILE_ARI_K", 5000L, as.integer)
success_factor_corr <- get_env("SUCCESS_FACTOR_CORR", 0.90, as.numeric)
success_parameter_corr <- get_env("SUCCESS_PARAMETER_CORR", 0.90, as.numeric)

make_equal_mixture_params <- function(H, G, sep, variance_mode = mixture_variance_mode) {
  variance_mode <- match.arg(variance_mode, c("unequal", "equal"))
  G <- normalize_G_counts(G, H)
  lapply(seq_len(H), function(h) {
    Gh <- G[h]
    if (Gh == 1L) {
      list(pi = 1, mu = 0, sd = 1)
    } else if (Gh == 2L) {
      sd <- if (variance_mode == "unequal") c(0.25, 0.65) else c(0.35, 0.35)
      list(pi = c(0.5, 0.5), mu = c(-sep, sep), sd = sd)
    } else if (Gh == 3L) {
      sd <- if (variance_mode == "unequal") g3_raw_sd else c(0.35, 0.35, 0.35)
      list(pi = c(0.30, 0.40, 0.30), mu = c(-sep, 0, sep), sd = sd)
    } else {
      stop("This comparison currently supports G_h in {1, 2, 3}.")
    }
  })
}

make_random_mixture_params <- function(H, G, sep, seed, variance_mode = mixture_variance_mode) {
  variance_mode <- match.arg(variance_mode, c("unequal", "equal"))
  set.seed(seed)
  G <- normalize_G_counts(G, H)

  lapply(seq_len(H), function(h) {
    Gh <- G[h]
    if (Gh == 1L) return(list(pi = 1, mu = 0, sd = 1))
    if (!Gh %in% c(2L, 3L)) {
      stop("This comparison currently supports random mixtures only for G_h in {1, 2, 3}.")
    }
    pi_raw <- rgamma(Gh, shape = runif(Gh, 1.5, 5.0), rate = 1)
    pi <- pi_raw / sum(pi_raw)
    base_mu <- if (Gh == 2L) c(-sep, sep) else c(-sep, 0, sep)
    jitter <- runif(Gh, -0.35 * sep, 0.35 * sep)
    mu <- sort(base_mu + jitter)
    if (Gh > 1L) {
      min_gap <- 0.45 * sep
      for (g in 2:Gh) {
        if (mu[g] - mu[g - 1L] < min_gap) mu[g] <- mu[g - 1L] + min_gap
      }
      mu <- mu - mean(mu)
    }
    sd <- if (variance_mode == "unequal") runif(Gh, 0.25, 0.70) else rep(0.35, Gh)
    list(pi = pi, mu = mu, sd = sd)
  })
}

make_viroli_smoke_mixture_params <- function(H, G, sep) {
  G <- normalize_G_counts(G, H)
  lapply(seq_len(H), function(h) {
    Gh <- G[h]
    if (Gh == 1L) {
      list(pi = 1, mu = 0, sd = 1)
    } else if (Gh == 2L) {
      list(pi = c(0.50, 0.50), mu = sep * c(-1, 1), sd = c(0.55, 0.85))
    } else if (Gh == 3L) {
      list(pi = c(0.30, 0.40, 0.30), mu = sep * c(-1.35, 0, 1.35), sd = c(0.45, 0.65, 0.45))
    } else {
      stop("The Viroli smoke DGP currently supports G_h in {1, 2, 3}.")
    }
  })
}

make_mixture_params <- function(H, G, sep, seed, mode = mixture_param_mode) {
  mode <- match.arg(mode, c("equal", "random", "viroli_smoke"))
  if (mode == "equal") {
    make_equal_mixture_params(H, G, sep, variance_mode = mixture_variance_mode)
  } else if (mode == "viroli_smoke") {
    make_viroli_smoke_mixture_params(H, G, sep)
  } else {
    make_random_mixture_params(H, G, sep, seed, variance_mode = mixture_variance_mode)
  }
}

make_item_intercepts <- function(p, H, block_id, seed, mode = intercept_mode) {
  make_sample_size_item_intercepts(
    p = p,
    H = H,
    block_id = block_id,
    seed = seed,
    mode = mode,
    intercept_sd = intercept_sd,
    intercept_block_span = intercept_block_span,
    intercept_clip = intercept_clip
  )
}

make_block_sizes <- function(p, H, mode = block_size_mode) {
  make_sample_size_block_sizes(p = p, H = H, mode = mode)
}

make_dgp_loadings <- function(design_name, p, H) {
  make_sample_size_loadings(
    design = design_name,
    p = p,
    H = H,
    block_size_mode = block_size_mode,
    loading_sign_mode = loading_sign_mode
  )
}

simulate_original_binary_probit <- function(n, p, H, G, sep, loading_design, seed) {
  set.seed(seed)
  mixture_params <- make_mixture_params(H, G, sep, seed = seed + 7919L)
  loading_out <- make_dgp_loadings(loading_design, p, H)
  alpha <- make_item_intercepts(
    p = p,
    H = H,
    block_id = loading_out$block_id,
    seed = seed + 3571L
  )
  F <- matrix(NA_real_, n, H)
  component <- matrix(NA_integer_, n, H)
  standardized_params <- vector("list", H)

  for (h in seq_len(H)) {
    draw_h <- sample_standardized_mixture(n, mixture_params[[h]])
    F[, h] <- draw_h$x
    component[, h] <- draw_h$component
    standardized_params[[h]] <- draw_h$parameters
  }

  Z_latent <- sweep(F %*% t(loading_out$Lambda), 2L, alpha, "+") +
    matrix(rnorm(n * p), n, p)
  X_binary <- 1L * (Z_latent > 0)
  list(
    X_binary = X_binary,
    Z_latent = Z_latent,
    F = F,
    Lambda = loading_out$Lambda,
    alpha = alpha,
    component = component,
    mixture_params = standardized_params,
    block_id = loading_out$block_id,
    block_sizes = loading_out$block_sizes,
    prevalence = colMeans(X_binary)
  )
}

align_factors <- function(F_true, F_est) {
  C <- cor(F_true, F_est)
  H_true <- ncol(F_true)
  H_est <- ncol(F_est)
  remaining_true <- seq_len(H_true)
  remaining_est <- seq_len(H_est)
  n_match <- min(H_true, H_est)
  est_index <- rep(NA_integer_, H_true)

  if (n_match > 0L) {
    for (step in seq_len(n_match)) {
      sub <- abs(C[remaining_true, remaining_est, drop = FALSE])
      pick <- which(sub == max(sub, na.rm = TRUE), arr.ind = TRUE)[1L, ]
      true_h <- remaining_true[pick[1L]]
      est_h <- remaining_est[pick[2L]]
      est_index[true_h] <- est_h
      remaining_true <- setdiff(remaining_true, true_h)
      remaining_est <- setdiff(remaining_est, est_h)
    }
  }

  signs <- rep(1, H_true)
  matched_abs_cor <- rep(0, H_true)
  matched <- which(!is.na(est_index))
  if (length(matched) > 0L) {
    signs[matched] <- sign(C[cbind(matched, est_index[matched])])
    signs[signs == 0] <- 1
    matched_abs_cor[matched] <- abs(C[cbind(matched, est_index[matched])])
  }

  list(
    est_index = est_index,
    unmatched_est_index = remaining_est,
    signs = signs,
    matched_abs_cor = matched_abs_cor,
    mean_abs_cor = mean(matched_abs_cor),
    n_matched = n_match,
    H_true = H_true,
    H_est = H_est
  )
}

align_loadings_to_truth <- function(Lambda_true, Lambda_est_raw, F_true = NULL, F_est = NULL) {
  H_true <- ncol(Lambda_true)
  H_est <- ncol(Lambda_est_raw)
  if (H_true > H_est) stop("Cannot align fewer estimated loading columns than true loading columns.")
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' is required for globally optimal loading-column alignment.")
  }

  dist_mat <- matrix(NA_real_, H_true, H_est)
  sign_mat <- matrix(1, H_true, H_est)
  for (h_true in seq_len(H_true)) {
    lambda_true_h <- Lambda_true[, h_true]
    for (h_est in seq_len(H_est)) {
      lambda_est_h <- Lambda_est_raw[, h_est]
      pos_dist <- sum((lambda_true_h - lambda_est_h)^2)
      neg_dist <- sum((lambda_true_h + lambda_est_h)^2)
      if (pos_dist <= neg_dist) {
        dist_mat[h_true, h_est] <- pos_dist
        sign_mat[h_true, h_est] <- 1
      } else {
        dist_mat[h_true, h_est] <- neg_dist
        sign_mat[h_true, h_est] <- -1
      }
    }
  }

  if (H_true <= H_est) {
    est_index <- as.integer(clue::solve_LSAP(dist_mat))
  } else {
    true_for_est <- as.integer(clue::solve_LSAP(t(dist_mat)))
    est_index <- rep(NA_integer_, H_true)
    est_index[true_for_est] <- seq_len(H_est)
  }
  signs <- sign_mat[cbind(seq_len(H_true), est_index)]

  loading_cor <- rep(NA_real_, H_true)
  for (h in seq_len(H_true)) {
    loading_cor[h] <- suppressWarnings(abs(cor(Lambda_true[, h], Lambda_est_raw[, est_index[h]])))
  }

  matched_abs_cor <- rep(NA_real_, H_true)
  if (!is.null(F_true) && !is.null(F_est)) {
    C <- cor(F_true, F_est)
    matched_abs_cor <- abs(C[cbind(seq_len(H_true), est_index)])
  }

  list(
    est_index = est_index,
    unmatched_est_index = setdiff(seq_len(H_est), est_index),
    signs = signs,
    matched_abs_cor = matched_abs_cor,
    mean_abs_cor = mean(matched_abs_cor, na.rm = TRUE),
    matched_loading_abs_cor = loading_cor,
    mean_loading_abs_cor = mean(loading_cor, na.rm = TRUE),
    n_matched = length(est_index),
    H_true = H_true,
    H_est = H_est
  )
}

align_lambda_to_truth <- function(Lambda_true, Lambda_est_raw, align) {
  p <- nrow(Lambda_true)
  H_true <- ncol(Lambda_true)
  out <- matrix(0, p, H_true)
  matched <- which(!is.na(align$est_index))
  if (length(matched) > 0L) {
    out[, matched] <- sweep(
      Lambda_est_raw[, align$est_index[matched], drop = FALSE],
      2L,
      align$signs[matched],
      "*"
    )
  }
  dimnames(out) <- dimnames(Lambda_true)
  out
}

choose_alignment <- function(Lambda_true, Lambda_est, F_true, F_est, mode = alignment_mode) {
  # One alignment helper keeps the reported factor, loading, and mixture
  # recovery metrics comparable.  The default aligns by factor-score
  # correlation; smoke-test reproduction can instead align by loading columns.
  if (identical(mode, "loadings")) {
    align_loadings_to_truth(Lambda_true, Lambda_est, F_true = F_true, F_est = F_est)
  } else {
    align_factors(F_true, F_est)
  }
}

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  n_tab <- sum(tab)
  choose2 <- function(x) x * (x - 1) / 2
  sum_comb <- sum(choose2(tab))
  row_comb <- sum(choose2(rowSums(tab)))
  col_comb <- sum(choose2(colSums(tab)))
  total_comb <- choose2(n_tab)
  expected <- row_comb * col_comb / total_comb
  max_index <- 0.5 * (row_comb + col_comb)
  if (abs(max_index - expected) < 1e-12) return(0)
  (sum_comb - expected) / (max_index - expected)
}

rbind_fill <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (!length(x)) return(data.frame())
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) d[[nm]] <- NA
    d[, all_names, drop = FALSE]
  })
  do.call(rbind, x)
}

maybe_write_parameter_recovery <- function(parameter_table, out_file, K_joint) {
  if (isTRUE(write_parameter_tables) && K_joint <= write_parameter_table_max_K) {
    write.csv(parameter_table, out_file, row.names = FALSE)
  }
}

joint_class_index <- function(component, G) {
  component <- as.matrix(component)
  H <- ncol(component)
  G <- normalize_G_counts(G, H)
  multipliers <- c(1L, cumprod(G[-length(G)]))
  as.integer(1L + (component - 1L) %*% multipliers)
}

joint_profile_grid <- function(H, G) {
  G <- normalize_G_counts(G, H)
  grid <- expand.grid(lapply(G, seq_len), KEEP.OUT.ATTRS = FALSE)
  names(grid) <- paste0("factor_", seq_len(H), "_component")
  grid
}

true_joint_mixture_parameters <- function(mixture_params) {
  H <- length(mixture_params)
  G <- vapply(mixture_params, function(z) length(z$pi), integer(1L))
  grid <- joint_profile_grid(H, G)
  K <- nrow(grid)
  weight <- rep(1, K)
  mu <- matrix(NA_real_, K, H)
  var <- matrix(NA_real_, K, H)
  for (h in seq_len(H)) {
    comp <- grid[[h]]
    weight <- weight * mixture_params[[h]]$pi[comp]
    mu[, h] <- mixture_params[[h]]$mu[comp]
    var[, h] <- mixture_params[[h]]$sd[comp]^2
  }
  out <- cbind(joint_class = seq_len(K), grid, true_weight = weight)
  for (h in seq_len(H)) out[[paste0("true_mu_", h)]] <- mu[, h]
  for (h in seq_len(H)) out[[paste0("true_var_", h)]] <- var[, h]
  out
}

product_joint_mixture_parameters_from_fits <- function(mixture_fits, G = NULL) {
  H <- length(mixture_fits)
  if (is.null(G)) {
    G <- vapply(mixture_fits, function(z) length(z$pi), integer(1L))
  }
  G <- normalize_G_counts(G, H)
  grid <- joint_profile_grid(H, G)
  K <- nrow(grid)
  weight <- rep(1, K)
  mu <- matrix(NA_real_, K, H)
  var <- matrix(NA_real_, K, H)
  for (h in seq_len(H)) {
    comp <- grid[[h]]
    fit_h <- mixture_fits[[h]]
    fit_var_h <- if (!is.null(fit_h$var)) fit_h$var else fit_h$sd^2
    weight <- weight * fit_h$pi[comp]
    mu[, h] <- fit_h$mu[comp]
    var[, h] <- fit_var_h[comp]
  }
  list(pi = weight / sum(weight), mu = mu, sig2 = var)
}

sort_mixture_fit_by_mean <- function(fit) {
  ord <- order(fit$mu)
  fit$pi <- fit$pi[ord]
  fit$mu <- fit$mu[ord]
  if (!is.null(fit$var)) fit$var <- fit$var[ord]
  if (!is.null(fit$sd)) fit$sd <- fit$sd[ord]
  fit
}

align_product_mixture_fits <- function(mixture_fits, align) {
  H <- align$H_true
  out <- vector("list", H)
  for (h in seq_len(H)) {
    est_h <- align$est_index[h]
    if (is.na(est_h)) {
      out[[h]] <- list(pi = 1, mu = 0, sd = 1)
    } else {
      fit_h <- mixture_fits[[est_h]]
      sign_h <- align$signs[h]
      out[[h]] <- list(
        pi = fit_h$pi,
        mu = sign_h * fit_h$mu,
        var = if (!is.null(fit_h$var)) fit_h$var else fit_h$sd^2
      )
    }
  }
  out
}

marginal_mixture_parameter_summary <- function(true_mixture_params, mixture_fits, G) {
  H <- length(true_mixture_params)
  G <- normalize_G_counts(G, H)
  rows <- vector("list", sum(G))
  idx <- 0L
  for (h in seq_len(H)) {
    true_h <- true_mixture_params[[h]]
    fit_h <- sort_mixture_fit_by_mean(mixture_fits[[h]])
    fit_var <- if (!is.null(fit_h$var)) fit_h$var else fit_h$sd^2
    for (g in seq_len(G[h])) {
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        factor = h,
        component = g,
        true_weight = true_h$pi[g],
        est_weight = fit_h$pi[g],
        true_mu = true_h$mu[g],
        est_mu = fit_h$mu[g],
        true_var = true_h$sd[g]^2,
        est_var = fit_var[g],
        stringsAsFactors = FALSE
      )
    }
  }
  table <- do.call(rbind, rows)
  table$mu_error <- table$est_mu - table$true_mu
  table$var_error <- table$est_var - table$true_var
  table$weight_error <- table$est_weight - table$true_weight
  summary <- data.frame(
    marginal_mu_rmse = sqrt(mean(table$mu_error^2)),
    marginal_var_rmse = sqrt(mean(table$var_error^2)),
    marginal_weight_rmse = sqrt(mean(table$weight_error^2)),
    marginal_mu_corr = safe_cor(table$true_mu, table$est_mu),
    marginal_log_var_corr = safe_cor(log(pmax(table$true_var, 1e-8)), log(pmax(table$est_var, 1e-8))),
    marginal_weight_corr = safe_cor(table$true_weight, table$est_weight),
    stringsAsFactors = FALSE
  )
  list(table = table, summary = summary)
}

empty_joint_parameter_summary <- function(reason = "not_computed") {
  data.frame(
    joint_mu_rmse = NA_real_,
    joint_var_rmse = NA_real_,
    joint_mu_corr = NA_real_,
    joint_var_corr = NA_real_,
    joint_weight_corr = NA_real_,
    joint_weight_rmse = NA_real_,
    joint_weight_l1 = NA_real_,
    joint_weight_max_abs_error = NA_real_,
    est_effective_classes_001 = NA_real_,
    est_effective_classes_01 = NA_real_,
    joint_alignment_method = reason,
    stringsAsFactors = FALSE
  )
}

dgp_diagnostics <- function(sim, G) {
  H <- ncol(sim$F)
  G <- normalize_G_counts(G, H)
  prevalence <- sim$prevalence
  component_counts <- unlist(lapply(seq_len(H), function(h) {
    tabulate(sim$component[, h], nbins = G[h])
  }), use.names = FALSE)
  unique_profiles <- unique(as.data.frame(sim$component))
  data.frame(
    dgp_n_items_all_zero = sum(prevalence <= 0),
    dgp_n_items_all_one = sum(prevalence >= 1),
    dgp_n_items_prev_lt_01 = sum(prevalence < 0.01),
    dgp_n_items_prev_gt_99 = sum(prevalence > 0.99),
    dgp_min_item_prevalence = min(prevalence),
    dgp_q05_item_prevalence = unname(quantile(prevalence, 0.05, names = FALSE)),
    dgp_median_item_prevalence = median(prevalence),
    dgp_q95_item_prevalence = unname(quantile(prevalence, 0.95, names = FALSE)),
    dgp_max_item_prevalence = max(prevalence),
    dgp_min_marginal_component_count = min(component_counts),
    dgp_median_marginal_component_count = median(component_counts),
    dgp_n_empty_marginal_components = sum(component_counts == 0L),
    dgp_observed_joint_profiles = nrow(unique_profiles),
    dgp_possible_joint_profiles = prod(as.numeric(G)),
    stringsAsFactors = FALSE
  )
}

align_joint_mfa_axis_parameters <- function(est_mu, est_sig2, align) {
  H <- align$H_true
  K <- nrow(est_mu)
  mu_aligned <- matrix(0, K, H)
  sig2_aligned <- matrix(1, K, H)
  for (h in seq_len(H)) {
    est_h <- align$est_index[h]
    if (!is.na(est_h)) {
      mu_aligned[, h] <- align$signs[h] * est_mu[, est_h]
      sig2_aligned[, h] <- est_sig2[, est_h]
    }
  }
  list(mu = mu_aligned, sig2 = sig2_aligned)
}

aligned_joint_mixture_parameter_summary <- function(true_params, est_pi, est_mu, est_sig2,
                                                    est_class = NULL,
                                                    alignment_method = "provided") {
  H <- ncol(est_mu)
  true_mu <- as.matrix(true_params[, paste0("true_mu_", seq_len(H)), drop = FALSE])
  true_var <- as.matrix(true_params[, paste0("true_var_", seq_len(H)), drop = FALSE])
  true_weight <- true_params$true_weight

  out <- true_params
  out$est_class <- if (is.null(est_class)) seq_len(nrow(est_mu)) else est_class
  out$est_weight <- est_pi
  for (h in seq_len(H)) out[[paste0("est_mu_", h)]] <- est_mu[, h]
  for (h in seq_len(H)) out[[paste0("est_var_", h)]] <- est_sig2[, h]
  out$component_mu_rmse <- sqrt(rowMeans((true_mu - est_mu)^2))
  out$component_var_rmse <- sqrt(rowMeans((true_var - est_sig2)^2))
  out$component_weight_abs_error <- abs(true_weight - est_pi)
  out$alignment_method <- alignment_method

  summary <- data.frame(
    joint_mu_rmse = sqrt(mean((true_mu - est_mu)^2)),
    joint_var_rmse = sqrt(mean((true_var - est_sig2)^2)),
    joint_mu_corr = suppressWarnings(cor(as.vector(true_mu), as.vector(est_mu))),
    joint_var_corr = suppressWarnings(cor(as.vector(true_var), as.vector(est_sig2))),
    joint_weight_corr = suppressWarnings(cor(true_weight, est_pi)),
    joint_weight_rmse = sqrt(mean((true_weight - est_pi)^2)),
    joint_weight_l1 = sum(abs(true_weight - est_pi)),
    joint_weight_max_abs_error = max(abs(true_weight - est_pi)),
    est_effective_classes_001 = sum(est_pi > 0.001),
    est_effective_classes_01 = sum(est_pi > 0.01),
    joint_alignment_method = alignment_method,
    stringsAsFactors = FALSE
  )

  list(table = out, summary = summary)
}

align_product_mixture_parameters_by_axis <- function(true_params, mixture_fits, G = NULL) {
  # Product-mixture fits already have independent coordinate-wise components.
  # After factor sign/permutation alignment, only within-coordinate component
  # labels remain.  Sorting each marginal by its mean is O(sum_h G_h), whereas
  # solving a dense assignment over the full G^H joint grid is O((G^H)^2) memory.
  mixture_fits <- lapply(mixture_fits, sort_mixture_fit_by_mean)
  est <- product_joint_mixture_parameters_from_fits(mixture_fits, G = G)
  aligned_joint_mixture_parameter_summary(
    true_params = true_params,
    est_pi = est$pi,
    est_mu = est$mu,
    est_sig2 = est$sig2,
    alignment_method = "axis_sorted_product"
  )
}

align_joint_mixture_parameters <- function(true_params, est_pi, est_mu, est_sig2) {
  H <- ncol(est_mu)
  K <- nrow(est_mu)
  max_K <- get0("max_joint_alignment_K", ifnotfound = 5000L)
  if (K > max_K) {
    stop(
      "Dense joint-mixture alignment requested with K = ", K,
      ", which would allocate a K-by-K distance matrix. ",
      "Use align_product_mixture_parameters_by_axis() for product-mixture fits ",
      "or raise MAX_JOINT_ALIGNMENT_K only for small diagnostic runs.",
      call. = FALSE
    )
  }
  true_mu <- as.matrix(true_params[, paste0("true_mu_", seq_len(H)), drop = FALSE])
  true_var <- as.matrix(true_params[, paste0("true_var_", seq_len(H)), drop = FALSE])

  scale_mu <- mean(apply(true_mu, 2L, sd), na.rm = TRUE)
  if (!is.finite(scale_mu) || scale_mu <= 1e-8) scale_mu <- 1
  log_true_var <- log(pmax(true_var, 1e-8))
  log_est_var <- log(pmax(est_sig2, 1e-8))

  dist_mat <- matrix(NA_real_, K, K)
  for (k_true in seq_len(K)) {
    for (k_est in seq_len(K)) {
      pieces <- c(true_mu[k_true, ], est_mu[k_est, ], log_true_var[k_true, ], log_est_var[k_est, ])
      if (any(!is.finite(pieces))) {
        dist_mat[k_true, k_est] <- 1e12
      } else {
        mu_dist <- mean(((true_mu[k_true, ] - est_mu[k_est, ]) / scale_mu)^2)
        var_dist <- mean((log_true_var[k_true, ] - log_est_var[k_est, ])^2)
        dist_mat[k_true, k_est] <- mu_dist + 0.25 * var_dist
      }
    }
  }
  dist_mat[!is.finite(dist_mat)] <- 1e12

  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' is required for globally optimal mixture-label alignment.")
  }
  est_index <- as.integer(clue::solve_LSAP(dist_mat))

  est_mu_aligned <- est_mu[est_index, , drop = FALSE]
  est_var_aligned <- est_sig2[est_index, , drop = FALSE]
  est_weight_aligned <- est_pi[est_index]

  aligned_joint_mixture_parameter_summary(
    true_params = true_params,
    est_pi = est_weight_aligned,
    est_mu = est_mu_aligned,
    est_sig2 = est_var_aligned,
    est_class = est_index,
    alignment_method = "dense_joint_lsap"
  )
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L || sd(x) <= 1e-12 || sd(y) <= 1e-12) return(NA_real_)
  suppressWarnings(cor(x, y))
}

safe_group_key <- function(data, group_cols) {
  group_frame <- data[, group_cols, drop = FALSE]
  group_frame[] <- lapply(group_frame, function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- "<NA>"
    x
  })
  interaction(group_frame, drop = TRUE, sep = " | ")
}

parameter_block_from_table <- function(parameter_table, true_prefix, est_prefix, transform = identity) {
  true_cols <- grep(paste0("^", true_prefix), names(parameter_table), value = TRUE)
  if (!length(true_cols)) return(NULL)
  suffix <- sub(paste0("^", true_prefix), "", true_cols)
  est_cols <- paste0(est_prefix, suffix)
  if (!all(est_cols %in% names(parameter_table))) return(NULL)
  list(
    truth = as.numeric(transform(as.matrix(parameter_table[, true_cols, drop = FALSE]))),
    estimate = as.numeric(transform(as.matrix(parameter_table[, est_cols, drop = FALSE])))
  )
}

combined_parameter_correlations <- function(Lambda_true, Lambda_est_aligned, parameter_table,
                                            alpha_true = NULL, alpha_est = NULL) {
  alpha_block <- if (!is.null(alpha_true) && !is.null(alpha_est)) {
    list(truth = as.numeric(alpha_true), estimate = as.numeric(alpha_est))
  } else {
    NULL
  }
  flat_blocks <- list(
    alpha = alpha_block,
    loading = list(truth = as.numeric(Lambda_true), estimate = as.numeric(Lambda_est_aligned)),
    joint_mu = parameter_block_from_table(parameter_table, "true_mu_", "est_mu_"),
    joint_var = parameter_block_from_table(parameter_table, "true_var_", "est_var_"),
    joint_weight = list(truth = parameter_table$true_weight, estimate = parameter_table$est_weight)
  )
  flat_blocks <- flat_blocks[!vapply(flat_blocks, is.null, logical(1L))]
  flat_truth <- unlist(lapply(flat_blocks, `[[`, "truth"), use.names = FALSE)
  flat_estimate <- unlist(lapply(flat_blocks, `[[`, "estimate"), use.names = FALSE)

  blocks <- list(
    alpha = alpha_block,
    loading = list(truth = as.numeric(Lambda_true), estimate = as.numeric(Lambda_est_aligned)),
    joint_mu = parameter_block_from_table(parameter_table, "true_mu_", "est_mu_"),
    joint_log_var = parameter_block_from_table(
      parameter_table,
      "true_var_",
      "est_var_",
      transform = function(x) log(pmax(x, 1e-8))
    ),
    joint_weight = list(truth = parameter_table$true_weight, estimate = parameter_table$est_weight)
  )
  blocks <- blocks[!vapply(blocks, is.null, logical(1L))]
  valid_blocks <- blocks[vapply(blocks, function(block) {
    is.finite(safe_cor(block$truth, block$estimate))
  }, logical(1L))]
  if (!length(valid_blocks)) {
    return(data.frame(
      flat_parameter_corr = safe_cor(flat_truth, flat_estimate),
      flat_parameter_blocks = paste(names(flat_blocks), collapse = "+"),
      all_parameter_corr = NA_real_,
      all_parameter_corr_raw = NA_real_,
      all_parameter_blocks = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  standardized <- lapply(valid_blocks, function(block) {
    truth_sd <- sd(block$truth, na.rm = TRUE)
    list(
      truth = (block$truth - mean(block$truth, na.rm = TRUE)) / truth_sd,
      estimate = (block$estimate - mean(block$estimate, na.rm = TRUE)) / truth_sd
    )
  })

  data.frame(
    flat_parameter_corr = safe_cor(flat_truth, flat_estimate),
    flat_parameter_blocks = paste(names(flat_blocks), collapse = "+"),
    all_parameter_corr = safe_cor(
      unlist(lapply(standardized, `[[`, "truth"), use.names = FALSE),
      unlist(lapply(standardized, `[[`, "estimate"), use.names = FALSE)
    ),
    all_parameter_corr_raw = safe_cor(
      unlist(lapply(valid_blocks, `[[`, "truth"), use.names = FALSE),
      unlist(lapply(valid_blocks, `[[`, "estimate"), use.names = FALSE)
    ),
    all_parameter_blocks = paste(names(valid_blocks), collapse = "+"),
    stringsAsFactors = FALSE
  )
}

class_map_from_mixtures_local <- function(F_hat, mixture_fits) {
  out <- sapply(seq_len(ncol(F_hat)), function(h) {
    max.col(mixture_responsibilities(F_hat[, h], mixture_fits[[h]]), ties.method = "first")
  })
  if (ncol(F_hat) == 1L) out <- matrix(out, ncol = 1L)
  out
}

sample_dirichlet <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

sample_inv_gamma <- function(shape, rate) {
  1 / rgamma(1L, shape = shape, rate = rate)
}

rmvnorm_chol <- function(mean, cov) {
  mean + as.numeric(t(chol(cov)) %*% rnorm(length(mean)))
}

initialize_joint_mfa_state <- function(X, H, K, seed) {
  set.seed(seed)
  n <- nrow(X)
  Z <- initialize_binary_Z(X, seed = seed)
  svd_out <- svd_scores_from_Z(Z, H = H, center_Z = TRUE)
  F <- scale(svd_out$S)
  if (any(!is.finite(F))) F <- matrix(rnorm(n * H), n, H)

  if (n >= K) {
    km <- try(kmeans(F, centers = K, nstart = 5, iter.max = 30), silent = TRUE)
    C <- if (inherits(km, "try-error")) sample.int(K, n, replace = TRUE) else km$cluster
  } else {
    C <- sample.int(K, n, replace = TRUE)
  }

  Lambda <- update_working_loadings_no_intercept(Z, F, loading_penalty = 0)$Lambda
  alpha <- colMeans(Z - F %*% t(Lambda))
  list(Z = Z, F = as.matrix(F), C = as.integer(C), alpha = as.numeric(alpha), Lambda = as.matrix(Lambda))
}

sample_binary_Z_joint_mfa <- function(X, F, Lambda, alpha) {
  mean_mat <- sweep(F %*% t(Lambda), 2L, alpha, "+")
  Z <- matrix(NA_real_, nrow(X), ncol(X))
  one <- X == 1
  zero <- !one
  Z[one] <- rtruncnorm_binary_vec(mean_mat[one], 1, 0, Inf)
  Z[zero] <- rtruncnorm_binary_vec(mean_mat[zero], 1, -Inf, 0)
  Z
}

sample_alpha_lambda_given_Z_F <- function(
    Z,
    F,
    tau_intercept = 5,
    tau_lambda = 1.5,
    parallel = FALSE,
    workers = NULL) {
  n <- nrow(Z)
  p <- ncol(Z)
  H <- ncol(F)
  W <- cbind(1, F)
  prior_prec <- diag(c(1 / tau_intercept^2, rep(1 / tau_lambda^2, H)), H + 1L)
  V <- solve(crossprod(W) + prior_prec)
  chol_V <- chol(V)
  WtZ <- crossprod(W, Z)

  rows <- parallel_lapply(seq_len(p), function(j) {
    m <- V %*% WtZ[, j]
    beta <- as.numeric(m + t(chol_V) %*% rnorm(H + 1L))
    c(alpha = beta[1L], lambda = beta[-1L])
  }, parallel = parallel, workers = workers)
  draw_mat <- do.call(rbind, rows)
  alpha <- draw_mat[, 1L]
  Lambda <- draw_mat[, -1L, drop = FALSE]

  list(alpha = alpha, Lambda = Lambda)
}

normalize_joint_mfa_scale <- function(F, alpha, Lambda, mu, sig2, pi_vec, min_scale = 1e-4) {
  # Fix the factor scale convention without changing the probit linear
  # predictor.  For each factor h, compute the current marginal mixture mean
  # m_h and sd s_h, set eta_h^*=(eta_h-m_h)/s_h, and transform
  # Lambda_h^*=s_h Lambda_h, alpha^*=alpha + Lambda m.
  H <- ncol(F)
  pi_vec <- pi_vec / sum(pi_vec)
  marginal_mean <- as.numeric(crossprod(pi_vec, mu))
  marginal_second <- as.numeric(crossprod(pi_vec, sig2 + mu^2))
  marginal_var <- pmax(marginal_second - marginal_mean^2, min_scale^2)
  marginal_scale <- sqrt(marginal_var)
  marginal_scale[!is.finite(marginal_scale) | marginal_scale < min_scale] <- min_scale

  alpha_new <- alpha + as.numeric(Lambda %*% marginal_mean)
  Lambda_new <- sweep(Lambda, 2L, marginal_scale, "*")
  F_new <- sweep(sweep(F, 2L, marginal_mean, "-"), 2L, marginal_scale, "/")
  mu_new <- sweep(sweep(mu, 2L, marginal_mean, "-"), 2L, marginal_scale, "/")
  sig2_new <- sweep(sig2, 2L, marginal_scale^2, "/")

  list(
    F = F_new,
    alpha = alpha_new,
    Lambda = Lambda_new,
    mu = mu_new,
    sig2 = sig2_new,
    marginal_mean_before = marginal_mean,
    marginal_scale_before = marginal_scale
  )
}

scalar_effective_sample_size <- function(x, max_lag = NULL) {
  # Simple positive-sequence univariate ESS.  We use this as a scalar mixing
  # diagnostic for Gibbs output after burn-in and thinning.
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3L || !is.finite(var(x)) || var(x) <= 0) return(NA_real_)
  if (is.null(max_lag)) max_lag <- min(n - 1L, 250L)
  ac <- as.numeric(acf(x, lag.max = max_lag, plot = FALSE, demean = TRUE)$acf[-1L])
  ac <- ac[is.finite(ac)]
  if (!length(ac)) return(n)
  first_nonpositive <- which(ac <= 0)[1L]
  use_lags <- if (is.na(first_nonpositive)) seq_along(ac) else seq_len(first_nonpositive - 1L)
  if (!length(use_lags)) return(n)
  tau <- 1 + 2 * sum(ac[use_lags])
  if (!is.finite(tau) || tau <= 0) return(NA_real_)
  min(n, n / tau)
}

compute_parameter_ess_table <- function(trace_list) {
  make_table <- function(group, mat) {
    if (is.null(mat) || !length(mat)) return(NULL)
    mat <- as.matrix(mat)
    data.frame(
      parameter_group = group,
      parameter = colnames(mat),
      ess = apply(mat, 2L, scalar_effective_sample_size),
      n_draws = nrow(mat),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, Filter(Negate(is.null), list(
    make_table("alpha", trace_list$alpha),
    make_table("lambda", trace_list$lambda),
    make_table("pi", trace_list$pi),
    make_table("mu", trace_list$mu),
    make_table("sig2", trace_list$sig2)
  )))
  rownames(out) <- NULL
  out
}

summarize_parameter_ess <- function(ess_table) {
  if (is.null(ess_table) || !nrow(ess_table)) {
    return(data.frame(
      gibbs_min_parameter_ess = NA_real_,
      gibbs_median_parameter_ess = NA_real_,
      gibbs_mean_parameter_ess = NA_real_,
      gibbs_min_alpha_ess = NA_real_,
      gibbs_median_alpha_ess = NA_real_,
      gibbs_min_lambda_ess = NA_real_,
      gibbs_median_lambda_ess = NA_real_,
      gibbs_min_pi_ess = NA_real_,
      gibbs_median_pi_ess = NA_real_,
      gibbs_min_mu_ess = NA_real_,
      gibbs_median_mu_ess = NA_real_,
      gibbs_min_sig2_ess = NA_real_,
      gibbs_median_sig2_ess = NA_real_
    ))
  }
  group_value <- function(group, FUN) {
    z <- ess_table$ess[ess_table$parameter_group == group]
    if (!length(z) || all(!is.finite(z))) return(NA_real_)
    FUN(z, na.rm = TRUE)
  }
  data.frame(
    gibbs_min_parameter_ess = if (all(!is.finite(ess_table$ess))) NA_real_ else min(ess_table$ess, na.rm = TRUE),
    gibbs_median_parameter_ess = median(ess_table$ess, na.rm = TRUE),
    gibbs_mean_parameter_ess = mean(ess_table$ess, na.rm = TRUE),
    gibbs_min_alpha_ess = group_value("alpha", min),
    gibbs_median_alpha_ess = group_value("alpha", median),
    gibbs_min_lambda_ess = group_value("lambda", min),
    gibbs_median_lambda_ess = group_value("lambda", median),
    gibbs_min_pi_ess = group_value("pi", min),
    gibbs_median_pi_ess = group_value("pi", median),
    gibbs_min_mu_ess = group_value("mu", min),
    gibbs_median_mu_ess = group_value("mu", median),
    gibbs_min_sig2_ess = group_value("sig2", min),
    gibbs_median_sig2_ess = group_value("sig2", median)
  )
}

summarize_gibbs_timing <- function(history, burn, thin) {
  if (is.null(history) || !nrow(history)) {
    return(data.frame())
  }
  keep <- history$iteration > burn & ((history$iteration - burn) %% thin == 0L)
  kept <- history[keep, , drop = FALSE]
  if (!nrow(kept)) kept <- history
  median_or_na <- function(x) {
    if (is.null(x) || !length(x) || all(!is.finite(x))) NA_real_ else median(x, na.rm = TRUE)
  }
  data.frame(
    gibbs_median_iteration_seconds = median_or_na(kept$iteration_seconds),
    gibbs_median_z_sample_seconds = median_or_na(kept$z_sample_seconds),
    gibbs_median_factor_sample_seconds = median_or_na(kept$factor_sample_seconds),
    gibbs_median_class_sample_seconds = median_or_na(kept$class_sample_seconds),
    gibbs_median_mixture_weight_seconds = median_or_na(kept$mixture_weight_seconds),
    gibbs_median_mixture_parameter_seconds = median_or_na(kept$mixture_parameter_seconds),
    gibbs_median_regression_seconds = median_or_na(kept$regression_seconds),
    gibbs_median_normalization_seconds = median_or_na(kept$normalization_seconds),
    gibbs_median_keep_draw_seconds = median_or_na(kept$keep_draw_seconds),
    gibbs_median_diagnostics_seconds = median_or_na(kept$diagnostics_seconds),
    stringsAsFactors = FALSE
  )
}

fit_joint_class_probit_mfa_gibbs <- function(
    X,
    H,
    K,
    n_iter = 200L,
    burn = floor(n_iter / 2),
    thin = 5L,
    tau_lambda = 1.5,
    tau_intercept = 5,
    alpha_dirichlet = 1 / K,
    mu0 = 0,
    kappa0 = 0.05,
    a0 = 3,
    b0 = 2,
    min_var = 0.05,
    stop_on_stability = TRUE,
    min_iter = max(50L, floor(n_iter / 3)),
    check_window = max(25L, floor(n_iter / 10)),
    lambda_stability_tol = 0.01,
    intercept_stability_tol = 0.01,
    occupied_stability_tol = 1L,
    normalize_scale = TRUE,
    min_scale = 1e-4,
    parallel = FALSE,
    workers = NULL,
    compute_parameter_ess = TRUE,
    seed = 1L,
    verbose = TRUE) {
  set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  state <- initialize_joint_mfa_state(X, H, K, seed = seed)
  Z <- state$Z
  F <- state$F
  C <- state$C
  alpha <- state$alpha
  Lambda <- state$Lambda

  pi_vec <- rep(1 / K, K)
  mu <- matrix(0, K, H)
  sig2 <- matrix(1, K, H)
  history <- data.frame(
    iteration = seq_len(n_iter),
    occupied_classes = NA_integer_,
    lambda_norm = NA_real_,
    intercept_norm = NA_real_,
    max_abs_mixture_mean_before_normalize = NA_real_,
    max_abs_log_mixture_scale_before_normalize = NA_real_,
    lambda_rel_range = NA_real_,
    intercept_rel_range = NA_real_,
    occupied_range = NA_integer_,
    stable = FALSE,
    z_sample_seconds = NA_real_,
    factor_sample_seconds = NA_real_,
    class_sample_seconds = NA_real_,
    mixture_weight_seconds = NA_real_,
    mixture_parameter_seconds = NA_real_,
    regression_seconds = NA_real_,
    normalization_seconds = NA_real_,
    keep_draw_seconds = NA_real_,
    diagnostics_seconds = NA_real_,
    iteration_seconds = NA_real_
  )

  keep_F <- matrix(0, n, H)
  keep_alpha <- numeric(p)
  keep_alpha_raw <- numeric(p)
  keep_Lambda <- matrix(0, p, H)
  keep_pi <- numeric(K)
  keep_mu <- matrix(0, K, H)
  keep_sig2 <- matrix(0, K, H)
  n_keep <- 0L
  keep_target <- max(0L, floor((n_iter - burn) / max(thin, 1L)))
  if (isTRUE(compute_parameter_ess) && keep_target > 0L) {
    alpha_trace <- matrix(NA_real_, keep_target, p)
    colnames(alpha_trace) <- paste0("alpha_", seq_len(p))
    lambda_trace <- matrix(NA_real_, keep_target, p * H)
    colnames(lambda_trace) <- as.vector(outer(paste0("item", seq_len(p)), paste0("factor", seq_len(H)), paste, sep = "_"))
    pi_trace <- matrix(NA_real_, keep_target, K)
    colnames(pi_trace) <- paste0("pi_", seq_len(K))
    mu_trace <- matrix(NA_real_, keep_target, K * H)
    colnames(mu_trace) <- as.vector(outer(paste0("class", seq_len(K)), paste0("factor", seq_len(H)), paste, sep = "_"))
    sig2_trace <- matrix(NA_real_, keep_target, K * H)
    colnames(sig2_trace) <- as.vector(outer(paste0("class", seq_len(K)), paste0("factor", seq_len(H)), paste, sep = "_"))
  } else {
    alpha_trace <- lambda_trace <- pi_trace <- mu_trace <- sig2_trace <- NULL
  }
  converged <- FALSE
  completed_iter <- n_iter

  t0 <- proc.time()[["elapsed"]]
  for (iter in seq_len(n_iter)) {
    iter_start <- Sys.time()
    z_start <- Sys.time()
    Z <- sample_binary_Z_joint_mfa(X = X, F = F, Lambda = Lambda, alpha = alpha)
    history$z_sample_seconds[iter] <- as.numeric(difftime(Sys.time(), z_start, units = "secs"))

    factor_start <- Sys.time()
    LtL <- crossprod(Lambda)
    Ltz <- sweep(Z, 2L, alpha, "-") %*% Lambda
    F_rows <- parallel_lapply(seq_len(n), function(i) {
      c_i <- C[i]
      prior_prec <- diag(1 / pmax(sig2[c_i, ], min_var), H)
      V <- solve(LtL + prior_prec)
      m <- V %*% (Ltz[i, ] + prior_prec %*% mu[c_i, ])
      rmvnorm_chol(as.numeric(m), V)
    }, parallel = parallel, workers = workers)
    F <- do.call(rbind, F_rows)
    history$factor_sample_seconds[iter] <- as.numeric(difftime(Sys.time(), factor_start, units = "secs"))

    class_start <- Sys.time()
    log_prob <- matrix(NA_real_, n, K)
    for (k in seq_len(K)) {
      lp <- log(pi_vec[k] + 1e-300)
      centered <- sweep(F, 2L, mu[k, ], "-")
      log_dens <- -0.5 * rowSums(sweep(centered^2, 2L, pmax(sig2[k, ], min_var), "/")) -
        0.5 * sum(log(2 * pi * pmax(sig2[k, ], min_var)))
      log_prob[, k] <- lp + log_dens
    }
    log_prob <- log_prob - apply(log_prob, 1L, max)
    prob <- exp(log_prob)
    prob <- prob / rowSums(prob)
    cumulative_prob <- t(apply(prob, 1L, cumsum))
    C <- 1L + rowSums(runif(n) > cumulative_prob)
    C[C < 1L] <- 1L
    C[C > K] <- K
    history$class_sample_seconds[iter] <- as.numeric(difftime(Sys.time(), class_start, units = "secs"))

    weight_start <- Sys.time()
    counts <- tabulate(C, nbins = K)
    pi_vec <- sample_dirichlet(counts + alpha_dirichlet)
    history$mixture_weight_seconds[iter] <- as.numeric(difftime(Sys.time(), weight_start, units = "secs"))

    mixture_parameter_start <- Sys.time()
    for (k in seq_len(K)) {
      idx <- which(C == k)
      n_k <- length(idx)
      for (h in seq_len(H)) {
        if (n_k > 0L) {
          x <- F[idx, h]
          xbar <- mean(x)
          ss <- sum((x - xbar)^2)
        } else {
          xbar <- 0
          ss <- 0
        }
        kappa_n <- kappa0 + n_k
        a_n <- a0 + n_k / 2
        b_n <- b0 + 0.5 * ss + (kappa0 * n_k * (xbar - mu0)^2) / (2 * kappa_n)
        sig2[k, h] <- max(sample_inv_gamma(a_n, b_n), min_var)
        mu_n <- (kappa0 * mu0 + n_k * xbar) / kappa_n
        mu[k, h] <- rnorm(1L, mu_n, sqrt(sig2[k, h] / kappa_n))
      }
    }
    history$mixture_parameter_seconds[iter] <- as.numeric(difftime(Sys.time(), mixture_parameter_start, units = "secs"))

    regression_start <- Sys.time()
    regression_draw <- sample_alpha_lambda_given_Z_F(
      Z = Z,
      F = F,
      tau_intercept = tau_intercept,
      tau_lambda = tau_lambda,
      parallel = parallel,
      workers = workers
    )
    alpha <- regression_draw$alpha
    Lambda <- regression_draw$Lambda
    alpha_raw <- alpha
    history$regression_seconds[iter] <- as.numeric(difftime(Sys.time(), regression_start, units = "secs"))

    normalization_start <- Sys.time()
    if (isTRUE(normalize_scale)) {
      norm_out <- normalize_joint_mfa_scale(
        F = F,
        alpha = alpha,
        Lambda = Lambda,
        mu = mu,
        sig2 = sig2,
        pi_vec = pi_vec,
        min_scale = min_scale
      )
      F <- norm_out$F
      alpha <- norm_out$alpha
      Lambda <- norm_out$Lambda
      mu <- norm_out$mu
      sig2 <- norm_out$sig2
      history$max_abs_mixture_mean_before_normalize[iter] <-
        max(abs(norm_out$marginal_mean_before), na.rm = TRUE)
      history$max_abs_log_mixture_scale_before_normalize[iter] <-
        max(abs(log(norm_out$marginal_scale_before)), na.rm = TRUE)
    }
    history$normalization_seconds[iter] <- as.numeric(difftime(Sys.time(), normalization_start, units = "secs"))

    diagnostics_start <- Sys.time()
    history$occupied_classes[iter] <- sum(counts > 0L)
    history$lambda_norm[iter] <- sqrt(sum(Lambda^2))
    history$intercept_norm[iter] <- sqrt(sum(alpha^2))

    if (iter >= check_window) {
      idx_window <- seq.int(iter - check_window + 1L, iter)
      lambda_window <- history$lambda_norm[idx_window]
      intercept_window <- history$intercept_norm[idx_window]
      occupied_window <- history$occupied_classes[idx_window]
      history$lambda_rel_range[iter] <- diff(range(lambda_window, na.rm = TRUE)) /
        max(mean(abs(lambda_window), na.rm = TRUE), 1e-8)
      history$intercept_rel_range[iter] <- diff(range(intercept_window, na.rm = TRUE)) /
        max(mean(abs(intercept_window), na.rm = TRUE), 1e-8)
      history$occupied_range[iter] <- diff(range(occupied_window, na.rm = TRUE))
      history$stable[iter] <- iter >= min_iter &&
        is.finite(history$lambda_rel_range[iter]) &&
        is.finite(history$intercept_rel_range[iter]) &&
        history$lambda_rel_range[iter] <= lambda_stability_tol &&
        history$intercept_rel_range[iter] <= intercept_stability_tol &&
        history$occupied_range[iter] <= occupied_stability_tol
    }
    history$diagnostics_seconds[iter] <- as.numeric(difftime(Sys.time(), diagnostics_start, units = "secs"))

    keep_start <- Sys.time()
    if (iter > burn && ((iter - burn) %% thin == 0L)) {
      keep_F <- keep_F + F
      keep_alpha <- keep_alpha + alpha
      keep_alpha_raw <- keep_alpha_raw + alpha_raw
      keep_Lambda <- keep_Lambda + Lambda
      keep_pi <- keep_pi + pi_vec
      keep_mu <- keep_mu + mu
      keep_sig2 <- keep_sig2 + sig2
      n_keep <- n_keep + 1L
      if (isTRUE(compute_parameter_ess) && n_keep <= keep_target) {
        alpha_trace[n_keep, ] <- alpha
        lambda_trace[n_keep, ] <- as.numeric(Lambda)
        pi_trace[n_keep, ] <- pi_vec
        mu_trace[n_keep, ] <- as.numeric(mu)
        sig2_trace[n_keep, ] <- as.numeric(sig2)
      }
    }
    history$keep_draw_seconds[iter] <- as.numeric(difftime(Sys.time(), keep_start, units = "secs"))
    history$iteration_seconds[iter] <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))

    if (isTRUE(verbose) && (iter == 1L || iter %% max(1L, floor(n_iter / 5)) == 0L)) {
      cat(sprintf(
        "joint-MFA Gibbs iter %d/%d | occupied=%d/%d | ||Lambda||=%.2f | lambda_range=%s\n",
        iter, n_iter, sum(counts > 0L), K, sqrt(sum(Lambda^2)),
        ifelse(is.na(history$lambda_rel_range[iter]), "NA", sprintf("%.3f", history$lambda_rel_range[iter]))
      ))
    }

    if (isTRUE(stop_on_stability) && isTRUE(history$stable[iter])) {
      converged <- TRUE
      completed_iter <- iter
      if (isTRUE(verbose)) {
        cat(sprintf(
          "joint-MFA Gibbs stability stop at iter %d | lambda_range=%.4f | intercept_range=%.4f | occupied_range=%d\n",
          iter,
          history$lambda_rel_range[iter],
          history$intercept_rel_range[iter],
          history$occupied_range[iter]
        ))
      }
      break
    }
  }

  history <- history[seq_len(completed_iter), , drop = FALSE]

  if (n_keep > 0L) {
    F_hat <- keep_F / n_keep
    alpha_hat <- keep_alpha / n_keep
    alpha_raw_hat <- keep_alpha_raw / n_keep
    Lambda_hat <- keep_Lambda / n_keep
    pi_hat <- keep_pi / n_keep
    mu_hat <- keep_mu / n_keep
    sig2_hat <- keep_sig2 / n_keep
  } else {
    F_hat <- F
    alpha_hat <- alpha
    alpha_raw_hat <- if (exists("alpha_raw", inherits = FALSE)) alpha_raw else alpha
    Lambda_hat <- Lambda
    pi_hat <- pi_vec
    mu_hat <- mu
    sig2_hat <- sig2
  }
  if (isTRUE(compute_parameter_ess) && n_keep > 0L) {
    trace_list <- list(
      alpha = alpha_trace[seq_len(n_keep), , drop = FALSE],
      lambda = lambda_trace[seq_len(n_keep), , drop = FALSE],
      pi = pi_trace[seq_len(n_keep), , drop = FALSE],
      mu = mu_trace[seq_len(n_keep), , drop = FALSE],
      sig2 = sig2_trace[seq_len(n_keep), , drop = FALSE]
    )
    ess_table <- compute_parameter_ess_table(trace_list)
    ess_summary <- summarize_parameter_ess(ess_table)
  } else {
    ess_table <- NULL
    ess_summary <- summarize_parameter_ess(NULL)
  }

  list(
    F_hat = F_hat,
    alpha_hat = alpha_hat,
    alpha_raw_hat = alpha_raw_hat,
    Lambda_hat = Lambda_hat,
    C = C,
    pi = pi_hat,
    mu = mu_hat,
    sig2 = sig2_hat,
    final_pi = pi_vec,
    final_mu = mu,
    final_sig2 = sig2,
    history = history,
    ess_table = ess_table,
    ess_summary = ess_summary,
    n_keep = n_keep,
    converged = converged,
    completed_iter = completed_iter,
    stability_criteria = list(
      stop_on_stability = stop_on_stability,
      min_iter = min_iter,
      check_window = check_window,
      lambda_stability_tol = lambda_stability_tol,
      intercept_stability_tol = intercept_stability_tol,
      occupied_stability_tol = occupied_stability_tol,
      normalize_scale = normalize_scale,
      min_scale = min_scale
    ),
    seconds = proc.time()[["elapsed"]] - t0,
    seconds_per_iter = (proc.time()[["elapsed"]] - t0) / completed_iter
  )
}

evaluate_fit <- function(
    method,
    sim,
    G,
    F_hat,
    Lambda_hat,
    seconds,
    alpha_hat = NULL,
    alpha_raw_hat = NULL,
    C_hat = NULL,
    mixture_fits = NULL) {
  align <- choose_alignment(sim$Lambda, Lambda_hat, sim$F, F_hat)
  Lambda_aligned <- align_lambda_to_truth(sim$Lambda, Lambda_hat, align)
  true_alpha <- if (!is.null(sim$alpha)) sim$alpha else rep(0, nrow(sim$Lambda))
  true_eta <- sweep(sim$F %*% t(sim$Lambda), 2L, true_alpha, "+")
  F_aligned <- sweep(F_hat[, align$est_index, drop = FALSE], 2L, align$signs, "*")
  standardize_scores <- function(F) {
    out <- scale(F)
    out[!is.finite(out)] <- 0
    out
  }
  factor_score_raw_rmse <- sqrt(mean((sim$F - F_aligned)^2))
  factor_score_rmse <- sqrt(mean((standardize_scores(sim$F) - standardize_scores(F_aligned))^2))
  if (is.null(alpha_hat)) alpha_hat <- rep(0, nrow(Lambda_aligned))
  if (is.null(alpha_raw_hat)) alpha_raw_hat <- rep(NA_real_, length(alpha_hat))
  eta_hat <- sweep(F_aligned %*% t(Lambda_aligned), 2L, alpha_hat, "+")
  alpha_raw_rmse <- if (all(!is.finite(alpha_raw_hat))) {
    NA_real_
  } else {
    sqrt(mean((true_alpha - alpha_raw_hat)^2, na.rm = TRUE))
  }
  prob_rmse <- sqrt(mean((pnorm(true_eta) - pnorm(eta_hat))^2))

  K_profile <- prod(as.numeric(G))
  if (is.finite(K_profile) && K_profile <= max_joint_profile_ari_K) {
    true_joint <- joint_class_index(sim$component, G)
    if (!is.null(C_hat)) {
      joint_ari <- adjusted_rand_index(true_joint, C_hat)
    } else if (!is.null(mixture_fits)) {
      class_hat <- class_map_from_mixtures_local(F_hat, mixture_fits)
      joint_ari <- adjusted_rand_index(true_joint, joint_class_index(class_hat, G))
    } else {
      joint_ari <- NA_real_
    }
  } else {
    joint_ari <- NA_real_
  }

  data.frame(
    method = method,
    alignment_mode = alignment_mode,
    mean_factor_abs_cor = align$mean_abs_cor,
    mean_loading_abs_cor = if (!is.null(align$mean_loading_abs_cor)) align$mean_loading_abs_cor else NA_real_,
    min_factor_abs_cor = min(align$matched_abs_cor),
    factor_score_rmse = factor_score_rmse,
    factor_score_raw_rmse = factor_score_raw_rmse,
    lambda_corr = suppressWarnings(cor(as.vector(sim$Lambda), as.vector(Lambda_aligned))),
    lambda_rmse = sqrt(mean((sim$Lambda - Lambda_aligned)^2)),
    alpha_corr = safe_cor(true_alpha, alpha_hat),
    alpha_rmse = sqrt(mean((true_alpha - alpha_hat)^2)),
    alpha_raw_corr = safe_cor(true_alpha, alpha_raw_hat),
    alpha_raw_rmse = alpha_raw_rmse,
    probability_rmse = prob_rmse,
    joint_profile_ari = joint_ari,
    seconds = seconds,
    stringsAsFactors = FALSE
  )
}

plot_joint_mfa_history <- function(history, out_file, title = "Joint-mixture Gibbs diagnostics") {
  if (is.null(history) || !nrow(history)) return(invisible(FALSE))
  plot_history_series <- function(x, y, type = "l", lwd = 2, col = "#377eb8",
                                  xlab = "iteration", ylab = "", main = "") {
    if (all(!is.finite(y))) {
      plot.new()
      title(main = main, xlab = xlab, ylab = ylab)
      text(0.5, 0.5, "not available")
    } else {
      plot(x, y, type = type, lwd = lwd, col = col, xlab = xlab, ylab = ylab, main = main)
    }
  }
  png(out_file, width = 1500, height = 1200, res = 170)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  plot_history_series(history$iteration, history$lambda_norm, col = "#377eb8",
                      ylab = "||Lambda||", main = "Loading norm")
  plot_history_series(history$iteration, history$intercept_norm, col = "#4daf4a",
                      ylab = "||alpha||", main = "Intercept norm")
  plot_history_series(history$iteration, history$occupied_classes, col = "#e41a1c",
                      ylab = "occupied classes", main = "Joint classes")
  plot_history_series(history$iteration, history$lambda_rel_range, col = "#984ea3",
                      ylab = "rolling relative range", main = "Lambda stability")
  mtext(title, outer = TRUE, cex = 1.05)
  invisible(TRUE)
}

write_ours_timing_history <- function(fit, out_file) {
  pieces <- list()
  if (!is.null(fit$pretrain_fit$em_history)) {
    d <- fit$pretrain_fit$em_history
    d$stage <- "lowrank_em_svd"
    pieces[[length(pieces) + 1L]] <- d
  } else if (!is.null(fit$pretrain_fit$history)) {
    d <- fit$pretrain_fit$history
    d$stage <- "sampled_z_pretraining"
    pieces[[length(pieces) + 1L]] <- d
  }
  if (!is.null(fit$pretrain_fit$rotation_history)) {
    d <- fit$pretrain_fit$rotation_history
    d$stage <- "mixture_rotation"
    pieces[[length(pieces) + 1L]] <- d
  }
  if (!is.null(fit$refine_fit$joint_refinement$history)) {
    d <- fit$refine_fit$joint_refinement$history
    d$stage <- "map_refinement"
    pieces[[length(pieces) + 1L]] <- d
  }
  timing <- rbind_fill(pieces)
  if (nrow(timing)) write.csv(timing, out_file, row.names = FALSE)
  invisible(timing)
}

summarize_ours_timing <- function(fit) {
  em_hist <- fit$pretrain_fit$em_history
  rot_hist <- fit$pretrain_fit$rotation_history
  ref_hist <- fit$refine_fit$joint_refinement$history
  nonzero_ref <- if (!is.null(ref_hist) && "iteration" %in% names(ref_hist)) {
    ref_hist[ref_hist$iteration > 0, , drop = FALSE]
  } else {
    data.frame()
  }
  median_or_na <- function(x) {
    if (is.null(x) || !length(x) || all(!is.finite(x))) NA_real_ else median(x, na.rm = TRUE)
  }
  data.frame(
    ours_median_em_iteration_seconds = median_or_na(em_hist$iteration_seconds),
    ours_median_em_e_step_seconds = median_or_na(em_hist$e_step_seconds),
    ours_median_em_projection_seconds = median_or_na(em_hist$projection_seconds),
    ours_median_em_objective_seconds = median_or_na(em_hist$objective_seconds),
    ours_median_rotation_iteration_seconds = median_or_na(rot_hist$iteration_seconds),
    ours_median_rotation_update_seconds = median_or_na(rot_hist$rotation_update_seconds),
    ours_median_rotation_mixture_refit_seconds = median_or_na(rot_hist$mixture_refit_seconds),
    ours_median_rotation_objective_seconds = median_or_na(rot_hist$objective_seconds),
    ours_median_refinement_iteration_seconds = median_or_na(nonzero_ref$iteration_seconds),
    ours_median_refinement_factor_update_seconds = median_or_na(nonzero_ref$factor_update_seconds),
    ours_median_refinement_lambda_update_seconds = median_or_na(nonzero_ref$lambda_update_seconds),
    ours_median_refinement_mixture_update_seconds = median_or_na(nonzero_ref$mixture_update_seconds),
    stringsAsFactors = FALSE
  )
}

fit_ours <- function(X, H, G, seed) {
  t0 <- proc.time()[["elapsed"]]
  if (identical(ours_pretraining_method, "em_svd")) {
    fit <- fit_binary_probit_em_svd_pretrain_then_refine(
      X = X,
      H = H,
      G_fixed = G,
      em_max_iter = em_svd_iter,
      em_tol_loglik = em_svd_tol_loglik,
      em_tol_L = em_svd_tol_L,
      em_init_method = em_svd_init_method,
      em_init_z = em_svd_init_z,
      em_random_starts = em_svd_random_starts,
      em_random_start_scale = em_svd_random_start_scale,
      rotation_random_starts = rotation_random_starts,
      rotation_loading_l1_penalty = rotation_loading_l1_penalty,
      rotation_max_outer = rotation_max_outer,
      n_mix_starts = rotation_n_mix_starts,
      mixture_max_iter = mixture_max_iter,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      refine_mu_prior_mean = refine_mu_prior_mean,
      refine_mu_prior_kappa = refine_mu_prior_kappa,
      refine_var_prior_shape = refine_var_prior_shape,
      refine_var_prior_scale = refine_var_prior_scale,
      refine_weight_prior_alpha = refine_weight_prior_alpha,
      min_mixture_var = min_mixture_var,
      grid_size = rotation_grid_size,
      rotation_sweep = rotation_sweep,
      rotation_objective_tolerance = rotation_objective_tolerance,
      rotation_min_outer = rotation_min_outer,
      require_mixture_convergence_for_rotation_stop = rotation_require_mixture_convergence,
      loading_penalty = pretrain_loading_penalty,
      n_refine_iter = refine_iter,
      factor_update = factor_update,
      lambda_l1_penalty = lambda_l1_penalty,
      lasso_backend = lasso_backend,
      glmnet_standardize = glmnet_standardize,
      objective_tolerance = refine_objective_tolerance,
      min_refine_iter = refine_min_iter,
      return_best_refinement_iteration = refine_return_best_iteration,
      refinement_selection_objective = refine_selection_objective,
      require_mixture_convergence_for_stop = refine_require_mixture_convergence,
      mixture_refit = mixture_refit,
      enforce_monotone_refinement = refine_enforce_monotone,
      monotone_tolerance = refine_monotone_tolerance,
      parallel = parallel_ours,
      workers = parallel_workers,
      seed = seed,
      verbose = FALSE
    )
  } else {
    fit <- fit_binary_probit_pretrain_then_refine(
      X = X,
      H = H,
      G_fixed = G,
      n_aug_iter = pretrain_aug_iter,
      z_update = pretrain_z_update,
      pretrain_objective = pretrain_objective,
      pretrain_objective_tolerance = pretrain_objective_tolerance,
      pretrain_objective_patience = pretrain_objective_patience,
      min_aug_iter = pretrain_min_aug_iter,
      center_Z_for_svd = center_Z_for_svd,
      return_best_iteration = pretrain_return_best_iteration,
      n_refine_iter = refine_iter,
      factor_update = factor_update,
      lambda_l1_penalty = lambda_l1_penalty,
      lasso_backend = lasso_backend,
      glmnet_standardize = glmnet_standardize,
      mixture_max_iter = mixture_max_iter,
      mixture_update = mixture_update,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      refine_mu_prior_mean = refine_mu_prior_mean,
      refine_mu_prior_kappa = refine_mu_prior_kappa,
      refine_var_prior_shape = refine_var_prior_shape,
      refine_var_prior_scale = refine_var_prior_scale,
      refine_weight_prior_alpha = refine_weight_prior_alpha,
      min_mixture_var = min_mixture_var,
      objective_tolerance = refine_objective_tolerance,
      min_refine_iter = refine_min_iter,
      stopping_objective = refine_stopping_objective,
      enforce_monotone_refinement = refine_enforce_monotone,
      monotone_tolerance = refine_monotone_tolerance,
      return_best_refinement_iteration = refine_return_best_iteration,
      refinement_selection_objective = refine_selection_objective,
      require_mixture_convergence_for_stop = refine_require_mixture_convergence,
      mixture_refit = mixture_refit,
      parallel = parallel_ours,
      workers = parallel_workers,
      seed = seed,
      verbose = FALSE
    )
  }
  fit$seconds <- proc.time()[["elapsed"]] - t0
  fit
}

if (!is.null(G_configs_input)) {
  design_chunks <- lapply(H_values, function(H_value) {
    valid_input <- Filter(function(config) length(config) %in% c(1L, H_value), G_configs_input)
    if (!length(valid_input)) {
      stop("No valid G_CONFIGS entry was supplied for H = ", H_value, ".")
    }
    valid_configs <- lapply(valid_input, normalize_G_counts, H = H_value)
    config_labels <- vapply(valid_configs, format_G_config, character(1L))
    unique_config_labels <- unique(config_labels)
    expand.grid(
      rep = rep_values,
      H_true = H_value,
      G_config = unique_config_labels,
      np_index = seq_len(nrow(np_settings)),
      separation = separations,
      loading_design = loading_designs,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  })
  design_grid <- do.call(rbind, design_chunks)
  rownames(design_grid) <- NULL
} else {
  design_grid <- expand.grid(
    rep = rep_values,
    H_true = H_values,
    G_config = as.character(G_values),
    np_index = seq_len(nrow(np_settings)),
    separation = separations,
    loading_design = loading_designs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

checkpoint_file <- file.path(out_dir, "comparison_results_checkpoint.csv")
final_results_file <- file.path(out_dir, "comparison_results.csv")
existing_results <- data.frame()
if (isTRUE(resume_existing)) {
  if (file.exists(checkpoint_file)) {
    existing_results <- read.csv(checkpoint_file)
  } else if (file.exists(final_results_file)) {
    existing_results <- read.csv(final_results_file)
  }
}

make_result_key <- function(scenario, method) paste(scenario, method, sep = " || ")
completed_keys <- character(0L)
if (nrow(existing_results) && all(c("scenario", "method") %in% names(existing_results))) {
  completed_keys <- make_result_key(existing_results$scenario, existing_results$method)
  cat(sprintf("Resuming from %d existing result rows in %s\n", nrow(existing_results), out_dir))
}

all_results <- vector("list", nrow(existing_results) + nrow(design_grid) * 3L)
result_idx <- 0L
if (nrow(existing_results)) {
  for (i in seq_len(nrow(existing_results))) {
    result_idx <- result_idx + 1L
    all_results[[result_idx]] <- existing_results[i, , drop = FALSE]
  }
}

for (row_idx in seq_len(nrow(design_grid))) {
  row <- design_grid[row_idx, ]
  np_row <- np_settings[row$np_index, ]
  H_scenario <- as.integer(row$H_true)
  G_scenario <- normalize_G_counts(parse_g_count_string(row$G_config), H_scenario)
  G_label <- format_G_config(G_scenario)
  K_joint <- prod(as.numeric(G_scenario))
  mfa_alpha_dirichlet <- if (is.finite(mfa_alpha_dirichlet_override)) {
    mfa_alpha_dirichlet_override
  } else {
    1 / max(1, K_joint)
  }
  seed <- seed_base + 100000L * H_scenario + 50000L * sum(G_scenario * seq_along(G_scenario)) +
    10000L * row$rep + 1000L * row_idx
  scenario <- sprintf(
    "rep%d_%s_n%d_p%d_H%d_G%s_sep%s",
    row$rep, row$loading_design, np_row$n, np_row$p, H_scenario, G_label, row$separation
  )
  cat("\nScenario:", scenario, "\n")

  need_ours <- isTRUE(run_ours) &&
    !(make_result_key(scenario, "independent_marginal_mixture") %in% completed_keys)
  need_joint_mfa <- isTRUE(run_joint_mfa) &&
    !(make_result_key(scenario, "joint_mixture_factor_gibbs") %in% completed_keys)
  need_viroli <- isTRUE(run_viroli) &&
    !(make_result_key(scenario, viroli_method_name) %in% completed_keys)

  if (!need_ours && !need_joint_mfa && !need_viroli) {
    cat("Skipping completed scenario.\n")
    next
  }

  sim <- simulate_original_binary_probit(
    n = np_row$n,
    p = np_row$p,
      H = H_scenario,
      G = G_scenario,
    sep = row$separation,
    loading_design = row$loading_design,
    seed = seed
  )
  dgp_summary <- dgp_diagnostics(sim, G_scenario)
  if (is.finite(K_joint) && K_joint <= max_joint_parameter_K) {
    true_joint_params <- true_joint_mixture_parameters(sim$mixture_params)
  } else {
    true_joint_params <- NULL
    if (need_joint_mfa) {
      cat(sprintf(
        "Skipping joint MFA for K=%s because dense joint parameter recovery is capped at MAX_JOINT_PARAMETER_K=%s.\n",
        format(K_joint, scientific = FALSE),
        max_joint_parameter_K
      ))
      need_joint_mfa <- FALSE
    }
  }

  if (need_ours) {
    cat("Fitting independent marginal-mixture probit factor model...\n")
    ours <- fit_ours(sim$X_binary, H = H_scenario, G = G_scenario, seed = seed + 11L)
    ours_eval <- evaluate_fit(
      method = "independent_marginal_mixture",
      sim = sim,
      G = G_scenario,
      F_hat = ours$refine_fit$F_hat,
      Lambda_hat = ours$refine_fit$Lambda_hat,
      seconds = ours$seconds,
      alpha_hat = ours$refine_fit$alpha_hat,
      mixture_fits = ours$refine_fit$mixture_fits
    )
    ours_align <- choose_alignment(sim$Lambda, ours$refine_fit$Lambda_hat, sim$F, ours$refine_fit$F_hat)
    ours_lambda_aligned <- align_lambda_to_truth(sim$Lambda, ours$refine_fit$Lambda_hat, ours_align)
    ours_mixture_fits_aligned <- align_product_mixture_fits(ours$refine_fit$mixture_fits, ours_align)
    ours_marginal_recovery <- marginal_mixture_parameter_summary(
      sim$mixture_params,
      ours_mixture_fits_aligned,
      G_scenario
    )
    if (!is.null(true_joint_params)) {
      ours_param_recovery <- align_product_mixture_parameters_by_axis(
        true_params = true_joint_params,
        mixture_fits = ours_mixture_fits_aligned,
        G = G_scenario
      )
      ours_all_param_recovery <- combined_parameter_correlations(
        Lambda_true = sim$Lambda,
        Lambda_est_aligned = ours_lambda_aligned,
        parameter_table = ours_param_recovery$table,
        alpha_true = sim$alpha,
        alpha_est = ours$refine_fit$alpha_hat
      )
      ours_param_table <- cbind(method = "independent_marginal_mixture", ours_param_recovery$table)
      maybe_write_parameter_recovery(
        ours_param_table,
        file.path(out_dir, paste0("joint_parameter_recovery_ours_", scenario, ".csv")),
        K_joint = K_joint
      )
    } else {
      ours_param_recovery <- list(table = data.frame(), summary = empty_joint_parameter_summary("skipped_large_K"))
      ours_all_param_recovery <- data.frame(
        flat_parameter_corr = NA_real_,
        flat_parameter_blocks = "skipped_large_K",
        all_parameter_corr = NA_real_,
        all_parameter_corr_raw = NA_real_,
        all_parameter_blocks = "skipped_large_K",
        stringsAsFactors = FALSE
      )
    }
    if (isTRUE(write_parameter_tables)) {
      write.csv(
        cbind(method = "independent_marginal_mixture", ours_marginal_recovery$table),
        file.path(out_dir, paste0("marginal_parameter_recovery_ours_", scenario, ".csv")),
        row.names = FALSE
      )
    }
    if (isTRUE(write_iteration_histories)) {
      write_ours_timing_history(
        ours,
        file.path(out_dir, paste0("ours_timing_history_", scenario, ".csv"))
      )
    }
    result_idx <- result_idx + 1L
    all_results[[result_idx]] <- cbind(
      scenario = scenario,
      rep = row$rep,
      np_regime = np_row$np_regime,
      n = np_row$n,
      p = np_row$p,
      H_true = H_scenario,
      G_true = G_label,
      G_config = G_label,
      K_joint = K_joint,
      separation = row$separation,
      mixture_param_mode = mixture_param_mode,
      mixture_variance_mode = mixture_variance_mode,
      intercept_mode = intercept_mode,
      intercept_sd = intercept_sd,
      intercept_block_span = intercept_block_span,
      intercept_clip = intercept_clip,
      loading_design = row$loading_design,
      block_size_mode = block_size_mode,
      loading_sign_mode = loading_sign_mode,
      alignment_mode = alignment_mode,
      dgp_summary,
      ours_pretraining_method = ours_pretraining_method,
      em_svd_init_method = if (identical(ours_pretraining_method, "em_svd")) em_svd_init_method else NA_character_,
      em_svd_init_z = if (identical(ours_pretraining_method, "em_svd")) em_svd_init_z else NA_character_,
      selected_em_start = if (identical(ours_pretraining_method, "em_svd") && !is.null(ours$pretrain_fit$selected_em_start)) {
        ours$pretrain_fit$selected_em_start
      } else {
        NA_character_
      },
      mixture_update = mixture_update,
      pretrain_loading_penalty = pretrain_loading_penalty,
      rotation_loading_l1_penalty = rotation_loading_l1_penalty,
      lambda_l1_penalty = lambda_l1_penalty,
      lasso_backend = lasso_backend,
      glmnet_standardize = glmnet_standardize,
      refine_enforce_monotone = refine_enforce_monotone,
      refine_monotone_tolerance = refine_monotone_tolerance,
      mu_prior_mean = mu_prior_mean,
      mu_prior_kappa = mu_prior_kappa,
      var_prior_shape = var_prior_shape,
      var_prior_scale = var_prior_scale,
      weight_prior_alpha = weight_prior_alpha,
      refine_mu_prior_mean = refine_mu_prior_mean,
      refine_mu_prior_kappa = refine_mu_prior_kappa,
      refine_var_prior_shape = refine_var_prior_shape,
      refine_var_prior_scale = refine_var_prior_scale,
      refine_weight_prior_alpha = refine_weight_prior_alpha,
      min_mixture_var = min_mixture_var,
      parallel_enabled = parallel_ours,
      parallel_workers = if (is.null(parallel_workers)) NA_integer_ else parallel_workers,
      ours_eval,
      ours_convergence_summary(ours),
      ours_marginal_recovery$summary,
      ours_param_recovery$summary,
      ours_all_param_recovery,
      summarize_ours_timing(ours)
    )
    completed_keys <- c(completed_keys, make_result_key(scenario, "independent_marginal_mixture"))
  }

  if (need_joint_mfa) {
    cat(sprintf("Fitting joint-mixture probit factor Gibbs baseline with K=%s...\n", format(K_joint, scientific = FALSE)))
    mfa <- fit_joint_class_probit_mfa_gibbs(
      X = sim$X_binary,
      H = H_scenario,
      K = K_joint,
      n_iter = mfa_iter,
      burn = mfa_burn,
      thin = mfa_thin,
      tau_lambda = mfa_tau_lambda,
      tau_intercept = mfa_tau_intercept,
      alpha_dirichlet = mfa_alpha_dirichlet,
      stop_on_stability = mfa_stop_on_stability,
      min_iter = mfa_min_iter,
      check_window = mfa_check_window,
      lambda_stability_tol = mfa_lambda_stability_tol,
      intercept_stability_tol = mfa_intercept_stability_tol,
      occupied_stability_tol = mfa_occupied_stability_tol,
      normalize_scale = mfa_normalize_scale,
      min_scale = mfa_min_scale,
      parallel = parallel_gibbs,
      workers = parallel_workers,
      compute_parameter_ess = mfa_compute_parameter_ess,
      seed = seed + 23L,
      verbose = mfa_verbose
    )
    mfa_eval <- evaluate_fit(
      method = "joint_mixture_factor_gibbs",
      sim = sim,
      G = G_scenario,
      F_hat = mfa$F_hat,
      Lambda_hat = mfa$Lambda_hat,
      seconds = mfa$seconds,
      alpha_hat = mfa$alpha_hat,
      alpha_raw_hat = mfa$alpha_raw_hat,
      C_hat = mfa$C
    )
    mfa_align <- choose_alignment(sim$Lambda, mfa$Lambda_hat, sim$F, mfa$F_hat)
    mfa_lambda_aligned <- align_lambda_to_truth(sim$Lambda, mfa$Lambda_hat, mfa_align)
    mfa_axis_params <- align_joint_mfa_axis_parameters(mfa$mu, mfa$sig2, mfa_align)
    mfa_param_recovery <- align_joint_mixture_parameters(
      true_params = true_joint_params,
      est_pi = mfa$pi,
      est_mu = mfa_axis_params$mu,
      est_sig2 = mfa_axis_params$sig2
    )
    mfa_all_param_recovery <- combined_parameter_correlations(
      Lambda_true = sim$Lambda,
      Lambda_est_aligned = mfa_lambda_aligned,
      parameter_table = mfa_param_recovery$table,
      alpha_true = sim$alpha,
      alpha_est = mfa$alpha_hat
    )
    mfa_param_table <- cbind(method = "joint_mixture_factor_gibbs", mfa_param_recovery$table)
    maybe_write_parameter_recovery(
      mfa_param_table,
      file.path(out_dir, paste0("joint_parameter_recovery_gibbs_", scenario, ".csv")),
      K_joint = K_joint
    )
    result_idx <- result_idx + 1L
    all_results[[result_idx]] <- cbind(
      scenario = scenario,
      rep = row$rep,
      np_regime = np_row$np_regime,
      n = np_row$n,
      p = np_row$p,
      H_true = H_scenario,
      G_true = G_label,
      G_config = G_label,
      K_joint = K_joint,
      separation = row$separation,
      mixture_param_mode = mixture_param_mode,
      mixture_variance_mode = mixture_variance_mode,
      intercept_mode = intercept_mode,
      intercept_sd = intercept_sd,
      intercept_block_span = intercept_block_span,
      intercept_clip = intercept_clip,
      loading_design = row$loading_design,
      block_size_mode = block_size_mode,
      loading_sign_mode = loading_sign_mode,
      alignment_mode = alignment_mode,
      dgp_summary,
      ours_pretraining_method = NA_character_,
      parallel_enabled = parallel_gibbs,
      parallel_workers = if (is.null(parallel_workers)) NA_integer_ else parallel_workers,
      mfa_eval,
      mfa_param_recovery$summary,
      mfa_all_param_recovery,
      mfa$ess_summary,
      summarize_gibbs_timing(mfa$history, burn = mfa_burn, thin = mfa_thin),
      occupied_classes = tail(mfa$history$occupied_classes, 1L),
      posterior_kept_draws = mfa$n_keep,
      mfa_converged = mfa$converged,
      mfa_completed_iter = mfa$completed_iter,
      mfa_final_lambda_rel_range = tail(mfa$history$lambda_rel_range, 1L),
      mfa_final_intercept_rel_range = tail(mfa$history$intercept_rel_range, 1L),
      mfa_final_occupied_range = tail(mfa$history$occupied_range, 1L)
    )
    completed_keys <- c(completed_keys, make_result_key(scenario, "joint_mixture_factor_gibbs"))

    if (isTRUE(write_iteration_histories)) {
      write.csv(
        mfa$history,
        file.path(out_dir, paste0("joint_mfa_history_", scenario, ".csv")),
        row.names = FALSE
      )
      plot_joint_mfa_history(
        mfa$history,
        file.path(out_dir, paste0("joint_mfa_history_", scenario, ".png")),
        title = scenario
      )
    }
    if (!is.null(mfa$ess_table) && nrow(mfa$ess_table)) {
      write.csv(
        mfa$ess_table,
        file.path(out_dir, paste0("joint_mfa_parameter_ess_", scenario, ".csv")),
        row.names = FALSE
      )
    }
  }

  if (need_viroli) {
    cat(sprintf(
      "Fitting Viroli-style independent-mixture probit factor Gibbs baseline with G=%s...\n",
      G_label
    ))
    viroli <- fit_viroli_probit_independent_gibbs(
      X = sim$X_binary,
      H = H_scenario,
      G = G_scenario,
      n_iter = viroli_iter,
      burn = viroli_burn,
      thin = viroli_thin,
      tau_lambda = viroli_tau_lambda,
      tau_intercept = viroli_tau_intercept,
      lambda_l1_penalty = viroli_lambda_l1_penalty,
      alpha_dirichlet = viroli_alpha_dirichlet,
      min_scale = viroli_min_scale,
      normalize_each_draw = viroli_normalize_each_draw,
      parallel = parallel_gibbs,
      workers = parallel_workers,
      compute_parameter_ess = viroli_compute_parameter_ess,
      seed = if (is.finite(viroli_seed_override)) as.integer(viroli_seed_override) else seed + 37L,
      verbose = viroli_verbose
    )
    viroli_eval <- evaluate_fit(
      method = viroli_method_name,
      sim = sim,
      G = G_scenario,
      F_hat = viroli$F_hat,
      Lambda_hat = viroli$Lambda_hat,
      seconds = viroli$seconds,
      alpha_hat = viroli$alpha_hat,
      mixture_fits = viroli$mixture_fits
    )
    viroli_align <- choose_alignment(sim$Lambda, viroli$Lambda_hat, sim$F, viroli$F_hat)
    viroli_lambda_aligned <- align_lambda_to_truth(sim$Lambda, viroli$Lambda_hat, viroli_align)
    viroli_mixture_fits_aligned <- align_product_mixture_fits(viroli$mixture_fits, viroli_align)
    viroli_marginal_recovery <- marginal_mixture_parameter_summary(
      sim$mixture_params,
      viroli_mixture_fits_aligned,
      G_scenario
    )
    if (!is.null(true_joint_params)) {
      viroli_param_recovery <- align_product_mixture_parameters_by_axis(
        true_params = true_joint_params,
        mixture_fits = viroli_mixture_fits_aligned,
        G = G_scenario
      )
      viroli_all_param_recovery <- combined_parameter_correlations(
        Lambda_true = sim$Lambda,
        Lambda_est_aligned = viroli_lambda_aligned,
        parameter_table = viroli_param_recovery$table,
        alpha_true = sim$alpha,
        alpha_est = viroli$alpha_hat
      )
      viroli_param_table <- cbind(method = viroli_method_name, viroli_param_recovery$table)
      maybe_write_parameter_recovery(
        viroli_param_table,
        file.path(out_dir, paste0("joint_parameter_recovery_", viroli_method_name, "_", scenario, ".csv")),
        K_joint = K_joint
      )
    } else {
      viroli_param_recovery <- list(table = data.frame(), summary = empty_joint_parameter_summary("skipped_large_K"))
      viroli_all_param_recovery <- data.frame(
        flat_parameter_corr = NA_real_,
        flat_parameter_blocks = "skipped_large_K",
        all_parameter_corr = NA_real_,
        all_parameter_corr_raw = NA_real_,
        all_parameter_blocks = "skipped_large_K",
        stringsAsFactors = FALSE
      )
    }
    if (isTRUE(write_parameter_tables)) {
      write.csv(
        cbind(method = viroli_method_name, viroli_marginal_recovery$table),
        file.path(out_dir, paste0("marginal_parameter_recovery_", viroli_method_name, "_", scenario, ".csv")),
        row.names = FALSE
      )
    }
    result_idx <- result_idx + 1L
    all_results[[result_idx]] <- cbind(
      scenario = scenario,
      rep = row$rep,
      np_regime = np_row$np_regime,
      n = np_row$n,
      p = np_row$p,
      H_true = H_scenario,
      G_true = G_label,
      G_config = G_label,
      K_joint = K_joint,
      separation = row$separation,
      mixture_param_mode = mixture_param_mode,
      mixture_variance_mode = mixture_variance_mode,
      intercept_mode = intercept_mode,
      intercept_sd = intercept_sd,
      intercept_block_span = intercept_block_span,
      intercept_clip = intercept_clip,
      loading_design = row$loading_design,
      block_size_mode = block_size_mode,
      loading_sign_mode = loading_sign_mode,
      alignment_mode = alignment_mode,
      dgp_summary,
      ours_pretraining_method = NA_character_,
      viroli_lambda_l1_penalty = viroli_lambda_l1_penalty,
      viroli_loading_prior = viroli$loading_prior,
      min_mixture_var = min_mixture_var,
      viroli_seed = if (is.finite(viroli_seed_override)) as.integer(viroli_seed_override) else seed + 37L,
      parallel_enabled = parallel_gibbs,
      parallel_workers = if (is.null(parallel_workers)) NA_integer_ else parallel_workers,
      viroli_eval,
      viroli_marginal_recovery$summary,
      viroli_param_recovery$summary,
      viroli_all_param_recovery,
      viroli$ess_summary,
      summarize_gibbs_timing(viroli$history, burn = viroli_burn, thin = viroli_thin),
      occupied_classes = tail(viroli$history$occupied_factor_classes, 1L),
      posterior_kept_draws = viroli$n_keep,
      mfa_converged = NA,
      mfa_completed_iter = nrow(viroli$history),
      mfa_final_lambda_rel_range = NA_real_,
      mfa_final_intercept_rel_range = NA_real_,
      mfa_final_occupied_range = NA_real_
    )
    completed_keys <- c(completed_keys, make_result_key(scenario, viroli_method_name))

    if (isTRUE(write_iteration_histories)) {
      write.csv(
        viroli$history,
        file.path(out_dir, paste0("viroli_history_", scenario, ".csv")),
        row.names = FALSE
      )
    }
    if (!is.null(viroli$ess_table) && nrow(viroli$ess_table)) {
      write.csv(
        viroli$ess_table,
        file.path(out_dir, paste0("viroli_parameter_ess_", scenario, ".csv")),
        row.names = FALSE
      )
    }
  }

  results_so_far <- rbind_fill(all_results[seq_len(result_idx)])
  write.csv(results_so_far, checkpoint_file, row.names = FALSE)
}

results <- rbind_fill(all_results[seq_len(result_idx)])
write.csv(results, final_results_file, row.names = FALSE)

summarize_results <- function(results) {
  if (!nrow(results)) return(data.frame())
  group_cols <- c(
    "method", "np_regime", "n", "p", "H_true", "G_true", "K_joint",
    "separation", "mixture_param_mode", "mixture_variance_mode", "loading_design",
    "block_size_mode", "loading_sign_mode", "alignment_mode",
    "pretrain_loading_penalty", "rotation_loading_l1_penalty", "lambda_l1_penalty",
    "viroli_lambda_l1_penalty", "viroli_loading_prior", "min_mixture_var"
  )
  metric_cols <- c(
    "mean_factor_abs_cor", "min_factor_abs_cor", "factor_score_rmse",
    "factor_score_raw_rmse", "lambda_corr", "lambda_rmse",
    "alpha_corr", "alpha_rmse", "alpha_raw_corr", "alpha_raw_rmse",
    "marginal_mu_rmse", "marginal_var_rmse", "marginal_weight_rmse",
    "marginal_mu_corr", "marginal_log_var_corr", "marginal_weight_corr",
    "probability_rmse", "joint_profile_ari", "joint_mu_rmse", "joint_var_rmse",
    "joint_mu_corr", "joint_var_corr", "joint_weight_corr", "joint_weight_rmse",
    "joint_weight_l1",
    "flat_parameter_corr", "all_parameter_corr", "all_parameter_corr_raw",
    "joint_weight_max_abs_error", "est_effective_classes_001",
    "est_effective_classes_01", "seconds"
  )
  character_cols <- intersect(c("flat_parameter_blocks", "all_parameter_blocks"), names(results))
  group_cols <- intersect(group_cols, names(results))
  metric_cols <- intersect(metric_cols, names(results))
  group_key <- safe_group_key(results, group_cols)
  out <- lapply(split(results, group_key), function(d) {
    base <- d[1L, group_cols, drop = FALSE]
    metrics <- as.data.frame(lapply(metric_cols, function(metric) {
      x <- d[[metric]]
      if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    }))
    names(metrics) <- metric_cols
    character_summaries <- if (length(character_cols)) {
      out <- as.data.frame(lapply(character_cols, function(metric) {
        value <- unique(d[[metric]][!is.na(d[[metric]]) & nzchar(d[[metric]])])
        if (length(value)) value[1L] else NA_character_
      }))
      names(out) <- character_cols
      out
    } else {
      data.frame(row.names = seq_len(nrow(base)))
    }
    cbind(base, metrics, character_summaries)
  })
  combined <- do.call(rbind, out)
  rownames(combined) <- NULL
  combined
}

summary <- summarize_results(results)
write.csv(summary, file.path(out_dir, "comparison_summary.csv"), row.names = FALSE)

first_n_reaching <- function(d, metric, threshold) {
  if (!metric %in% names(d)) return(NA_integer_)
  ok <- is.finite(d[[metric]]) & d[[metric]] >= threshold
  if (!any(ok)) return(NA_integer_)
  min(d$n[ok], na.rm = TRUE)
}

make_sample_size_threshold_summary <- function(summary) {
  if (!nrow(summary)) return(data.frame())
  group_cols <- c(
    "method", "H_true", "G_true", "K_joint", "mixture_param_mode",
    "mixture_variance_mode", "loading_design", "loading_sign_mode", "alignment_mode"
  )
  group_key <- safe_group_key(summary, group_cols)
  out <- lapply(split(summary, group_key), function(d) {
    d <- d[order(d$n, d$p), , drop = FALSE]
    base <- d[1L, group_cols, drop = FALSE]
    data.frame(
      base,
      n_for_factor_corr = first_n_reaching(d, "mean_factor_abs_cor", success_factor_corr),
      n_for_alpha_corr = first_n_reaching(d, "alpha_corr", success_parameter_corr),
      n_for_lambda_corr = first_n_reaching(d, "lambda_corr", success_factor_corr),
      n_for_joint_mu_corr = first_n_reaching(d, "joint_mu_corr", success_parameter_corr),
      n_for_joint_var_corr = first_n_reaching(d, "joint_var_corr", success_parameter_corr),
      n_for_joint_weight_corr = first_n_reaching(d, "joint_weight_corr", success_parameter_corr),
      factor_corr_threshold = success_factor_corr,
      parameter_corr_threshold = success_parameter_corr,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

threshold_summary <- make_sample_size_threshold_summary(summary)
write.csv(
  threshold_summary,
  file.path(out_dir, "sample_size_threshold_summary.csv"),
  row.names = FALSE
)

plot_metric_bars <- function(summary, out_file) {
  metrics <- c(
    "mean_factor_abs_cor",
    "alpha_corr",
    "alpha_rmse",
    "alpha_raw_corr",
    "alpha_raw_rmse",
    "lambda_corr",
    "joint_mu_rmse",
    "joint_mu_corr",
    "joint_var_rmse",
    "joint_weight_l1",
    "probability_rmse",
    "joint_profile_ari",
    "seconds"
  )
  metrics <- metrics[metrics %in% names(summary)]
  if (!length(metrics) || !nrow(summary)) return(invisible(FALSE))

  labels <- paste(summary$method, summary$loading_design, sep = "\n")
  png(out_file, width = 1800, height = 1100, res = 170)
  op <- par(mfrow = c(2, 3), mar = c(8, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  for (metric in metrics) {
    vals <- summary[[metric]]
    ylim <- range(c(0, vals), na.rm = TRUE)
    if (metric %in% c("mean_factor_abs_cor", "alpha_corr", "alpha_raw_corr", "lambda_corr", "joint_profile_ari")) {
      ylim <- c(0, 1)
    }
    metric_bar_colors <- c(
      independent_marginal_mixture = "#377eb8",
      joint_mixture_factor_gibbs = "#e41a1c",
      viroli_independent_factor_gibbs = "#4daf4a",
      viroli_laplace_gibbs = "#4daf4a",
      viroli_gaussian_gibbs = "#984ea3"
    )
    cols <- metric_bar_colors[summary$method]
    cols[is.na(cols)] <- "#777777"
    barplot(
      vals,
      names.arg = labels,
      las = 2,
      cex.names = 0.62,
      col = cols,
      border = NA,
      ylim = ylim,
      main = metric,
      ylab = metric
    )
    box()
  }
  mtext("Original simulation: independent product mixture vs joint-mixture factor Gibbs", outer = TRUE, cex = 1.05)
  invisible(TRUE)
}

plot_sample_size_correlation_lines <- function(summary, out_file) {
  corr_metrics <- c(
    flat_parameter_corr = "flat all-parameter",
    all_parameter_corr = "all parameters",
    alpha_corr = "intercept",
    alpha_raw_corr = "raw Gibbs intercept",
    mean_factor_abs_cor = "factor score",
    lambda_corr = "loading",
    joint_mu_corr = "joint mixture mean",
    joint_var_corr = "joint mixture variance",
    joint_weight_corr = "joint mixture weight"
  )
  corr_metrics <- corr_metrics[names(corr_metrics) %in% names(summary)]
  finite_metric <- vapply(names(corr_metrics), function(metric) {
    any(is.finite(summary[[metric]]))
  }, logical(1L))
  corr_metrics <- corr_metrics[finite_metric]
  if (!length(corr_metrics) || !nrow(summary)) return(invisible(FALSE))

  H_values_plot <- sort(unique(summary$H_true))
  n_panels <- length(corr_metrics) * length(H_values_plot)
  n_col <- length(corr_metrics)
  n_row <- length(H_values_plot)
  png(out_file, width = max(1500, 460 * n_col), height = max(700, 430 * n_row), res = 170)
  op <- par(
    mfrow = c(n_row, n_col),
    mar = c(4.2, 4.2, 3.2, 1),
    oma = c(0, 0, 3, 0)
  )
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  method_colors <- c(
    independent_marginal_mixture = "#377eb8",
    joint_mixture_factor_gibbs = "#e41a1c",
    viroli_independent_factor_gibbs = "#4daf4a",
    viroli_laplace_gibbs = "#4daf4a",
    viroli_gaussian_gibbs = "#984ea3"
  )
  method_labels <- c(
    independent_marginal_mixture = "product mixture",
    joint_mixture_factor_gibbs = "joint Gibbs",
    viroli_independent_factor_gibbs = "Viroli Gibbs",
    viroli_laplace_gibbs = "Viroli Laplace",
    viroli_gaussian_gibbs = "Viroli Gaussian"
  )
  p_values <- sort(unique(summary$p))
  p_ltys <- setNames(rep_len(c(1, 2, 3, 4, 5, 6), length(p_values)), as.character(p_values))
  p_pchs <- setNames(rep_len(c(16, 17, 15, 18, 8, 3), length(p_values)), as.character(p_values))

  legend_drawn <- FALSE
  for (H_value in H_values_plot) {
    for (metric in names(corr_metrics)) {
      d_metric <- summary[summary$H_true == H_value & is.finite(summary[[metric]]), , drop = FALSE]
      if (!nrow(d_metric)) {
        plot.new()
        title(main = sprintf("%s | H=%s", corr_metrics[[metric]], H_value))
        next
      }

      xlim <- range(d_metric$n, na.rm = TRUE)
      if (diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5)
      plot(
        NA,
        xlim = xlim,
        ylim = c(-1, 1),
        xlab = "n",
        ylab = "correlation",
        main = sprintf("%s | H=%s", corr_metrics[[metric]], H_value)
      )
      abline(h = c(0, 0.9), col = c("#c7c7c7", "#777777"), lty = c(1, 3), lwd = c(1, 1.2))
      box()

      for (method in unique(d_metric$method)) {
        for (p_value in p_values) {
          d_line <- d_metric[d_metric$method == method & d_metric$p == p_value, , drop = FALSE]
          if (!nrow(d_line)) next
          d_line <- d_line[order(d_line$n), , drop = FALSE]
          lines(
            d_line$n,
            d_line[[metric]],
            col = method_colors[[method]],
            lty = p_ltys[[as.character(p_value)]],
            lwd = 2
          )
          points(
            d_line$n,
            d_line[[metric]],
            col = method_colors[[method]],
            pch = p_pchs[[as.character(p_value)]],
            cex = 0.8
          )
        }
      }

      if (!legend_drawn) {
        method_legend <- unique(d_metric$method)
        legend(
          "bottomright",
          legend = c(
            method_labels[method_legend],
            paste0("p=", p_values)
          ),
          col = c(method_colors[method_legend], rep("#333333", length(p_values))),
          lty = c(rep(1, length(method_legend)), p_ltys[as.character(p_values)]),
          pch = c(rep(NA, length(method_legend)), p_pchs[as.character(p_values)]),
          lwd = c(rep(2, length(method_legend)), rep(1.5, length(p_values))),
          bty = "n",
          cex = 0.75
        )
        legend_drawn <- TRUE
      }
    }
  }

  mtext("Sample-size correlation recovery by method and p", outer = TRUE, cex = 1.1)
  invisible(TRUE)
}

summarize_flat_parameter_reps <- function(results) {
  if (!nrow(results) || !"flat_parameter_corr" %in% names(results)) return(data.frame())
  group_cols <- c(
    "method", "n", "p", "H_true", "G_true", "K_joint", "separation",
    "mixture_param_mode", "mixture_variance_mode", "loading_design",
    "loading_sign_mode", "alignment_mode"
  )
  group_cols <- intersect(group_cols, names(results))
  group_key <- safe_group_key(results, group_cols)
  out <- lapply(split(results, group_key), function(d) {
    x <- d$flat_parameter_corr[is.finite(d$flat_parameter_corr)]
    base <- d[1L, group_cols, drop = FALSE]
    data.frame(
      base,
      mean_flat_parameter_corr = if (length(x)) mean(x) else NA_real_,
      sd_flat_parameter_corr = if (length(x) > 1L) sd(x) else NA_real_,
      n_reps = length(x),
      lower_2sd = if (length(x) > 1L) mean(x) - 2 * sd(x) else NA_real_,
      upper_2sd = if (length(x) > 1L) mean(x) + 2 * sd(x) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  combined <- do.call(rbind, out)
  rownames(combined) <- NULL
  combined
}

plot_flat_parameter_mean_sd <- function(flat_summary, out_file) {
  if (!nrow(flat_summary)) return(invisible(FALSE))
  panel_grid <- unique(flat_summary[, c("H_true", "G_true"), drop = FALSE])
  panel_grid <- panel_grid[order(panel_grid$H_true, panel_grid$G_true), , drop = FALSE]
  p_values <- sort(unique(flat_summary$p))
  n_col <- min(2L, nrow(panel_grid))
  n_row <- ceiling(nrow(panel_grid) / n_col)
  png(out_file, width = max(1300, 680 * n_col), height = max(900, 620 * n_row), res = 170)
  op <- par(mfrow = c(n_row, n_col), mar = c(4.4, 4.4, 3.2, 1), oma = c(4, 0, 3, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  method_colors <- c(
    independent_marginal_mixture = "#377eb8",
    joint_mixture_factor_gibbs = "#e41a1c",
    viroli_independent_factor_gibbs = "#4daf4a",
    viroli_laplace_gibbs = "#4daf4a",
    viroli_gaussian_gibbs = "#984ea3"
  )
  method_labels <- c(
    independent_marginal_mixture = "product mixture",
    joint_mixture_factor_gibbs = "joint Gibbs",
    viroli_independent_factor_gibbs = "Viroli Gibbs",
    viroli_laplace_gibbs = "Viroli Laplace",
    viroli_gaussian_gibbs = "Viroli Gaussian"
  )
  p_ltys <- setNames(rep_len(c(1, 2, 3, 4, 5, 6), length(p_values)), as.character(p_values))

  for (panel_idx in seq_len(nrow(panel_grid))) {
    H_value <- panel_grid$H_true[panel_idx]
    G_value <- as.character(panel_grid$G_true[panel_idx])
    d_H <- flat_summary[
      flat_summary$H_true == H_value & as.character(flat_summary$G_true) == G_value,
      ,
      drop = FALSE
    ]
    K_label <- if ("K_joint" %in% names(d_H)) d_H$K_joint[1L] else NA_integer_
    xlim <- range(d_H$n, na.rm = TRUE)
    if (diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5)
    ylim <- range(c(d_H$lower_2sd, d_H$upper_2sd, d_H$mean_flat_parameter_corr, 0, 1), na.rm = TRUE)
    ylim <- c(max(-1, ylim[1L]), min(1, ylim[2L]))
    plot(
      NA,
      xlim = xlim,
      ylim = ylim,
      xlab = "n",
      ylab = "flat parameter correlation",
      main = sprintf("H=%s, G=%s (K=%s)", H_value, G_value, K_label),
      log = "x",
      xaxt = "n"
    )
    axis(1, at = sort(unique(d_H$n)), labels = sort(unique(d_H$n)))
    abline(h = c(0, 0.9), col = c("#c7c7c7", "#777777"), lty = c(1, 3), lwd = c(1, 1.2))
    box()

    for (method in unique(d_H$method)) {
      for (p_value in p_values) {
        d_line <- d_H[d_H$method == method & d_H$p == p_value, , drop = FALSE]
        if (!nrow(d_line)) next
        d_line <- d_line[order(d_line$n), , drop = FALSE]
        col <- method_colors[[method]]
        lty <- p_ltys[[as.character(p_value)]]
        lines(d_line$n, d_line$mean_flat_parameter_corr, col = col, lty = lty, lwd = 2)
        points(d_line$n, d_line$mean_flat_parameter_corr, col = col, pch = 16, cex = 0.85)
        has_sd <- is.finite(d_line$lower_2sd) & is.finite(d_line$upper_2sd)
        if (any(has_sd)) {
          arrows(
            x0 = d_line$n[has_sd],
            y0 = pmax(-1, d_line$lower_2sd[has_sd]),
            x1 = d_line$n[has_sd],
            y1 = pmin(1, d_line$upper_2sd[has_sd]),
            angle = 90,
            code = 3,
            length = 0.04,
            col = col,
            lwd = 1.3
          )
        }
      }
    }
  }

  legend_methods <- unique(flat_summary$method)
  legend(
    "bottom",
    inset = c(0, -0.16),
    xpd = NA,
    horiz = TRUE,
    legend = c(method_labels[legend_methods], paste0("p=", p_values)),
    col = c(method_colors[legend_methods], rep("#333333", length(p_values))),
    lty = c(rep(1, length(legend_methods)), p_ltys[as.character(p_values)]),
    pch = c(rep(16, length(legend_methods)), rep(NA, length(p_values))),
    lwd = c(rep(2, length(legend_methods)), rep(1.5, length(p_values))),
    bty = "n",
    cex = 0.85
  )
  mtext("Flattened posterior-mean parameter correlation: mean +/- 2 sd across reps", outer = TRUE, cex = 1.05)
  invisible(TRUE)
}

flat_parameter_summary <- summarize_flat_parameter_reps(results)
write.csv(
  flat_parameter_summary,
  file.path(out_dir, "flat_parameter_correlation_summary.csv"),
  row.names = FALSE
)

plot_metric_bars(summary, file.path(out_dir, "comparison_metric_bars.png"))
plot_sample_size_correlation_lines(summary, file.path(out_dir, "sample_size_correlation_lines.png"))
plot_flat_parameter_mean_sd(
  flat_parameter_summary,
  file.path(out_dir, "flat_parameter_correlation_mean_2sd.png")
)

cat("\nComparison summary:\n")
print(summary)
cat("\nWrote results to:\n", out_dir, "\n")
