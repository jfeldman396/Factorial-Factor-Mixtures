#!/usr/bin/env Rscript

# Final sample-size simulation launcher for the product-mixture probit factor
# model paper experiments.
#
# The launcher deliberately keeps the scientific design in one readable place.
# It calls compare_original_simulation_joint_mfa_gibbs.R in three phases, all
# with the same seeds and scenario labels:
#
#   1. Product MAP for the full p grid.
#   2. Viroli-style independent-mixture Gibbs with a Laplace loading prior for
#      the smaller p grid.
#   3. Viroli-style independent-mixture Gibbs with a diffuse Gaussian loading
#      prior for the smaller p grid.
#
# The shared checkpoint file lets the run resume cleanly and ensures that all
# methods are compared against the same simulated data in a given scenario.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
driver <- file.path(script_dir, "compare_original_simulation_joint_mfa_gibbs.R")

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

split_csv <- function(x) {
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

np_grid_from_values <- function(n_values, p_values) {
  pieces <- character(0)
  for (p in p_values) {
    for (n in n_values) {
      pieces <- c(pieces, sprintf("n%dp%d:%d:%d", n, p, n, p))
    }
  }
  paste(pieces, collapse = ",")
}

run_driver <- function(env) {
  env_vec <- sprintf("%s=%s", names(env), shQuote(env, type = "sh"))
  cmd <- paste(
    "env",
    paste(env_vec, collapse = " "),
    shQuote(file.path(R.home("bin"), "Rscript"), type = "sh"),
    shQuote(driver, type = "sh")
  )
  cat("\nRunning driver phase with settings:\n")
  cat(paste(sprintf("  %s=%s", names(env), env), collapse = "\n"), "\n")
  status <- system(cmd)
  if (!identical(status, 0L)) stop("Simulation driver failed with status ", status)
}

override_env <- function(base, updates) {
  out <- base
  for (nm in names(updates)) out[[nm]] <- updates[[nm]]
  out
}

n_values <- as.integer(split_csv(get_env("N_VALUES", "100,200,300,400")))
p_values_product <- as.integer(split_csv(get_env("P_VALUES_PRODUCT", "500,1000,2000,4000")))
p_values_viroli <- as.integer(split_csv(get_env("P_VALUES_VIROLI", "500,1000")))
block_modes <- split_csv(get_env("BLOCK_SIZE_MODES", "moderate_ifeval_like"))
parallel_workers <- get_env("PARALLEL_WORKERS", "18")
run_label <- get_env("RUN_LABEL", "final_cross_ifeval_product_viroli")
out_dir <- get_env(
  "OUT_DIR",
  file.path(repo_root, "results", "full", run_label)
)

common_env <- c(
  OUT_DIR = out_dir,
  SEED = get_env("SEED", "20260731"),
  RESUME_EXISTING = get_env("RESUME_EXISTING", "TRUE"),
  REP_VALUES = get_env("REP_VALUES", paste(seq_len(25L), collapse = ",")),
  H_VALUES = get_env("H_VALUES", "5,10,15,20"),
  G_VALUES = get_env("G_VALUES", "2,3"),
  SEPARATIONS = get_env("SEPARATIONS", "1,2"),
  NP_GRID = np_grid_from_values(n_values, p_values_product),
  BLOCK_SIZE_MODE = block_modes[1L],
  LOADING_DESIGNS = get_env("LOADING_DESIGNS", "balanced_moderate_dense_signed_cross"),
  LOADING_SIGN_MODE = get_env("LOADING_SIGN_MODE", "block"),
  ALIGNMENT_MODE = get_env("ALIGNMENT_MODE", "loadings"),
  MIXTURE_PARAM_MODE = get_env("MIXTURE_PARAM_MODE", "viroli_smoke"),
  MIXTURE_VARIANCE_MODE = get_env("MIXTURE_VARIANCE_MODE", "unequal"),
  INTERCEPT_MODE = get_env("INTERCEPT_MODE", "ifeval_like"),
  INTERCEPT_SD = get_env("INTERCEPT_SD", "0.45"),
  INTERCEPT_BLOCK_SPAN = get_env("INTERCEPT_BLOCK_SPAN", "1.6"),
  INTERCEPT_CLIP = get_env("INTERCEPT_CLIP", "1.75"),
  OURS_PRETRAINING_METHOD = "em_svd",
  EM_SVD_INIT = get_env("EM_SVD_INIT", "both"),
  EM_SVD_INIT_Z = get_env("EM_SVD_INIT_Z", "expectation"),
  EM_SVD_ITER = get_env("EM_SVD_ITER", "50"),
  EM_SVD_TOL_LOGLIK = get_env("EM_SVD_TOL_LOGLIK", "1e-5"),
  EM_SVD_TOL_L = get_env("EM_SVD_TOL_L", "1e-4"),
  EM_SVD_RANDOM_STARTS = get_env("EM_SVD_RANDOM_STARTS", "0"),
  ROTATION_ITER = get_env("ROTATION_ITER", "20"),
  ROTATION_MIN_OUTER = get_env("ROTATION_MIN_OUTER", "2"),
  ROTATION_OBJECTIVE_TOLERANCE = get_env("ROTATION_OBJECTIVE_TOLERANCE", "1e-3"),
  ROTATION_REQUIRE_MIXTURE_CONVERGENCE = get_env("ROTATION_REQUIRE_MIXTURE_CONVERGENCE", "TRUE"),
  ROTATION_RANDOM_STARTS = get_env("ROTATION_RANDOM_STARTS", "1"),
  ROTATION_N_MIX_STARTS = get_env("ROTATION_N_MIX_STARTS", "3"),
  ROTATION_GRID_SIZE = get_env("ROTATION_GRID_SIZE", "21"),
  ROTATION_SWEEP = get_env("ROTATION_SWEEP", "full"),
  PRETRAIN_LOADING_PENALTY = get_env("PRETRAIN_LOADING_PENALTY", "10"),
  ROTATION_LOADING_L1_PENALTY = get_env("ROTATION_LOADING_L1_PENALTY", "10"),
  REFINE_ITER = get_env("REFINE_ITER", "50"),
  REFINE_MIN_ITER = get_env("REFINE_MIN_ITER", "3"),
  REFINE_OBJECTIVE_TOLERANCE = get_env("REFINE_OBJECTIVE_TOLERANCE", "1e-3"),
  REFINE_REQUIRE_MIXTURE_CONVERGENCE = get_env("REFINE_REQUIRE_MIXTURE_CONVERGENCE", "TRUE"),
  REFINE_RETURN_BEST_ITERATION = get_env("REFINE_RETURN_BEST_ITERATION", "TRUE"),
  REFINE_SELECTION_OBJECTIVE = get_env("REFINE_SELECTION_OBJECTIVE", "posterior_objective"),
  REFINE_ENFORCE_MONOTONE = get_env("REFINE_ENFORCE_MONOTONE", "TRUE"),
  LAMBDA_L1_PENALTY = get_env("LAMBDA_L1_PENALTY", "10"),
  LASSO_BACKEND = get_env("LASSO_BACKEND", "glmnet"),
  GLMNET_STANDARDIZE = get_env("GLMNET_STANDARDIZE", "FALSE"),
  FACTOR_UPDATE = get_env("FACTOR_UPDATE", "marginal"),
  MIXTURE_UPDATE = get_env("MIXTURE_UPDATE", "map"),
  MIXTURE_REFIT = get_env("MIXTURE_REFIT", "em"),
  MIXTURE_MAX_ITER = get_env("MIXTURE_MAX_ITER", "100"),
  MU_PRIOR_MEAN = get_env("MU_PRIOR_MEAN", "0"),
  MU_PRIOR_KAPPA = get_env("MU_PRIOR_KAPPA", "0.05"),
  VAR_PRIOR_SHAPE = get_env("VAR_PRIOR_SHAPE", "3"),
  VAR_PRIOR_SCALE = get_env("VAR_PRIOR_SCALE", "2"),
  WEIGHT_PRIOR_ALPHA = get_env("WEIGHT_PRIOR_ALPHA", "1"),
  REFINE_MU_PRIOR_MEAN = get_env("REFINE_MU_PRIOR_MEAN", "0"),
  REFINE_MU_PRIOR_KAPPA = get_env("REFINE_MU_PRIOR_KAPPA", "0.05"),
  REFINE_VAR_PRIOR_SHAPE = get_env("REFINE_VAR_PRIOR_SHAPE", "3"),
  REFINE_VAR_PRIOR_SCALE = get_env("REFINE_VAR_PRIOR_SCALE", "2"),
  REFINE_WEIGHT_PRIOR_ALPHA = get_env("REFINE_WEIGHT_PRIOR_ALPHA", "1"),
  MIN_MIXTURE_VAR = get_env("MIN_MIXTURE_VAR", "0.05"),
  PARALLEL_OURS = get_env("PARALLEL_OURS", "TRUE"),
  PARALLEL_GIBBS = get_env("PARALLEL_GIBBS", "TRUE"),
  PARALLEL_WORKERS = parallel_workers,
  MAX_JOINT_PARAMETER_K = get_env("MAX_JOINT_PARAMETER_K", "5000"),
  MAX_JOINT_PROFILE_ARI_K = get_env("MAX_JOINT_PROFILE_ARI_K", "5000"),
  WRITE_PARAMETER_TABLES = get_env("WRITE_PARAMETER_TABLES", "TRUE"),
  WRITE_ITERATION_HISTORIES = get_env("WRITE_ITERATION_HISTORIES", "FALSE"),
  MFA_ITER = get_env("MFA_ITER", "2000"),
  MFA_BURN = get_env("MFA_BURN", "1000"),
  MFA_THIN = get_env("MFA_THIN", "1"),
  VIROLI_ITER = get_env("VIROLI_ITER", "2000"),
  VIROLI_BURN = get_env("VIROLI_BURN", "1000"),
  VIROLI_THIN = get_env("VIROLI_THIN", "1"),
  VIROLI_COMPUTE_PARAMETER_ESS = get_env("VIROLI_COMPUTE_PARAMETER_ESS", "TRUE"),
  VIROLI_NORMALIZE_EACH_DRAW = get_env("VIROLI_NORMALIZE_EACH_DRAW", "TRUE"),
  VIROLI_MIN_SCALE = get_env("VIROLI_MIN_SCALE", "1e-4"),
  VIROLI_VERBOSE = get_env("VIROLI_VERBOSE", "FALSE")
)

cat("Final product/Viroli simulation launcher\n")
cat("Output directory:", out_dir, "\n")
cat("n grid:", paste(n_values, collapse = ", "), "\n")
cat("product p grid:", paste(p_values_product, collapse = ", "), "\n")
cat("Viroli p grid:", paste(p_values_viroli, collapse = ", "), "\n")
cat("H grid:", common_env[["H_VALUES"]], "\n")
cat("G grid:", common_env[["G_VALUES"]], "\n")
cat("separation grid:", common_env[["SEPARATIONS"]], "\n")
cat("loading design:", common_env[["LOADING_DESIGNS"]], "\n")
cat("block mode(s):", paste(block_modes, collapse = ", "), "\n")

for (block_mode in block_modes) {
  phase_env <- override_env(common_env, c(BLOCK_SIZE_MODE = block_mode))
  run_driver(override_env(
    phase_env,
    c(
    NP_GRID = np_grid_from_values(n_values, p_values_product),
    RUN_OURS = "TRUE",
    RUN_JOINT_MFA = "FALSE",
    RUN_VIROLI = "FALSE"
    )
  ))

  if (length(p_values_viroli)) {
    run_driver(override_env(
      phase_env,
      c(
      NP_GRID = np_grid_from_values(n_values, p_values_viroli),
      RUN_OURS = "FALSE",
      RUN_JOINT_MFA = "FALSE",
      RUN_VIROLI = "TRUE",
      VIROLI_METHOD_NAME = "viroli_laplace_gibbs",
      VIROLI_LAMBDA_L1_PENALTY = get_env("VIROLI_LAPLACE_L1_PENALTY", "10")
      )
    ))
    run_driver(override_env(
      phase_env,
      c(
      NP_GRID = np_grid_from_values(n_values, p_values_viroli),
      RUN_OURS = "FALSE",
      RUN_JOINT_MFA = "FALSE",
      RUN_VIROLI = "TRUE",
      VIROLI_METHOD_NAME = "viroli_gaussian_gibbs",
      VIROLI_LAMBDA_L1_PENALTY = "0"
      )
    ))
  }
}

cat("\nFinal simulation launcher completed. Checkpoint:\n")
cat(file.path(out_dir, "comparison_results_checkpoint.csv"), "\n")
