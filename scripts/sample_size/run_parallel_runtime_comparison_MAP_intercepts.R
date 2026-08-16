#!/usr/bin/env Rscript

# Runtime comparison with a smaller grid run twice:
#   1. product-mixture and joint-Gibbs fits in serial,
#   2. product-mixture and joint-Gibbs fits with safe within-fit parallelism.
#
# The statistical model and simulated data-generating settings match the
# moderate cross-loading simulations; only the parallel flags differ.

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))

run_one_mode <- function(mode, parallel_enabled) {
  Sys.setenv(
    OUT_DIR = file.path(repo_root, "results", "full", paste0("parallel_runtime_", mode, "_p500_MAP_intercepts")),
    INTERCEPT_MODE = "ifeval_like",
    INTERCEPT_SD = "0.45",
    INTERCEPT_BLOCK_SPAN = "1.6",
    INTERCEPT_CLIP = "1.75",
    MIXTURE_PARAM_MODE = "equal",
    MIXTURE_VARIANCE_MODE = "unequal",
    H_VALUES = "3,4",
    G_CONFIGS = paste("3,3,1", "3,3,1,3", sep = ";"),
    NP_GRID = paste(
      "n100p500:100:500",
      "n500p500:500:500",
      "n1000p500:1000:500",
      sep = ","
    ),
    LOADING_DESIGNS = paste(
      "balanced_moderate_few_positive_cross",
      "balanced_moderate_dense_signed_cross",
      sep = ","
    ),
    REP_VALUES = paste(seq_len(5L), collapse = ","),
    SEPARATIONS = "1.0",
    PRETRAIN_AUG_ITER = "10",
    CENTER_Z_FOR_SVD = "TRUE",
    PRETRAIN_MIN_AUG_ITER = "3",
    PRETRAIN_OBJECTIVE = "full_data_loglik",
    PRETRAIN_OBJECTIVE_TOLERANCE = "1e-3",
    PRETRAIN_OBJECTIVE_PATIENCE = "2",
    PRETRAIN_RETURN_BEST_ITERATION = "TRUE",
    REFINE_ITER = "10",
    REFINE_MIN_ITER = "3",
    REFINE_OBJECTIVE_TOLERANCE = "1e-3",
    REFINE_STOPPING_OBJECTIVE = "posterior_objective",
    REFINE_RETURN_BEST_ITERATION = "TRUE",
    REFINE_SELECTION_OBJECTIVE = "posterior_objective",
    MIXTURE_UPDATE = "map",
    MU_PRIOR_MEAN = "0",
    MU_PRIOR_KAPPA = "0.01",
    VAR_PRIOR_SHAPE = "2",
    VAR_PRIOR_SCALE = "1.5",
    WEIGHT_PRIOR_ALPHA = "1",
    MFA_ITER = "2000",
    MFA_BURN = "1000",
    MFA_THIN = "1",
    MFA_STOP_ON_STABILITY = "FALSE",
    MFA_NORMALIZE_SCALE = "TRUE",
    MFA_VERBOSE = "FALSE",
    PARALLEL_OURS = if (isTRUE(parallel_enabled)) "TRUE" else "FALSE",
    PARALLEL_GIBBS = if (isTRUE(parallel_enabled)) "TRUE" else "FALSE",
    PARALLEL_WORKERS = Sys.getenv("PARALLEL_WORKERS", unset = "4"),
    RESUME_EXISTING = "TRUE"
  )

  source(file.path(script_dir, "compare_original_simulation_joint_mfa_gibbs.R"), local = TRUE)
  source(file.path(script_dir, "plot_sample_size_rmse_panels.R"), local = TRUE)
  source(file.path(script_dir, "plot_sample_size_timing_lines.R"), local = TRUE)
}

run_one_mode("serial", FALSE)
run_one_mode("parallel", TRUE)
