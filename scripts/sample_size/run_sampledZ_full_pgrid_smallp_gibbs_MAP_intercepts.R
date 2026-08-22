#!/usr/bin/env Rscript

# Full sampled-Z sample-size study.
#
# This launcher reruns the main simulation after changing pretraining to sample
# augmented probit Z values by default.  It covers the full
# (n, p, H, G, loading design, block-size design) grid, but only runs the
# joint-mixture Gibbs baseline for the two smallest p values.  The product
# mixture estimator is run everywhere.
#
# Default grid:
#   n in {100, 500, 1000, 2000}
#   p in {250, 500, 1000, 2000}
#   H in {3, 4}
#   G in {2, 3}
#   loading designs: Sparse and Cross
#   block-size designs: balanced and moderately IFEval-like unbalanced
#
# Outputs are written separately by block-size design:
#   results/full/sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts
#   results/full/sampledZ_pgrid_unbalanced_crossloading_smallp_gibbs_MAP_intercepts
#
# Parallelization:
#   The product-mixture estimator uses within-fit parallelism for item
#   regressions and factor-score updates.  The joint Gibbs baseline is kept
#   serial by default because it was slower in the focused timing runs.  For
#   the large-p phase, where Gibbs is disabled, LARGE_P_PARALLEL_OURS and
#   LARGE_P_PARALLEL_WORKERS can be set separately from the small-p phase.

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
  status <- system(cmd)
  if (!identical(status, 0L)) {
    stop("Simulation driver failed with status ", status)
  }
}

block_modes <- split_csv(get_env("BLOCK_SIZE_MODES", "balanced,moderate_ifeval_like"))
n_values <- as.integer(split_csv(get_env("N_VALUES", "100,500,1000,2000")))
p_small <- as.integer(split_csv(get_env("P_VALUES_GIBBS", "250,500")))
p_large <- as.integer(split_csv(get_env("P_VALUES_OURS_ONLY", "1000,2000")))
small_p_parallel_ours <- get_env("PARALLEL_OURS", "TRUE")
large_p_parallel_ours <- get_env("LARGE_P_PARALLEL_OURS", small_p_parallel_ours)
small_p_parallel_workers <- get_env("PARALLEL_WORKERS", "18")
large_p_parallel_workers <- get_env("LARGE_P_PARALLEL_WORKERS", small_p_parallel_workers)

common_env <- c(
  INTERCEPT_MODE = "ifeval_like",
  INTERCEPT_SD = "0.45",
  INTERCEPT_BLOCK_SPAN = "1.6",
  INTERCEPT_CLIP = "1.75",
  MIXTURE_PARAM_MODE = "equal",
  MIXTURE_VARIANCE_MODE = "unequal",
  H_VALUES = get_env("H_VALUES", "3,4"),
  G_VALUES = get_env("G_VALUES", "2,3"),
  LOADING_DESIGNS = get_env(
    "LOADING_DESIGNS",
    paste(
      "balanced_moderate_few_positive_cross",
      "balanced_moderate_dense_signed_cross",
      sep = ","
    )
  ),
  REP_VALUES = get_env("REP_VALUES", paste(seq_len(25L), collapse = ",")),
  SEPARATIONS = "1.0",
  PRETRAIN_Z_UPDATE = "sample",
  PRETRAIN_AUG_ITER = get_env("PRETRAIN_AUG_ITER", "10"),
  CENTER_Z_FOR_SVD = "TRUE",
  PRETRAIN_MIN_AUG_ITER = "3",
  PRETRAIN_OBJECTIVE = "full_data_loglik",
  PRETRAIN_OBJECTIVE_TOLERANCE = "1e-3",
  PRETRAIN_OBJECTIVE_PATIENCE = "2",
  PRETRAIN_RETURN_BEST_ITERATION = "TRUE",
  REFINE_ITER = get_env("REFINE_ITER", "10"),
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
  RUN_OURS = "TRUE",
  PARALLEL_GIBBS = get_env("PARALLEL_GIBBS", "FALSE"),
  MFA_ITER = get_env("MFA_ITER", "2000"),
  MFA_BURN = get_env("MFA_BURN", "1000"),
  MFA_THIN = "1",
  MFA_STOP_ON_STABILITY = "FALSE",
  MFA_NORMALIZE_SCALE = "TRUE",
  MFA_VERBOSE = "FALSE",
  RESUME_EXISTING = get_env("RESUME_EXISTING", "TRUE")
)

block_labels <- c(
  balanced = "balanced",
  moderate_ifeval_like = "unbalanced"
)

for (block_mode in block_modes) {
  block_label <- block_labels[[block_mode]]
  if (is.null(block_label)) block_label <- block_mode
  out_dir <- file.path(
    repo_root,
    "results",
    "full",
    sprintf("sampledZ_pgrid_%s_crossloading_smallp_gibbs_MAP_intercepts", block_label)
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  base_env <- c(
    OUT_DIR = out_dir,
    BLOCK_SIZE_MODE = block_mode,
    common_env
  )

  cat("\nBlock mode:", block_mode, "\n")
  cat("Output:", out_dir, "\n")

  if (length(p_small) > 0L) {
    cat("Running small-p grid with product mixture and joint Gibbs. p = ",
        paste(p_small, collapse = ", "), "\n", sep = "")
    cat("Product-mixture parallelism: ",
        small_p_parallel_ours, " with workers=", small_p_parallel_workers, "\n", sep = "")
    run_driver(c(
      base_env,
      PARALLEL_OURS = small_p_parallel_ours,
      PARALLEL_WORKERS = small_p_parallel_workers,
      NP_GRID = np_grid_from_values(n_values, p_small),
      RUN_JOINT_MFA = "TRUE"
    ))
  }

  if (length(p_large) > 0L) {
    cat("Running large-p grid with product mixture only. p = ",
        paste(p_large, collapse = ", "), "\n", sep = "")
    cat("Large-p product-mixture parallelism: ",
        large_p_parallel_ours, " with workers=", large_p_parallel_workers, "\n", sep = "")
    run_driver(c(
      base_env,
      PARALLEL_OURS = large_p_parallel_ours,
      PARALLEL_WORKERS = large_p_parallel_workers,
      NP_GRID = np_grid_from_values(n_values, p_large),
      RUN_JOINT_MFA = "FALSE"
    ))
  }
}

cat("\nCompleted sampled-Z full p-grid launcher.\n")
