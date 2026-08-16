#!/usr/bin/env Rscript

# Focused check for two practical claims:
#   1. the proposed product-mixture probit factor estimator benefits from
#      within-fit parallelization;
#   2. parameter recovery remains accurate when one factor coordinate is
#      standard Gaussian, i.e. one coordinate has G_h = 1.
#
# The script runs matched serial and parallel fits of the proposed method only.
# The joint-mixture Gibbs sampler is intentionally disabled because the focused
# timing benchmark showed that its current parallel path is slower at these
# sizes.

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
launcher_script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
launcher_repo_root <- normalizePath(file.path(launcher_script_dir, "..", ".."))

source_isolated <- function(path) {
  source(path, local = new.env(parent = globalenv()))
}

run_one_mode <- function(mode, parallel_enabled) {
  Sys.setenv(
    OUT_DIR = file.path(
      launcher_repo_root,
      "results",
      "full",
      paste0("gaussian_coordinate_", mode, "_ours_H4_G2221_G3331_n1000_p500_25reps")
    ),
    RUN_OURS = "TRUE",
    RUN_JOINT_MFA = "FALSE",
    INTERCEPT_MODE = "ifeval_like",
    INTERCEPT_SD = "0.45",
    INTERCEPT_BLOCK_SPAN = "1.6",
    INTERCEPT_CLIP = "1.75",
    MIXTURE_PARAM_MODE = "equal",
    MIXTURE_VARIANCE_MODE = "unequal",
    H_VALUES = "4",
    G_CONFIGS = paste("2,2,2,1", "3,3,3,1", sep = ";"),
    NP_GRID = "n1000p500:1000:500",
    LOADING_DESIGNS = "balanced_moderate_few_positive_cross",
    REP_VALUES = paste(seq_len(25L), collapse = ","),
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
    PARALLEL_OURS = if (isTRUE(parallel_enabled)) "TRUE" else "FALSE",
    PARALLEL_GIBBS = "FALSE",
    PARALLEL_WORKERS = Sys.getenv("PARALLEL_WORKERS", unset = "18"),
    RESUME_EXISTING = "TRUE"
  )

  source_isolated(file.path(launcher_script_dir, "compare_original_simulation_joint_mfa_gibbs.R"))
  source_isolated(file.path(launcher_script_dir, "plot_sample_size_rmse_panels.R"))
  source_isolated(file.path(launcher_script_dir, "plot_sample_size_timing_lines.R"))
}

run_one_mode("parallel18", TRUE)
run_one_mode("serial", FALSE)
