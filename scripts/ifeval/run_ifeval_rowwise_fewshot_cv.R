#!/usr/bin/env Rscript

# Command-line entry point for few-shot row-wise P-IFA cross-validation.
# Configuration is supplied through environment variables so runs are easy to
# reproduce on a workstation or scheduler.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  script_path <- sub("^--file=", "", file_arg[1L])
  script_path <- gsub("~\\+~", " ", script_path)
  dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source(file.path(repo_root, "R", "pifa_rowwise_cv.R"))

parse_int_grid_local <- function(x, default) {
  if (!nzchar(x)) return(default)
  x <- gsub("[[:space:]]+", "", x)
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    bounds <- as.integer(strsplit(x, ":", fixed = TRUE)[[1L]])
    return(seq(bounds[1L], bounds[2L]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
}

parse_num_grid_local <- function(x, default) {
  if (!nzchar(x)) return(default)
  as.numeric(strsplit(gsub("[[:space:]]+", "", x), ",", fixed = TRUE)[[1L]])
}

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, if (default) "TRUE" else "FALSE")
  tolower(value) %in% c("true", "1", "yes")
}

read_binary_matrix_with_missing <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- raw[[1L]]
  usable <- vapply(seq_len(ncol(X)), function(j) {
    observed <- X[!is.na(X[, j]), j]
    length(observed) >= 2L && all(observed %in% c(0, 1)) && length(unique(observed)) == 2L
  }, logical(1L))
  X[, usable, drop = FALSE]
}

matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  file.path(repo_root, "data", "ifeval", "openeval_ifeval_only_binary_matrix.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(repo_root, "results", "full", "ifeval_rowwise_fewshot_cv")
)
X <- read_binary_matrix_with_missing(matrix_path)

result <- cv_pifa_rowwise(
  X = X,
  H_grid = parse_int_grid_local(Sys.getenv("H_GRID", "2:8"), 2:8),
  G_values = parse_int_grid_local(Sys.getenv("G_VALUES", "1,2,3"), 1:3),
  lambda_grid = parse_num_grid_local(
    Sys.getenv("LAMBDA_L1_GRID", "0,1,2,4,8,12"),
    c(0, 1, 2, 4, 8, 12)
  ),
  K = as.integer(Sys.getenv("K_FOLDS", "5")),
  reveal_fraction = as.numeric(Sys.getenv("REVEAL_FRACTION", "0.20")),
  n_reveal_repeats = as.integer(Sys.getenv("N_REVEAL_REPEATS", "3")),
  seed = as.integer(Sys.getenv("SEED", "1")),
  parallel = env_flag("PARALLEL", TRUE),
  workers = as.integer(Sys.getenv("WORKERS", "18")),
  fit_control = list(
    n_aug_iter = as.integer(Sys.getenv("PRETRAIN_AUG_ITER", "200")),
    n_refine_iter = as.integer(Sys.getenv("REFINE_ITER", "100")),
    n_ica_starts = as.integer(Sys.getenv("ROTATION_ICA_STARTS", "1")),
    max_outer = as.integer(Sys.getenv("MAX_OUTER", "20")),
    n_mix_starts = as.integer(Sys.getenv("N_MIX_STARTS", "2")),
    mixture_max_iter = as.integer(Sys.getenv("MIXTURE_MAX_ITER", "200")),
    refinement_objective_tolerance = as.numeric(
      Sys.getenv("REFINE_OBJECTIVE_TOLERANCE", "1e-3")
    ),
    refinement_require_mixture_convergence = env_flag(
      "REFINE_REQUIRE_MIXTURE_CONVERGENCE", TRUE
    ),
    factor_score_bound = as.numeric(Sys.getenv("FACTOR_SCORE_BOUND", "5")),
    heldout_maxit = as.integer(Sys.getenv("HELDOUT_FACTOR_MAXIT", "100")),
    heldout_max_component_starts = as.integer(
      Sys.getenv("HELDOUT_MAX_COMPONENT_STARTS", "2")
    )
  ),
  fit_selected_full_data = env_flag("FIT_SELECTED_FULL_DATA", TRUE),
  save_training_fits = env_flag("SAVE_TRAINING_FITS", FALSE),
  resume = env_flag("RESUME_EXISTING", TRUE),
  task_batch_size = as.integer(Sys.getenv("TASK_BATCH_SIZE", Sys.getenv("WORKERS", "18"))),
  out_dir = out_dir
)

print(result)
cat("Outputs saved in: ", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
