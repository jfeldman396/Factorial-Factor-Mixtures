#!/usr/bin/env Rscript

# Benchmark Product MAP inner parallelization.  This intentionally runs through
# the main simulation driver so the benchmark uses the same DGP, algorithm
# settings, convergence checks, and evaluation metrics as the simulation study.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else sys.frame(1)$ofile
repo_root <- normalizePath(file.path(dirname(normalizePath(this_file)), "..", ".."))
driver <- file.path(repo_root, "scripts", "sample_size", "compare_original_simulation_joint_mfa_gibbs.R")

get_env <- function(name, default, cast = identity) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  cast(value)
}

parse_int_vec <- function(x) {
  as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

bench_n <- get_env("BENCH_N", 100L, as.integer)
bench_p <- get_env("BENCH_P", 1000L, as.integer)
bench_H <- get_env("BENCH_H", 10L, as.integer)
bench_G <- get_env("BENCH_G", 2L, as.integer)
bench_rep <- get_env("BENCH_REP", 1L, as.integer)
bench_sep <- get_env("BENCH_SEP", 1, as.numeric)
bench_workers <- get_env("BENCH_WORKERS", c(1L, 4L, 8L, 18L), parse_int_vec)
bench_penalty <- get_env("BENCH_LAPLACE_PENALTY", 10, as.numeric)
bench_seed <- get_env("BENCH_SEED", 20260731L, as.integer)

bench_loading_design <- get_env("BENCH_LOADING_DESIGN", "balanced_moderate_dense_signed_cross", as.character)
bench_block_mode <- get_env("BENCH_BLOCK_SIZE_MODE", "ifeval_like", as.character)
bench_strength <- get_env("BENCH_LOADING_STRENGTH", "weak", as.character)
bench_cross_prob <- get_env("BENCH_CROSS_LOADING_PROB", 0.2, as.numeric)
bench_primary_range <- get_env("BENCH_PRIMARY_LOADING_RANGE", "1.25,1.75", as.character)
bench_cross_range <- get_env("BENCH_CROSS_LOADING_RANGE", "1.25,1.75", as.character)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
stable_dir <- file.path(repo_root, "results", "diagnostics", "parallel_worker_benchmarks")
out_base <- get_env(
  "BENCH_OUT_DIR",
  file.path(stable_dir, paste0("product_map_n", bench_n, "_p", bench_p, "_H", bench_H, "_G", bench_G, "_", run_id)),
  as.character
)
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

base_env <- c(
  SEED = as.character(bench_seed),
  RESUME_EXISTING = "FALSE",
  LOADING_DESIGNS = bench_loading_design,
  LOADING_SIGN_MODE = "block",
  CROSS_SIGN_MODE = "random",
  ALIGNMENT_MODE = "loadings",
  MIXTURE_PARAM_MODE = "viroli_smoke",
  MIXTURE_VARIANCE_MODE = "unequal",
  INTERCEPT_MODE = "ifeval_like",
  INTERCEPT_SD = "0.45",
  INTERCEPT_BLOCK_SPAN = "1.6",
  INTERCEPT_CLIP = "1.75",
  SEPARATIONS = as.character(bench_sep),
  OURS_PRETRAINING_METHOD = "em_svd",
  EM_SVD_INIT = "both",
  EM_SVD_INIT_Z = "expectation",
  EM_SVD_ITER = "50",
  EM_SVD_TOL_LOGLIK = "1e-5",
  EM_SVD_TOL_L = "1e-4",
  EM_SVD_RANDOM_STARTS = "0",
  PRETRAIN_LOADING_PENALTY = as.character(bench_penalty),
  ROTATION_OPTIMIZER = "riemannian",
  ROTATION_ITER = "20",
  ROTATION_MIN_OUTER = "2",
  ROTATION_OBJECTIVE_TOLERANCE = "1e-3",
  ROTATION_REQUIRE_MIXTURE_CONVERGENCE = "TRUE",
  ROTATION_RANDOM_STARTS = "1",
  ROTATION_N_MIX_STARTS = "3",
  ROTATION_GRID_SIZE = "21",
  RIEMANNIAN_ROTATION_STEPS = "10",
  RIEMANNIAN_ETA0 = "1",
  RIEMANNIAN_BETA = "0.5",
  RIEMANNIAN_MIN_ETA = "1e-8",
  RIEMANNIAN_GRAD_TOL = "1e-6",
  RIEMANNIAN_UPDATE = "cayley",
  ROTATION_LOADING_L1_PENALTY = as.character(bench_penalty),
  REFINE_ITER = "50",
  REFINE_MIN_ITER = "3",
  REFINE_OBJECTIVE_TOLERANCE = "1e-3",
  REFINE_REQUIRE_MIXTURE_CONVERGENCE = "TRUE",
  REFINE_RETURN_BEST_ITERATION = "TRUE",
  REFINE_SELECTION_OBJECTIVE = "posterior_objective",
  REFINE_ENFORCE_MONOTONE = "TRUE",
  LAMBDA_L1_PENALTY = as.character(bench_penalty),
  LASSO_BACKEND = "glmnet",
  GLMNET_STANDARDIZE = "FALSE",
  FACTOR_UPDATE = "marginal",
  MIXTURE_UPDATE = "map",
  MIXTURE_REFIT = "em",
  MIXTURE_MAX_ITER = "100",
  MU_PRIOR_MEAN = "0",
  MU_PRIOR_KAPPA = "0.05",
  VAR_PRIOR_SHAPE = "3",
  VAR_PRIOR_SCALE = "2",
  WEIGHT_PRIOR_ALPHA = "1",
  REFINE_MU_PRIOR_MEAN = "0",
  REFINE_MU_PRIOR_KAPPA = "0.05",
  REFINE_VAR_PRIOR_SHAPE = "3",
  REFINE_VAR_PRIOR_SCALE = "2",
  REFINE_WEIGHT_PRIOR_ALPHA = "1",
  MIN_MIXTURE_VAR = "0.05",
  WRITE_PARAMETER_TABLES = "FALSE",
  WRITE_ITERATION_HISTORIES = "TRUE",
  MAX_JOINT_PARAMETER_K = "5000",
  MAX_JOINT_PROFILE_ARI_K = "5000",
  REP_VALUES = as.character(bench_rep),
  H_VALUES = as.character(bench_H),
  G_VALUES = as.character(bench_G),
  NP_GRID = sprintf("n%sp%s:%s:%s", bench_n, bench_p, bench_n, bench_p),
  BLOCK_SIZE_MODE = bench_block_mode,
  LOADING_STRENGTH = bench_strength,
  PRIMARY_LOADING_RANGE = bench_primary_range,
  CROSS_LOADING_RANGE = bench_cross_range,
  CROSS_LOADING_PROB = as.character(bench_cross_prob),
  RUN_OURS = "TRUE",
  RUN_JOINT_MFA = "FALSE",
  RUN_VIROLI = "FALSE",
  PARALLEL_GIBBS = "FALSE"
)

read_one_result <- function(case_dir, workers, parallel_enabled) {
  result_file <- file.path(case_dir, "comparison_results.csv")
  if (!file.exists(result_file)) {
    stop("Missing result file: ", result_file, call. = FALSE)
  }
  result <- read.csv(result_file, check.names = FALSE)
  result$benchmark_workers <- workers
  result$benchmark_parallel <- parallel_enabled

  history_file <- list.files(case_dir, pattern = "^ours_timing_history_.*[.]csv$", full.names = TRUE)
  if (!length(history_file)) {
    stop("Missing Product MAP timing history in: ", case_dir, call. = FALSE)
  }
  history <- read.csv(history_file[[1]], check.names = FALSE)
  history$benchmark_workers <- workers
  history$benchmark_parallel <- parallel_enabled

  stage_split <- split(history, history$stage)
  stage_timing <- do.call(rbind, lapply(names(stage_split), function(stage_name) {
    d <- stage_split[[stage_name]]
    data.frame(
      benchmark_workers = workers,
      benchmark_parallel = parallel_enabled,
      stage = stage_name,
      n_iter = nrow(d),
      mean_iteration_seconds = mean(d$iteration_seconds, na.rm = TRUE),
      median_iteration_seconds = median(d$iteration_seconds, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  list(result = result, stage_timing = stage_timing, history = history)
}

all_results <- list()
all_stage_timings <- list()
all_histories <- list()

cat("Product MAP inner parallel benchmark\n")
cat(sprintf(
  "Setting: n=%d, p=%d, H=%d, G=%d, sep=%s, block=%s, loadings=%s, strength=%s, cross_prob=%s\n",
  bench_n, bench_p, bench_H, bench_G, bench_sep, bench_block_mode,
  bench_loading_design, bench_strength, bench_cross_prob
))
cat(sprintf(
  "Laplace/loading penalty: %s, workers tested: %s\n",
  bench_penalty, paste(bench_workers, collapse = ", ")
))
cat("Output directory:\n", out_base, "\n", sep = "")

for (workers in bench_workers) {
  parallel_enabled <- workers > 1L
  label <- if (parallel_enabled) paste0("parallel_w", workers) else "serial"
  case_dir <- file.path(out_base, label)
  dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)

  run_env <- c(
    base_env,
    OUT_DIR = case_dir,
    PARALLEL_OURS = if (parallel_enabled) "TRUE" else "FALSE",
    PARALLEL_WORKERS = as.character(workers)
  )

  cat("\nRunning ", label, "...\n", sep = "")
  log_file <- file.path(case_dir, "run.log")
  old_env <- Sys.getenv(names(run_env), unset = NA_character_)
  names(old_env) <- names(run_env)
  do.call(Sys.setenv, as.list(run_env))
  status <- system2("Rscript", args = shQuote(driver), stdout = log_file, stderr = log_file)
  restore_env <- as.list(old_env[!is.na(old_env)])
  unset_env <- names(old_env)[is.na(old_env)]
  if (length(restore_env)) do.call(Sys.setenv, restore_env)
  if (length(unset_env)) Sys.unsetenv(unset_env)
  if (!identical(status, 0L)) {
    stop("Benchmark case failed: ", label, ". See ", log_file, call. = FALSE)
  }

  one <- read_one_result(case_dir, workers, parallel_enabled)
  all_results[[label]] <- one$result
  all_stage_timings[[label]] <- one$stage_timing
  all_histories[[label]] <- one$history
  cat(sprintf("%s complete: wall %.2fs\n", label, one$result$seconds[[1]]))
}

results <- do.call(rbind, all_results)
stage_timings <- do.call(rbind, all_stage_timings)
histories <- do.call(rbind, all_histories)

result_stub <- sprintf("product_map_inner_parallel_n%d_p%d_H%d_G%d", bench_n, bench_p, bench_H, bench_G)
results_file <- file.path(out_base, paste0(result_stub, "_results.csv"))
stage_timing_file <- file.path(out_base, paste0(result_stub, "_stage_timing.csv"))
history_file <- file.path(out_base, paste0(result_stub, "_all_histories.csv"))
write.csv(results, results_file, row.names = FALSE)
write.csv(stage_timings, stage_timing_file, row.names = FALSE)
write.csv(histories, history_file, row.names = FALSE)

dir.create(stable_dir, recursive = TRUE, showWarnings = FALSE)
stable_results_file <- file.path(stable_dir, paste0(result_stub, "_latest_results.csv"))
stable_stage_timing_file <- file.path(stable_dir, paste0(result_stub, "_latest_stage_timing.csv"))
stable_history_file <- file.path(stable_dir, paste0(result_stub, "_latest_all_histories.csv"))
write.csv(results, stable_results_file, row.names = FALSE)
write.csv(stage_timings, stable_stage_timing_file, row.names = FALSE)
write.csv(histories, stable_history_file, row.names = FALSE)

plot_file <- file.path(out_base, paste0(result_stub, "_timing.png"))
stable_plot_file <- file.path(stable_dir, paste0(result_stub, "_latest_timing.png"))
labels <- ifelse(results$benchmark_parallel, paste0(results$benchmark_workers, " workers"), "serial")
png(plot_file, width = 1700, height = 850, res = 150)
op <- par(mfrow = c(1, 2), mar = c(7, 5, 4, 1), oma = c(0, 0, 2, 0))
barplot(
  results$seconds,
  names.arg = labels,
  las = 2,
  col = "#377eb8",
  ylab = "seconds",
  main = "End-to-end"
)
stage_wide <- reshape(
  stage_timings[, c("benchmark_workers", "stage", "median_iteration_seconds")],
  idvar = "benchmark_workers",
  timevar = "stage",
  direction = "wide"
)
stage_mat <- as.matrix(stage_wide[, grep("^median_iteration_seconds[.]", names(stage_wide)), drop = FALSE])
rownames(stage_mat) <- ifelse(stage_wide$benchmark_workers > 1L, paste0(stage_wide$benchmark_workers, " workers"), "serial")
colnames(stage_mat) <- sub("^median_iteration_seconds[.]", "", colnames(stage_mat))
barplot(
  t(stage_mat),
  beside = TRUE,
  las = 2,
  col = c("#4daf4a", "#984ea3", "#ff7f00")[seq_len(ncol(stage_mat))],
  ylab = "median seconds",
  main = "Median iteration time by stage"
)
legend("topright", legend = colnames(stage_mat), fill = c("#4daf4a", "#984ea3", "#ff7f00")[seq_len(ncol(stage_mat))], bty = "n", cex = 0.8)
mtext(
  sprintf("Product MAP inner parallel benchmark: n=%d, p=%d, H=%d, G=%d", bench_n, bench_p, bench_H, bench_G),
  outer = TRUE,
  cex = 1.05,
  font = 2
)
par(op)
dev.off()
file.copy(plot_file, stable_plot_file, overwrite = TRUE)

cat("\nSaved benchmark summary:\n")
cat(results_file, "\n", sep = "")
cat(stage_timing_file, "\n", sep = "")
cat(history_file, "\n", sep = "")
cat(stable_results_file, "\n", sep = "")
cat(stable_stage_timing_file, "\n", sep = "")
cat(stable_history_file, "\n", sep = "")
cat(stable_plot_file, "\n", sep = "")

display_cols <- c(
  "benchmark_workers", "benchmark_parallel", "seconds",
  "factor_score_rmse", "lambda_rmse", "alpha_rmse", "marginal_mu_rmse",
  "marginal_var_rmse", "marginal_weight_rmse", "converged",
  "em_svd_completed_iter", "refinement_completed_iter"
)
print(results[, intersect(display_cols, names(results)), drop = FALSE], row.names = FALSE)

print(stage_timings, row.names = FALSE)
