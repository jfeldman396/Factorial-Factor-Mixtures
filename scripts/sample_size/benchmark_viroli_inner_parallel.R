#!/usr/bin/env Rscript

# Benchmark inner-chain parallelization for the Viroli-style probit Gibbs
# sampler.  This script intentionally calls the main simulation driver so the
# benchmark uses the same DGP, evaluation code, priors, normalization, and ESS
# path as the simulation study.

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
bench_p <- get_env("BENCH_P", 500L, as.integer)
bench_H <- get_env("BENCH_H", 10L, as.integer)
bench_G <- get_env("BENCH_G", 2L, as.integer)
bench_rep <- get_env("BENCH_REP", 1L, as.integer)
bench_sep <- get_env("BENCH_SEP", 1, as.numeric)
bench_workers <- get_env("BENCH_WORKERS", c(1L, 4L, 8L, 18L), parse_int_vec)
bench_iter <- get_env("BENCH_ITER", 2000L, as.integer)
bench_burn <- get_env("BENCH_BURN", 1000L, as.integer)
bench_penalty <- get_env("BENCH_LAPLACE_PENALTY", 10, as.numeric)
bench_seed <- get_env("BENCH_SEED", 20260731L, as.integer)

bench_loading_design <- get_env("BENCH_LOADING_DESIGN", "balanced_moderate_dense_signed_cross", as.character)
bench_block_mode <- get_env("BENCH_BLOCK_SIZE_MODE", "ifeval_like", as.character)
bench_strength <- get_env("BENCH_LOADING_STRENGTH", "weak", as.character)
bench_cross_prob <- get_env("BENCH_CROSS_LOADING_PROB", 0.2, as.numeric)
bench_primary_range <- get_env("BENCH_PRIMARY_LOADING_RANGE", "1.25,1.75", as.character)
bench_cross_range <- get_env("BENCH_CROSS_LOADING_RANGE", "1.25,1.75", as.character)

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_base <- get_env(
  "BENCH_OUT_DIR",
  file.path(repo_root, "results", "diagnostics", "viroli_inner_parallel_benchmark", run_id),
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
  WRITE_PARAMETER_TABLES = "FALSE",
  WRITE_ITERATION_HISTORIES = "TRUE",
  MAX_JOINT_PARAMETER_K = "5000",
  MAX_JOINT_PROFILE_ARI_K = "5000",
  MFA_ITER = as.character(bench_iter),
  MFA_BURN = as.character(bench_burn),
  MFA_THIN = "1",
  VIROLI_ITER = as.character(bench_iter),
  VIROLI_BURN = as.character(bench_burn),
  VIROLI_THIN = "1",
  VIROLI_COMPUTE_PARAMETER_ESS = "TRUE",
  VIROLI_NORMALIZE_EACH_DRAW = "TRUE",
  VIROLI_MIN_SCALE = "1e-4",
  VIROLI_VERBOSE = "FALSE",
  REP_VALUES = as.character(bench_rep),
  H_VALUES = as.character(bench_H),
  G_VALUES = as.character(bench_G),
  NP_GRID = sprintf("n%sp%s:%s:%s", bench_n, bench_p, bench_n, bench_p),
  BLOCK_SIZE_MODE = bench_block_mode,
  LOADING_STRENGTH = bench_strength,
  PRIMARY_LOADING_RANGE = bench_primary_range,
  CROSS_LOADING_RANGE = bench_cross_range,
  CROSS_LOADING_PROB = as.character(bench_cross_prob),
  RUN_OURS = "FALSE",
  RUN_JOINT_MFA = "FALSE",
  RUN_VIROLI = "TRUE",
  VIROLI_METHOD_NAME = "viroli_laplace_gibbs",
  VIROLI_LAMBDA_L1_PENALTY = as.character(bench_penalty),
  PARALLEL_OURS = "FALSE"
)

read_one_result <- function(case_dir, workers, parallel_enabled) {
  result_file <- file.path(case_dir, "comparison_results.csv")
  if (!file.exists(result_file)) {
    stop("Missing result file: ", result_file, call. = FALSE)
  }
  result <- read.csv(result_file, check.names = FALSE)
  result$benchmark_workers <- workers
  result$benchmark_parallel <- parallel_enabled

  history_file <- list.files(case_dir, pattern = "^viroli_history_.*[.]csv$", full.names = TRUE)
  if (!length(history_file)) {
    stop("Missing Viroli history in: ", case_dir, call. = FALSE)
  }
  history <- read.csv(history_file[[1]], check.names = FALSE)
  timing_cols <- grep("_seconds$", names(history), value = TRUE)
  timing <- data.frame(
    benchmark_workers = workers,
    benchmark_parallel = parallel_enabled,
    n_iter = nrow(history),
    mean_iteration_seconds = mean(history$iteration_seconds, na.rm = TRUE),
    median_iteration_seconds = median(history$iteration_seconds, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  for (nm in timing_cols) {
    timing[[paste0("mean_", nm)]] <- mean(history[[nm]], na.rm = TRUE)
    timing[[paste0("median_", nm)]] <- median(history[[nm]], na.rm = TRUE)
  }

  list(result = result, timing = timing)
}

all_results <- list()
all_timings <- list()

cat("Viroli inner parallel benchmark\n")
cat(sprintf(
  "Setting: n=%d, p=%d, H=%d, G=%d, sep=%s, block=%s, loadings=%s, strength=%s, cross_prob=%s\n",
  bench_n, bench_p, bench_H, bench_G, bench_sep, bench_block_mode,
  bench_loading_design, bench_strength, bench_cross_prob
))
cat(sprintf(
  "Iterations: %d, burn: %d, Laplace penalty: %s, workers tested: %s\n",
  bench_iter, bench_burn, bench_penalty, paste(bench_workers, collapse = ", ")
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
    PARALLEL_GIBBS = if (parallel_enabled) "TRUE" else "FALSE",
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
  all_timings[[label]] <- one$timing
  cat(sprintf(
    "%s complete: wall %.2fs, mean iter %.4fs\n",
    label, one$result$seconds[[1]], one$timing$mean_iteration_seconds[[1]]
  ))
}

results <- do.call(rbind, all_results)
timings <- do.call(rbind, all_timings)

results_file <- file.path(out_base, "viroli_inner_parallel_benchmark_results.csv")
timings_file <- file.path(out_base, "viroli_inner_parallel_iteration_timing.csv")
write.csv(results, results_file, row.names = FALSE)
write.csv(timings, timings_file, row.names = FALSE)

plot_file <- file.path(out_base, "viroli_inner_parallel_timing.png")
png(plot_file, width = 1600, height = 750, res = 150)
op <- par(mfrow = c(1, 2), mar = c(6, 5, 4, 1), oma = c(0, 0, 2, 0))
labels <- ifelse(timings$benchmark_parallel, paste0(timings$benchmark_workers, " workers"), "serial")
barplot(
  results$seconds,
  names.arg = labels,
  las = 2,
  col = "#377eb8",
  ylab = "seconds",
  main = "End-to-end"
)
barplot(
  timings$mean_iteration_seconds,
  names.arg = labels,
  las = 2,
  col = "#4daf4a",
  ylab = "seconds",
  main = "Mean time per Gibbs iteration"
)
mtext(
  sprintf("Viroli Laplace inner parallel benchmark: n=%d, p=%d, H=%d, G=%d", bench_n, bench_p, bench_H, bench_G),
  outer = TRUE,
  cex = 1.05,
  font = 2
)
par(op)
dev.off()

cat("\nSaved benchmark summary:\n")
cat(results_file, "\n", sep = "")
cat(timings_file, "\n", sep = "")
cat(plot_file, "\n", sep = "")

display_cols <- c(
  "benchmark_workers", "benchmark_parallel", "seconds", "seconds_per_iter",
  "factor_score_rmse", "lambda_rmse", "alpha_rmse", "marginal_mu_rmse",
  "marginal_var_rmse", "marginal_weight_rmse"
)
print(results[, intersect(display_cols, names(results)), drop = FALSE], row.names = FALSE)

timing_cols <- c(
  "benchmark_workers", "mean_iteration_seconds", "mean_z_sample_seconds",
  "mean_factor_sample_seconds", "mean_class_sample_seconds",
  "mean_mixture_parameter_seconds", "mean_regression_seconds",
  "mean_normalization_seconds", "mean_keep_draw_seconds"
)
print(timings[, intersect(timing_cols, names(timings)), drop = FALSE], row.names = FALSE)
