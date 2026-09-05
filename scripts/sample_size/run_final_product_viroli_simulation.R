#!/usr/bin/env Rscript

# Final high-dimensional sample-size launcher.
#
# This is the public entry point for the simulation study in the GitHub repo.
# It varies sample size, item count, rank, component count, loading strength,
# cross-loading density, and block balance.  Product MAP is run everywhere.
# The two Viroli-style Gibbs baselines are run only for p in {500, 1000} by
# default because those are the direct timing/accuracy comparison cells.
#
# Replications are run serially at the launcher level.  Within each replication,
# the methods use their pre-declared optimized internal worker settings:
# Product MAP uses 18 workers for all p, while Viroli/Gibbs uses one serial
# chain worker for p=500 and four internal workers for p=1000.

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

parse_ints <- function(x) as.integer(split_csv(x))
parse_nums <- function(x) as.numeric(split_csv(x))

safe_token <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "value")
}

format_num <- function(x) format(x, scientific = FALSE, trim = TRUE)

np_grid_from_values <- function(n_values, p_values) {
  pieces <- character(0)
  for (p in p_values) {
    for (n in n_values) {
      pieces <- c(pieces, sprintf("n%dp%d:%d:%d", n, p, n, p))
    }
  }
  paste(pieces, collapse = ",")
}

rbind_fill <- function(x) {
  x <- x[!vapply(x, is.null, logical(1L))]
  if (!length(x)) return(data.frame())
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(d) {
    missing <- setdiff(all_names, names(d))
    for (nm in missing) d[[nm]] <- NA
    d[, all_names, drop = FALSE]
  })
  do.call(rbind, x)
}

loading_ranges <- function(strength) {
  strength <- match.arg(strength, c("weak", "strong"))
  switch(
    strength,
    weak = list(primary = c(1.25, 1.75), cross = c(1.25, 1.75)),
    strong = list(primary = c(2.50, 3.00), cross = c(2.50, 3.00))
  )
}

combine_results <- function(out_dir, chunk_dir) {
  task_dirs <- unique(dirname(list.files(
    chunk_dir,
    pattern = "comparison_results(_checkpoint)?\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )))
  if (!length(task_dirs)) return(data.frame())
  files <- vapply(task_dirs, function(d) {
    final <- file.path(d, "comparison_results.csv")
    checkpoint <- file.path(d, "comparison_results_checkpoint.csv")
    if (file.exists(final)) final else checkpoint
  }, character(1L))
  files <- files[file.exists(files)]
  pieces <- lapply(files, function(path) {
    d <- read.csv(path, check.names = FALSE)
    d$source_chunk <- basename(dirname(path))
    d
  })
  results <- rbind_fill(pieces)
  if (nrow(results)) {
    write.csv(results, file.path(out_dir, "comparison_results.csv"), row.names = FALSE)
    metric_cols <- intersect(
      c(
        "factor_score_rmse", "mean_factor_abs_cor", "lambda_rmse", "alpha_rmse",
        "marginal_mu_rmse", "marginal_var_rmse", "marginal_weight_rmse",
        "stage1_signal_relative_frobenius_error", "stage1_sinTheta_op",
        "seconds", "ess_min", "ess_median", "ess_mean"
      ),
      names(results)
    )
    group_cols <- intersect(
      c(
        "method", "n", "p", "H_true", "G_true", "loading_strength",
        "cross_loading_prob", "block_size_mode", "loading_design",
        "viroli_loading_prior"
      ),
      names(results)
    )
    if (length(group_cols) && length(metric_cols)) {
      group_frame <- as.data.frame(
        lapply(results[, group_cols, drop = FALSE], as.character),
        stringsAsFactors = FALSE
      )
      group_key <- do.call(interaction, c(group_frame, list(drop = TRUE, sep = " | ")))
      summary <- lapply(split(results, group_key), function(d) {
        base <- d[1L, group_cols, drop = FALSE]
        metric_summary <- do.call(cbind, lapply(metric_cols, function(metric) {
          x <- d[[metric]]
          data.frame(
            setNames(list(sum(is.finite(x))), paste0(metric, "_n")),
            setNames(list(if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)), paste0(metric, "_mean")),
            setNames(list(if (sum(is.finite(x)) <= 1L) NA_real_ else sd(x, na.rm = TRUE)), paste0(metric, "_sd")),
            check.names = FALSE
          )
        }))
        cbind(base, metric_summary)
      })
      summary <- do.call(rbind, summary)
      rownames(summary) <- NULL
      write.csv(summary, file.path(out_dir, "comparison_summary.csv"), row.names = FALSE)
    }
  }
  results
}

run_one_task <- function(task) {
  task_dir <- task$out_dir
  final_file <- file.path(task_dir, "comparison_results.csv")
  checkpoint_file <- file.path(task_dir, "comparison_results_checkpoint.csv")
  if (file.exists(final_file) || file.exists(checkpoint_file)) {
    return(data.frame(task = task$name, status = "skipped_existing", seconds = 0))
  }
  dir.create(task_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(task_dir, "run.log")
  env <- task$env
  env_vec <- sprintf("%s=%s", names(env), shQuote(unname(env), type = "sh"))
  t0 <- proc.time()[["elapsed"]]
  status <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = shQuote(driver, type = "sh"),
    env = env_vec,
    stdout = log_file,
    stderr = log_file,
    wait = TRUE
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (!identical(status, 0L)) {
    return(data.frame(task = task$name, status = paste0("failed_", status), seconds = elapsed))
  }
  data.frame(task = task$name, status = "completed", seconds = elapsed)
}

run_task_pool <- function(tasks, workers, out_dir, chunk_dir, phase_label) {
  if (!length(tasks)) return(invisible(data.frame()))
  workers <- max(1L, min(as.integer(workers), length(tasks)))
  cat(sprintf("\nPhase %s: %d task chunks, %d worker(s)\n", phase_label, length(tasks), workers))
  runner <- function(task) {
    cat(sprintf("  [%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), task$name))
    run_one_task(task)
  }
  task_status <- if (workers > 1L) {
    parallel::mclapply(tasks, runner, mc.cores = workers, mc.preschedule = FALSE)
  } else {
    lapply(tasks, runner)
  }
  status <- do.call(rbind, task_status)
  status_file <- file.path(out_dir, paste0("task_status_", safe_token(phase_label), ".csv"))
  write.csv(status, status_file, row.names = FALSE)
  combine_results(out_dir, chunk_dir)
  cat(sprintf(
    "Phase %s complete: %d completed, %d skipped, %d failed\n",
    phase_label,
    sum(status$status == "completed"),
    sum(status$status == "skipped_existing"),
    sum(!status$status %in% c("completed", "skipped_existing"))
  ))
  invisible(status)
}

n_values <- parse_ints(get_env("N_VALUES", "100,200"))
p_values_product <- parse_ints(get_env("P_VALUES_PRODUCT", "500,1000,1500,2000"))
p_values_gibbs <- parse_ints(get_env("P_VALUES_GIBBS", "500,1000"))
h_values <- parse_ints(get_env("H_VALUES", "5,10,15,20"))
g_values <- parse_ints(get_env("G_VALUES", "2,3"))
loading_strengths <- split_csv(get_env("LOADING_STRENGTHS", "weak,strong"))
cross_loading_probs <- parse_nums(get_env("CROSS_LOADING_PROBS", "0.075,0.2"))
block_modes <- split_csv(get_env("BLOCK_SIZE_MODES", "balanced,ifeval_like"))
rep_values <- parse_ints(get_env("REP_VALUES", paste(seq_len(25L), collapse = ",")))
separations <- get_env("SEPARATIONS", "1")
run_label <- get_env("RUN_LABEL", "signal_support_grid_optimized_workers")
out_dir <- get_env("OUT_DIR", file.path(repo_root, "results", "full", run_label))
chunk_dir <- file.path(out_dir, "chunks")
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

product_task_workers <- as.integer(get_env("TASK_WORKERS_PRODUCT", "1"))
product_internal_workers <- get_env("PRODUCT_INTERNAL_WORKERS", "18")
product_internal_workers_small_p <- get_env(
  "PRODUCT_INTERNAL_WORKERS_SMALL_P",
  product_internal_workers
)
product_internal_workers_large_p <- get_env(
  "PRODUCT_INTERNAL_WORKERS_LARGE_P",
  product_internal_workers
)
gibbs_parallel_p_min <- as.integer(get_env("GIBBS_PARALLEL_P_MIN", "1000"))
gibbs_task_workers <- as.integer(get_env("TASK_WORKERS_GIBBS", "1"))
gibbs_internal_workers_serial <- get_env("GIBBS_INTERNAL_WORKERS_SERIAL", "1")
gibbs_internal_workers_parallel <- get_env("GIBBS_INTERNAL_WORKERS_PARALLEL", "4")

common_env <- c(
  SEED = get_env("SEED", "20260731"),
  RESUME_EXISTING = "TRUE",
  LOADING_DESIGNS = get_env("LOADING_DESIGNS", "balanced_moderate_dense_signed_cross"),
  LOADING_SIGN_MODE = get_env("LOADING_SIGN_MODE", "block"),
  CROSS_SIGN_MODE = get_env("CROSS_SIGN_MODE", "random"),
  ALIGNMENT_MODE = get_env("ALIGNMENT_MODE", "loadings"),
  MIXTURE_PARAM_MODE = get_env("MIXTURE_PARAM_MODE", "viroli_smoke"),
  MIXTURE_VARIANCE_MODE = get_env("MIXTURE_VARIANCE_MODE", "unequal"),
  INTERCEPT_MODE = get_env("INTERCEPT_MODE", "ifeval_like"),
  INTERCEPT_SD = get_env("INTERCEPT_SD", "0.45"),
  INTERCEPT_BLOCK_SPAN = get_env("INTERCEPT_BLOCK_SPAN", "1.6"),
  INTERCEPT_CLIP = get_env("INTERCEPT_CLIP", "1.75"),
  SEPARATIONS = separations,
  OURS_PRETRAINING_METHOD = get_env("OURS_PRETRAINING_METHOD", "em_svd"),
  EM_SVD_INIT = get_env("EM_SVD_INIT", "both"),
  EM_SVD_INIT_Z = get_env("EM_SVD_INIT_Z", "expectation"),
  EM_SVD_ITER = get_env("EM_SVD_ITER", "50"),
  EM_SVD_TOL_LOGLIK = get_env("EM_SVD_TOL_LOGLIK", "1e-5"),
  EM_SVD_TOL_L = get_env("EM_SVD_TOL_L", "1e-4"),
  EM_SVD_RANDOM_STARTS = get_env("EM_SVD_RANDOM_STARTS", "0"),
  PRETRAIN_LOADING_PENALTY = get_env("PRETRAIN_LOADING_PENALTY", "10"),
  ROTATION_OPTIMIZER = get_env("ROTATION_OPTIMIZER", "riemannian"),
  ROTATION_ITER = get_env("ROTATION_ITER", "20"),
  ROTATION_MIN_OUTER = get_env("ROTATION_MIN_OUTER", "2"),
  ROTATION_OBJECTIVE_TOLERANCE = get_env("ROTATION_OBJECTIVE_TOLERANCE", "1e-3"),
  ROTATION_REQUIRE_MIXTURE_CONVERGENCE = get_env("ROTATION_REQUIRE_MIXTURE_CONVERGENCE", "TRUE"),
  ROTATION_RANDOM_STARTS = get_env("ROTATION_RANDOM_STARTS", "1"),
  ROTATION_N_MIX_STARTS = get_env("ROTATION_N_MIX_STARTS", "3"),
  ROTATION_GRID_SIZE = get_env("ROTATION_GRID_SIZE", "21"),
  RIEMANNIAN_ROTATION_STEPS = get_env("RIEMANNIAN_ROTATION_STEPS", "10"),
  RIEMANNIAN_ETA0 = get_env("RIEMANNIAN_ETA0", "1"),
  RIEMANNIAN_BETA = get_env("RIEMANNIAN_BETA", "0.5"),
  RIEMANNIAN_MIN_ETA = get_env("RIEMANNIAN_MIN_ETA", "1e-8"),
  RIEMANNIAN_GRAD_TOL = get_env("RIEMANNIAN_GRAD_TOL", "1e-6"),
  RIEMANNIAN_UPDATE = get_env("RIEMANNIAN_UPDATE", "cayley"),
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
  WRITE_PARAMETER_TABLES = get_env("WRITE_PARAMETER_TABLES", "FALSE"),
  WRITE_ITERATION_HISTORIES = get_env("WRITE_ITERATION_HISTORIES", "FALSE"),
  MAX_JOINT_PARAMETER_K = get_env("MAX_JOINT_PARAMETER_K", "5000"),
  MAX_JOINT_PROFILE_ARI_K = get_env("MAX_JOINT_PROFILE_ARI_K", "5000"),
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

make_tasks <- function(
    phase,
    p_values,
    method,
    product_inner_workers = product_internal_workers_small_p,
    gibbs_parallel = FALSE,
    gibbs_inner_workers = gibbs_internal_workers_serial) {
  tasks <- list()
  for (block_mode in block_modes) {
    for (strength in loading_strengths) {
      ranges <- loading_ranges(strength)
      for (cross_prob in cross_loading_probs) {
        for (rep in rep_values) {
          for (H in h_values) {
            for (G in g_values) {
              name <- paste(
                phase,
                method,
                paste0("block", safe_token(block_mode)),
                paste0("strength", safe_token(strength)),
                paste0("cp", safe_token(format_num(cross_prob))),
                paste0("rep", rep),
                paste0("H", H),
                paste0("G", G),
                sep = "_"
              )
              task_dir <- file.path(chunk_dir, name)
              env <- c(
                common_env,
                OUT_DIR = task_dir,
                REP_VALUES = as.character(rep),
                H_VALUES = as.character(H),
                G_VALUES = as.character(G),
                NP_GRID = np_grid_from_values(n_values, p_values),
                BLOCK_SIZE_MODE = block_mode,
                LOADING_STRENGTH = strength,
                PRIMARY_LOADING_RANGE = paste(format_num(ranges$primary), collapse = ","),
                CROSS_LOADING_RANGE = paste(format_num(ranges$cross), collapse = ","),
                CROSS_LOADING_PROB = format_num(cross_prob)
              )
              if (method == "product_map") {
                env <- c(
                  env,
                  RUN_OURS = "TRUE",
                  RUN_JOINT_MFA = "FALSE",
                  RUN_VIROLI = "FALSE",
                  PARALLEL_OURS = "TRUE",
                  PARALLEL_GIBBS = "FALSE",
                  PARALLEL_WORKERS = product_inner_workers
                )
              } else if (method == "viroli_laplace") {
                env <- c(
                  env,
                  RUN_OURS = "FALSE",
                  RUN_JOINT_MFA = "FALSE",
                  RUN_VIROLI = "TRUE",
                  VIROLI_METHOD_NAME = "viroli_laplace_gibbs",
                  VIROLI_LAMBDA_L1_PENALTY = get_env("VIROLI_LAPLACE_L1_PENALTY", "10"),
                  PARALLEL_OURS = "FALSE",
                  PARALLEL_GIBBS = if (isTRUE(gibbs_parallel)) "TRUE" else "FALSE",
                  PARALLEL_WORKERS = gibbs_inner_workers
                )
              } else if (method == "viroli_gaussian") {
                env <- c(
                  env,
                  RUN_OURS = "FALSE",
                  RUN_JOINT_MFA = "FALSE",
                  RUN_VIROLI = "TRUE",
                  VIROLI_METHOD_NAME = "viroli_gaussian_gibbs",
                  VIROLI_LAMBDA_L1_PENALTY = "0",
                  PARALLEL_OURS = "FALSE",
                  PARALLEL_GIBBS = if (isTRUE(gibbs_parallel)) "TRUE" else "FALSE",
                  PARALLEL_WORKERS = gibbs_inner_workers
                )
              } else {
                stop("Unknown method: ", method, call. = FALSE)
              }
              tasks[[length(tasks) + 1L]] <- list(name = name, out_dir = task_dir, env = env)
            }
          }
        }
      }
    }
  }
  tasks
}

cat("Final signal-support simulation launcher\n")
cat("Output directory:", out_dir, "\n")
cat("n grid:", paste(n_values, collapse = ", "), "\n")
cat("Product p grid:", paste(p_values_product, collapse = ", "), "\n")
cat("Gibbs p grid:", paste(p_values_gibbs, collapse = ", "), "\n")
cat("H grid:", paste(h_values, collapse = ", "), "\n")
cat("G grid:", paste(g_values, collapse = ", "), "\n")
cat("Loading strengths:", paste(loading_strengths, collapse = ", "), "\n")
cat("Cross-loading probabilities:", paste(cross_loading_probs, collapse = ", "), "\n")
cat("Block modes:", paste(block_modes, collapse = ", "), "\n")
cat("Reps per cell:", length(rep_values), "\n")
cat("Product task workers:", product_task_workers, "\n")
cat("Product internal workers:", product_internal_workers, "\n")
cat("Gibbs/Viroli task workers:", gibbs_task_workers, "\n")
cat("Gibbs/Viroli serial-inner p values:", paste(p_values_gibbs[p_values_gibbs < gibbs_parallel_p_min], collapse = ", "), "\n")
cat("Gibbs/Viroli parallel-inner p values:", paste(p_values_gibbs[p_values_gibbs >= gibbs_parallel_p_min], collapse = ", "), "\n")
cat("Gibbs/Viroli serial-inner workers:", gibbs_internal_workers_serial, "\n")
cat("Gibbs/Viroli parallel-inner workers:", gibbs_internal_workers_parallel, "\n")

small_product_p <- p_values_product[p_values_product <= 1000]
large_product_p <- p_values_product[p_values_product > 1000]
gibbs_serial_p <- p_values_gibbs[p_values_gibbs < gibbs_parallel_p_min]
gibbs_parallel_p <- p_values_gibbs[p_values_gibbs >= gibbs_parallel_p_min]

if (length(small_product_p)) {
  run_task_pool(
    make_tasks(
      "smallp",
      small_product_p,
      "product_map",
      product_inner_workers = product_internal_workers_small_p
    ),
    product_task_workers,
    out_dir,
    chunk_dir,
    "product_map_small_p"
  )
  if (length(gibbs_serial_p)) {
    run_task_pool(
      make_tasks(
        "smallp_serial",
        gibbs_serial_p,
        "viroli_laplace",
        gibbs_parallel = FALSE,
        gibbs_inner_workers = gibbs_internal_workers_serial
      ),
      gibbs_task_workers,
      out_dir,
      chunk_dir,
      "viroli_laplace_serial_inner_p"
    )
    run_task_pool(
      make_tasks(
        "smallp_serial",
        gibbs_serial_p,
        "viroli_gaussian",
        gibbs_parallel = FALSE,
        gibbs_inner_workers = gibbs_internal_workers_serial
      ),
      gibbs_task_workers,
      out_dir,
      chunk_dir,
      "viroli_gaussian_serial_inner_p"
    )
  }
  if (length(gibbs_parallel_p)) {
    run_task_pool(
      make_tasks(
        "smallp_parallel",
        gibbs_parallel_p,
        "viroli_laplace",
        gibbs_parallel = TRUE,
        gibbs_inner_workers = gibbs_internal_workers_parallel
      ),
      gibbs_task_workers,
      out_dir,
      chunk_dir,
      "viroli_laplace_parallel_inner_p"
    )
    run_task_pool(
      make_tasks(
        "smallp_parallel",
        gibbs_parallel_p,
        "viroli_gaussian",
        gibbs_parallel = TRUE,
        gibbs_inner_workers = gibbs_internal_workers_parallel
      ),
      gibbs_task_workers,
      out_dir,
      chunk_dir,
      "viroli_gaussian_parallel_inner_p"
    )
  }
}

if (length(large_product_p)) {
  run_task_pool(
    make_tasks(
      "largep",
      large_product_p,
      "product_map",
      product_inner_workers = product_internal_workers_large_p
    ),
    product_task_workers,
    out_dir,
    chunk_dir,
    "product_map_large_p"
  )
}

combined <- combine_results(out_dir, chunk_dir)
cat("\nLauncher complete. Current combined result rows:", nrow(combined), "\n")
cat("Raw results:", file.path(out_dir, "comparison_results.csv"), "\n")
cat("Summary:", file.path(out_dir, "comparison_summary.csv"), "\n")
