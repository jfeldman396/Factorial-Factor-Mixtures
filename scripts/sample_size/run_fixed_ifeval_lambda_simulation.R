#!/usr/bin/env Rscript

# Fixed-DGP IFEval-like simulation launcher.
#
# Scientific target:
#   Hold the population model fixed within each design cell, then simulate new
#   datasets across replications.  This isolates estimator behavior from random
#   changes in the loading matrix.  Smaller p settings are nested block-wise
#   subsets of one master IFEval-like loading matrix.
#
# Default design:
#   - IFEval-like unbalanced blocks with at least 30 primary items per factor.
#   - Nonzero loadings have magnitudes Uniform(2, 3).
#   - Cross-loadings are signed and moderately dense.
#   - Product MAP is run on all p values.
#   - Viroli Laplace and diffuse Gaussian Gibbs are run on p <= 1000.
#   - Mixture complexity is varied through all-2 versus all-3 component
#     settings; mixed/alternating G settings are intentionally excluded.

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

safe_token <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "value")
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

rbind_fill <- function(x) {
  x <- x[!vapply(x, is.null, logical(1L))]
  if (!length(x)) return(data.frame())
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(d) {
    for (nm in setdiff(all_names, names(d))) d[[nm]] <- NA
    d[, all_names, drop = FALSE]
  })
  do.call(rbind, x)
}

g_config_for_type <- function(type, H) {
  type <- match.arg(type, c("all2", "all3"))
  if (type == "all2") return(rep(2L, H))
  if (type == "all3") return(rep(3L, H))
}

format_g_config <- function(G) paste(as.integer(G), collapse = ",")
label_g_config <- function(G) paste(as.integer(G), collapse = "-")

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
  write.csv(results, file.path(out_dir, "comparison_results.csv"), row.names = FALSE)
  results
}

run_one_task <- function(task) {
  final_file <- file.path(task$out_dir, "comparison_results.csv")
  checkpoint_file <- file.path(task$out_dir, "comparison_results_checkpoint.csv")
  if (file.exists(final_file) || file.exists(checkpoint_file)) {
    return(data.frame(task = task$name, status = "skipped_existing", seconds = 0))
  }
  dir.create(task$out_dir, recursive = TRUE, showWarnings = FALSE)
  env_vec <- sprintf("%s=%s", names(task$env), shQuote(unname(task$env), type = "sh"))
  t0 <- proc.time()[["elapsed"]]
  status <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = shQuote(driver, type = "sh"),
    env = env_vec,
    stdout = file.path(task$out_dir, "run.log"),
    stderr = file.path(task$out_dir, "run.log"),
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
  cat(sprintf("\nPhase %s: %d task chunks, %d launcher worker(s)\n", phase_label, length(tasks), workers))
  runner <- function(task) {
    cat(sprintf("  [%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), task$name))
    run_one_task(task)
  }
  status <- if (workers > 1L) {
    parallel::mclapply(tasks, runner, mc.cores = workers, mc.preschedule = FALSE)
  } else {
    lapply(tasks, runner)
  }
  status <- do.call(rbind, status)
  write.csv(status, file.path(out_dir, paste0("task_status_", safe_token(phase_label), ".csv")), row.names = FALSE)
  combined <- combine_results(out_dir, chunk_dir)
  cat(sprintf(
    "Phase %s checkpoint: %d completed, %d skipped, %d failed, %d result rows\n",
    phase_label,
    sum(status$status == "completed"),
    sum(status$status == "skipped_existing"),
    sum(!status$status %in% c("completed", "skipped_existing")),
    nrow(combined)
  ))
  invisible(status)
}

n_values <- parse_ints(get_env("N_VALUES", "100,200,400"))
p_values_product <- parse_ints(get_env("P_VALUES_PRODUCT", "500,1000,1500,2000"))
p_values_gibbs <- parse_ints(get_env("P_VALUES_GIBBS", "500,1000"))
h_values <- parse_ints(get_env("H_VALUES", "5,10"))
g_config_types <- split_csv(get_env("G_CONFIG_TYPES", "all2,all3"))
rep_values <- parse_ints(get_env("REP_VALUES", paste(seq_len(25L), collapse = ",")))
run_label <- get_env("RUN_LABEL", "fixed_ifeval_lambda_min30_u2_3_cp0_05_h5_h10")
out_dir <- get_env("OUT_DIR", file.path(repo_root, "results", "full", run_label))
chunk_dir <- file.path(out_dir, "chunks")
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

product_task_workers <- as.integer(get_env("TASK_WORKERS_PRODUCT", "1"))
gibbs_task_workers <- as.integer(get_env("TASK_WORKERS_GIBBS", "1"))
product_internal_workers <- get_env("PRODUCT_INTERNAL_WORKERS", "18")
gibbs_parallel_p_min <- as.integer(get_env("GIBBS_PARALLEL_P_MIN", "0"))
gibbs_internal_workers_serial <- get_env("GIBBS_INTERNAL_WORKERS_SERIAL", "4")
gibbs_internal_workers_parallel <- get_env("GIBBS_INTERNAL_WORKERS_PARALLEL", "4")
dgp_p_max <- as.character(max(c(p_values_product, p_values_gibbs)))

common_env <- c(
  SEED = get_env("SEED", "20260731"),
  RESUME_EXISTING = "TRUE",
  FIX_DGP_PARAMETERS = "TRUE",
  DGP_P_MAX = dgp_p_max,
  LOADING_DESIGNS = "balanced_moderate_dense_signed_cross",
  BLOCK_SIZE_MODE = "ifeval_min30",
  LOADING_SIGN_MODE = "block",
  LOADING_STRENGTH = "strong",
  PRIMARY_LOADING_RANGE = get_env("PRIMARY_LOADING_RANGE", "2,3"),
  CROSS_LOADING_RANGE = get_env("CROSS_LOADING_RANGE", "2,3"),
  CROSS_LOADING_PROB = get_env("CROSS_LOADING_PROB", "0.05"),
  CROSS_SIGN_MODE = get_env("CROSS_SIGN_MODE", "random"),
  ALIGNMENT_MODE = get_env("ALIGNMENT_MODE", "loadings"),
  MIXTURE_PARAM_MODE = get_env("MIXTURE_PARAM_MODE", "viroli_smoke"),
  MIXTURE_VARIANCE_MODE = get_env("MIXTURE_VARIANCE_MODE", "unequal"),
  INTERCEPT_MODE = get_env("INTERCEPT_MODE", "ifeval_like"),
  INTERCEPT_SD = get_env("INTERCEPT_SD", "0.45"),
  INTERCEPT_BLOCK_SPAN = get_env("INTERCEPT_BLOCK_SPAN", "1.6"),
  INTERCEPT_CLIP = get_env("INTERCEPT_CLIP", "1.75"),
  SEPARATIONS = get_env("SEPARATIONS", "1"),
  OURS_PRETRAINING_METHOD = get_env("OURS_PRETRAINING_METHOD", "em_svd"),
  EM_SVD_INIT = get_env("EM_SVD_INIT", "both"),
  EM_SVD_INIT_Z = get_env("EM_SVD_INIT_Z", "expectation"),
  EM_SVD_ITER = get_env("EM_SVD_ITER", "50"),
  EM_SVD_TOL_LOGLIK = get_env("EM_SVD_TOL_LOGLIK", "1e-5"),
  EM_SVD_TOL_L = get_env("EM_SVD_TOL_L", "1e-4"),
  PRETRAIN_LOADING_PENALTY = get_env("PRETRAIN_LOADING_PENALTY", "10"),
  ROTATION_OPTIMIZER = get_env("ROTATION_OPTIMIZER", "riemannian"),
  ROTATION_ITER = get_env("ROTATION_ITER", "20"),
  ROTATION_REQUIRE_MIXTURE_CONVERGENCE = get_env("ROTATION_REQUIRE_MIXTURE_CONVERGENCE", "TRUE"),
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
  MIN_MIXTURE_VAR = get_env("MIN_MIXTURE_VAR", "0.05"),
  WRITE_PARAMETER_TABLES = get_env("WRITE_PARAMETER_TABLES", "FALSE"),
  WRITE_ITERATION_HISTORIES = get_env("WRITE_ITERATION_HISTORIES", "FALSE"),
  MAX_JOINT_PARAMETER_K = get_env("MAX_JOINT_PARAMETER_K", "5000"),
  MAX_JOINT_PROFILE_ARI_K = get_env("MAX_JOINT_PROFILE_ARI_K", "5000"),
  VIROLI_ITER = get_env("VIROLI_ITER", "2000"),
  VIROLI_BURN = get_env("VIROLI_BURN", "1000"),
  VIROLI_THIN = get_env("VIROLI_THIN", "1"),
  VIROLI_COMPUTE_PARAMETER_ESS = get_env("VIROLI_COMPUTE_PARAMETER_ESS", "TRUE"),
  VIROLI_NORMALIZE_EACH_DRAW = get_env("VIROLI_NORMALIZE_EACH_DRAW", "TRUE"),
  VIROLI_VERBOSE = get_env("VIROLI_VERBOSE", "FALSE")
)

make_tasks <- function(phase, p_values, method, gibbs_parallel = FALSE, gibbs_inner_workers = "1") {
  tasks <- list()
  for (rep in rep_values) {
    for (H in h_values) {
      for (g_type in g_config_types) {
        G <- g_config_for_type(g_type, H)
        g_label <- label_g_config(G)
        name <- paste(phase, method, paste0("rep", rep), paste0("H", H), paste0("G", safe_token(g_label)), sep = "_")
        task_dir <- file.path(chunk_dir, name)
        env <- c(
          common_env,
          OUT_DIR = task_dir,
          REP_VALUES = as.character(rep),
          H_VALUES = as.character(H),
          G_CONFIGS = format_g_config(G),
          NP_GRID = np_grid_from_values(n_values, p_values)
        )
        if (method == "product_map") {
          env <- c(
            env,
            RUN_OURS = "TRUE",
            RUN_JOINT_MFA = "FALSE",
            RUN_VIROLI = "FALSE",
            PARALLEL_OURS = "TRUE",
            PARALLEL_GIBBS = "FALSE",
            PARALLEL_WORKERS = product_internal_workers
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
  tasks
}

cat("Fixed-DGP IFEval-like Lambda simulation\n")
cat("Output directory:", out_dir, "\n")
cat("n grid:", paste(n_values, collapse = ", "), "\n")
cat("Product p grid:", paste(p_values_product, collapse = ", "), "\n")
cat("Gibbs p grid:", paste(p_values_gibbs, collapse = ", "), "\n")
cat("H grid:", paste(h_values, collapse = ", "), "\n")
cat("G configurations:", paste(g_config_types, collapse = ", "), "\n")
cat("Block mode: ifeval_min30; nonzero loading magnitudes: Uniform(2, 3)\n")
cat("Cross-loading probability:", common_env[["CROSS_LOADING_PROB"]], "\n")
cat("Replications:", length(rep_values), "\n")

run_task_pool(
  make_tasks("fixed_ifeval", p_values_product, "product_map"),
  product_task_workers,
  out_dir,
  chunk_dir,
  "product_map_all_p"
)

gibbs_serial_p <- p_values_gibbs[p_values_gibbs < gibbs_parallel_p_min]
gibbs_parallel_p <- p_values_gibbs[p_values_gibbs >= gibbs_parallel_p_min]

if (length(gibbs_serial_p)) {
  run_task_pool(
    make_tasks("fixed_ifeval_serial", gibbs_serial_p, "viroli_laplace", gibbs_parallel = FALSE, gibbs_inner_workers = gibbs_internal_workers_serial),
    gibbs_task_workers,
    out_dir,
    chunk_dir,
    "viroli_laplace_serial_p"
  )
  run_task_pool(
    make_tasks("fixed_ifeval_serial", gibbs_serial_p, "viroli_gaussian", gibbs_parallel = FALSE, gibbs_inner_workers = gibbs_internal_workers_serial),
    gibbs_task_workers,
    out_dir,
    chunk_dir,
    "viroli_gaussian_serial_p"
  )
}

if (length(gibbs_parallel_p)) {
  run_task_pool(
    make_tasks("fixed_ifeval_parallel", gibbs_parallel_p, "viroli_laplace", gibbs_parallel = TRUE, gibbs_inner_workers = gibbs_internal_workers_parallel),
    gibbs_task_workers,
    out_dir,
    chunk_dir,
    "viroli_laplace_parallel_p"
  )
  run_task_pool(
    make_tasks("fixed_ifeval_parallel", gibbs_parallel_p, "viroli_gaussian", gibbs_parallel = TRUE, gibbs_inner_workers = gibbs_internal_workers_parallel),
    gibbs_task_workers,
    out_dir,
    chunk_dir,
    "viroli_gaussian_parallel_p"
  )
}

combined <- combine_results(out_dir, chunk_dir)
cat("\nLauncher complete. Current combined result rows:", nrow(combined), "\n")
cat("Raw results:", file.path(out_dir, "comparison_results.csv"), "\n")
