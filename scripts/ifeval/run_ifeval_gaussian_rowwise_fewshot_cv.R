#!/usr/bin/env Rscript

# Rank-H Gaussian probit factor-model comparison under the same row-wise
# few-shot protocol used by the P-IFA model-selection study. Each fold excludes
# complete model rows during training. For every held-out row, the shared
# reveal plan is used to estimate a Gaussian-prior MAP factor score, and the
# unrevealed responses are scored by Bernoulli log score, accuracy, and Brier.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

source(file.path(repo_root, "R", "pifa_rowwise_cv.R"))
source(file.path(repo_root, "R", "viroli_probit_independent_gibbs.R"))

env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, if (default) "TRUE" else "FALSE")
  tolower(value) %in% c("true", "1", "yes")
}

read_binary_matrix_complete <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- raw[[1L]]
  usable <- vapply(seq_len(ncol(X)), function(j) {
    observed <- X[!is.na(X[, j]), j]
    length(observed) >= 2L && all(observed %in% c(0, 1)) &&
      length(unique(observed)) == 2L && !anyNA(X[, j])
  }, logical(1L))
  X[, usable, drop = FALSE]
}

gaussian_mixture_prior <- function(H) {
  lapply(seq_len(H), function(h) {
    list(pi = 1, mu = 0, var = 1, G = 1L, converged = TRUE)
  })
}

score_gaussian_heldout_rows <- function(X,
                                        fold,
                                        reveal_plan,
                                        fit,
                                        heldout_maxit = 100L,
                                        factor_score_bound = 5,
                                        probability_epsilon = 1e-10) {
  units <- Filter(function(z) identical(z$fold, as.integer(fold)), reveal_plan)
  H <- ncol(fit$Lambda_hat)
  prior <- gaussian_mixture_prior(H)
  score_rows <- vector("list", length(units))
  diagnostic_rows <- vector("list", length(units))

  for (u in seq_along(units)) {
    unit <- units[[u]]
    reveal_mask <- seq_len(ncol(X)) %in% unit$reveal_index
    factor_fit <- update_one_factor_score_mixture_missing(
      x_i = X[unit$row_index, ],
      obs_i = reveal_mask,
      f_init = rep(0, H),
      Lambda = fit$Lambda_hat,
      alpha = fit$alpha_hat,
      mixture_fits = prior,
      mixture_prior_weight = 1,
      maxit = heldout_maxit,
      factor_score_bound = factor_score_bound,
      return_diagnostics = TRUE
    )
    converged <- identical(factor_fit$convergence, 0L) &&
      is.finite(factor_fit$value) && all(is.finite(factor_fit$par))
    eta_columns <- as.list(stats::setNames(factor_fit$par, paste0("eta_", seq_len(H))))
    diagnostic_rows[[u]] <- cbind(
      data.frame(
        fold = fold,
        row_index = unit$row_index,
        row_id = unit$row_id,
        reveal_repeat = unit$reveal_repeat,
        converged = converged,
        convergence_code = factor_fit$convergence,
        objective = -factor_fit$value,
        optimizer_iterations = unname(factor_fit$counts[["function"]] %||% NA_integer_),
        n_revealed = unit$n_revealed,
        n_predicted = unit$n_predicted,
        message = factor_fit$message,
        stringsAsFactors = FALSE
      ),
      as.data.frame(eta_columns, optional = TRUE)
    )
    if (!converged) next

    prediction_index <- unit$prediction_index
    linear_predictor <- fit$alpha_hat[prediction_index] +
      as.numeric(fit$Lambda_hat[prediction_index, , drop = FALSE] %*% factor_fit$par)
    probability <- pmin(
      pmax(pnorm(linear_predictor), probability_epsilon),
      1 - probability_epsilon
    )
    outcome <- as.numeric(X[unit$row_index, prediction_index])
    loglik <- outcome * log(probability) + (1 - outcome) * log1p(-probability)
    score_rows[[u]] <- data.frame(
      fold = fold,
      row_index = unit$row_index,
      row_id = unit$row_id,
      reveal_repeat = unit$reveal_repeat,
      total_loglik = sum(loglik),
      n_test_entries = length(outcome),
      mean_loglik = mean(loglik),
      total_squared_error = sum((probability - outcome)^2),
      total_correct = sum(as.numeric(probability >= 0.5) == outcome),
      stringsAsFactors = FALSE
    )
  }

  score_rows <- Filter(Negate(is.null), score_rows)
  list(
    unit_scores = if (length(score_rows)) do.call(rbind, score_rows) else data.frame(),
    diagnostics = do.call(rbind, diagnostic_rows)
  )
}

matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  file.path(repo_root, "data", "ifeval_threshold_1", "openeval_ifeval_only_binary_matrix.csv")
)
pifa_cv_dir <- Sys.getenv(
  "PIFA_CV_DIR",
  file.path(repo_root, "results", "full", "ifeval_rowwise_fewshot_cv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(pifa_cv_dir, "gaussian_H4_lambda4_rowwise_fewshot_cv")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

H <- as.integer(Sys.getenv("H_FIXED", "4"))
lambda_l1_penalty <- as.numeric(Sys.getenv("LAMBDA_L1_PENALTY", "4"))
n_iter <- as.integer(Sys.getenv("GIBBS_ITER", "2000"))
burn <- as.integer(Sys.getenv("GIBBS_BURN", "1000"))
thin <- as.integer(Sys.getenv("GIBBS_THIN", "1"))
tau_lambda <- as.numeric(Sys.getenv("TAU_LAMBDA", "1.5"))
tau_intercept <- as.numeric(Sys.getenv("TAU_INTERCEPT", "5"))
heldout_maxit <- as.integer(Sys.getenv("HELDOUT_FACTOR_MAXIT", "100"))
factor_score_bound <- as.numeric(Sys.getenv("FACTOR_SCORE_BOUND", "5"))
seed <- as.integer(Sys.getenv("SEED", "1"))
fold_workers <- as.integer(Sys.getenv("FOLD_WORKERS", "5"))
save_fold_fits <- env_flag("SAVE_FOLD_FITS", TRUE)

if (burn >= n_iter) stop("GIBBS_BURN must be smaller than GIBBS_ITER.")
X <- read_binary_matrix_complete(matrix_path)
fold_assignments <- readRDS(file.path(pifa_cv_dir, "rowwise_cv_fold_assignments.rds"))
reveal_plan <- readRDS(file.path(pifa_cv_dir, "rowwise_cv_reveal_masks.rds"))
K <- max(fold_assignments)
if (length(fold_assignments) != nrow(X)) {
  stop("Saved row folds do not match the supplied matrix.")
}
plan_index <- pifa_reveal_plan_index(reveal_plan)
if (!identical(sort(unique(plan_index$row_index)), seq_len(nrow(X)))) {
  stop("Saved reveal plan does not cover every supplied row.")
}

settings <- data.frame(
  H = H,
  prior = "standard_normal",
  lambda_l1_penalty = lambda_l1_penalty,
  n_iter = n_iter,
  burn = burn,
  thin = thin,
  tau_lambda = tau_lambda,
  tau_intercept = tau_intercept,
  K_folds = K,
  reveal_fraction = unique(plan_index$n_revealed / plan_index$n_observed)[1L],
  reveal_repeats = max(plan_index$reveal_repeat),
  factor_score_bound = factor_score_bound,
  prediction_rule = "posterior_mean_parameters_plus_gaussian_MAP_factor",
  n_models = nrow(X),
  n_items = ncol(X)
)
write.csv(settings, file.path(out_dir, "gaussian_rowwise_settings.csv"), row.names = FALSE)

fit_fold <- function(fold) {
  train_rows <- which(fold_assignments != fold)
  start_time <- proc.time()[["elapsed"]]
  fit <- fit_viroli_probit_independent_gibbs(
    X = X[train_rows, , drop = FALSE],
    H = H,
    G = rep(1L, H),
    n_iter = n_iter,
    burn = burn,
    thin = thin,
    tau_lambda = tau_lambda,
    tau_intercept = tau_intercept,
    normalize_each_draw = TRUE,
    lambda_l1_penalty = lambda_l1_penalty,
    parallel = FALSE,
    workers = 1L,
    compute_parameter_ess = FALSE,
    retain_loading_second_moment = FALSE,
    align_retained_draws = TRUE,
    seed = seed + 100000L * fold,
    verbose = FALSE
  )
  fit_seconds <- proc.time()[["elapsed"]] - start_time
  evaluation_start <- proc.time()[["elapsed"]]
  scored <- score_gaussian_heldout_rows(
    X = X,
    fold = fold,
    reveal_plan = reveal_plan,
    fit = fit,
    heldout_maxit = heldout_maxit,
    factor_score_bound = factor_score_bound
  )
  evaluation_seconds <- proc.time()[["elapsed"]] - evaluation_start
  unit_scores <- scored$unit_scores
  fold_score <- data.frame(
    method = "Gaussian probit Gibbs",
    H = H,
    lambda_l1_penalty = lambda_l1_penalty,
    fold = fold,
    total_loglik = sum(unit_scores$total_loglik),
    n_test_entries = sum(unit_scores$n_test_entries),
    mean_loglik = sum(unit_scores$total_loglik) / sum(unit_scores$n_test_entries),
    accuracy = sum(unit_scores$total_correct) / sum(unit_scores$n_test_entries),
    brier = sum(unit_scores$total_squared_error) / sum(unit_scores$n_test_entries),
    heldout_factor_failures = sum(!scored$diagnostics$converged),
    fit_seconds = fit_seconds,
    evaluation_seconds = evaluation_seconds,
    n_kept_draws = fit$n_keep,
    alignment_mean_cost = fit$draw_alignment$mean_cost
  )
  if (save_fold_fits) {
    saveRDS(fit, file.path(out_dir, sprintf("gaussian_rowwise_fold%02d_fit.rds", fold)))
  }
  list(fold_score = fold_score, unit_scores = unit_scores, diagnostics = scored$diagnostics)
}

folds <- seq_len(K)
results <- if (.Platform$OS.type == "unix" && fold_workers > 1L) {
  parallel::mclapply(folds, fit_fold, mc.cores = min(fold_workers, K))
} else {
  lapply(folds, fit_fold)
}
fold_scores <- do.call(rbind, lapply(results, `[[`, "fold_score"))
unit_scores <- do.call(rbind, lapply(results, `[[`, "unit_scores"))
diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
write.csv(fold_scores, file.path(out_dir, "gaussian_rowwise_fold_scores.csv"), row.names = FALSE)
write.csv(unit_scores, file.path(out_dir, "gaussian_rowwise_unit_scores.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out_dir, "gaussian_rowwise_factor_diagnostics.csv"), row.names = FALSE)

summary <- data.frame(
  method = "Gaussian probit Gibbs",
  H = H,
  lambda_l1_penalty = lambda_l1_penalty,
  mean_loglik = sum(fold_scores$total_loglik) / sum(fold_scores$n_test_entries),
  se_loglik = sd(fold_scores$mean_loglik) / sqrt(nrow(fold_scores)),
  mean_accuracy = weighted.mean(fold_scores$accuracy, fold_scores$n_test_entries),
  mean_brier = weighted.mean(fold_scores$brier, fold_scores$n_test_entries),
  n_completed_folds = nrow(fold_scores),
  heldout_factor_failures = sum(fold_scores$heldout_factor_failures),
  mean_fit_seconds = mean(fold_scores$fit_seconds),
  total_fit_seconds = sum(fold_scores$fit_seconds),
  wall_time_proxy_seconds = max(fold_scores$fit_seconds)
)
write.csv(summary, file.path(out_dir, "gaussian_rowwise_summary.csv"), row.names = FALSE)

pifa_results <- read.csv(file.path(pifa_cv_dir, "rowwise_cv_results.csv"), stringsAsFactors = FALSE)
pifa_comparison <- pifa_results[
  pifa_results$H == H & pifa_results$lambda_l1_penalty == lambda_l1_penalty &
    pifa_results$G_config %in% c("1,2,2,3", "1,3,3,3"),
  c("candidate_id", "H", "G_config", "lambda_l1_penalty", "mean_loglik",
    "se_loglik", "mean_accuracy", "mean_brier")
]
comparison <- rbind(
  data.frame(
    method = paste0("P-IFA G=", pifa_comparison$G_config),
    H = pifa_comparison$H,
    lambda_l1_penalty = pifa_comparison$lambda_l1_penalty,
    mean_loglik = pifa_comparison$mean_loglik,
    se_loglik = pifa_comparison$se_loglik,
    mean_accuracy = pifa_comparison$mean_accuracy,
    mean_brier = pifa_comparison$mean_brier
  ),
  summary[, c("method", "H", "lambda_l1_penalty", "mean_loglik", "se_loglik",
              "mean_accuracy", "mean_brier")]
)
comparison$delta_loglik_from_gaussian <- comparison$mean_loglik - summary$mean_loglik
write.csv(comparison, file.path(out_dir, "pifa_gaussian_rowwise_comparison.csv"), row.names = FALSE)

pifa_fold_scores <- read.csv(
  file.path(pifa_cv_dir, "rowwise_cv_fold_scores.csv"),
  stringsAsFactors = FALSE
)
paired_folds <- do.call(rbind, lapply(seq_len(nrow(pifa_comparison)), function(i) {
  candidate_fold <- pifa_fold_scores[
    pifa_fold_scores$H == H &
      pifa_fold_scores$G_config == pifa_comparison$G_config[i] &
      pifa_fold_scores$lambda_l1_penalty == lambda_l1_penalty,
    ,
    drop = FALSE
  ]
  candidate_fold <- candidate_fold[order(candidate_fold$fold), , drop = FALSE]
  gaussian_fold <- fold_scores[order(fold_scores$fold), , drop = FALSE]
  data.frame(
    G_config = pifa_comparison$G_config[i],
    fold = candidate_fold$fold,
    pifa_loglik = candidate_fold$mean_loglik,
    gaussian_loglik = gaussian_fold$mean_loglik,
    loglik_difference = candidate_fold$mean_loglik - gaussian_fold$mean_loglik,
    pifa_accuracy = candidate_fold$accuracy,
    gaussian_accuracy = gaussian_fold$accuracy,
    pifa_brier = candidate_fold$brier,
    gaussian_brier = gaussian_fold$brier
  )
}))
write.csv(
  paired_folds,
  file.path(out_dir, "pifa_gaussian_paired_fold_scores.csv"),
  row.names = FALSE
)

print(summary)
print(comparison)
cat("Gaussian row-wise outputs saved in: ", normalizePath(out_dir), "\n", sep = "")
