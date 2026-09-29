# Few-shot row-wise cross-validation for the P-IFA MAP estimator.
#
# The training unit is an LLM (a row of the binary response matrix). Each
# training-fold model is fitted without the held-out rows. A small revealed
# subset of each held-out row is then used only to estimate that row's factor
# vector; the remaining observed responses are used only for prediction.

.pifa_rowwise_cv_root <- local({
  files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1L))
  files <- files[!is.na(files)]
  if (length(files)) {
    normalizePath(file.path(dirname(tail(files, 1L)), ".."), mustWork = FALSE)
  } else {
    normalizePath(getwd(), mustWork = FALSE)
  }
})

.old_functions_only <- getOption("pifa.cv.functions_only")
options(pifa.cv.functions_only = TRUE)
source(file.path(
  .pifa_rowwise_cv_root,
  "scripts", "ifeval", "cv_ifeval_rank_lambda_models.R"
))
options(pifa.cv.functions_only = .old_functions_only)
rm(.old_functions_only)

`%||%` <- function(x, y) if (is.null(x)) y else x

rbind_fill <- function(rows) {
  rows <- Filter(function(x) !is.null(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(columns, names(x))
    for (column in missing) x[[column]] <- NA
    x[columns]
  })
  do.call(rbind, rows)
}

valid_pifa_G_configs <- function(H, G_values = 1:3) {
  H <- as.integer(H)
  G_values <- sort(unique(as.integer(G_values)))
  if (length(H) != 1L || is.na(H) || H < 1L) stop("H must be a positive integer.")
  if (!length(G_values) || any(is.na(G_values)) || any(G_values < 1L)) {
    stop("G_values must contain positive integers.")
  }
  configs <- expand_columnwise_G_configs(
    H = H,
    component_values = G_values,
    max_gaussian_coords = 1L,
    unique_up_to_permutation = TRUE
  )
  stopifnot(
    all(vapply(configs, function(G) sum(G == 1L) <= 1L, logical(1L))),
    !anyDuplicated(vapply(configs, G_config_label, character(1L)))
  )
  configs
}

make_pifa_candidate_grid <- function(H_grid = 2:8,
                                     G_values = 1:3,
                                     lambda_grid = c(0, 1, 2, 4, 8, 12)) {
  rows <- list()
  k <- 0L
  for (H in as.integer(H_grid)) {
    for (G in valid_pifa_G_configs(H, G_values)) {
      for (lambda in as.numeric(lambda_grid)) {
        k <- k + 1L
        rows[[k]] <- data.frame(
          candidate_id = sprintf(
            "H%d_G%s_lambda%s",
            H,
            gsub(",", "-", G_config_label(G)),
            gsub("\\.", "p", as.character(lambda))
          ),
          H = H,
          G_config = G_config_label(G),
          lambda_l1_penalty = lambda,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

make_pifa_row_folds <- function(n, K = 5L, seed = 1L) {
  n <- as.integer(n)
  K <- as.integer(K)
  if (n < 2L || K < 2L || K > n) stop("Require 2 <= K <= n.")
  set.seed(seed)
  sample(rep(seq_len(K), length.out = n), n, replace = FALSE)
}

make_pifa_reveal_plan <- function(X,
                                  fold_assignments,
                                  reveal_fraction = 0.20,
                                  n_reveal_repeats = 3L,
                                  seed = 1L) {
  X <- as.matrix(X)
  if (length(fold_assignments) != nrow(X)) {
    stop("fold_assignments must have one entry per row of X.")
  }
  if (!is.finite(reveal_fraction) || reveal_fraction <= 0 || reveal_fraction >= 1) {
    stop("reveal_fraction must lie strictly between zero and one.")
  }
  n_reveal_repeats <- as.integer(n_reveal_repeats)
  if (n_reveal_repeats < 1L) stop("n_reveal_repeats must be positive.")

  set.seed(seed)
  out <- vector("list", nrow(X) * n_reveal_repeats)
  cursor <- 0L
  for (i in seq_len(nrow(X))) {
    observed <- which(!is.na(X[i, ]))
    if (length(observed) < 2L) {
      stop("Every row must have at least two observed responses; row ", i, " does not.")
    }
    n_reveal <- min(
      length(observed) - 1L,
      max(1L, as.integer(round(reveal_fraction * length(observed))))
    )
    for (repeat_id in seq_len(n_reveal_repeats)) {
      cursor <- cursor + 1L
      reveal <- sort(sample(observed, n_reveal, replace = FALSE))
      predict <- setdiff(observed, reveal)
      out[[cursor]] <- list(
        fold = as.integer(fold_assignments[i]),
        row_index = i,
        row_id = if (is.null(rownames(X))) as.character(i) else rownames(X)[i],
        reveal_repeat = repeat_id,
        reveal_index = reveal,
        prediction_index = predict,
        n_observed = length(observed),
        n_revealed = length(reveal),
        n_predicted = length(predict)
      )
    }
  }
  out
}

pifa_reveal_plan_index <- function(reveal_plan) {
  do.call(rbind, lapply(reveal_plan, function(z) {
    data.frame(
      fold = z$fold,
      row_index = z$row_index,
      row_id = z$row_id,
      reveal_repeat = z$reveal_repeat,
      n_observed = z$n_observed,
      n_revealed = z$n_revealed,
      n_predicted = z$n_predicted,
      reveal_index = paste(z$reveal_index, collapse = "|"),
      prediction_index = paste(z$prediction_index, collapse = "|"),
      stringsAsFactors = FALSE
    )
  }))
}

component_profile_starts <- function(mixture_fits, max_profiles = 2L) {
  H <- length(mixture_fits)
  marginal_mean <- vapply(mixture_fits, function(fit) {
    sum(fit$pi * fit$mu)
  }, numeric(1L))
  starts <- list(rep(0, H), marginal_mean)

  max_profiles <- max(0L, as.integer(max_profiles))
  if (max_profiles > 0L) {
    component_grid <- expand.grid(lapply(mixture_fits, function(fit) seq_along(fit$pi)))
    profile_weight <- apply(component_grid, 1L, function(index) {
      prod(vapply(seq_len(H), function(h) {
        mixture_fits[[h]]$pi[index[h]]
      }, numeric(1L)))
    })
    ord <- order(profile_weight, decreasing = TRUE)
    for (index in ord[seq_len(min(max_profiles, length(ord)))]) {
      component_index <- as.integer(component_grid[index, ])
      starts[[length(starts) + 1L]] <- vapply(seq_len(H), function(h) {
        mixture_fits[[h]]$mu[component_index[h]]
      }, numeric(1L))
    }
  }
  starts[!duplicated(vapply(starts, function(z) paste(signif(z, 12), collapse = ","), character(1L)))]
}

estimate_heldout_factor_map <- function(x_i,
                                        reveal_index,
                                        fit,
                                        mixture_prior_weight = 1,
                                        maxit = 100L,
                                        factor_score_bound = 5,
                                        max_component_starts = 2L) {
  x_i <- as.numeric(x_i)
  reveal_index <- as.integer(reveal_index)
  if (!length(reveal_index) || anyNA(x_i[reveal_index])) {
    stop("The reveal set must contain at least one non-missing response.")
  }
  reveal_mask <- seq_along(x_i) %in% reveal_index
  starts <- component_profile_starts(fit$mixture_fits, max_component_starts)
  attempts <- lapply(starts, function(start) {
    update_one_factor_score_mixture_missing(
      x_i = x_i,
      obs_i = reveal_mask,
      f_init = start,
      Lambda = fit$Lambda,
      alpha = fit$alpha,
      mixture_fits = fit$mixture_fits,
      mixture_prior_weight = mixture_prior_weight,
      maxit = maxit,
      factor_score_bound = factor_score_bound,
      return_diagnostics = TRUE
    )
  })
  successful <- vapply(attempts, function(z) {
    identical(z$convergence, 0L) && is.finite(z$value) && all(is.finite(z$par))
  }, logical(1L))
  if (!any(successful)) {
    return(list(
      eta_hat = rep(NA_real_, ncol(fit$Lambda)),
      converged = FALSE,
      convergence_code = 100L,
      objective = NA_real_,
      iterations = NA_integer_,
      n_starts = length(attempts),
      selected_start = NA_integer_,
      message = "all held-out factor starts failed"
    ))
  }
  valid <- which(successful)
  selected <- valid[which.min(vapply(attempts[valid], `[[`, numeric(1L), "value"))]
  best <- attempts[[selected]]
  list(
    eta_hat = best$par,
    converged = TRUE,
    convergence_code = best$convergence,
    objective = -best$value,
    iterations = unname(best$counts[["function"]] %||% NA_integer_),
    n_starts = length(attempts),
    selected_start = selected,
    message = best$message
  )
}

score_pifa_heldout_rows <- function(X,
                                    fold,
                                    reveal_plan,
                                    fit,
                                    candidate_id,
                                    mixture_prior_weight = 1,
                                    heldout_maxit = 100L,
                                    factor_score_bound = 5,
                                    max_component_starts = 2L,
                                    probability_epsilon = 1e-10) {
  units <- Filter(function(z) identical(z$fold, as.integer(fold)), reveal_plan)
  H <- ncol(fit$Lambda)
  score_rows <- vector("list", length(units))
  diagnostic_rows <- vector("list", length(units))

  for (u in seq_along(units)) {
    unit <- units[[u]]
    factor_fit <- estimate_heldout_factor_map(
      x_i = X[unit$row_index, ],
      reveal_index = unit$reveal_index,
      fit = fit,
      mixture_prior_weight = mixture_prior_weight,
      maxit = heldout_maxit,
      factor_score_bound = factor_score_bound,
      max_component_starts = max_component_starts
    )
    eta_columns <- as.list(stats::setNames(factor_fit$eta_hat, paste0("eta_", seq_len(H))))
    diagnostic_rows[[u]] <- cbind(
      data.frame(
        candidate_id = candidate_id,
        fold = fold,
        row_index = unit$row_index,
        row_id = unit$row_id,
        reveal_repeat = unit$reveal_repeat,
        converged = factor_fit$converged,
        convergence_code = factor_fit$convergence_code,
        objective = factor_fit$objective,
        optimizer_iterations = factor_fit$iterations,
        n_starts = factor_fit$n_starts,
        selected_start = factor_fit$selected_start,
        n_revealed = unit$n_revealed,
        n_predicted = unit$n_predicted,
        message = factor_fit$message,
        stringsAsFactors = FALSE
      ),
      as.data.frame(eta_columns, optional = TRUE)
    )

    if (!factor_fit$converged) next
    prediction_index <- unit$prediction_index
    linear_predictor <- fit$alpha[prediction_index] +
      as.numeric(fit$Lambda[prediction_index, , drop = FALSE] %*% factor_fit$eta_hat)
    probability <- pmin(
      pmax(pnorm(linear_predictor), probability_epsilon),
      1 - probability_epsilon
    )
    outcome <- as.numeric(X[unit$row_index, prediction_index])
    loglik <- outcome * log(probability) + (1 - outcome) * log1p(-probability)
    score_rows[[u]] <- data.frame(
      candidate_id = candidate_id,
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

  scores <- Filter(Negate(is.null), score_rows)
  list(
    unit_scores = if (length(scores)) do.call(rbind, scores) else data.frame(),
    diagnostics = do.call(rbind, diagnostic_rows)
  )
}

pifa_fit_control_defaults <- function() {
  list(
    n_aug_iter = 200L,
    n_refine_iter = 100L,
    z_update = "expectation",
    n_random_starts = 1L,
    n_ica_starts = 1L,
    rotation_optimizer = "riemannian",
    max_outer = 20L,
    n_mix_starts = 2L,
    mixture_max_iter = 200L,
    mixture_prior_weight = 1,
    maxit_per_subject = 45L,
    min_mixture_var = 0.05,
    refinement_objective_tolerance = 1e-3,
    refinement_min_iter = 2L,
    refinement_require_mixture_convergence = TRUE,
    refinement_monotone_tolerance = 1e-8,
    normalize_factor_scale = TRUE,
    factor_score_bound = 5
  )
}

fit_pifa_rowwise_training_task <- function(task,
                                           X,
                                           fold_assignments,
                                           reveal_plan,
                                           fit_control,
                                           seed,
                                           fit_workers = 1L,
                                           save_fit = FALSE) {
  candidate_id <- task$candidate_id
  fold <- as.integer(task$fold)
  H <- as.integer(task$H)
  G <- as.integer(strsplit(task$G_config, ",", fixed = TRUE)[[1L]])
  lambda <- as.numeric(task$lambda_l1_penalty)
  train_rows <- which(fold_assignments != fold)
  test_rows <- which(fold_assignments == fold)
  if (length(intersect(train_rows, test_rows))) stop("Row-fold leakage detected.")

  X_train <- X[train_rows, , drop = FALSE]
  W_train <- !is.na(X_train)
  X_train_work <- X_train
  X_train_work[!W_train] <- 0
  start_time <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    fit_mixture_missing_probit(
      X = X_train_work,
      W = W_train,
      H = H,
      G = G,
      lambda_l1_penalty = lambda,
      fold = fold,
      workers = fit_workers,
      seed = seed,
      n_aug_iter = fit_control$n_aug_iter,
      n_refine_iter = fit_control$n_refine_iter,
      z_update = fit_control$z_update,
      loading_penalty = lambda,
      n_random_starts = fit_control$n_random_starts,
      n_ica_starts = fit_control$n_ica_starts,
      rotation_optimizer = fit_control$rotation_optimizer,
      max_outer = fit_control$max_outer,
      n_mix_starts = fit_control$n_mix_starts,
      mixture_max_iter = fit_control$mixture_max_iter,
      mixture_prior_weight = fit_control$mixture_prior_weight,
      maxit_per_subject = fit_control$maxit_per_subject,
      min_mixture_var = fit_control$min_mixture_var,
      refinement_objective_tolerance = fit_control$refinement_objective_tolerance,
      refinement_min_iter = fit_control$refinement_min_iter,
      refinement_require_mixture_convergence = fit_control$refinement_require_mixture_convergence,
      refinement_monotone_tolerance = fit_control$refinement_monotone_tolerance,
      normalize_factor_scale = fit_control$normalize_factor_scale,
      factor_score_bound = fit_control$factor_score_bound
    ),
    error = function(e) e
  )
  fit_seconds <- proc.time()[["elapsed"]] - start_time
  if (inherits(fit, "error")) {
    return(list(
      fold_score = data.frame(
        candidate_id = candidate_id,
        H = H,
        G_config = task$G_config,
        lambda_l1_penalty = lambda,
        fold = fold,
        total_loglik = NA_real_,
        n_test_entries = 0L,
        mean_loglik = NA_real_,
        accuracy = NA_real_,
        brier = NA_real_,
        training_fit_converged = FALSE,
        heldout_factor_failures = length(test_rows),
        fit_seconds = fit_seconds,
        evaluation_seconds = 0,
        error = conditionMessage(fit),
        stringsAsFactors = FALSE
      ),
      unit_scores = data.frame(),
      diagnostics = data.frame(),
      fit = NULL
    ))
  }

  training_fit_converged <- isTRUE(fit$pretraining_converged) &&
    isTRUE(fit$rotation_converged) &&
    isTRUE(fit$refinement_converged) &&
    isTRUE(fit$rotation_all_mixtures_converged) &&
    (!isTRUE(fit_control$refinement_require_mixture_convergence) ||
       isTRUE(fit$refinement_all_mixtures_converged))
  evaluation_start <- proc.time()[["elapsed"]]
  scored <- score_pifa_heldout_rows(
    X = X,
    fold = fold,
    reveal_plan = reveal_plan,
    fit = fit,
    candidate_id = candidate_id,
    mixture_prior_weight = fit_control$mixture_prior_weight,
    heldout_maxit = fit_control$heldout_maxit %||% 100L,
    factor_score_bound = fit_control$factor_score_bound,
    max_component_starts = fit_control$heldout_max_component_starts %||% 2L
  )
  evaluation_seconds <- proc.time()[["elapsed"]] - evaluation_start
  unit_scores <- scored$unit_scores
  n_test_entries <- if (nrow(unit_scores)) sum(unit_scores$n_test_entries) else 0L
  total_loglik <- if (nrow(unit_scores)) sum(unit_scores$total_loglik) else NA_real_
  total_squared_error <- if (nrow(unit_scores)) sum(unit_scores$total_squared_error) else NA_real_
  total_correct <- if (nrow(unit_scores)) sum(unit_scores$total_correct) else NA_real_
  heldout_failures <- sum(!scored$diagnostics$converged)

  list(
    fold_score = data.frame(
      candidate_id = candidate_id,
      H = H,
      G_config = task$G_config,
      lambda_l1_penalty = lambda,
      fold = fold,
      total_loglik = total_loglik,
      n_test_entries = n_test_entries,
      mean_loglik = total_loglik / n_test_entries,
      accuracy = total_correct / n_test_entries,
      brier = total_squared_error / n_test_entries,
      training_fit_converged = training_fit_converged,
      heldout_factor_failures = heldout_failures,
      fit_seconds = fit_seconds,
      evaluation_seconds = evaluation_seconds,
      error = "",
      stringsAsFactors = FALSE
    ),
    unit_scores = unit_scores,
    diagnostics = scored$diagnostics,
    fit = if (isTRUE(save_fit)) fit else NULL
  )
}

summarize_pifa_rowwise_scores <- function(fold_scores, K) {
  candidate_ids <- unique(fold_scores$candidate_id)
  out <- lapply(candidate_ids, function(candidate_id) {
    d <- fold_scores[fold_scores$candidate_id == candidate_id, , drop = FALSE]
    valid <- is.finite(d$mean_loglik) & d$n_test_entries > 0L
    dv <- d[valid, , drop = FALSE]
    data.frame(
      candidate_id = candidate_id,
      H = d$H[1L],
      G_config = d$G_config[1L],
      lambda_l1_penalty = d$lambda_l1_penalty[1L],
      mean_loglik = if (nrow(dv)) sum(dv$total_loglik) / sum(dv$n_test_entries) else NA_real_,
      total_loglik = if (nrow(dv)) sum(dv$total_loglik) else NA_real_,
      n_test_entries = if (nrow(dv)) sum(dv$n_test_entries) else 0L,
      se_loglik = if (nrow(dv) > 1L) stats::sd(dv$mean_loglik) / sqrt(nrow(dv)) else NA_real_,
      mean_accuracy = if (nrow(dv)) weighted.mean(dv$accuracy, dv$n_test_entries) else NA_real_,
      mean_brier = if (nrow(dv)) weighted.mean(dv$brier, dv$n_test_entries) else NA_real_,
      n_completed_folds = nrow(dv),
      training_convergence_failures = sum(!d$training_fit_converged),
      heldout_factor_failures = sum(d$heldout_factor_failures),
      convergence_failures = sum(!d$training_fit_converged) + sum(d$heldout_factor_failures),
      elapsed_time = sum(d$fit_seconds + d$evaluation_seconds),
      eligible = nrow(dv) == K && all(d$training_fit_converged) &&
        sum(d$heldout_factor_failures) == 0L,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out[order(out$mean_loglik, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
}

cv_pifa_rowwise <- function(
    X,
    H_grid = 2:8,
    G_values = 1:3,
    lambda_grid = c(0, 1, 2, 4, 8, 12),
    K = 5,
    reveal_fraction = 0.20,
    n_reveal_repeats = 3,
    seed = 1,
    parallel = TRUE,
    workers = NULL,
    fit_control = list(),
    fit_selected_full_data = TRUE,
    save_training_fits = FALSE,
    resume = TRUE,
    task_batch_size = NULL,
    out_dir = NULL) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  if (!all(is.na(X) | X %in% c(0, 1))) stop("X must contain only 0, 1, or NA.")
  H_grid <- as.integer(H_grid)
  if (any(H_grid < 1L | H_grid >= nrow(X) | H_grid > ncol(X))) {
    stop("Each H must be positive and smaller than nrow(X) and no larger than ncol(X).")
  }
  workers <- resolve_workers(workers)
  fit_control <- utils::modifyList(pifa_fit_control_defaults(), fit_control)
  candidates <- make_pifa_candidate_grid(H_grid, G_values, lambda_grid)
  fold_assignments <- make_pifa_row_folds(nrow(X), K, seed)
  reveal_plan <- make_pifa_reveal_plan(
    X,
    fold_assignments,
    reveal_fraction,
    n_reveal_repeats,
    seed = seed + 1L
  )
  tasks <- merge(candidates, data.frame(fold = seq_len(K)), all = TRUE)
  tasks <- tasks[order(tasks$candidate_id, tasks$fold), , drop = FALSE]
  task_keys <- paste(tasks$candidate_id, tasks$fold, sep = "::fold=")
  checkpoint_path <- if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    file.path(out_dir, "rowwise_cv_task_checkpoint.rds")
  } else {
    NULL
  }
  checkpoint_meta <- list(
    dimensions = dim(X),
    observed_entries = sum(!is.na(X)),
    candidate_ids = candidates$candidate_id,
    K = K,
    reveal_fraction = reveal_fraction,
    n_reveal_repeats = n_reveal_repeats,
    seed = seed,
    fold_assignments = fold_assignments,
    fit_control = fit_control
  )
  task_results_by_key <- list()
  if (isTRUE(resume) && !is.null(checkpoint_path) && file.exists(checkpoint_path)) {
    checkpoint <- readRDS(checkpoint_path)
    if (!identical(checkpoint$meta, checkpoint_meta)) {
      stop("Existing row-wise CV checkpoint does not match the requested run settings.")
    }
    task_results_by_key <- checkpoint$results
    message("Resuming ", length(task_results_by_key), " completed tasks from checkpoint.")
  }
  missing_indices <- which(!task_keys %in% names(task_results_by_key))
  if (is.null(task_batch_size)) task_batch_size <- max(1L, workers)
  task_batch_size <- max(1L, as.integer(task_batch_size))

  if (length(missing_indices)) {
    batches <- split(
      missing_indices,
      ceiling(seq_along(missing_indices) / task_batch_size)
    )
    for (batch_number in seq_along(batches)) {
      batch_indices <- batches[[batch_number]]
      batch_results <- parallel_lapply(
        batch_indices,
        function(index) {
          fit_pifa_rowwise_training_task(
            task = tasks[index, , drop = FALSE],
            X = X,
            fold_assignments = fold_assignments,
            reveal_plan = reveal_plan,
            fit_control = fit_control,
            seed = seed + 100000L * index,
            fit_workers = 1L,
            save_fit = save_training_fits
          )
        },
        parallel = isTRUE(parallel),
        workers = min(workers, length(batch_indices))
      )
      names(batch_results) <- task_keys[batch_indices]
      task_results_by_key[names(batch_results)] <- batch_results
      if (!is.null(checkpoint_path)) {
        saveRDS(
          list(meta = checkpoint_meta, results = task_results_by_key),
          checkpoint_path
        )
      }
      message(
        "Completed row-wise CV batch ", batch_number, "/", length(batches),
        " (", length(task_results_by_key), "/", nrow(tasks), " tasks)."
      )
    }
  }
  task_results <- unname(task_results_by_key[task_keys])
  fold_scores <- rbind_fill(lapply(task_results, `[[`, "fold_score"))
  unit_scores_list <- lapply(task_results, `[[`, "unit_scores")
  unit_scores_list <- Filter(function(z) nrow(z) > 0L, unit_scores_list)
  unit_scores <- rbind_fill(unit_scores_list)
  diagnostics_list <- lapply(task_results, `[[`, "diagnostics")
  diagnostics_list <- Filter(function(z) nrow(z) > 0L, diagnostics_list)
  diagnostics <- rbind_fill(diagnostics_list)
  cv_results <- summarize_pifa_rowwise_scores(fold_scores, K)
  eligible <- cv_results[cv_results$eligible %in% TRUE, , drop = FALSE]
  if (!nrow(eligible)) {
    stop("No candidate completed all folds without convergence failures.")
  }
  best_model_spec <- eligible[which.max(eligible$mean_loglik), c(
    "candidate_id", "H", "G_config", "lambda_l1_penalty", "mean_loglik",
    "se_loglik", "mean_accuracy", "mean_brier"
  ), drop = FALSE]

  final_full_data_fit <- NULL
  if (isTRUE(fit_selected_full_data)) {
    W_full <- !is.na(X)
    X_work <- X
    X_work[!W_full] <- 0
    selected_G <- as.integer(strsplit(best_model_spec$G_config, ",", fixed = TRUE)[[1L]])
    final_full_data_fit <- fit_mixture_missing_probit(
      X = X_work,
      W = W_full,
      H = best_model_spec$H,
      G = selected_G,
      lambda_l1_penalty = best_model_spec$lambda_l1_penalty,
      fold = K + 1L,
      workers = if (isTRUE(parallel)) workers else 1L,
      seed = seed + 900000000L,
      n_aug_iter = fit_control$n_aug_iter,
      n_refine_iter = fit_control$n_refine_iter,
      z_update = fit_control$z_update,
      loading_penalty = best_model_spec$lambda_l1_penalty,
      n_random_starts = fit_control$n_random_starts,
      n_ica_starts = fit_control$n_ica_starts,
      rotation_optimizer = fit_control$rotation_optimizer,
      max_outer = fit_control$max_outer,
      n_mix_starts = fit_control$n_mix_starts,
      mixture_max_iter = fit_control$mixture_max_iter,
      mixture_prior_weight = fit_control$mixture_prior_weight,
      maxit_per_subject = fit_control$maxit_per_subject,
      min_mixture_var = fit_control$min_mixture_var,
      refinement_objective_tolerance = fit_control$refinement_objective_tolerance,
      refinement_min_iter = fit_control$refinement_min_iter,
      refinement_require_mixture_convergence = fit_control$refinement_require_mixture_convergence,
      refinement_monotone_tolerance = fit_control$refinement_monotone_tolerance,
      normalize_factor_scale = fit_control$normalize_factor_scale,
      factor_score_bound = fit_control$factor_score_bound
    )
  }

  out <- list(
    cv_results = cv_results,
    best_model_spec = best_model_spec,
    fold_assignments = fold_assignments,
    reveal_masks = reveal_plan,
    fold_level_scores = fold_scores,
    heldout_unit_scores = unit_scores,
    heldout_factor_diagnostics = diagnostics,
    candidate_grid = candidates,
    training_fits = if (isTRUE(save_training_fits)) lapply(task_results, `[[`, "fit") else NULL,
    final_full_data_fit = final_full_data_fit,
    settings = list(
      K = K,
      reveal_fraction = reveal_fraction,
      n_reveal_repeats = n_reveal_repeats,
      seed = seed,
      fit_control = fit_control
    )
  )
  class(out) <- c("pifa_rowwise_cv", "list")

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(cv_results, file.path(out_dir, "rowwise_cv_results.csv"), row.names = FALSE)
    utils::write.csv(best_model_spec, file.path(out_dir, "rowwise_cv_selected_model.csv"), row.names = FALSE)
    utils::write.csv(fold_scores, file.path(out_dir, "rowwise_cv_fold_scores.csv"), row.names = FALSE)
    utils::write.csv(diagnostics, file.path(out_dir, "rowwise_cv_factor_diagnostics.csv"), row.names = FALSE)
    utils::write.csv(pifa_reveal_plan_index(reveal_plan), file.path(out_dir, "rowwise_cv_reveal_plan_index.csv"), row.names = FALSE)
    utils::write.csv(
      data.frame(
        row_index = seq_along(fold_assignments),
        row_id = if (is.null(rownames(X))) as.character(seq_len(nrow(X))) else rownames(X),
        fold = fold_assignments
      ),
      file.path(out_dir, "rowwise_cv_fold_assignments.csv"),
      row.names = FALSE
    )
    saveRDS(fold_assignments, file.path(out_dir, "rowwise_cv_fold_assignments.rds"))
    saveRDS(reveal_plan, file.path(out_dir, "rowwise_cv_reveal_masks.rds"))
    if (!is.null(final_full_data_fit)) {
      saveRDS(final_full_data_fit, file.path(out_dir, "rowwise_cv_selected_full_data_fit.rds"))
    }
    saveRDS(out, file.path(out_dir, "rowwise_cv_object.rds"))
  }
  out
}

print.pifa_rowwise_cv <- function(x, ...) {
  cat("Few-shot row-wise P-IFA cross-validation\n")
  cat("Candidates:", nrow(x$cv_results), "\n")
  cat("Selected model:\n")
  print(x$best_model_spec, row.names = FALSE)
  invisible(x)
}
