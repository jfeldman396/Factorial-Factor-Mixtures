#!/usr/bin/env Rscript

# Product-mixture pretraining wrapper using the hybrid soft-threshold EM-SVD
# signal estimator.
#
# The only new ingredient relative to R/probit_ifa_em_svd_pretraining.R is the
# first-stage estimate of the centered probit signal L.  After that point, the
# pipeline is intentionally identical: take the rank-H SVD of L, rotate the
# scores toward independent marginal mixtures, estimate Lambda from the rotated
# scores, and pass the result to the same MAP refinement routine.

fit_lowrank_probit_em_svd_soft <- function(
    X,
    H,
    shrinkage = NULL,
    shrinkage_ratio = 0.5,
    shrinkage_reference = c("dHplus1", "dH", "sqrt_np"),
    rank_cap = NULL,
    max_iter = 50L,
    tol_loglik = 1e-5,
    tol_L = 1e-4,
    seed = 1L,
    verbose = FALSE) {
  shrinkage_reference <- match.arg(shrinkage_reference)
  set.seed(seed)
  if (is.null(shrinkage)) {
    ref <- initial_soft_svd_shrinkage_reference(
      X = X,
      H = H,
      reference = shrinkage_reference
    )
    shrinkage <- as.numeric(shrinkage_ratio) * ref
  } else {
    ref <- NA_real_
  }
  fit <- fit_probit_signal_em_svd_soft_one(
    X = X,
    H = H,
    shrinkage = shrinkage,
    rank_cap = rank_cap,
    max_iter = max_iter,
    tol_loglik = tol_loglik,
    tol_L = tol_L,
    verbose = verbose
  )
  list(
    alpha = fit$alpha_hat,
    L = fit$L_hat,
    probit_loglik = fit$loglik,
    history = fit$history,
    converged = fit$converged,
    n_completed = fit$iterations,
    shrinkage = shrinkage,
    shrinkage_reference = shrinkage_reference,
    shrinkage_reference_value = ref,
    shrinkage_ratio = if (is.finite(ref)) shrinkage / max(ref, .Machine$double.eps) else NA_real_,
    rank_cap = rank_cap
  )
}

pretrain_probit_ifa_em_svd_soft <- function(
    X,
    H,
    G_fixed,
    em_max_iter = 50L,
    em_tol_loglik = 1e-5,
    em_tol_L = 1e-4,
    soft_shrinkage = NULL,
    soft_shrinkage_ratio = 0.5,
    soft_shrinkage_reference = c("dHplus1", "dH", "sqrt_np"),
    soft_rank_cap = NULL,
    rotation_loading_l1_penalty = 0,
    rotation_random_starts = 3L,
    rotation_ica_starts = 0L,
    rotation_ica_functions = c("logcosh", "exp"),
    rotation_ica_max_iter = 200L,
    rotation_ica_tol = 1e-4,
    rotation_max_outer = 5L,
    n_mix_starts = 5L,
    mixture_max_iter = 50L,
    mixture_update = c("map", "mle"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 1.5,
    weight_prior_alpha = 1,
    grid_size = 31L,
    rotation_objective_tolerance = 1e-4,
    rotation_min_outer = 2L,
    require_mixture_convergence_for_rotation_stop = TRUE,
    rotation_optimizer = c("riemannian", "givens"),
    rotation_sweep = c("full", "hybrid", "multi_disjoint", "disjoint", "promising"),
    riemannian_rotation_steps = 10L,
    riemannian_eta0 = 1,
    riemannian_beta = 0.5,
    riemannian_min_eta = 1e-8,
    riemannian_grad_tol = 1e-6,
    riemannian_update = c("cayley", "expm"),
    loading_penalty = 0,
    seed = 1L,
    parallel = FALSE,
    workers = NULL,
    verbose = FALSE) {
  X <- as.matrix(X)
  H <- as.integer(H)
  G_fixed <- normalize_G_fixed(G_fixed, H)
  mixture_update <- match.arg(mixture_update)
  rotation_optimizer <- match.arg(rotation_optimizer)
  rotation_sweep <- match.arg(rotation_sweep)
  riemannian_update <- match.arg(riemannian_update)
  soft_shrinkage_reference <- match.arg(soft_shrinkage_reference)

  lowrank <- fit_lowrank_probit_em_svd_soft(
    X = X,
    H = H,
    shrinkage = soft_shrinkage,
    shrinkage_ratio = soft_shrinkage_ratio,
    shrinkage_reference = soft_shrinkage_reference,
    rank_cap = soft_rank_cap,
    max_iter = em_max_iter,
    tol_loglik = em_tol_loglik,
    tol_L = em_tol_L,
    seed = seed,
    verbose = verbose
  )

  spectral <- spectral_scores_from_lowrank_signal(lowrank$L, H)
  rotation <- rotate_em_svd_scores_with_mixtures(
    S = spectral$S,
    G_fixed = G_fixed,
    loading_basis = spectral$B,
    rotation_loading_l1_penalty = rotation_loading_l1_penalty,
    n_random_starts = rotation_random_starts,
    n_ica_starts = rotation_ica_starts,
    ica_functions = rotation_ica_functions,
    ica_max_iter = rotation_ica_max_iter,
    ica_tol = rotation_ica_tol,
    max_outer = rotation_max_outer,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha,
    grid_size = grid_size,
    rotation_objective_tolerance = rotation_objective_tolerance,
    rotation_min_outer = rotation_min_outer,
    require_mixture_convergence_for_rotation_stop = require_mixture_convergence_for_rotation_stop,
    rotation_optimizer = rotation_optimizer,
    rotation_sweep = rotation_sweep,
    riemannian_rotation_steps = riemannian_rotation_steps,
    riemannian_eta0 = riemannian_eta0,
    riemannian_beta = riemannian_beta,
    riemannian_min_eta = riemannian_min_eta,
    riemannian_grad_tol = riemannian_grad_tol,
    riemannian_update = riemannian_update,
    seed = seed + 10000L,
    parallel = parallel,
    workers = workers,
    verbose = verbose
  )

  F_hat <- rotation$F_hat
  Lambda_hat <- spectral$B %*% rotation$R
  if (loading_penalty > 0) {
    Lambda_hat <- soft_threshold(Lambda_hat, loading_penalty)
  }
  rownames(Lambda_hat) <- colnames(X)
  colnames(Lambda_hat) <- paste0("factor_", seq_len(H))

  responsibilities <- lapply(seq_len(H), function(h) {
    mixture_responsibilities(F_hat[, h], rotation$fits[[h]])
  })
  class_map <- sapply(responsibilities, max.col, ties.method = "first")
  if (H == 1L) class_map <- matrix(class_map, ncol = 1L)
  colnames(class_map) <- paste0("factor_", seq_len(H))

  fitted <- sweep(F_hat %*% t(Lambda_hat), 2L, lowrank$alpha, "+")
  residual <- probit_ifa_truncated_mean(X, lowrank$alpha, lowrank$L) - fitted

  history <- data.frame(
    stage = c("lowrank_em_svd_soft", "mixture_rotation"),
    iteration = c(lowrank$n_completed, rotation$rotation_completed_outer),
    probit_loglik = c(lowrank$probit_loglik, lowrank$probit_loglik),
    mixture_loglik = c(NA_real_, rotation$loglik),
    objective = c(lowrank$probit_loglik, lowrank$probit_loglik + rotation$loglik)
  )

  list(
    model = "probit_ifa_em_svd_soft_spectral_mixture_pretraining",
    X = X,
    H = H,
    H_selection = NULL,
    H_selection_strategy = "fixed_supplied_H",
    S = spectral$S,
    svd_fit = spectral,
    lowrank_fit = lowrank,
    alpha_hat = lowrank$alpha,
    L_hat = lowrank$L,
    R = rotation$R,
    F_hat = F_hat,
    Lambda_hat = Lambda_hat,
    Lambda_ls = Lambda_hat,
    loading_update = "em_svd_soft_lowrank_signal",
    loading_penalty = loading_penalty,
    mixture_fits = rotation$fits,
    G_hat = rotation$G_hat,
    G_fixed = G_fixed,
    pretrain_G_selection = "fixed",
    mixture_update = mixture_update,
    rotation_optimizer = rotation_optimizer,
    rotation_fit = rotation,
    responsibilities = responsibilities,
    class_map = class_map,
    profile_id = binary_profile_id(class_map),
    fitted = fitted,
    residual = residual,
    Psi_hat = diag(1, ncol(X)),
    Psi_hat_unconstrained = diag(pmax(colMeans(residual^2), 1e-4), ncol(X)),
    history = history,
    em_history = lowrank$history,
    pretraining_converged = isTRUE(lowrank$converged),
    pretraining_completed_iter = lowrank$n_completed,
    selected_pretraining_iteration = lowrank$n_completed,
    em_init_method = "soft_threshold_svd",
    em_init_z = "expectation",
    selected_em_start = "soft_threshold_svd",
    all_em_start_loglik = lowrank$probit_loglik,
    soft_shrinkage = lowrank$shrinkage,
    soft_shrinkage_ratio = lowrank$shrinkage_ratio,
    soft_shrinkage_reference = lowrank$shrinkage_reference,
    soft_rank_cap = lowrank$rank_cap,
    rotation_completed_outer = rotation$rotation_completed_outer,
    rotation_converged = isTRUE(rotation$rotation_converged),
    rotation_history = rotation$rotation_history,
    rotation_step_history = rotation$rotation_step_history,
    z_update = "none_em_svd_soft",
    fix_psi_identity = TRUE,
    estimate_intercept = TRUE,
    call = match.call()
  )
}

fit_binary_probit_em_svd_soft_pretrain_then_refine <- function(
    X,
    H,
    G_fixed,
    n_refine_iter = 5L,
    n_mix_starts = 3L,
    mixture_max_iter = 20L,
    em_max_iter = 50L,
    em_tol_loglik = 1e-5,
    em_tol_L = 1e-4,
    mixture_update = c("map", "mle"),
    mixture_refit = c("em", "fixed_responsibility_mstep", "posterior_moment_mstep"),
    mu_prior_mean = 0,
    mu_prior_kappa = 0.01,
    var_prior_shape = 2,
    var_prior_scale = 1.5,
    weight_prior_alpha = 1,
    refine_mu_prior_mean = NULL,
    refine_mu_prior_kappa = NULL,
    refine_var_prior_shape = NULL,
    refine_var_prior_scale = NULL,
    refine_weight_prior_alpha = NULL,
    factor_update = c("marginal", "conditional_soft", "conditional_hard"),
    min_mixture_var = 1e-3,
    lambda_l1_penalty = 0,
    lasso_backend = c("proximal", "glmnet"),
    glmnet_standardize = FALSE,
    objective_tolerance = 1e-5,
    min_refine_iter = 1L,
    enforce_monotone_refinement = TRUE,
    monotone_tolerance = 1e-8,
    return_best_refinement_iteration = FALSE,
    refinement_selection_objective = c("posterior_objective", "joint_objective", "binary_loglik", "mixture_loglik"),
    require_mixture_convergence_for_stop = FALSE,
    soft_shrinkage = NULL,
    soft_shrinkage_ratio = 0.5,
    soft_shrinkage_reference = c("dHplus1", "dH", "sqrt_np"),
    soft_rank_cap = NULL,
    seed = 1L,
    parallel = FALSE,
    workers = NULL,
    verbose = FALSE,
    ...) {
  mixture_update <- match.arg(mixture_update)
  mixture_refit <- match.arg(mixture_refit)
  factor_update <- match.arg(factor_update)
  refinement_selection_objective <- match.arg(refinement_selection_objective)
  lasso_backend <- match.arg(lasso_backend)
  soft_shrinkage_reference <- match.arg(soft_shrinkage_reference)
  if (is.null(refine_mu_prior_mean)) refine_mu_prior_mean <- mu_prior_mean
  if (is.null(refine_mu_prior_kappa)) refine_mu_prior_kappa <- mu_prior_kappa
  if (is.null(refine_var_prior_shape)) refine_var_prior_shape <- var_prior_shape
  if (is.null(refine_var_prior_scale)) refine_var_prior_scale <- var_prior_scale
  if (is.null(refine_weight_prior_alpha)) refine_weight_prior_alpha <- weight_prior_alpha

  pretrain_fit <- pretrain_probit_ifa_em_svd_soft(
    X = X,
    H = H,
    G_fixed = G_fixed,
    em_max_iter = em_max_iter,
    em_tol_loglik = em_tol_loglik,
    em_tol_L = em_tol_L,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mu_prior_mean = mu_prior_mean,
    mu_prior_kappa = mu_prior_kappa,
    var_prior_shape = var_prior_shape,
    var_prior_scale = var_prior_scale,
    weight_prior_alpha = weight_prior_alpha,
    soft_shrinkage = soft_shrinkage,
    soft_shrinkage_ratio = soft_shrinkage_ratio,
    soft_shrinkage_reference = soft_shrinkage_reference,
    soft_rank_cap = soft_rank_cap,
    seed = seed,
    parallel = parallel,
    workers = workers,
    verbose = verbose,
    ...
  )

  refine_fit <- fit_binary_probit_refinement(
    X = X,
    pretrain_fit = pretrain_fit,
    n_refine_iter = n_refine_iter,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mixture_refit = mixture_refit,
    mu_prior_mean = refine_mu_prior_mean,
    mu_prior_kappa = refine_mu_prior_kappa,
    var_prior_shape = refine_var_prior_shape,
    var_prior_scale = refine_var_prior_scale,
    weight_prior_alpha = refine_weight_prior_alpha,
    factor_update = factor_update,
    min_mixture_var = min_mixture_var,
    lambda_l1_penalty = lambda_l1_penalty,
    lasso_backend = lasso_backend,
    glmnet_standardize = glmnet_standardize,
    objective_tolerance = objective_tolerance,
    min_refine_iter = min_refine_iter,
    stopping_objective = "posterior_objective",
    enforce_monotone_refinement = enforce_monotone_refinement,
    monotone_tolerance = monotone_tolerance,
    return_best_refinement_iteration = return_best_refinement_iteration,
    refinement_selection_objective = refinement_selection_objective,
    require_mixture_convergence_for_stop = require_mixture_convergence_for_stop,
    parallel = parallel,
    workers = workers,
    verbose = verbose,
    ...
  )

  list(pretrain_fit = pretrain_fit, refine_fit = refine_fit)
}
