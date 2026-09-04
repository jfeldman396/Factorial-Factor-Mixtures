#!/usr/bin/env Rscript

# Held-out cell CV for the IFEval independent-mixture probit factor model.
#
# This script tunes rank H, column-specific component counts G_h, and the
# entrywise Lambda penalty by held-out predictive log likelihood.
#
# The CV split removes individual model-by-item cells.  Held-out cells are not
# used in the augmentation, loading updates, or factor-score updates; they are
# scored only after fitting.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
ifeval_script_dir <- script_dir
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

source(file.path(script_dir, "ifeval_imfm_helpers.R"))

parse_int_grid <- function(x, default) {
  if (is.na(x) || !nzchar(x)) return(default)
  x <- gsub("[[:space:]]+", "", x)
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    bounds <- as.integer(strsplit(x, ":", fixed = TRUE)[[1L]])
    return(seq(bounds[1L], bounds[2L]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
}

parse_num_grid <- function(x, default) {
  if (is.na(x) || !nzchar(x)) return(default)
  as.numeric(strsplit(gsub("[[:space:]]+", "", x), ",", fixed = TRUE)[[1L]])
}

expand_columnwise_G_configs <- function(H,
                                        component_values = c(2L, 3L),
                                        max_gaussian_coords = 0L,
                                        unique_up_to_permutation = FALSE) {
  # Enumerate all per-factor component-count vectors.  For H = 4 and
  # component_values = c(2, 3), this returns the 16 vectors in {2,3}^4.
  # If component_values includes 1, max_gaussian_coords limits how many
  # coordinates are allowed to be Gaussian rather than non-Gaussian mixtures.
  # When unique_up_to_permutation is TRUE, retain only sorted component-count
  # multisets.  This avoids refitting likelihood-equivalent models that differ
  # only by factor-coordinate labels, e.g. (3,3,1,3) versus (3,1,3,3).
  grid <- expand.grid(rep(list(as.integer(component_values)), H))
  configs <- lapply(seq_len(nrow(grid)), function(i) as.integer(grid[i, ]))
  max_gaussian_coords <- as.integer(max_gaussian_coords)
  if (any(as.integer(component_values) == 1L) && is.finite(max_gaussian_coords)) {
    configs <- Filter(function(G) sum(G == 1L) <= max_gaussian_coords, configs)
  }
  if (isTRUE(unique_up_to_permutation)) {
    sorted_labels <- vapply(configs, function(G) paste(sort(G), collapse = ","), character(1L))
    configs <- configs[!duplicated(sorted_labels)]
    configs <- lapply(configs, sort)
  }
  configs
}

G_config_label <- function(G) {
  paste(as.integer(G), collapse = ",")
}

G_config_seed <- function(G) {
  sum(as.integer(G) * seq_along(G))
}

clip01_obs <- function(p, n) {
  eps <- 0.5 / max(1L, n)
  pmin(pmax(p, eps), 1 - eps)
}

read_binary_matrix <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- raw[[1L]]
  complete <- colSums(is.na(X)) == 0L
  X <- X[, complete, drop = FALSE]
  nonconstant <- colSums(X == 1) > 0L & colSums(X == 0) > 0L
  X <- X[, nonconstant, drop = FALSE]
  if (is.null(colnames(X))) colnames(X) <- paste0("item_", seq_len(ncol(X)))
  X
}

initialize_alpha_missing <- function(X, W) {
  n_obs <- colSums(W)
  p_obs <- colSums(X * W) / pmax(n_obs, 1L)
  qnorm(clip01_obs(p_obs, n_obs))
}

expected_Z_missing <- function(X, W, F_hat = NULL, Lambda = NULL, alpha = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(F_hat)) {
    alpha <- initialize_alpha_missing(X, W)
    eta <- matrix(rep(alpha, each = n), n, p, dimnames = dimnames(X))
  } else {
    eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  }

  Z <- eta
  for (j in seq_len(p)) {
    obs <- W[, j]
    if (!any(obs)) next
    mu <- eta[obs, j]
    y <- X[obs, j]
    lower <- ifelse(y == 1, 0, -Inf)
    upper <- ifelse(y == 1, Inf, 0)
    Z[obs, j] <- truncnorm_binary_moments_vec(mu, 1, lower, upper)$mean
  }
  dimnames(Z) <- dimnames(X)
  Z
}

sample_Z_missing <- function(X, W, F_hat = NULL, Lambda = NULL, alpha = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(F_hat)) {
    alpha <- initialize_alpha_missing(X, W)
    eta <- matrix(rep(alpha, each = n), n, p, dimnames = dimnames(X))
  } else {
    eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  }

  Z <- eta
  for (j in seq_len(p)) {
    obs <- W[, j]
    if (!any(obs)) next
    mu <- eta[obs, j]
    y <- X[obs, j]
    lower <- ifelse(y == 1, 0, -Inf)
    upper <- ifelse(y == 1, Inf, 0)
    Z[obs, j] <- rtruncnorm_binary_vec(mu, 1, lower, upper)
  }
  dimnames(Z) <- dimnames(X)
  Z
}

update_Z_missing <- function(X, W, F_hat = NULL, Lambda = NULL, alpha = NULL,
                             z_update = c("sample", "expectation"), seed = NULL) {
  z_update <- match.arg(z_update)
  if (z_update == "sample") {
    sample_Z_missing(X, W, F_hat, Lambda, alpha, seed = seed)
  } else {
    expected_Z_missing(X, W, F_hat, Lambda, alpha)
  }
}

binary_loglik_masked_alpha <- function(X, W, F_hat, Lambda, alpha) {
  eta <- sweep(F_hat %*% t(Lambda), 2L, alpha, "+")
  p1 <- pmax(pnorm(eta), 1e-12)
  p0 <- pmax(pnorm(-eta), 1e-12)
  ll <- X * log(p1) + (1 - X) * log(p0)
  sum(ll[W])
}

score_cells <- function(X, heldout, fit) {
  eta <- sweep(fit$F_hat %*% t(fit$Lambda), 2L, fit$alpha, "+")
  prob <- pmin(pmax(pnorm(eta), 1e-12), 1 - 1e-12)
  y <- X[heldout]
  phat <- prob[heldout]
  ll <- y * log(phat) + (1 - y) * log(1 - phat)
  data.frame(
    heldout_total_loglik = sum(ll),
    heldout_loglik_per_response = mean(ll),
    heldout_accuracy = mean(as.numeric(phat >= 0.5) == y),
    heldout_brier = mean((phat - y)^2),
    n_heldout = length(y)
  )
}

update_working_loadings_missing <- function(Z, W, F_hat, loading_penalty = 0.05) {
  p <- ncol(Z)
  H <- ncol(F_hat)
  alpha <- numeric(p)
  Lambda <- matrix(0, p, H, dimnames = list(colnames(Z), paste0("factor_", seq_len(H))))

  for (j in seq_len(p)) {
    obs <- W[, j]
    if (sum(obs) <= H + 1L) {
      alpha[j] <- if (any(obs)) mean(Z[obs, j]) else 0
      next
    }
    D <- cbind(1, F_hat[obs, , drop = FALSE])
    y <- Z[obs, j]
    ridge <- diag(c(0, rep(1e-5, H)), H + 1L)
    coef <- tryCatch(solve(crossprod(D) + ridge, crossprod(D, y)),
                     error = function(e) rep(0, H + 1L))
    alpha[j] <- coef[1L]
    Lambda[j, ] <- soft_threshold(coef[-1L], loading_penalty)
  }

  list(alpha = alpha, Lambda = Lambda)
}

update_loadings_probit_missing <- function(X, W, F_hat, Lambda_init, alpha_init,
                                           lambda_l1_penalty, workers) {
  p <- ncol(X)
  H <- ncol(F_hat)
  rows <- parallel_lapply(seq_len(p), function(j) {
    obs <- W[, j]
    y <- X[obs, j]
    if (length(unique(y)) < 2L) {
      a <- qnorm(clip01_obs(mean(y), length(y)))
      return(c(a, rep(0, H)))
    }
    fit_probit_lasso_item_alpha(
      y = y,
      F_hat = F_hat[obs, , drop = FALSE],
      theta_init = c(alpha_init[j], Lambda_init[j, ]),
      lambda_l1_penalty = lambda_l1_penalty,
      maxit = 160L,
      tol = 1e-5
    )
  }, parallel = workers > 1L, workers = workers)

  mat <- do.call(rbind, rows)
  list(
    alpha = as.numeric(mat[, 1L]),
    Lambda = matrix(
      mat[, -1L, drop = FALSE],
      nrow = p,
      ncol = H,
      dimnames = list(colnames(X), paste0("factor_", seq_len(H)))
    )
  )
}

initialize_mixture_from_Z <- function(Z, W, H, G, loading_penalty, workers,
                                      n_random_starts, max_outer,
                                      n_mix_starts, mixture_max_iter,
                                      mixture_update, seed) {
  svd_out <- svd_scores_from_Z(Z, H = H, center_Z = TRUE)
  G_fixed <- if (length(G) == 1L) rep(as.integer(G), H) else as.integer(G)
  if (length(G_fixed) != H) {
    stop("G must be either scalar or a length-H component-count vector.")
  }
  rot <- estimate_mixture_ica_unknown_G(
    S = svd_out$S,
    G_fixed = G_fixed,
    n_random_starts = n_random_starts,
    max_outer = max_outer,
    n_mix_starts = n_mix_starts,
    mixture_max_iter = mixture_max_iter,
    mixture_update = mixture_update,
    mu_prior_kappa = 0.05,
    var_prior_shape = 4,
    var_prior_scale = 0.35,
    weight_prior_alpha = 1.2,
    grid_size = 17L,
    seed = seed,
    parallel = workers > 1L,
    workers = workers,
    verbose = FALSE
  )
  load <- update_working_loadings_missing(Z, W, rot$F_hat, loading_penalty)
  list(F_hat = rot$F_hat, Lambda = load$Lambda, alpha = load$alpha,
       mixture_fits = rot$fits, R = rot$R)
}

update_one_factor_score_mixture_missing <- function(x_i, obs_i, f_init, Lambda,
                                                    alpha, mixture_fits,
                                                    mixture_prior_weight,
                                                    maxit = 50L) {
  y <- as.numeric(x_i[obs_i])
  L <- Lambda[obs_i, , drop = FALSE]
  a <- alpha[obs_i]
  H <- length(f_init)

  objective <- function(f) {
    eta <- a + as.numeric(L %*% f)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    ll <- sum(y * log(p1) + (1 - y) * log(p0))
    lp <- 0
    for (h in seq_len(H)) {
      lp <- lp + single_factor_log_prior_grad(f[h], mixture_fits[[h]])$log_density
    }
    -(ll + mixture_prior_weight * lp)
  }

  gradient <- function(f) {
    eta <- a + as.numeric(L %*% f)
    phi <- dnorm(eta)
    p1 <- pmax(pnorm(eta), 1e-12)
    p0 <- pmax(pnorm(-eta), 1e-12)
    score_eta <- y * phi / p1 - (1 - y) * phi / p0
    grad <- as.numeric(crossprod(L, score_eta))
    for (h in seq_len(H)) {
      grad[h] <- grad[h] + mixture_prior_weight *
        single_factor_log_prior_grad(f[h], mixture_fits[[h]])$grad
    }
    -grad
  }

  opt <- tryCatch(
    optim(
      par = as.numeric(f_init),
      fn = objective,
      gr = gradient,
      method = "L-BFGS-B",
      control = list(maxit = maxit, factr = 1e7)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || !all(is.finite(opt$par))) f_init else opt$par
}

update_factor_scores_mixture_missing <- function(X, W, F_hat, Lambda, alpha,
                                                 mixture_fits,
                                                 mixture_prior_weight,
                                                 maxit_per_subject,
                                                 workers) {
  H <- ncol(F_hat)
  rows <- parallel_lapply(seq_len(nrow(X)), function(i) {
    update_one_factor_score_mixture_missing(
      x_i = X[i, ],
      obs_i = W[i, ],
      f_init = F_hat[i, ],
      Lambda = Lambda,
      alpha = alpha,
      mixture_fits = mixture_fits,
      mixture_prior_weight = mixture_prior_weight,
      maxit = maxit_per_subject
    )
  }, parallel = workers > 1L, workers = workers)
  out <- matrix(unlist(rows), ncol = H, byrow = TRUE)
  rownames(out) <- rownames(F_hat)
  colnames(out) <- colnames(F_hat)
  out
}

fit_mixture_missing_probit <- function(X, W, H, G, lambda_l1_penalty, fold,
                                       workers, seed,
                                       n_aug_iter, n_refine_iter,
                                       z_update,
                                       loading_penalty,
                                       n_random_starts, max_outer,
                                       n_mix_starts, mixture_max_iter,
                                       mixture_prior_weight,
                                       maxit_per_subject,
                                       min_mixture_var) {
  G_fixed <- if (length(G) == 1L) rep(as.integer(G), H) else as.integer(G)
  if (length(G_fixed) != H) {
    stop("G must be either scalar or a length-H component-count vector.")
  }
  G_config <- G_config_label(G_fixed)
  set.seed(seed + 100000L * fold + 1000L * H +
             10L * G_config_seed(G_fixed) + round(lambda_l1_penalty))
  z_update <- match.arg(z_update, c("sample", "expectation"))
  Z <- update_Z_missing(
    X,
    W,
    z_update = z_update,
    seed = seed + 100000L * fold + 1000L * H +
      10L * G_config_seed(G_fixed) + round(lambda_l1_penalty)
  )
  history <- data.frame()
  current <- NULL

  for (iter in seq_len(n_aug_iter)) {
    current <- initialize_mixture_from_Z(
      Z = Z,
      W = W,
      H = H,
      G = G,
      loading_penalty = loading_penalty,
      workers = workers,
      n_random_starts = n_random_starts,
      max_outer = max_outer,
      n_mix_starts = n_mix_starts,
      mixture_max_iter = mixture_max_iter,
      mixture_update = "map",
      seed = seed + iter
    )
    train_ll <- binary_loglik_masked_alpha(X, W, current$F_hat, current$Lambda, current$alpha)
    history <- rbind(history, data.frame(
      stage = "pretrain",
      iteration = iter,
      train_loglik_per_response = train_ll / sum(W),
      mixture_loglik = mixture_prior_loglik(current$F_hat, current$mixture_fits)
    ))
    Z <- update_Z_missing(
      X,
      W,
      current$F_hat,
      current$Lambda,
      current$alpha,
      z_update = z_update,
      seed = seed + 100000L * fold + 1000L * H +
        10L * G_config_seed(G_fixed) + round(lambda_l1_penalty) + iter
    )
  }

  F_hat <- current$F_hat
  Lambda <- current$Lambda
  alpha <- current$alpha
  mixture_fits <- current$mixture_fits

  load <- update_loadings_probit_missing(
    X, W, F_hat, Lambda, alpha, lambda_l1_penalty, workers
  )
  Lambda <- load$Lambda
  alpha <- load$alpha

  prev_score <- -Inf
  for (iter in seq_len(n_refine_iter)) {
    F_hat <- update_factor_scores_mixture_missing(
      X = X,
      W = W,
      F_hat = F_hat,
      Lambda = Lambda,
      alpha = alpha,
      mixture_fits = mixture_fits,
      mixture_prior_weight = mixture_prior_weight,
      maxit_per_subject = maxit_per_subject,
      workers = workers
    )
    scaled <- normalize_factor_scale_alpha(
      F_hat = F_hat,
      Lambda = Lambda,
      mixture_fits = mixture_fits,
      target_scale = 1
    )
    F_hat <- scaled$F_hat
    Lambda <- scaled$Lambda
    mixture_fits <- scaled$mixture_fits

    load <- update_loadings_probit_missing(
      X, W, F_hat, Lambda, alpha, lambda_l1_penalty, workers
    )
    Lambda <- load$Lambda
    alpha <- load$alpha

    mixture_fits <- update_mixture_fits_refinement(
      F_hat = F_hat,
      mixture_fits = mixture_fits,
      n_starts = n_mix_starts,
      max_iter = mixture_max_iter,
      min_var = min_mixture_var,
      mixture_update = "map",
      mu_prior_kappa = 0.05,
      var_prior_shape = 4,
      var_prior_scale = 0.35,
      weight_prior_alpha = 1.2,
      parallel = workers > 1L,
      workers = workers
    )

    train_ll <- binary_loglik_masked_alpha(X, W, F_hat, Lambda, alpha)
    history <- rbind(history, data.frame(
      stage = "refine",
      iteration = iter,
      train_loglik_per_response = train_ll / sum(W),
      mixture_loglik = mixture_prior_loglik(F_hat, mixture_fits)
    ))
    if (iter >= 3L && is.finite(prev_score) &&
        abs(train_ll / sum(W) - prev_score) < 1e-4) {
      break
    }
    prev_score <- train_ll / sum(W)
  }

  ord <- ordered_component_labels(F_hat, mixture_fits)
  list(
    method = "independent_mixture_probit",
    H = H,
    G = max(G_fixed),
    G_config = G_config,
    F_hat = F_hat,
    Lambda = Lambda,
    alpha = alpha,
    mixture_fits = mixture_fits,
    class_map = ord$class_map,
    responsibilities = ord$responsibilities,
    profile_id = profile_id_from_class_map(ord$class_map),
    history = history
  )
}

append_csv <- function(row, path) {
  write.table(
    row,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(path),
    append = file.exists(path)
  )
}

make_summary <- function(scores) {
  if (!"G_config" %in% names(scores)) {
    scores$G_config <- ifelse(is.na(scores$G), NA_character_, as.character(scores$G))
  }
  split_key <- paste(
    scores$G_config,
    scores$H,
    scores$lambda_l1_penalty,
    sep = "|"
  )
  rows <- lapply(split(scores, split_key), function(d) {
    data.frame(
      G = d$G[1L],
      G_config = d$G_config[1L],
      H = d$H[1L],
      lambda_l1_penalty = d$lambda_l1_penalty[1L],
      n_completed_folds = nrow(d),
      mean_heldout_loglik_per_response = mean(d$heldout_loglik_per_response),
      sd_heldout_loglik_per_response = if (nrow(d) > 1L) sd(d$heldout_loglik_per_response) else NA_real_,
      se_heldout_loglik_per_response = if (nrow(d) > 1L) sd(d$heldout_loglik_per_response) / sqrt(nrow(d)) else NA_real_,
      mean_heldout_accuracy = mean(d$heldout_accuracy),
      mean_heldout_brier = mean(d$heldout_brier),
      mean_train_loglik_per_response = mean(d$train_loglik_per_response),
      mean_fit_seconds = mean(d$fit_seconds)
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$G, out$G_config, out$H, out$lambda_l1_penalty), ]
  rownames(out) <- NULL
  out
}

write_cv_plots <- function(summary_scores, out_dir) {
  if (!nrow(summary_scores)) return(invisible(NULL))

  panels <- split(summary_scores, summary_scores$G_config)
  np <- length(panels)
  if (np > 20L) {
    top <- summary_scores[order(summary_scores$mean_heldout_loglik_per_response, decreasing = TRUE), ]
    top <- head(top, 30L)
    write.csv(top, file.path(out_dir, "rank_lambda_top_candidates_by_heldout_ll.csv"),
              row.names = FALSE)

    png(file.path(out_dir, "rank_lambda_top_candidates_by_heldout_ll.png"),
        width = 2400, height = 1700, res = 170)
    op <- par(mar = c(5, 10, 3, 1))
    on.exit({
      par(op)
      dev.off()
    }, add = TRUE)
    labels <- sprintf(
      "H=%s; G=[%s]; lambda=%s; folds=%s",
      top$H,
      ifelse(is.na(top$G_config), "NA", top$G_config),
      top$lambda_l1_penalty,
      top$n_completed_folds
    )
    dotchart(
      rev(top$mean_heldout_loglik_per_response),
      labels = rev(labels),
      pch = 19,
      xlab = "held-out log likelihood / response",
      main = "Top IFEval rank / column-G / penalty candidates"
    )
    abline(v = max(top$mean_heldout_loglik_per_response, na.rm = TRUE),
           col = "#B23A48", lty = 2)
    return(invisible(NULL))
  }

  png(file.path(out_dir, "rank_lambda_heldout_loglik_lines.png"),
      width = 2200, height = 1200, res = 170)
  nr <- ceiling(np / 2)
  op <- par(mfrow = c(nr, 2), mar = c(4, 4.2, 3, 1), oma = c(0, 0, 2, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  cols <- c("#355C9A", "#B23A48", "#2A9D8F", "#6D597A", "#E9A03F", "#4A5568")
  for (nm in names(panels)) {
    d <- panels[[nm]]
    lambda_vals <- sort(unique(d$lambda_l1_penalty))
    ylim <- range(d$mean_heldout_loglik_per_response, na.rm = TRUE)
    plot(NA, xlim = range(d$H), ylim = ylim, xlab = "rank H",
         ylab = "held-out log likelihood / response", main = nm)
    for (k in seq_along(lambda_vals)) {
      dd <- d[d$lambda_l1_penalty == lambda_vals[k], ]
      dd <- dd[order(dd$H), ]
      lines(dd$H, dd$mean_heldout_loglik_per_response, type = "b",
            pch = 19, col = cols[(k - 1L) %% length(cols) + 1L])
    }
    legend("bottomright", legend = paste0("lambda=", lambda_vals),
           col = cols[seq_along(lambda_vals)], lty = 1, pch = 19, cex = 0.7,
           bty = "n")
  }
  title("IFEval rank/penalty CV; larger held-out log likelihood is better", outer = TRUE)
}

refresh_outputs <- function(scores_path, out_dir) {
  if (!file.exists(scores_path)) return(invisible(NULL))
  scores <- read.csv(scores_path, stringsAsFactors = FALSE)
  if (!nrow(scores)) return(invisible(NULL))
  summary_scores <- make_summary(scores)
  write.csv(summary_scores, file.path(out_dir, "ifeval_rank_lambda_cv_summary_partial.csv"),
            row.names = FALSE)

  complete <- summary_scores[summary_scores$n_completed_folds == max(summary_scores$n_completed_folds), ]
  selected <- complete[which.max(complete$mean_heldout_loglik_per_response), , drop = FALSE]
  rownames(selected) <- NULL
  write.csv(selected, file.path(out_dir, "ifeval_rank_lambda_selected_by_heldout_ll_partial.csv"),
            row.names = FALSE)
  if (tolower(Sys.getenv("REFRESH_PLOTS", "TRUE")) %in% c("false", "0", "no")) {
    return(invisible(NULL))
  }
  tryCatch(
    write_cv_plots(summary_scores, out_dir),
    error = function(e) message("Skipping plot refresh: ", conditionMessage(e))
  )
}

run_rscript_with_env <- function(script_path, env_values) {
  env_names <- names(env_values)
  old_values <- Sys.getenv(env_names, unset = NA_character_)
  on.exit({
    for (k in seq_along(env_names)) {
      if (is.na(old_values[k])) {
        Sys.unsetenv(env_names[k])
      } else {
        do.call(Sys.setenv, stats::setNames(as.list(old_values[k]), env_names[k]))
      }
    }
  }, add = TRUE)
  do.call(Sys.setenv, as.list(env_values))
  status <- system2("Rscript", args = shQuote(script_path))
  if (!identical(status, 0L)) {
    warning("Selected-model refit returned nonzero status: ", status)
  }
  invisible(status)
}

matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  file.path(repo_root, "data", "ifeval", "openeval_ifeval_only_binary_matrix.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(repo_root, "results", "ifeval_rank_lambda_cv_map")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

X_full <- read_binary_matrix(matrix_path)
max_feasible_H <- min(nrow(X_full) - 1L, ncol(X_full))
H_grid <- parse_int_grid(Sys.getenv("H_GRID"), default = 1:min(6L, max_feasible_H))
H_grid <- H_grid[H_grid >= 1L & H_grid <= max_feasible_H]
G_grid <- parse_int_grid(Sys.getenv("G_GRID"), default = c(2L, 3L))
G_mode <- Sys.getenv("G_MODE", "fixed")
G_mode <- match.arg(G_mode, c("fixed", "column_grid"))
G_component_values <- parse_int_grid(
  Sys.getenv("G_COMPONENT_VALUES"),
  default = c(2L, 3L)
)
max_gaussian_coords <- as.integer(Sys.getenv("MAX_GAUSSIAN_COORDS", "0"))
unique_G_up_to_permutation <- isTRUE(
  tolower(Sys.getenv("UNIQUE_G_CONFIGS_UP_TO_PERMUTATION", "FALSE")) %in%
    c("true", "1", "yes")
)
lambda_grid <- parse_num_grid(Sys.getenv("LAMBDA_L1_GRID"), default = c(0, 1, 2, 4, 8, 12))
K_folds <- as.integer(Sys.getenv("K_FOLDS", "3"))
workers <- resolve_workers(as.integer(Sys.getenv("WORKERS", "6")))
seed <- as.integer(Sys.getenv("SEED", "20260812"))
resume_existing <- isTRUE(tolower(Sys.getenv("RESUME_EXISTING", "TRUE")) %in% c("true", "1", "yes"))

n_aug_iter <- as.integer(Sys.getenv("PRETRAIN_AUG_ITER", "5"))
n_refine_iter <- as.integer(Sys.getenv("REFINE_ITER", "6"))
pretrain_z_update <- Sys.getenv("PRETRAIN_Z_UPDATE", "sample")
if (!pretrain_z_update %in% c("sample", "expectation")) {
  stop("PRETRAIN_Z_UPDATE must be either 'sample' or 'expectation'.")
}
loading_penalty <- as.numeric(Sys.getenv("PRETRAIN_LOADING_PENALTY", "0.05"))
n_random_starts <- as.integer(Sys.getenv("N_RANDOM_STARTS", "1"))
max_outer <- as.integer(Sys.getenv("MAX_OUTER", "3"))
n_mix_starts <- as.integer(Sys.getenv("N_MIX_STARTS", "2"))
mixture_max_iter <- as.integer(Sys.getenv("MIXTURE_MAX_ITER", "15"))
mixture_prior_weight <- as.numeric(Sys.getenv("MIXTURE_PRIOR_WEIGHT", "0.35"))
maxit_per_subject <- as.integer(Sys.getenv("MAXIT_PER_SUBJECT", "45"))
min_mixture_var <- as.numeric(Sys.getenv("MIN_MIXTURE_VAR", "0.05"))
save_fits <- isTRUE(tolower(Sys.getenv("SAVE_FITS", "FALSE")) %in% c("true", "1", "yes"))
fit_selected_after_cv <- isTRUE(tolower(Sys.getenv("FIT_SELECTED_AFTER_CV", "TRUE")) %in% c("true", "1", "yes"))
item_metadata_path <- Sys.getenv(
  "ITEM_METADATA_PATH",
  file.path(repo_root, "data", "ifeval", "openeval_item_metadata.csv")
)

set.seed(seed)
fold_id <- matrix(
  sample.int(K_folds, length(X_full), replace = TRUE),
  nrow = nrow(X_full),
  ncol = ncol(X_full),
  dimnames = dimnames(X_full)
)

scores_path <- file.path(out_dir, "ifeval_rank_lambda_cv_fold_scores.csv")
history_path <- file.path(out_dir, "ifeval_rank_lambda_cv_histories.csv")
existing <- if (resume_existing && file.exists(scores_path)) {
  read.csv(scores_path, stringsAsFactors = FALSE)
} else {
  data.frame()
}

already_done <- function(method, G, G_config, H, lambda, fold) {
  if (!resume_existing || !nrow(existing)) return(FALSE)
  if (!"G_config" %in% names(existing)) {
    existing$G_config <- ifelse(is.na(existing$G), NA_character_, as.character(existing$G))
  }
  same_method <- existing$method == method
  same_G <- if (is.na(G)) is.na(existing$G) else existing$G == G
  same_G_config <- if (is.na(G_config)) {
    is.na(existing$G_config)
  } else {
    existing$G_config == G_config
  }
  same_H <- existing$H == H
  same_lambda <- abs(existing$lambda_l1_penalty - lambda) < 1e-12
  same_fold <- existing$fold == fold
  any(same_method & same_G & same_G_config & same_H & same_lambda & same_fold)
}

message("IFEval matrix after filtering: n=", nrow(X_full), ", p=", ncol(X_full))
message("Method: independent_mixture_probit")
message("H grid: ", paste(H_grid, collapse = ", "))
message("G grid: ", paste(G_grid, collapse = ", "))
message("G mode: ", G_mode)
if (G_mode == "column_grid") {
  message("Columnwise component values: ", paste(G_component_values, collapse = ", "))
  message("Max Gaussian coordinates: ", max_gaussian_coords)
  message("Unique G configs up to permutation: ", unique_G_up_to_permutation)
}
message("lambda_l1 grid: ", paste(lambda_grid, collapse = ", "))
message("Folds: ", K_folds, "; MAP refinement; workers=", workers)
message("Pretraining Z update: ", pretrain_z_update)
message("Output directory: ", normalizePath(out_dir, mustWork = FALSE))

method_name <- "independent_mixture_probit"
for (G in G_grid) {
  for (H in H_grid) {
    G_configs <- if (G_mode == "column_grid") {
      expand_columnwise_G_configs(
        H,
        component_values = G_component_values,
        max_gaussian_coords = max_gaussian_coords,
        unique_up_to_permutation = unique_G_up_to_permutation
      )
    } else {
      list(rep(as.integer(G), H))
    }
    for (G_candidate in G_configs) {
      if (G_mode == "column_grid" && max(G_candidate) != G) {
        next
      }
      G_config <- G_config_label(G_candidate)
      for (lambda_l1_penalty in lambda_grid) {
        for (fold in seq_len(K_folds)) {
          if (already_done(method_name, G, G_config, H, lambda_l1_penalty, fold)) {
            message("Skipping completed: method=", method_name, ", G=", G,
                    ", G_config=", G_config, ", H=", H,
                    ", lambda=", lambda_l1_penalty, ", fold=", fold)
            next
          }

          message("CV fit: method=", method_name, ", G=", G,
                  ", G_config=", G_config,
                  ", H=", H, ", lambda=", lambda_l1_penalty,
                  ", fold=", fold, "/", K_folds)
          heldout <- fold_id == fold
          W <- !heldout
          start_time <- proc.time()[["elapsed"]]

          fit <- fit_mixture_missing_probit(
            X = X_full,
            W = W,
            H = H,
            G = G_candidate,
            lambda_l1_penalty = lambda_l1_penalty,
            fold = fold,
            workers = workers,
            seed = seed,
            n_aug_iter = n_aug_iter,
            n_refine_iter = n_refine_iter,
            z_update = pretrain_z_update,
            loading_penalty = loading_penalty,
            n_random_starts = n_random_starts,
            max_outer = max_outer,
            n_mix_starts = n_mix_starts,
            mixture_max_iter = mixture_max_iter,
            mixture_prior_weight = mixture_prior_weight,
            maxit_per_subject = maxit_per_subject,
            min_mixture_var = min_mixture_var
          )

          elapsed <- proc.time()[["elapsed"]] - start_time
          sc <- score_cells(X_full, heldout, fit)
          train_total_ll <- binary_loglik_masked_alpha(
            X_full, W, fit$F_hat, fit$Lambda, fit$alpha
          )
          sc$train_total_loglik <- train_total_ll
          sc$train_loglik_per_response <- train_total_ll / sum(W)
          sc$n_train <- sum(W)
          sc$fit_seconds <- elapsed

          row <- cbind(
            method = method_name,
            G = G,
            G_config = G_config,
            H = H,
            lambda_l1_penalty = lambda_l1_penalty,
            fold = fold,
            sc
          )
          append_csv(row, scores_path)

          hist <- cbind(
            method = method_name,
            G = G,
            G_config = G_config,
            H = H,
            lambda_l1_penalty = lambda_l1_penalty,
            fold = fold,
            pretrain_z_update = pretrain_z_update,
            fit$history
          )
          append_csv(hist, history_path)

          if (save_fits) {
            fit_file <- sprintf("%s_G%s_config%s_H%02d_lambda%s_fold%02d.rds",
                                method_name,
                                G,
                                gsub(",", "-", G_config),
                                H,
                                gsub("\\.", "p", as.character(lambda_l1_penalty)),
                                fold)
            saveRDS(fit, file.path(out_dir, fit_file))
          }

          refresh_outputs(scores_path, out_dir)
          message("  heldout ll/resp=", signif(sc$heldout_loglik_per_response, 4),
                  ", train ll/resp=", signif(sc$train_loglik_per_response, 4),
                  ", seconds=", signif(elapsed, 4))
        }
      }
      }
    }
}

scores <- read.csv(scores_path, stringsAsFactors = FALSE)
summary_scores <- make_summary(scores)
write.csv(summary_scores, file.path(out_dir, "ifeval_rank_lambda_cv_summary.csv"),
          row.names = FALSE)
selected <- summary_scores[which.max(summary_scores$mean_heldout_loglik_per_response), , drop = FALSE]
rownames(selected) <- NULL
write.csv(selected, file.path(out_dir, "ifeval_rank_lambda_selected_by_heldout_ll.csv"),
          row.names = FALSE)
if (!(tolower(Sys.getenv("REFRESH_PLOTS", "TRUE")) %in% c("false", "0", "no"))) {
  tryCatch(
    write_cv_plots(summary_scores, out_dir),
    error = function(e) message("Skipping final CV plots: ", conditionMessage(e))
  )
}

if (fit_selected_after_cv) {
  mixture_G_config <- if (!is.na(selected$G_config) && nzchar(selected$G_config)) {
    selected$G_config
  } else {
    paste(rep(as.integer(selected$G), as.integer(selected$H)), collapse = ",")
  }
  mixture_G_label <- gsub(",", "-", mixture_G_config)
  mixture_out <- file.path(
    out_dir,
    sprintf(
      "selected_mixture_H%d_Gconfig%s_lambda%s",
      selected$H,
      mixture_G_label,
      gsub("\\.", "p", as.character(selected$lambda_l1_penalty))
    )
  )
  dir.create(mixture_out, recursive = TRUE, showWarnings = FALSE)
  message("Fitting selected full-data mixture model in: ", mixture_out)
  run_rscript_with_env(
    file.path(ifeval_script_dir, "fit_interpret_ifeval_mixture.R"),
    c(
      MATRIX_PATH = matrix_path,
      ITEM_METADATA_PATH = item_metadata_path,
      OUT_DIR = mixture_out,
      H_FIXED = as.character(selected$H),
      G_FIXED = mixture_G_config,
      WORKERS = as.character(workers),
      REFINEMENT_LAMBDA_L1_PENALTY = as.character(selected$lambda_l1_penalty),
      PRETRAIN_LOADING_PENALTY = as.character(loading_penalty),
      PRETRAIN_Z_UPDATE = pretrain_z_update,
      PRETRAIN_AUG_ITER = as.character(n_aug_iter),
      REFINE_ITER = as.character(n_refine_iter),
      MIXTURE_MAX_ITER = as.character(mixture_max_iter),
      REQUIRE_MIXTURE_CONVERGENCE = "TRUE"
    )
  )
}

cat("\nSelected models by held-out predictive log likelihood:\n")
print(selected)
cat("\nOutputs saved in: ", normalizePath(out_dir), "\n", sep = "")
