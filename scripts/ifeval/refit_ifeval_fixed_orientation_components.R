#!/usr/bin/env Rscript

# Refit a selected IFEval P-IFA model with component counts attached to the
# selected model's displayed factor coordinates. Unlike a fresh rotation fit,
# this diagnostic does not allow the requested G_h slots to migrate across
# factors during rotation.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
source(file.path(script_dir, "ifeval_imfm_helpers.R"))

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

get_fit_field <- function(fit, primary, fallback) {
  value <- fit[[primary]]
  if (is.null(value)) value <- fit[[fallback]]
  if (is.null(value)) stop("Base fit has neither `", primary, "` nor `", fallback, "`.")
  value
}

plot_factor_marginals <- function(fit, path) {
  H <- ncol(fit$F_hat)
  width <- max(8, 4 * H)
  png(path, width = width, height = 4.2, units = "in", res = 220)
  old <- par(mfrow = c(1, H), mar = c(4.2, 4.2, 3, 1), oma = c(0, 0, 2.5, 0))
  on.exit({
    par(old)
    dev.off()
  }, add = TRUE)

  colors <- c("#2f6fb8", "#d43f3a", "#2e8b57", "#8e63b6")
  for (h in seq_len(H)) {
    scores <- fit$F_hat[, h]
    mixture <- fit$mixture_fits[[h]]
    limits <- range(c(scores, mixture$mu + 4 * sqrt(mixture$var),
                      mixture$mu - 4 * sqrt(mixture$var)), finite = TRUE)
    hist(scores, probability = TRUE, breaks = "FD", col = "#dce7f2",
         border = "white", xlim = limits, main = paste0("Factor ", h),
         xlab = "Factor score", ylab = "Density")
    grid <- seq(limits[1L], limits[2L], length.out = 500L)
    component_density <- vapply(seq_along(mixture$pi), function(g) {
      mixture$pi[g] * dnorm(grid, mixture$mu[g], sqrt(mixture$var[g]))
    }, numeric(length(grid)))
    if (is.null(dim(component_density))) {
      component_density <- matrix(component_density, ncol = 1L)
    }
    lines(grid, rowSums(component_density), lwd = 2.5, col = "black")
    for (g in seq_len(ncol(component_density))) {
      lines(grid, component_density[, g], lwd = 2,
            col = colors[1L + ((g - 1L) %% length(colors))])
      abline(v = mixture$mu[g], lty = 3,
             col = colors[1L + ((g - 1L) %% length(colors))])
    }
  }
  mtext(sprintf("Fixed-orientation P-IFA: G = (%s)", paste(fit$G_hat, collapse = ",")),
        outer = TRUE, cex = 1.2, font = 2)
}

base_fit_path <- Sys.getenv("BASE_FIT")
if (!nzchar(base_fit_path)) stop("Set BASE_FIT to the selected full-data fit RDS.")
matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  file.path(repo_root, "data", "ifeval_threshold_1", "openeval_ifeval_only_binary_matrix.csv")
)
item_metadata_path <- Sys.getenv(
  "ITEM_METADATA_PATH",
  file.path(repo_root, "data", "ifeval_threshold_1", "openeval_item_metadata.csv")
)
out_dir <- Sys.getenv("OUT_DIR", file.path(repo_root, "results", "diagnostics", "ifeval_fixed_orientation"))
G_fixed <- as.integer(strsplit(gsub("[[:space:]]+", "", Sys.getenv("G_FIXED", "2,2,3,1")),
                              ",", fixed = TRUE)[[1L]])
lambda <- as.numeric(Sys.getenv("LAMBDA_L1_PENALTY", "4"))
workers <- as.integer(Sys.getenv("WORKERS", "18"))
seed <- as.integer(Sys.getenv("SEED", "20260929"))

base_fit <- readRDS(base_fit_path)
X <- read_binary_matrix_complete(matrix_path)
F_start <- as.matrix(base_fit$F_hat)
Lambda_start <- as.matrix(get_fit_field(base_fit, "Lambda", "Lambda_hat"))
alpha_start <- as.numeric(get_fit_field(base_fit, "alpha", "alpha_hat"))
H <- ncol(F_start)
if (length(G_fixed) != H) stop("G_FIXED must have one value per factor in BASE_FIT.")
if (!identical(dim(X), c(nrow(F_start), nrow(Lambda_start)))) {
  stop("The data matrix dimensions do not match the base fit.")
}

set.seed(seed)
mixture_fits <- lapply(seq_len(H), function(h) {
  previous <- base_fit$mixture_fits[[h]]
  fit_gmm_1d(
    x = F_start[, h],
    G = G_fixed[h],
    n_starts = 2L,
    max_iter = 200L,
    min_var = 0.05,
    mixture_update = "map",
    mu_prior_mean = 0,
    mu_prior_kappa = 0.05,
    var_prior_shape = 3,
    var_prior_scale = 2,
    weight_prior_alpha = 1,
    init = if (length(previous$pi) == G_fixed[h]) previous else NULL
  )
})

pretrain_fit <- list(
  model = "fixed_orientation_component_refit",
  F_hat = F_start,
  Lambda_hat = Lambda_start,
  alpha_hat = alpha_start,
  mixture_fits = mixture_fits,
  G_hat = G_fixed,
  base_fit_path = normalizePath(base_fit_path, mustWork = TRUE)
)

start_time <- proc.time()[["elapsed"]]
fit <- fit_binary_probit_refinement(
  X = X,
  pretrain_fit = pretrain_fit,
  n_refine_iter = 100L,
  maxit_per_subject = 45L,
  n_mix_starts = 2L,
  mixture_max_iter = 200L,
  factor_update = "marginal",
  min_mixture_var = 0.05,
  mixture_update = "map",
  mixture_refit = "em",
  mu_prior_mean = 0,
  mu_prior_kappa = 0.05,
  var_prior_shape = 3,
  var_prior_scale = 2,
  weight_prior_alpha = 1,
  mixture_prior_weight = 1,
  factor_score_bound = 5,
  estimate_intercept = TRUE,
  preestimate_loadings = TRUE,
  lambda_l1_penalty = lambda,
  lasso_maxit = 160L,
  lasso_tol = 1e-5,
  lasso_backend = "proximal",
  normalize_factor_scale = TRUE,
  normalize_factor_location = TRUE,
  factor_scale_target = 1,
  factor_scale_method = "sd",
  objective_tolerance = 1e-3,
  min_refine_iter = 2L,
  stopping_objective = "posterior_objective",
  enforce_monotone_refinement = TRUE,
  monotone_tolerance = 1e-8,
  require_mixture_convergence_for_stop = TRUE,
  parallel = env_flag("PARALLEL", TRUE),
  workers = workers,
  verbose = TRUE
)
fit$elapsed_seconds <- proc.time()[["elapsed"]] - start_time

# Normalize location, scale, and signs for reporting while retaining F1--F4.
fit <- canonical_normalize_refined_fit(
  fit,
  sign_rule = "largest_loading_positive",
  order_rule = "none"
)
fit$X <- X
fit$requested_G_fixed <- G_fixed
fit$lambda_l1_penalty <- lambda
fit$matrix_path <- normalizePath(matrix_path, mustWork = TRUE)
fit$factor_accuracy_correlation <- vapply(seq_len(H), function(h) {
  cor(fit$F_hat[, h], rowMeans(X))
}, numeric(1L))

ord <- ordered_component_labels(fit$F_hat, fit$mixture_fits)
fit$class_map <- ord$class_map
fit$responsibilities <- ord$responsibilities
fit$profile_id <- profile_id_from_class_map(ord$class_map)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(fit, file.path(out_dir, "fixed_orientation_fit.rds"))

mixture_summary <- do.call(rbind, lapply(seq_len(H), function(h) {
  mixture <- fit$mixture_fits[[h]]
  data.frame(
    factor = h,
    component = seq_along(mixture$pi),
    weight = mixture$pi,
    mean = mixture$mu,
    sd = sqrt(mixture$var),
    converged = isTRUE(mixture$converged)
  )
}))
write.csv(mixture_summary, file.path(out_dir, "fixed_orientation_mixture_parameters.csv"),
          row.names = FALSE)
write.csv(
  data.frame(
    factor = seq_len(H),
    G = fit$G_hat,
    accuracy_correlation = fit$factor_accuracy_correlation,
    loading_l2 = sqrt(colSums(fit$Lambda_hat^2)),
    loading_nonzero = colSums(abs(fit$Lambda_hat) > 1e-8)
  ),
  file.path(out_dir, "fixed_orientation_factor_summary.csv"),
  row.names = FALSE
)

if (file.exists(item_metadata_path)) {
  metadata <- read.csv(
    item_metadata_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  item_index <- match(rownames(fit$Lambda_hat), metadata$item_id)
  if (anyNA(item_index)) stop("Item metadata does not cover every fitted loading row.")
  item_loadings <- metadata[item_index, c(
    "benchmark", "item_id", "prompt", "instruction_ids",
    "instruction_families", "n_instructions"
  )]
  for (h in seq_len(H)) item_loadings[[paste0("F", h)]] <- fit$Lambda_hat[, h]
  write.csv(
    item_loadings,
    file.path(out_dir, "fixed_orientation_item_loadings_metadata.csv"),
    row.names = FALSE
  )

  summarize_labels <- function(field) {
    expanded <- do.call(rbind, lapply(seq_len(nrow(item_loadings)), function(i) {
      labels <- strsplit(item_loadings[[field]][i], "|", fixed = TRUE)[[1L]]
      values <- as.data.frame(
        as.list(item_loadings[i, paste0("F", seq_len(H)), drop = FALSE]),
        stringsAsFactors = FALSE
      )
      values <- values[rep(1L, length(labels)), , drop = FALSE]
      data.frame(label = labels, values, check.names = FALSE)
    }))
    means <- aggregate(
      expanded[, paste0("F", seq_len(H)), drop = FALSE],
      list(label = expanded$label),
      mean
    )
    counts <- as.data.frame(table(expanded$label), stringsAsFactors = FALSE)
    names(counts) <- c("label", "n")
    merge(means, counts, by = "label", sort = FALSE)
  }

  write.csv(
    summarize_labels("instruction_families"),
    file.path(out_dir, "fixed_orientation_family_loading_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    summarize_labels("instruction_ids"),
    file.path(out_dir, "fixed_orientation_instruction_loading_summary.csv"),
    row.names = FALSE
  )
}
plot_factor_marginals(fit, file.path(out_dir, "fixed_orientation_factor_marginals.png"))

cat("Requested G:", paste(G_fixed, collapse = ","), "\n")
cat("Returned G:", paste(fit$G_hat, collapse = ","), "\n")
cat("Refinement converged:", fit$joint_refinement$converged, "\n")
cat("Completed iterations:", fit$joint_refinement$n_completed, "\n")
cat("Elapsed seconds:", signif(fit$elapsed_seconds, 5), "\n")
print(mixture_summary)
cat("Outputs saved in:", normalizePath(out_dir), "\n")
