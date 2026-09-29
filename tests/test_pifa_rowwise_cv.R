#!/usr/bin/env Rscript

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
test_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
root <- normalizePath(file.path(test_dir, ".."), mustWork = FALSE)
source(file.path(root, "R", "pifa_rowwise_cv.R"))

set.seed(11)
X <- matrix(rbinom(12L * 20L, 1L, 0.5), 12L, 20L)
rownames(X) <- paste0("model_", seq_len(nrow(X)))
X[1, c(3, 8)] <- NA
X[7, 15] <- NA

folds <- make_pifa_row_folds(nrow(X), K = 5L, seed = 12L)
stopifnot(
  setequal(unique(folds), 1:5),
  length(folds) == nrow(X),
  all(tabulate(folds, nbins = 5L) > 0L)
)

plan_a <- make_pifa_reveal_plan(X, folds, 0.20, 3L, seed = 13L)
plan_b <- make_pifa_reveal_plan(X, folds, 0.20, 3L, seed = 13L)
stopifnot(identical(plan_a, plan_b))
plan_index <- pifa_reveal_plan_index(plan_a)
stopifnot(
  all(nzchar(plan_index$reveal_index)),
  all(nzchar(plan_index$prediction_index))
)
for (unit in plan_a) {
  observed <- which(!is.na(X[unit$row_index, ]))
  stopifnot(
    !length(intersect(unit$reveal_index, unit$prediction_index)),
    setequal(c(unit$reveal_index, unit$prediction_index), observed),
    length(unit$reveal_index) == max(1L, round(0.20 * length(observed))),
    !anyNA(X[unit$row_index, unit$reveal_index]),
    !anyNA(X[unit$row_index, unit$prediction_index])
  )
}

for (H in c(2L, 8L)) {
  configs <- valid_pifa_G_configs(H, 1:3)
  labels <- vapply(configs, G_config_label, character(1L))
  stopifnot(
    !anyDuplicated(labels),
    all(vapply(configs, function(G) all(diff(G) >= 0L), logical(1L))),
    all(vapply(configs, function(G) sum(G == 1L) <= 1L, logical(1L)))
  )
}
grid <- make_pifa_candidate_grid(c(2L, 8L), 1:3, c(0, 4))
stopifnot(all(c(2L, 8L) %in% grid$H), 0 %in% grid$lambda_l1_penalty)

# Changing hidden prediction responses must not alter the held-out MAP factor.
fake_fit <- list(
  alpha = rep(0, ncol(X)),
  Lambda = matrix(rnorm(ncol(X) * 2L, sd = 0.2), ncol(X), 2L),
  mixture_fits = list(
    list(pi = c(0.5, 0.5), mu = c(-1, 1), var = c(0.5, 0.5), converged = TRUE),
    list(pi = c(0.5, 0.5), mu = c(-1, 1), var = c(0.5, 0.5), converged = TRUE)
  )
)
unit <- plan_a[[1L]]
x_original <- X[unit$row_index, ]
x_changed <- x_original
x_changed[unit$prediction_index] <- 1 - x_changed[unit$prediction_index]
eta_a <- estimate_heldout_factor_map(x_original, unit$reveal_index, fake_fit)
eta_b <- estimate_heldout_factor_map(x_changed, unit$reveal_index, fake_fit)
stopifnot(
  eta_a$converged,
  eta_b$converged,
  isTRUE(all.equal(eta_a$eta_hat, eta_b$eta_hat, tolerance = 1e-10))
)

# Rank-specific factor diagnostics must combine without losing coordinates.
combined <- rbind_fill(list(
  data.frame(candidate_id = "H2", eta_1 = 1, eta_2 = 2),
  data.frame(candidate_id = "H4", eta_1 = 3, eta_2 = 4, eta_3 = 5, eta_4 = 6)
))
stopifnot(
  identical(names(combined), c("candidate_id", "eta_1", "eta_2", "eta_3", "eta_4")),
  is.na(combined$eta_3[1L]),
  combined$eta_4[2L] == 6
)

cat("All few-shot row-wise CV tests passed.\n")
