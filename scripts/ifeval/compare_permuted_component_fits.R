#!/usr/bin/env Rscript

# Compare two full-data IFEval fits whose component-count vectors are
# permutations of one another, for example G=(3,3,1,3) versus G=(3,1,3,3).
# The comparison aligns the second fit to the first using the best
# permutation/sign match of the loading columns, then reports loading and
# factor-score correlations.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

base_dir <- Sys.getenv(
  "BASE_DIR",
  file.path(repo_root, "results", "full", "ifeval_threshold_sensitivity_2thirds_20260825", "threshold_0p5")
)
fit_a_dir <- Sys.getenv(
  "FIT_A_DIR",
  file.path(base_dir, "top3_full_data_refits", "rank3-H4-G3-3-1-3-lambda4")
)
fit_b_dir <- Sys.getenv(
  "FIT_B_DIR",
  file.path(base_dir, "permuted_component_refits", "H4_G3-1-3-3_lambda4")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(base_dir, "permuted_component_refits", "comparisons", "G3313_vs_G3133")
)
score_path <- Sys.getenv(
  "SCORE_PATH",
  file.path(base_dir, "componentwise_cv", "ifeval_rank_lambda_cv_fold_scores.csv")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

find_one <- function(dir, pattern) {
  hit <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(hit) != 1L) stop("Expected one file matching ", pattern, " in ", dir)
  hit
}

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  out <- lapply(seq_along(x), function(i) {
    cbind(x[i], permutations(x[-i]))
  })
  do.call(rbind, out)
}

best_permutation_alignment <- function(corr_mat) {
  H <- nrow(corr_mat)
  perms <- permutations(seq_len(H))
  best <- NULL
  for (i in seq_len(nrow(perms))) {
    perm <- perms[i, ]
    diag_corr <- corr_mat[cbind(seq_len(H), perm)]
    score <- sum(abs(diag_corr))
    if (is.null(best) || score > best$score) {
      signs <- sign(diag_corr)
      signs[signs == 0] <- 1
      best <- list(score = score, perm = perm, signs = signs, corrs = diag_corr)
    }
  }
  best
}

read_lambda <- function(dir) {
  x <- read.csv(find_one(dir, "openeval_item_intercepts_loadings_metadata[.]csv$"), check.names = FALSE)
  load_cols <- grep("^loading_factor_", names(x), value = TRUE)
  list(item = x, Lambda = as.matrix(x[, load_cols, drop = FALSE]), load_cols = load_cols)
}

fit_a <- readRDS(find_one(fit_a_dir, "_fit[.]rds$"))
fit_b <- readRDS(find_one(fit_b_dir, "_fit[.]rds$"))
lam_a <- read_lambda(fit_a_dir)
lam_b <- read_lambda(fit_b_dir)

if (ncol(lam_a$Lambda) != ncol(lam_b$Lambda)) {
  stop("The two fits have different H values.")
}
H <- ncol(lam_a$Lambda)

loading_corr_raw <- cor(lam_a$Lambda, lam_b$Lambda)
score_corr_raw <- cor(fit_a$F_hat, fit_b$F_hat)
align <- best_permutation_alignment(loading_corr_raw)

Lambda_b_aligned <- sweep(lam_b$Lambda[, align$perm, drop = FALSE], 2L, align$signs, "*")
F_b_aligned <- sweep(fit_b$F_hat[, align$perm, drop = FALSE], 2L, align$signs, "*")
loading_corr_aligned <- cor(lam_a$Lambda, Lambda_b_aligned)
score_corr_aligned <- cor(fit_a$F_hat, F_b_aligned)

fit_a_summary <- read.csv(find_one(fit_a_dir, "_fit_summary[.]csv$"), check.names = FALSE)
fit_b_summary <- read.csv(find_one(fit_b_dir, "_fit_summary[.]csv$"), check.names = FALSE)

cv_summary <- data.frame()
if (file.exists(score_path)) {
  scores <- read.csv(score_path, check.names = FALSE)
  target <- rbind(
    fit_a_summary[, c("H", "G_config", "lambda_l1_penalty")],
    fit_b_summary[, c("H", "G_config", "lambda_l1_penalty")]
  )
  for (i in seq_len(nrow(target))) {
    rows <- scores[
      scores$H == target$H[i] &
        scores$G_config == target$G_config[i] &
        scores$lambda_l1_penalty == target$lambda_l1_penalty[i],
      ,
      drop = FALSE
    ]
    cv_summary <- rbind(
      cv_summary,
      data.frame(
        fit = c("A", "B")[i],
        H = target$H[i],
        G_config = target$G_config[i],
        lambda_l1_penalty = target$lambda_l1_penalty[i],
        n_folds = nrow(rows),
        mean_heldout_loglik_per_response = mean(rows$heldout_loglik_per_response),
        sd_heldout_loglik_per_response = stats::sd(rows$heldout_loglik_per_response),
        stringsAsFactors = FALSE
      )
    )
  }
}

alignment_df <- data.frame(
  reference_factor = paste0("A_F", seq_len(H)),
  matched_fit_b_factor = paste0("B_F", align$perm),
  sign = align$signs,
  loading_correlation = align$corrs,
  aligned_loading_correlation = diag(loading_corr_aligned),
  aligned_factor_score_correlation = diag(score_corr_aligned),
  reference_G = strsplit(fit_a_summary$G_config[1], ",", fixed = TRUE)[[1L]],
  matched_fit_b_G = strsplit(fit_b_summary$G_config[1], ",", fixed = TRUE)[[1L]][align$perm],
  stringsAsFactors = FALSE
)

write.csv(alignment_df, file.path(out_dir, "permuted_fit_alignment_summary.csv"), row.names = FALSE)
write.csv(loading_corr_raw, file.path(out_dir, "raw_loading_correlation_matrix.csv"), row.names = TRUE)
write.csv(score_corr_raw, file.path(out_dir, "raw_factor_score_correlation_matrix.csv"), row.names = TRUE)
write.csv(loading_corr_aligned, file.path(out_dir, "aligned_loading_correlation_matrix.csv"), row.names = TRUE)
write.csv(score_corr_aligned, file.path(out_dir, "aligned_factor_score_correlation_matrix.csv"), row.names = TRUE)
write.csv(cv_summary, file.path(out_dir, "permuted_fit_cv_summary.csv"), row.names = FALSE)

comparison_summary <- data.frame(
  fit = c("A", "B"),
  directory = c(fit_a_dir, fit_b_dir),
  H = c(fit_a_summary$H[1], fit_b_summary$H[1]),
  G_config = c(fit_a_summary$G_config[1], fit_b_summary$G_config[1]),
  lambda_l1_penalty = c(fit_a_summary$lambda_l1_penalty[1], fit_b_summary$lambda_l1_penalty[1]),
  full_data_loglik_per_response = c(
    fit_a_summary$binary_loglik_per_response[1],
    fit_b_summary$binary_loglik_per_response[1]
  ),
  elapsed_sec = c(fit_a_summary$elapsed_sec[1], fit_b_summary$elapsed_sec[1]),
  stringsAsFactors = FALSE
)
write.csv(comparison_summary, file.path(out_dir, "permuted_fit_comparison_summary.csv"), row.names = FALSE)

png(file.path(out_dir, "raw_loading_correlation_heatmap.png"), width = 1100, height = 900, res = 160)
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)
image(
  x = seq_len(H), y = seq_len(H), z = t(loading_corr_raw[H:1, ]),
  col = colorRampPalette(c("#2F6DAE", "white", "#C43C4A"))(101),
  zlim = c(-1, 1), axes = FALSE,
  xlab = "fit B factors", ylab = "fit A factors",
  main = "Raw loading correlations"
)
axis(1, at = seq_len(H), labels = paste0("B_F", seq_len(H)))
axis(2, at = seq_len(H), labels = paste0("A_F", H:1), las = 1)
for (i in seq_len(H)) {
  for (j in seq_len(H)) {
    text(j, H - i + 1, sprintf("%.2f", loading_corr_raw[i, j]), cex = 0.9)
  }
}

message("Wrote comparison outputs to: ", normalizePath(out_dir, mustWork = FALSE))
message("Best alignment:")
print(alignment_df)
message("CV summary:")
print(cv_summary)
