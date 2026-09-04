#!/usr/bin/env Rscript

# Compare the item loadings on the Gaussian factor coordinates from two
# full-data IFEval refits:
#   rank 2: H=4, G=(2,3,3,1), whose Gaussian coordinate is F4
#   rank 3: H=4, G=(3,3,1,3), whose Gaussian coordinate is F3
#
# The script writes:
#   1. a row-per-item comparison CSV,
#   2. a compact summary CSV,
#   3. a scatter plot comparing the two Gaussian-factor loadings.

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
  "TOP_REFIT_DIR",
  file.path(
    repo_root,
    "results", "full", "ifeval_threshold_sensitivity_2thirds_20260825",
    "threshold_0p5", "top3_full_data_refits"
  )
)

out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(base_dir, "gaussian_factor_item_comparison")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

threshold <- as.numeric(Sys.getenv("LOADING_THRESHOLD", "0.5"))

rank2_path <- file.path(
  base_dir,
  "rank2-H4-G2-3-3-1-lambda4",
  "openeval_item_intercepts_loadings_metadata.csv"
)
rank3_path <- file.path(
  base_dir,
  "rank3-H4-G3-3-1-3-lambda4",
  "openeval_item_intercepts_loadings_metadata.csv"
)

if (!file.exists(rank2_path)) stop("Missing rank-2 item metadata: ", rank2_path)
if (!file.exists(rank3_path)) stop("Missing rank-3 item metadata: ", rank3_path)

add_primary_factor <- function(x) {
  load_cols <- grep("^loading_factor_", names(x), value = TRUE)
  abs_load <- abs(as.matrix(x[, load_cols, drop = FALSE]))
  x$primary_factor <- sub(
    "loading_factor_",
    "F",
    load_cols[max.col(abs_load, ties.method = "first")]
  )
  x$max_abs_loading <- apply(abs_load, 1L, max)
  x
}

rank2 <- add_primary_factor(read.csv(rank2_path, check.names = FALSE))
rank3 <- add_primary_factor(read.csv(rank3_path, check.names = FALSE))

keep_cols <- c(
  "item_id",
  "empirical_accuracy",
  "alpha",
  "question_snippet",
  "primary_semantic",
  "question_length_words",
  "numeric_tokens",
  "negation_exception",
  "long_context"
)

rank2_small <- rank2[, intersect(keep_cols, names(rank2)), drop = FALSE]
names(rank2_small)[names(rank2_small) == "alpha"] <- "rank2_alpha"
rank2_small$rank2_gaussian_factor <- "F4"
rank2_small$rank2_gaussian_loading <- rank2$loading_factor_4
rank2_small$rank2_abs_gaussian_loading <- abs(rank2$loading_factor_4)
rank2_small$rank2_primary_factor <- rank2$primary_factor
rank2_small$rank2_gaussian_is_primary <- rank2$primary_factor == "F4"

rank3_small <- rank3[, c("item_id", "alpha"), drop = FALSE]
names(rank3_small)[names(rank3_small) == "alpha"] <- "rank3_alpha"
rank3_small$rank3_gaussian_factor <- "F3"
rank3_small$rank3_gaussian_loading <- rank3$loading_factor_3
rank3_small$rank3_abs_gaussian_loading <- abs(rank3$loading_factor_3)
rank3_small$rank3_primary_factor <- rank3$primary_factor
rank3_small$rank3_gaussian_is_primary <- rank3$primary_factor == "F3"

comparison <- merge(rank2_small, rank3_small, by = "item_id", all = FALSE)
comparison$rank2_above_threshold <- comparison$rank2_abs_gaussian_loading >= threshold
comparison$rank3_above_threshold <- comparison$rank3_abs_gaussian_loading >= threshold
comparison$gaussian_overlap_class <- ifelse(
  comparison$rank2_above_threshold & comparison$rank3_above_threshold,
  "both",
  ifelse(
    comparison$rank2_above_threshold,
    "rank2_only",
    ifelse(comparison$rank3_above_threshold, "rank3_only", "neither")
  )
)
comparison$max_abs_gaussian_loading <- pmax(
  comparison$rank2_abs_gaussian_loading,
  comparison$rank3_abs_gaussian_loading
)
comparison <- comparison[order(comparison$gaussian_overlap_class, -comparison$max_abs_gaussian_loading), ]

summary_df <- data.frame(
  quantity = c(
    "items",
    "rank2_nonzero_gaussian_loadings",
    "rank3_nonzero_gaussian_loadings",
    "rank2_abs_ge_threshold",
    "rank3_abs_ge_threshold",
    "rank2_gaussian_primary",
    "rank3_gaussian_primary",
    "abs_ge_threshold_overlap",
    "rank2_only_abs_ge_threshold",
    "rank3_only_abs_ge_threshold",
    "loading_threshold"
  ),
  value = c(
    nrow(comparison),
    sum(abs(comparison$rank2_gaussian_loading) > 1e-8),
    sum(abs(comparison$rank3_gaussian_loading) > 1e-8),
    sum(comparison$rank2_above_threshold),
    sum(comparison$rank3_above_threshold),
    sum(comparison$rank2_gaussian_is_primary),
    sum(comparison$rank3_gaussian_is_primary),
    sum(comparison$rank2_above_threshold & comparison$rank3_above_threshold),
    sum(comparison$rank2_above_threshold & !comparison$rank3_above_threshold),
    sum(!comparison$rank2_above_threshold & comparison$rank3_above_threshold),
    threshold
  )
)

write.csv(
  comparison,
  file.path(out_dir, "rank2_rank3_gaussian_factor_item_comparison.csv"),
  row.names = FALSE
)
write.csv(
  summary_df,
  file.path(out_dir, "rank2_rank3_gaussian_factor_item_summary.csv"),
  row.names = FALSE
)

cols <- c(
  both = "#5B6770",
  rank2_only = "#2F6DAE",
  rank3_only = "#C43C4A",
  neither = "#D1D5DB"
)
plot_file <- file.path(out_dir, "rank2_rank3_gaussian_factor_loading_scatter.png")
png(plot_file, width = 1200, height = 1000, res = 170)
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)
plot(
  comparison$rank2_gaussian_loading,
  comparison$rank3_gaussian_loading,
  pch = 19,
  col = cols[comparison$gaussian_overlap_class],
  xlab = "rank 2 Gaussian loading: F4 in G=(2,3,3,1)",
  ylab = "rank 3 Gaussian loading: F3 in G=(3,3,1,3)",
  main = "IFEval Gaussian-factor item loadings",
  cex = 0.75
)
abline(h = c(-threshold, threshold), lty = 3, col = "#6B7280")
abline(v = c(-threshold, threshold), lty = 3, col = "#6B7280")
abline(h = 0, v = 0, lty = 1, col = "#E5E7EB")
legend(
  "topleft",
  legend = c("both", "rank 2 only", "rank 3 only", "neither"),
  col = cols[c("both", "rank2_only", "rank3_only", "neither")],
  pch = 19,
  bty = "n"
)

message("Wrote comparison CSV: ", file.path(out_dir, "rank2_rank3_gaussian_factor_item_comparison.csv"))
message("Wrote summary CSV: ", file.path(out_dir, "rank2_rank3_gaussian_factor_item_summary.csv"))
message("Wrote scatter plot: ", plot_file)
