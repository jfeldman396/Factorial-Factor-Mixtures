#!/usr/bin/env Rscript

# Validate row coverage, unique scientific keys, and matched DGP seeds for the
# authoritative three-arm recovery study.

options(stringsAsFactors = FALSE)

get_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

split_csv <- function(x) {
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

parse_ints <- function(x) as.integer(split_csv(x))

as_flag <- function(x) toupper(x) %in% c("TRUE", "T", "1", "YES", "Y")

key_string <- function(d, cols) {
  if (!nrow(d)) return(character(0))
  do.call(paste, c(lapply(d[cols], as.character), sep = "|"))
}

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1L]])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

prefix <- get_env(
  "RUN_LABEL_PREFIX",
  "fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full"
)
arms <- intersect(
  split_csv(get_env("ARM_FILTER", "separated,asymmetric_pi,moderate_overlap")),
  c("separated", "asymmetric_pi", "moderate_overlap")
)
n_values <- parse_ints(get_env("N_VALUES", "100,200,400"))
p_product <- parse_ints(get_env("P_VALUES_PRODUCT", "500,1000,1500,2000"))
p_gibbs <- parse_ints(get_env("P_VALUES_GIBBS", "500,1000"))
h_values <- parse_ints(get_env("H_VALUES", "5,10"))
g_values <- vapply(split_csv(get_env("G_CONFIG_TYPES", "all2,all3")), function(x) {
  as.integer(sub("^all", "", x))
}, integer(1L))
rep_values <- parse_ints(get_env("REP_VALUES", paste(seq_len(25L), collapse = ",")))
strict <- as_flag(get_env("STRICT_VALIDATION", "TRUE"))

method_specs <- list()
if (as_flag(get_env("RUN_PRODUCT_MAP", "TRUE"))) {
  method_specs[["independent_marginal_mixture"]] <- p_product
}
if (as_flag(get_env("RUN_VIROLI_LAPLACE", "TRUE"))) {
  method_specs[["viroli_laplace_gibbs"]] <- p_gibbs
}
if (as_flag(get_env("RUN_VIROLI_GAUSSIAN", "FALSE"))) {
  method_specs[["viroli_gaussian_gibbs"]] <- p_gibbs
}

expected_keys <- function(p_values) {
  d <- expand.grid(
    n = n_values,
    p = p_values,
    H_true = h_values,
    G_base = g_values,
    rep = rep_values,
    KEEP.OUT.ATTRS = FALSE
  )
  d$G_config <- mapply(
    function(g, h) paste(rep(g, h), collapse = "-"),
    d$G_base,
    d$H_true,
    USE.NAMES = FALSE
  )
  d[c("n", "p", "H_true", "G_config", "rep")]
}

integrity <- list()
seed_checks <- list()

for (arm in arms) {
  run_label <- paste(prefix, arm, sep = "_")
  result_file <- file.path(repo_root, "results", "full", run_label, "comparison_results.csv")
  if (!file.exists(result_file)) {
    for (method in names(method_specs)) {
      integrity[[length(integrity) + 1L]] <- data.frame(
        arm = arm,
        run_label = run_label,
        method = method,
        expected_rows = nrow(expected_keys(method_specs[[method]])),
        observed_rows = 0L,
        unique_keys = 0L,
        duplicate_keys = 0L,
        missing_keys = nrow(expected_keys(method_specs[[method]])),
        extra_keys = 0L,
        complete = FALSE
      )
    }
    next
  }

  results <- read.csv(result_file, check.names = FALSE)
  key_cols <- c("n", "p", "H_true", "G_config", "rep")
  absent <- setdiff(c("method", key_cols), names(results))
  if (length(absent)) stop("Missing required columns in ", result_file, ": ", paste(absent, collapse = ", "))

  for (method in names(method_specs)) {
    observed <- results[results$method == method, , drop = FALSE]
    expected <- expected_keys(method_specs[[method]])
    observed_key <- key_string(observed, key_cols)
    expected_key <- key_string(expected, key_cols)
    duplicate_count <- sum(duplicated(observed_key))
    missing_count <- length(setdiff(expected_key, observed_key))
    extra_count <- length(setdiff(observed_key, expected_key))
    integrity[[length(integrity) + 1L]] <- data.frame(
      arm = arm,
      run_label = run_label,
      method = method,
      expected_rows = length(expected_key),
      observed_rows = nrow(observed),
      unique_keys = length(unique(observed_key)),
      duplicate_keys = duplicate_count,
      missing_keys = missing_count,
      extra_keys = extra_count,
      complete = duplicate_count == 0L && missing_count == 0L && extra_count == 0L
    )
  }

  product <- results[results$method == "independent_marginal_mixture", , drop = FALSE]
  gibbs <- results[results$method == "viroli_laplace_gibbs", , drop = FALSE]
  seed_cols <- intersect(
    c("data_seed", "loading_parameter_seed", "mixture_parameter_seed"),
    intersect(names(product), names(gibbs))
  )
  if (nrow(product) && nrow(gibbs) && length(seed_cols)) {
    matched <- merge(
      product[c(key_cols, seed_cols)],
      gibbs[c(key_cols, seed_cols)],
      by = key_cols,
      suffixes = c("_product", "_gibbs")
    )
    mismatch <- rep(FALSE, nrow(matched))
    for (seed in seed_cols) {
      mismatch <- mismatch |
        as.character(matched[[paste0(seed, "_product")]]) !=
        as.character(matched[[paste0(seed, "_gibbs")]])
    }
    seed_checks[[length(seed_checks) + 1L]] <- data.frame(
      arm = arm,
      matched_rows = nrow(matched),
      seed_columns_checked = paste(seed_cols, collapse = ","),
      seed_mismatch_rows = sum(mismatch, na.rm = TRUE),
      seeds_match = !any(mismatch, na.rm = TRUE)
    )
  }
}

integrity <- if (length(integrity)) do.call(rbind, integrity) else data.frame()
seed_checks <- if (length(seed_checks)) do.call(rbind, seed_checks) else data.frame()
table_dir <- file.path(repo_root, "results", "selected_tables", "sample_size", "three_mixture_settings")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(integrity, file.path(table_dir, "three_arm_recovery_integrity.csv"), row.names = FALSE)
write.csv(seed_checks, file.path(table_dir, "three_arm_recovery_matched_seed_check.csv"), row.names = FALSE)

print(integrity, row.names = FALSE)
if (nrow(seed_checks)) print(seed_checks, row.names = FALSE)

failed <- !nrow(integrity) || any(!integrity$complete)
seed_failed <- nrow(seed_checks) && any(!seed_checks$seeds_match)
if (strict && (failed || seed_failed)) {
  stop("Three-arm recovery validation failed; inspect the integrity and seed-check CSV files.")
}
