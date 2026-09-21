#!/usr/bin/env Rscript

# Assemble the authoritative recovery artifacts from the completed separated
# study and the targeted asymmetric/overlap robustness runs. Every row retains
# its source artifact so differing replication coverage cannot be hidden.

options(stringsAsFactors = FALSE)

rbind_fill <- function(x) {
  x <- x[vapply(x, nrow, integer(1L)) > 0L]
  if (!length(x)) return(data.frame())
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(d) {
    for (nm in setdiff(all_names, names(d))) d[[nm]] <- NA
    d[, all_names, drop = FALSE]
  })
  do.call(rbind, x)
}

read_source <- function(repo_root, run_label, source_scope) {
  path <- file.path(repo_root, "results", "full", run_label, "comparison_results.csv")
  if (!file.exists(path)) stop("Missing source result: ", path)
  d <- read.csv(path, check.names = FALSE)
  d$source_run <- run_label
  d$source_scope <- source_scope
  d
}

key_string <- function(d) {
  do.call(paste, c(lapply(
    d[c("method", "n", "p", "H_true", "G_config", "rep")],
    as.character
  ), sep = "|"))
}

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1L]])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

output_prefix <- Sys.getenv(
  "RUN_LABEL_PREFIX",
  unset = "fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full"
)
table_dir <- file.path(
  repo_root,
  "results", "selected_tables", "sample_size", "three_mixture_settings"
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

source_labels <- list(
  separated = list(
    complete = "fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_h5_h10_final_combined"
  ),
  asymmetric_pi = list(
    n100 = "fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda3_n100_asymmetric_pi",
    product_completion = "fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_pmax2000_product_completion_asymmetric_pi",
    gibbs_check = "fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_asymmetric_pi_check"
  ),
  moderate_overlap = list(
    n100 = "fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda3_n100_moderate_overlap",
    product_completion = "fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_pmax2000_product_completion_moderate_overlap",
    gibbs_check = "fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_moderate_overlap_check"
  )
)

assemble_arm <- function(arm) {
  labels <- source_labels[[arm]]
  if (arm == "separated") {
    d <- read_source(repo_root, labels$complete, "full_25_rep_study")
    d <- d[d$method %in% c("independent_marginal_mixture", "viroli_laplace_gibbs"), , drop = FALSE]
  } else {
    n100 <- read_source(repo_root, labels$n100, "full_25_rep_n100")
    completion <- read_source(
      repo_root,
      labels$product_completion,
      "full_25_rep_product_n200_n400"
    )
    gibbs <- read_source(
      repo_root,
      labels$gibbs_check,
      "five_rep_gibbs_robustness_n200_n400"
    )
    d <- rbind_fill(list(
      n100[n100$n == 100 & n100$method %in% c(
        "independent_marginal_mixture", "viroli_laplace_gibbs"
      ), , drop = FALSE],
      completion[completion$n %in% c(200, 400) &
        completion$method == "independent_marginal_mixture", , drop = FALSE],
      gibbs[gibbs$n %in% c(200, 400) &
        gibbs$method == "viroli_laplace_gibbs", , drop = FALSE]
    ))
  }

  d$mixture_setting <- arm
  keys <- key_string(d)
  if (anyDuplicated(keys)) {
    stop("Duplicate scientific keys in assembled ", arm, " results.")
  }

  # Product MAP is complete at 25 reps across all p. Gibbs is complete at 25
  # reps for the separated study and n=100 robustness cells; the n=200/400
  # robustness checks intentionally contain five reps at p=500/1000.
  expected_reps <- ifelse(
    d$method == "independent_marginal_mixture" |
      arm == "separated" |
      d$n == 100,
    25L,
    5L
  )
  cell_id <- with(d, paste(method, n, p, H_true, G_config, sep = "|"))
  observed_by_cell <- table(cell_id)
  expected_by_cell <- tapply(expected_reps, cell_id, unique)
  if (any(lengths(expected_by_cell) != 1L)) {
    stop("Inconsistent expected replication counts in ", arm, ".")
  }
  expected_by_cell <- vapply(expected_by_cell, `[[`, integer(1L), 1L)
  observed_counts <- as.integer(observed_by_cell[names(expected_by_cell)])
  expected_counts <- unname(as.integer(expected_by_cell))
  if (any(is.na(observed_counts)) || any(observed_counts != expected_counts)) {
    bad <- names(expected_by_cell)[
      is.na(observed_counts) | observed_counts != expected_counts
    ]
    stop("Incomplete cells in ", arm, ": ", paste(bad, collapse = ", "))
  }

  output_label <- paste(output_prefix, arm, sep = "_")
  output_dir <- file.path(repo_root, "results", "full", output_label)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(d, file.path(output_dir, "comparison_results.csv"), row.names = FALSE)

  coverage <- aggregate(
    rep ~ mixture_setting + source_run + source_scope + method + n + p + H_true + G_config,
    data = d,
    FUN = function(x) length(unique(x))
  )
  names(coverage)[names(coverage) == "rep"] <- "completed_reps"
  coverage$expected_reps <- ifelse(
    coverage$method == "independent_marginal_mixture" |
      arm == "separated" |
      coverage$n == 100,
    25L,
    5L
  )
  coverage$complete <- coverage$completed_reps == coverage$expected_reps
  list(data = d, coverage = coverage, output_label = output_label)
}

assembled <- lapply(names(source_labels), assemble_arm)
names(assembled) <- names(source_labels)

coverage <- do.call(rbind, lapply(assembled, `[[`, "coverage"))
write.csv(
  coverage,
  file.path(table_dir, "three_arm_recovery_coverage_and_provenance.csv"),
  row.names = FALSE
)

arm_summary <- do.call(rbind, lapply(names(assembled), function(arm) {
  d <- assembled[[arm]]$data
  data.frame(
    mixture_setting = arm,
    method = names(table(d$method)),
    rows = as.integer(table(d$method)),
    unique_cells = as.integer(tapply(
      key_string(d),
      d$method,
      function(x) length(unique(sub("\\|[^|]+$", "", x)))
    )[names(table(d$method))])
  )
}))
write.csv(
  arm_summary,
  file.path(table_dir, "three_arm_recovery_row_summary.csv"),
  row.names = FALSE
)

cat("Assembled authoritative recovery artifacts:\n")
for (arm in names(assembled)) {
  d <- assembled[[arm]]$data
  cat("  ", arm, ": ", nrow(d), " rows -> results/full/",
      assembled[[arm]]$output_label, "/comparison_results.csv\n", sep = "")
}
cat("Coverage table: results/selected_tables/sample_size/three_mixture_settings/",
    "three_arm_recovery_coverage_and_provenance.csv\n", sep = "")
