#!/usr/bin/env Rscript

# Assemble the authoritative paper-facing result for each mixture arm. The
# matched p=500/1000 rows come from the completed subtype rerun; Product MAP
# p=1500/2000 rows come from the full recovery run.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1L]])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
saved_root <- file.path(repo_root, "results", "saved", "simulation")
arms <- c("separated", "asymmetric_pi", "moderate_overlap")
product_method <- "independent_marginal_mixture"
matched_methods <- c(product_method, "viroli_laplace_gibbs")

rbind_fill <- function(xs) {
  all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    for (name in setdiff(all_names, names(x))) x[[name]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, xs)
}

for (arm in arms) {
  arm_dir <- file.path(saved_root, arm)
  recovery_file <- file.path(arm_dir, "recovery_results.csv")
  subtype_file <- file.path(arm_dir, "subtype_results.csv")
  if (!file.exists(recovery_file) || !file.exists(subtype_file)) {
    stop("Missing saved recovery or subtype result for arm: ", arm, call. = FALSE)
  }

  recovery <- read.csv(recovery_file, check.names = FALSE)
  subtype <- read.csv(subtype_file, check.names = FALSE)
  matched <- subtype[
    subtype$method %in% matched_methods & subtype$p %in% c(500, 1000),
    , drop = FALSE
  ]
  product_extension <- recovery[
    recovery$method == product_method & recovery$p %in% c(1500, 2000),
    , drop = FALSE
  ]
  release <- rbind_fill(list(matched, product_extension))
  release$mixture_setting <- arm
  key <- do.call(
    paste,
    c(release[c("method", "n", "p", "H_true", "G_config", "rep")], sep = "|")
  )
  if (anyDuplicated(key)) stop("Duplicate release rows for arm: ", arm, call. = FALSE)
  write.csv(release, file.path(arm_dir, "analysis_results.csv"), row.names = FALSE)
  cat(arm, ": ", nrow(release), " rows\n", sep = "")
}
