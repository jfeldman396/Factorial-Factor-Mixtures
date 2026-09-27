#!/usr/bin/env Rscript

# Export compact IFEval results and manuscript artifacts from a completed
# threshold-sensitivity run. Fitted R objects, worker logs, and histories stay
# in results/full and are not part of the release.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[[1L]])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

split_fields <- function(x) {
  fields <- trimws(strsplit(x, "[ ,]+")[[1L]])
  fields[nzchar(fields)]
}

copy_files <- function(files, destination) {
  files <- files[file.exists(files)]
  if (!length(files)) return(character(0))
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  targets <- file.path(destination, basename(files))
  ok <- file.copy(files, targets, overwrite = TRUE)
  if (!all(ok)) stop("Could not copy: ", paste(files[!ok], collapse = ", "))
  targets
}

full_root <- get_env(
  "IFEVAL_FULL_ROOT",
  file.path(repo_root, "results", "full", "ifeval_threshold_sensitivity_current")
)
saved_root <- get_env(
  "IFEVAL_SAVED_ROOT",
  file.path(repo_root, "results", "saved", "ifeval")
)
plot_root <- get_env(
  "IFEVAL_PLOT_ROOT",
  file.path(repo_root, "results", "selected_plots", "ifeval")
)
table_root <- get_env(
  "IFEVAL_TABLE_ROOT",
  file.path(repo_root, "results", "selected_tables", "ifeval")
)
thresholds <- split_fields(get_env("THRESHOLD_LABELS", "0p5 0p67 1"))

manifest <- list()
for (threshold in thresholds) {
  threshold_root <- file.path(full_root, paste0("threshold_", threshold))
  cv_dir <- file.path(threshold_root, "componentwise_cv")
  required <- file.path(
    cv_dir,
    c(
      "ifeval_rank_lambda_cv_fold_scores.csv",
      "ifeval_rank_lambda_cv_summary.csv",
      "ifeval_rank_lambda_selected_by_heldout_ll.csv"
    )
  )
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop(
      "Threshold ", threshold, " is incomplete; missing: ",
      paste(basename(missing), collapse = ", "),
      call. = FALSE
    )
  }

  saved_threshold <- file.path(saved_root, paste0("threshold_", threshold))
  copied <- copy_files(required, file.path(saved_threshold, "model_selection"))

  selected <- read.csv(required[[3L]], check.names = FALSE)
  if (nrow(selected) != 1L) stop("Expected one selected model for threshold ", threshold)
  selected_dir <- file.path(
    cv_dir,
    sprintf(
      "selected_mixture_H%d_Gconfig%s_lambda%s",
      selected$H[[1L]],
      gsub(",", "-", selected$G_config[[1L]], fixed = TRUE),
      gsub("\\.", "p", as.character(selected$lambda_l1_penalty[[1L]]))
    )
  )
  if (!dir.exists(selected_dir)) stop("Missing selected-model directory: ", selected_dir)

  selected_tables <- list.files(
    selected_dir, pattern = "\\.(csv|json)$", full.names = TRUE, ignore.case = TRUE
  )
  copied <- c(copied, copy_files(selected_tables, file.path(saved_threshold, "selected_model")))

  cv_plots <- list.files(
    cv_dir, pattern = "\\.(png|pdf)$", full.names = TRUE, ignore.case = TRUE
  )
  selected_plots <- list.files(
    selected_dir, pattern = "\\.(png|pdf)$", full.names = TRUE, ignore.case = TRUE
  )
  copied <- c(
    copied,
    copy_files(cv_plots, file.path(plot_root, paste0("threshold_", threshold), "model_selection")),
    copy_files(selected_plots, file.path(plot_root, paste0("threshold_", threshold), "selected_model"))
  )

  visualization_dir <- file.path(threshold_root, "selected_visualizations")
  if (dir.exists(visualization_dir)) {
    viz_files <- list.files(visualization_dir, recursive = TRUE, full.names = TRUE)
    viz_plots <- viz_files[grepl("\\.(png|pdf)$", viz_files, ignore.case = TRUE)]
    viz_tables <- viz_files[grepl("\\.(csv|tex)$", viz_files, ignore.case = TRUE)]
    copied <- c(
      copied,
      copy_files(viz_plots, file.path(plot_root, paste0("threshold_", threshold), "selected_model")),
      copy_files(viz_tables, file.path(table_root, paste0("threshold_", threshold)))
    )
  }

  manifest[[length(manifest) + 1L]] <- data.frame(
    threshold = threshold,
    path = normalizePath(copied, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest)
dir.create(saved_root, recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, file.path(saved_root, "artifact_manifest.csv"), row.names = FALSE)
cat("IFEval release exported from: ", normalizePath(full_root), "\n", sep = "")
