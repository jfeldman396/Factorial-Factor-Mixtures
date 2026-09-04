#!/usr/bin/env Rscript

# Compare fitted IFEval factor loadings with the original IFEval instruction
# metadata.  For each full-data refit, each item is assigned to the factor on
# which it has the largest absolute loading, provided that loading exceeds a
# threshold.  The script then summarizes which IFEval instruction families and
# instruction IDs are concentrated in each factor.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

metadata_long_path <- Sys.getenv(
  "INSTRUCTION_LONG_PATH",
  file.path(repo_root, "data", "ifeval_threshold_0p5", "openeval_item_instruction_metadata_long.csv")
)
metadata_wide_path <- Sys.getenv(
  "ITEM_METADATA_PATH",
  file.path(repo_root, "data", "ifeval_threshold_0p5", "openeval_item_metadata.csv")
)
top_refit_dir <- Sys.getenv(
  "TOP_REFIT_DIR",
  file.path(
    repo_root,
    "results", "full", "ifeval_threshold_sensitivity_2thirds_20260825",
    "threshold_0p5", "top3_full_data_refits"
  )
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(top_refit_dir, "instruction_metadata_factor_comparison")
)
loading_threshold <- as.numeric(Sys.getenv("LOADING_THRESHOLD", "0.5"))

if (!file.exists(metadata_long_path)) {
  stop("Missing instruction metadata: ", metadata_long_path)
}
if (!file.exists(metadata_wide_path)) {
  stop("Missing item metadata: ", metadata_wide_path)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_dirs <- c(
  "rank1-H4-G3-3-3-3-lambda4",
  "rank2-H4-G2-3-3-1-lambda4",
  "rank3-H4-G3-3-1-3-lambda4"
)
model_labels <- c(
  "rank 1: H=4, G=(3,3,3,3)",
  "rank 2: H=4, G=(2,3,3,1)",
  "rank 3: H=4, G=(3,3,1,3)"
)

instruction_long <- read.csv(metadata_long_path, check.names = FALSE)
instruction_wide <- read.csv(metadata_wide_path, check.names = FALSE)
item_ids_all <- unique(instruction_long$item_id)
background_family <- table(unique(instruction_long[, c("item_id", "instruction_family")])$instruction_family)
background_id <- table(unique(instruction_long[, c("item_id", "instruction_id")])$instruction_id)

collapse_top <- function(tab, background, n_items, total_items, top_n = 5L) {
  if (length(tab) == 0L || n_items == 0L) return("")
  ord_names <- names(sort(tab, decreasing = TRUE))
  ord_names <- ord_names[seq_len(min(top_n, length(ord_names)))]
  vals <- vapply(ord_names, function(nm) {
    observed <- as.numeric(tab[[nm]])
    expected_rate <- as.numeric(background[[nm]]) / total_items
    observed_rate <- observed / n_items
    lift <- observed_rate / expected_rate
    sprintf("%s: %d, lift %.2f", nm, observed, lift)
  }, character(1L))
  paste(vals, collapse = "; ")
}

factor_rows <- list()
family_rows <- list()
instruction_rows <- list()
example_rows <- list()

for (m in seq_along(model_dirs)) {
  model_dir <- file.path(top_refit_dir, model_dirs[m])
  item_path <- file.path(model_dir, "openeval_item_intercepts_loadings_metadata.csv")
  if (!file.exists(item_path)) stop("Missing item/loadings file: ", item_path)

  item_df <- read.csv(item_path, check.names = FALSE)
  item_df <- merge(
    item_df,
    instruction_wide[, c("item_id", "instruction_ids", "instruction_families"), drop = FALSE],
    by = "item_id",
    all.x = TRUE
  )
  load_cols <- grep("^loading_factor_", names(item_df), value = TRUE)
  abs_loadings <- abs(as.matrix(item_df[, load_cols, drop = FALSE]))
  primary_col <- load_cols[max.col(abs_loadings, ties.method = "first")]

  item_df$primary_factor <- sub("loading_factor_", "F", primary_col)
  item_df$max_abs_loading <- apply(abs_loadings, 1L, max)

  for (h in seq_along(load_cols)) {
    factor_name <- paste0("F", h)
    selected <- item_df[
      item_df$primary_factor == factor_name & item_df$max_abs_loading >= loading_threshold,
      ,
      drop = FALSE
    ]
    selected_ids <- selected$item_id

    family_item <- unique(
      instruction_long[instruction_long$item_id %in% selected_ids, c("item_id", "instruction_family")]
    )
    instruction_item <- unique(
      instruction_long[instruction_long$item_id %in% selected_ids, c("item_id", "instruction_id", "instruction_family")]
    )
    family_tab <- table(family_item$instruction_family)
    instruction_tab <- table(instruction_item$instruction_id)

    factor_rows[[length(factor_rows) + 1L]] <- data.frame(
      model = model_labels[m],
      model_dir = model_dirs[m],
      factor = factor_name,
      n_primary_loading_items = nrow(selected),
      mean_signed_loading = if (nrow(selected)) mean(selected[[load_cols[h]]]) else NA_real_,
      mean_abs_loading = if (nrow(selected)) mean(abs(selected[[load_cols[h]]])) else NA_real_,
      top_instruction_families = collapse_top(
        family_tab, background_family, max(1L, nrow(selected)), length(item_ids_all), 5L
      ),
      top_instruction_ids = collapse_top(
        instruction_tab, background_id, max(1L, nrow(selected)), length(item_ids_all), 6L
      ),
      stringsAsFactors = FALSE
    )

    if (length(family_tab)) {
      for (nm in names(family_tab)) {
        family_rows[[length(family_rows) + 1L]] <- data.frame(
          model = model_labels[m],
          factor = factor_name,
          instruction_family = nm,
          n_items_with_family = as.integer(family_tab[[nm]]),
          selected_item_count = nrow(selected),
          background_item_count = as.integer(background_family[[nm]]),
          lift = (as.numeric(family_tab[[nm]]) / max(1L, nrow(selected))) /
            (as.numeric(background_family[[nm]]) / length(item_ids_all)),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(instruction_tab)) {
      for (nm in names(instruction_tab)) {
        fam <- instruction_item$instruction_family[match(nm, instruction_item$instruction_id)]
        instruction_rows[[length(instruction_rows) + 1L]] <- data.frame(
          model = model_labels[m],
          factor = factor_name,
          instruction_id = nm,
          instruction_family = fam,
          n_items_with_instruction = as.integer(instruction_tab[[nm]]),
          selected_item_count = nrow(selected),
          background_item_count = as.integer(background_id[[nm]]),
          lift = (as.numeric(instruction_tab[[nm]]) / max(1L, nrow(selected))) /
            (as.numeric(background_id[[nm]]) / length(item_ids_all)),
          stringsAsFactors = FALSE
        )
      }
    }

    if (nrow(selected)) {
      selected$signed_loading <- selected[[load_cols[h]]]
      selected <- selected[order(-abs(selected$signed_loading)), ]
      examples <- head(
        selected[, c("item_id", "empirical_accuracy", "signed_loading", "instruction_ids",
                     "instruction_families", "question_snippet")],
        5L
      )
      examples$model <- model_labels[m]
      examples$factor <- factor_name
      example_rows[[length(example_rows) + 1L]] <- examples[
        , c("model", "factor", "item_id", "empirical_accuracy", "signed_loading",
            "instruction_ids", "instruction_families", "question_snippet")
      ]
    }
  }
}

factor_summary <- do.call(rbind, factor_rows)
family_summary <- do.call(rbind, family_rows)
instruction_summary <- do.call(rbind, instruction_rows)
example_summary <- do.call(rbind, example_rows)

family_summary <- family_summary[order(family_summary$model, family_summary$factor, -family_summary$lift), ]
instruction_summary <- instruction_summary[
  order(instruction_summary$model, instruction_summary$factor, -instruction_summary$lift),
]

write.csv(
  factor_summary,
  file.path(out_dir, "factor_instruction_metadata_summary.csv"),
  row.names = FALSE
)
write.csv(
  family_summary,
  file.path(out_dir, "factor_instruction_family_enrichment.csv"),
  row.names = FALSE
)
write.csv(
  instruction_summary,
  file.path(out_dir, "factor_instruction_id_enrichment.csv"),
  row.names = FALSE
)
write.csv(
  example_summary,
  file.path(out_dir, "factor_instruction_examples.csv"),
  row.names = FALSE
)

message("Wrote factor summary: ", file.path(out_dir, "factor_instruction_metadata_summary.csv"))
message("Wrote family enrichment: ", file.path(out_dir, "factor_instruction_family_enrichment.csv"))
message("Wrote instruction enrichment: ", file.path(out_dir, "factor_instruction_id_enrichment.csv"))
message("Wrote examples: ", file.path(out_dir, "factor_instruction_examples.csv"))
