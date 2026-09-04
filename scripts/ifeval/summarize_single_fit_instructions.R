#!/usr/bin/env Rscript

# Summarize how original IFEval instruction categories align with one fitted
# independent-mixture probit factor model.  Each item is assigned to the factor
# with the largest absolute loading; optionally, items with multiple loadings
# above LOADING_THRESHOLD are also recorded as cross-loading memberships.

options(stringsAsFactors = FALSE)

fit_dir <- Sys.getenv("FIT_DIR", "")
instruction_long_path <- Sys.getenv(
  "INSTRUCTION_LONG_PATH",
  "data/ifeval_threshold_0p67/openeval_item_instruction_metadata_long.csv"
)
out_dir <- Sys.getenv("OUT_DIR", file.path(fit_dir, "instruction_factor_summary"))
loading_threshold <- as.numeric(Sys.getenv("LOADING_THRESHOLD", "0.5"))

if (!nzchar(fit_dir)) {
  stop("Set FIT_DIR to a directory containing openeval_item_intercepts_loadings_metadata.csv")
}
loadings_path <- file.path(fit_dir, "openeval_item_intercepts_loadings_metadata.csv")
if (!file.exists(loadings_path)) stop("Missing loadings file: ", loadings_path)
if (!file.exists(instruction_long_path)) {
  stop("Missing instruction metadata: ", instruction_long_path)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_df <- read.csv(loadings_path, check.names = FALSE)
inst <- read.csv(instruction_long_path, check.names = FALSE)
load_cols <- grep("^loading_factor_", names(load_df), value = TRUE)
if (!length(load_cols)) stop("No loading_factor_* columns found.")

H <- length(load_cols)
factor_labels <- paste0("F", seq_len(H))
abs_load <- abs(as.matrix(load_df[, load_cols, drop = FALSE]))
primary_index <- max.col(abs_load, ties.method = "first")
max_abs <- abs_load[cbind(seq_len(nrow(abs_load)), primary_index)]
second_abs <- apply(abs_load, 1L, function(z) {
  if (length(z) < 2L) NA_real_ else sort(z, decreasing = TRUE)[2L]
})

load_df$primary_factor <- factor_labels[primary_index]
load_df$max_abs_loading <- max_abs
load_df$second_abs_loading <- second_abs
load_df$cross_loading <- is.finite(second_abs) & second_abs >= loading_threshold

inst_join <- merge(
  inst,
  load_df[, c("item_id", "primary_factor", "max_abs_loading",
              "second_abs_loading", "cross_loading", load_cols), drop = FALSE],
  by = "item_id",
  all.x = FALSE,
  all.y = FALSE
)

unique_item_inst <- unique(inst_join[, c("item_id", "instruction_id",
                                         "instruction_family", "primary_factor")])
background_id <- table(unique(inst[, c("item_id", "instruction_id")])$instruction_id)
background_family <- table(unique(inst[, c("item_id", "instruction_family")])$instruction_family)
total_items <- length(unique(inst$item_id))

make_enrichment <- function(df, key_col, background) {
  keys <- sort(unique(df[[key_col]]))
  out <- list()
  for (key in keys) {
    for (f in factor_labels) {
      n_key_factor <- sum(df[[key_col]] == key & df$primary_factor == f)
      n_factor <- length(unique(df$item_id[df$primary_factor == f]))
      n_key <- as.integer(background[[key]])
      observed_rate <- n_key_factor / max(1L, n_factor)
      background_rate <- n_key / total_items
      out[[length(out) + 1L]] <- data.frame(
        variable = key_col,
        level = key,
        factor = f,
        n_items_with_level_on_factor = n_key_factor,
        n_primary_items_on_factor = n_factor,
        n_background_items_with_level = n_key,
        total_items = total_items,
        proportion_within_factor = observed_rate,
        background_proportion = background_rate,
        lift = observed_rate / background_rate,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

id_enrichment <- make_enrichment(unique_item_inst, "instruction_id", background_id)
family_enrichment <- make_enrichment(unique_item_inst, "instruction_family", background_family)

item_assignments <- load_df[
  order(load_df$primary_factor, -load_df$max_abs_loading),
  c("item_id", "empirical_accuracy", "alpha", "primary_factor",
    "max_abs_loading", "second_abs_loading", "cross_loading",
    load_cols, "question_snippet"),
  drop = FALSE
]

cross_memberships <- list()
for (h in seq_len(H)) {
  keep <- abs(load_df[[load_cols[h]]]) >= loading_threshold
  if (!any(keep)) next
  tmp <- load_df[keep, c("item_id", "primary_factor", "max_abs_loading",
                         "second_abs_loading", "cross_loading",
                         "question_snippet"), drop = FALSE]
  tmp$factor <- factor_labels[h]
  tmp$signed_loading <- load_df[[load_cols[h]]][keep]
  cross_memberships[[length(cross_memberships) + 1L]] <- tmp
}
cross_memberships <- if (length(cross_memberships)) {
  do.call(rbind, cross_memberships)
} else {
  data.frame()
}

write.csv(item_assignments, file.path(out_dir, "item_primary_factor_assignments.csv"), row.names = FALSE)
write.csv(inst_join, file.path(out_dir, "instruction_rows_with_primary_factor.csv"), row.names = FALSE)
write.csv(id_enrichment, file.path(out_dir, "instruction_id_factor_enrichment.csv"), row.names = FALSE)
write.csv(family_enrichment, file.path(out_dir, "instruction_family_factor_enrichment.csv"), row.names = FALSE)
write.csv(cross_memberships, file.path(out_dir, "loading_threshold_factor_memberships.csv"), row.names = FALSE)

plot_heatmap <- function(enrichment, variable_name, value_col, file_name, title) {
  dd <- enrichment[enrichment$variable == variable_name, , drop = FALSE]
  levels_order <- unique(dd$level)
  mat <- matrix(0, nrow = length(levels_order), ncol = H,
                dimnames = list(levels_order, factor_labels))
  for (i in seq_len(nrow(dd))) {
    mat[dd$level[i], dd$factor[i]] <- dd[[value_col]][i]
  }

  pal <- colorRampPalette(c("#f7fbff", "#9ecae1", "#3182bd", "#08519c"))(101)
  max_val <- max(mat, na.rm = TRUE)
  if (!is.finite(max_val) || max_val <= 0) max_val <- 1

  png(file.path(out_dir, file_name), width = 1500, height = 1800, res = 180)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mar = c(5, 16, 4, 3))
  image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1L, , drop = FALSE]),
    col = pal,
    breaks = seq(0, max_val, length.out = length(pal) + 1L),
    axes = FALSE,
    xlab = "primary loading factor",
    ylab = "",
    main = title
  )
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
  grid(nx = ncol(mat), ny = nrow(mat), col = adjustcolor("white", 0.7), lwd = 0.8)
  box()
  dev.off()
}

plot_heatmap(
  id_enrichment,
  "instruction_id",
  "proportion_within_factor",
  "instruction_id_by_primary_factor_proportion.png",
  "IFEval instruction IDs by primary loading factor"
)
plot_heatmap(
  family_enrichment,
  "instruction_family",
  "proportion_within_factor",
  "instruction_family_by_primary_factor_proportion.png",
  "IFEval instruction families by primary loading factor"
)
plot_heatmap(
  id_enrichment,
  "instruction_id",
  "lift",
  "instruction_id_by_primary_factor_lift.png",
  "IFEval instruction ID enrichment by primary loading factor"
)

cat("Wrote instruction/factor summaries to: ", normalizePath(out_dir), "\n", sep = "")
cat("Primary item counts by factor:\n")
print(table(item_assignments$primary_factor))
cat("\nTop instruction-ID enrichment per factor:\n")
for (f in factor_labels) {
  dd <- id_enrichment[id_enrichment$factor == f & id_enrichment$n_items_with_level_on_factor > 0, ]
  dd <- dd[order(-dd$lift, -dd$n_items_with_level_on_factor), ]
  print(head(dd[, c("level", "factor", "n_items_with_level_on_factor",
                   "n_primary_items_on_factor", "proportion_within_factor", "lift")], 8))
}
