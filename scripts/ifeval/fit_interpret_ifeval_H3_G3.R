#!/usr/bin/env Rscript

# Fit and interpret the selected OpenEval independent-mixture probit factor model.
#
# Input:
#   openeval_binary_matrix_complete_nonconstant.csv
#   openeval_item_metadata.csv
#
# Model:
#   X_ij | f_i ~ Bernoulli(Phi(alpha_j + lambda_j' f_i))
#   f_ih      ~ independent G-component Gaussian mixture
#
# This reproducible copy defaults to the IFEval H=3, G=3 analysis.

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
source(file.path(script_dir, "math500_intercept_imfm_fit.R"))

matrix_path <- Sys.getenv(
  "MATRIX_PATH",
  file.path(bundle_root, "data", "openeval_ifeval_only_binary_matrix.csv")
)
item_metadata_path <- Sys.getenv(
  "ITEM_METADATA_PATH",
  file.path(bundle_root, "data", "openeval_item_metadata.csv")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  file.path(bundle_root, "results", "reproduced_openeval_ifeval_H3_G3_interpretation")
)
parse_int_vector <- function(x, default) {
  if (is.na(x) || !nzchar(x)) return(as.integer(default))
  as.integer(strsplit(gsub("[[:space:]]+", "", x), ",", fixed = TRUE)[[1L]])
}

H_fixed <- as.integer(Sys.getenv("H_FIXED", "3"))
G_fixed <- parse_int_vector(Sys.getenv("G_FIXED", "3"), 3L)
if (!(length(G_fixed) %in% c(1L, H_fixed))) {
  stop("G_FIXED must be scalar or a comma-separated length-H_FIXED vector.")
}
G_label <- paste(G_fixed, collapse = "-")
G_report <- paste(G_fixed, collapse = ",")
workers <- as.integer(Sys.getenv("WORKERS", "8"))
top_n <- as.integer(Sys.getenv("TOP_N", "75"))
pretrain_loading_penalty <- as.numeric(Sys.getenv("PRETRAIN_LOADING_PENALTY", "0.05"))
refinement_lambda_l1_penalty <- as.numeric(Sys.getenv("REFINEMENT_LAMBDA_L1_PENALTY", "2"))
require_mixture_convergence <- toupper(Sys.getenv("REQUIRE_MIXTURE_CONVERGENCE", "FALSE")) %in%
  c("TRUE", "T", "1", "YES", "Y")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_binary_matrix <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  X <- as.matrix(raw[, -1L, drop = FALSE])
  storage.mode(X) <- "numeric"
  rownames(X) <- raw[[1L]]
  X <- X[, colSums(is.na(X)) == 0L, drop = FALSE]
  X <- X[, colSums(X == 1) > 0L & colSums(X == 0) > 0L, drop = FALSE]
  X
}

safe_json <- function(x) {
  tryCatch(jsonlite::fromJSON(x, simplifyVector = FALSE), error = function(e) NULL)
}

collapse_ws <- function(x) {
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

extract_first_string <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.character(x) && length(x) > 0L) return(x[[1L]])
  if (is.list(x)) {
    for (elt in x) {
      val <- extract_first_string(elt)
      if (!is.na(val)) return(val)
    }
  }
  NA_character_
}

extract_item_fields <- function(item_content) {
  outer <- safe_json(item_content)
  input_text <- NA_character_
  reference_text <- NA_character_
  category <- NA_character_
  source <- NA_character_
  prompt_type <- NA_character_

  if (!is.null(outer)) {
    if (!is.null(outer$input)) {
      input_text <- extract_first_string(outer$input)
    }
    if (!is.null(outer$references)) {
      reference_text <- extract_first_string(outer$references)
    }
  }

  inner <- if (!is.na(input_text)) safe_json(input_text) else NULL
  if (!is.null(inner)) {
    if (!is.null(inner$question)) {
      question <- inner$question
      prompt_type <- "question"
    } else if (!is.null(inner$Question)) {
      question <- inner$Question
      prompt_type <- "question"
    } else if (!is.null(inner$prompt)) {
      question <- inner$prompt
      prompt_type <- "instruction_prompt"
    } else if (!is.null(inner$`Pre-Revision Question`)) {
      question <- inner$`Pre-Revision Question`
      prompt_type <- "question"
    } else {
      question <- input_text
      prompt_type <- "raw_input"
    }
    if (!is.null(inner$category)) category <- as.character(inner$category)
    if (!is.null(inner$src)) source <- as.character(inner$src)
  } else {
    question <- input_text
    prompt_type <- "raw_input"
  }

  data.frame(
    prompt_type = prompt_type,
    question_text = collapse_ws(as.character(question)),
    reference_text = collapse_ws(as.character(reference_text)),
    category = category,
    source = source,
    stringsAsFactors = FALSE
  )
}

count_pat <- function(x, pattern) {
  vapply(gregexpr(pattern, x, ignore.case = TRUE, perl = TRUE), function(m) {
    if (length(m) == 1L && m[1L] < 0L) 0L else length(m)
  }, integer(1))
}

has_pat <- function(x, pattern) {
  grepl(pattern, x, ignore.case = TRUE, perl = TRUE)
}

semantic_tags <- function(question_text, benchmark) {
  q <- ifelse(is.na(question_text), "", question_text)
  words <- count_pat(q, "\\S+")
  nums <- count_pat(q, "\\b[0-9]+(?:\\.[0-9]+)?\\b")
  tags <- data.frame(
    question_length_words = words,
    numeric_tokens = nums,
    is_instruction = benchmark == "ifeval",
    is_science = benchmark == "gpqa",
    is_mmlu_pro = benchmark == "mmlu_pro",
    quantitative_symbolic = nums >= 2L | has_pat(q, "[=<>^]|sqrt|\\blog\\b|\\bsin\\b|\\bcos\\b|\\bprobability\\b|\\bpercent\\b|\\bratio\\b|\\bderive\\b|binomial|volatility|option"),
    negation_exception = has_pat(q, "\\b(not|except|least|false|incorrect|does not|cannot|which.*not)\\b"),
    mechanism_causal = has_pat(q, "\\b(mechanism|process|effect|pathway|cause|causes|function|role|responsible|therapy|gene|protein|cell|molecule)\\b"),
    long_context = words >= stats::quantile(words, 0.75, na.rm = TRUE),
    short_direct = words <= stats::quantile(words, 0.25, na.rm = TRUE),
    instruction_constraints = benchmark == "ifeval" | has_pat(q, "\\b(write|respond|include|do not|must|exactly|words|paragraphs|bullets|format|markdown)\\b"),
    stringsAsFactors = FALSE
  )

  tags$primary_semantic <- "general knowledge/reasoning"
  tags$primary_semantic[tags$short_direct] <- "short direct QA"
  tags$primary_semantic[tags$long_context] <- "long-context reasoning"
  tags$primary_semantic[tags$negation_exception] <- "negation/exception discrimination"
  tags$primary_semantic[tags$mechanism_causal] <- "mechanistic/scientific reasoning"
  tags$primary_semantic[tags$quantitative_symbolic] <- "quantitative/symbolic reasoning"
  tags$primary_semantic[tags$instruction_constraints] <- "instruction/constraint following"
  tags
}

summarize_mixture_profiles_ordered <- function(fit) {
  do.call(rbind, lapply(seq_len(fit$H), function(h) {
    cls <- fit$class_map[, h]
    tab <- tabulate(cls, nbins = fit$G_hat[h])
    data.frame(
      factor = paste0("F", h),
      group = seq_len(fit$G_hat[h]),
      n_models = tab,
      prop_models = tab / length(cls),
      mean = fit$mixture_fits[[h]]$mu,
      sd = sqrt(fit$mixture_fits[[h]]$var),
      weight = fit$mixture_fits[[h]]$pi
    )
  }))
}

plot_factor_scatter <- function(scores, out_dir) {
  png(file.path(out_dir, "openeval_factor_scatter_F1_F2.png"), width = 1300, height = 950, res = 160)
  op <- par(mar = c(5, 5, 3, 2))
  on.exit(par(op), add = TRUE)
  acc <- scores$accuracy
  pal <- colorRampPalette(c("#355C9A", "#F2C14E", "#B23A48"))(100)
  idx <- pmax(1, pmin(100, as.integer(cut(acc, breaks = 100, labels = FALSE))))
  plot(
    scores$factor_1,
    scores$factor_2,
    pch = 19,
    col = pal[idx],
    xlab = "factor 1 score",
    ylab = "factor 2 score",
    main = "OpenEval LLMs in refined factor space; color = accuracy"
  )
  elite <- head(order(scores$accuracy, decreasing = TRUE), 12)
  text(scores$factor_1[elite], scores$factor_2[elite], labels = scores$model_id[elite], pos = 3, cex = 0.45)
  box()
  dev.off()
}

plot_factor_scores_heatmap <- function(scores, H, out_dir) {
  png(file.path(out_dir, "openeval_factor_scores_by_llm.png"), width = 3200, height = 850, res = 180)
  op <- par(mar = c(12, 5, 4, 2))
  on.exit(par(op), add = TRUE)
  ord <- order(scores$accuracy, decreasing = TRUE)
  mat <- t(as.matrix(scores[ord, paste0("factor_", seq_len(H)), drop = FALSE]))
  colnames(mat) <- scores$model_id[ord]
  rownames(mat) <- paste0("F", seq_len(H))
  max_abs <- max(abs(mat), na.rm = TRUE)
  pal <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat),
        col = pal, breaks = seq(-max_abs, max_abs, length.out = 102),
        axes = FALSE, xlab = "", ylab = "factor",
        main = "Refined factor scores by LLM, ordered by accuracy")
  axis(2, at = seq_len(H), labels = rownames(mat), las = 1)
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.28, tick = FALSE)
  box()
  dev.off()
}

plot_group_heatmap <- function(scores, H, out_dir) {
  png(file.path(out_dir, "openeval_factor_groups_by_llm.png"), width = 3200, height = 850, res = 180)
  op <- par(mar = c(12, 5, 4, 2))
  on.exit(par(op), add = TRUE)
  ord <- order(scores$accuracy, decreasing = TRUE)
  mat <- t(as.matrix(scores[ord, paste0("group_factor_", seq_len(H)), drop = FALSE]))
  colnames(mat) <- scores$model_id[ord]
  rownames(mat) <- paste0("F", seq_len(H))
  G_max <- max(mat, na.rm = TRUE)
  group_cols <- colorRampPalette(c("#3B6EA8", "#F2C14E", "#9E2F44"))(G_max)
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat),
        col = group_cols, breaks = seq(0.5, G_max + 0.5, by = 1),
        axes = FALSE, xlab = "", ylab = "factor",
        main = "MAP mixture groups by factor and LLM, ordered by accuracy")
  axis(2, at = seq_len(H), labels = rownames(mat), las = 1)
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.28, tick = FALSE)
  legend("topright", fill = group_cols, legend = paste0("group ", seq_len(G_max)), bty = "n")
  box()
  dev.off()
}

plot_top_loading_heatmap <- function(loadings, H, out_dir, top_n = 140L) {
  Lambda <- as.matrix(loadings[, paste0("loading_factor_", seq_len(H)), drop = FALSE])
  rownames(Lambda) <- loadings$item_id
  keep <- head(order(rowSums(abs(Lambda)), decreasing = TRUE), min(top_n, nrow(Lambda)))
  mat <- Lambda[keep, , drop = FALSE]
  png(file.path(out_dir, "openeval_top_item_loading_heatmap.png"), width = 1250, height = 1650, res = 160)
  op <- par(mar = c(4, 9, 3, 2))
  on.exit(par(op), add = TRUE)
  max_abs <- max(abs(mat), na.rm = TRUE)
  pal <- colorRampPalette(c("#355C9A", "#F7F7F7", "#B23A48"))(101)
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1L, , drop = FALSE]),
        col = pal, breaks = seq(-max_abs, max_abs, length.out = 102),
        axes = FALSE, xlab = "factor", ylab = "top-loading items",
        main = "Largest refined item loadings")
  axis(1, at = seq_len(H), labels = paste0("F", seq_len(H)))
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.38)
  box()
  dev.off()
}

X <- read_binary_matrix(matrix_path)
model_id <- rownames(X)
n_aug_iter <- as.integer(Sys.getenv("PRETRAIN_AUG_ITER", "8"))
n_refine_iter <- as.integer(Sys.getenv("REFINE_ITER", "8"))
mixture_max_iter <- as.integer(Sys.getenv("MIXTURE_MAX_ITER", "20"))

t0 <- proc.time()[["elapsed"]]
pre <- fit_binary_probit_pretraining_intercept(
  X = X,
  H = H_fixed,
  G_fixed = G_fixed,
  n_aug_iter = n_aug_iter,
  z_update = "expectation",
  n_random_starts = 1L,
  max_outer = 4L,
  n_mix_starts = 3L,
  mixture_update = "map",
  mu_prior_kappa = 0.05,
  var_prior_shape = 4,
  var_prior_scale = 0.35,
  weight_prior_alpha = 1.2,
  loading_penalty = pretrain_loading_penalty,
  objective_tolerance = 5e-4,
  objective_tolerance_scale = "per_response",
  min_aug_iter = 5L,
  parallel = TRUE,
  workers = workers,
  seed = 20260724L,
  verbose = FALSE
)

ref <- fit_binary_probit_refinement_intercept(
  X = X,
  pretrain_fit = pre,
  n_refine_iter = n_refine_iter,
  maxit_per_subject = 60L,
  n_mix_starts = 3L,
  mixture_max_iter = mixture_max_iter,
  min_mixture_var = 0.05,
  mixture_update = "map",
  mu_prior_kappa = 0.05,
  var_prior_shape = 4,
  var_prior_scale = 0.35,
  weight_prior_alpha = 1.2,
  mixture_prior_weight = 0.2,
  lambda_l1_penalty = refinement_lambda_l1_penalty,
  objective_tolerance = 2e-4,
  objective_tolerance_scale = "relative_total",
  min_refine_iter = 4L,
  require_mixture_convergence = require_mixture_convergence,
  keep_best_binary_iterate = TRUE,
  parallel = TRUE,
  workers = workers,
  verbose = FALSE
)
elapsed <- proc.time()[["elapsed"]] - t0

ref$H <- ncol(ref$F_hat)
ref <- orient_factors_by_accuracy(ref)

scores <- data.frame(
  model_id = model_id,
  accuracy = rowMeans(X),
  ref$F_hat,
  ref$class_map,
  profile_id = ref$profile_id,
  check.names = FALSE
)
names(scores)[3:(2 + ref$H)] <- paste0("factor_", seq_len(ref$H))
names(scores)[(3 + ref$H):(2 + 2 * ref$H)] <- paste0("group_factor_", seq_len(ref$H))

loadings <- data.frame(
  item_id = colnames(X),
  empirical_accuracy = colMeans(X),
  alpha = ref$alpha_hat,
  ref$Lambda_hat,
  check.names = FALSE
)
names(loadings)[4:ncol(loadings)] <- paste0("loading_factor_", seq_len(ref$H))

item_meta <- read.csv(item_metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
parsed <- do.call(rbind, lapply(item_meta$item_content, extract_item_fields))
item_meta <- cbind(item_meta[, c("benchmark", "item_id", "schema_version")], parsed)
item_meta$question_snippet <- substr(item_meta$question_text, 1L, 360L)
tags <- semantic_tags(item_meta$question_text, item_meta$benchmark)
item_meta <- cbind(item_meta, tags)

load_meta <- merge(loadings, item_meta, by = "item_id", all.x = TRUE, sort = FALSE)
write.csv(scores, file.path(out_dir, "openeval_model_factor_scores_profiles.csv"), row.names = FALSE)
write.csv(load_meta, file.path(out_dir, "openeval_item_intercepts_loadings_metadata.csv"), row.names = FALSE)
write.csv(summarize_mixture_profiles_ordered(ref), file.path(out_dir, "openeval_factor_mixture_groups.csv"), row.names = FALSE)
write.csv(ref$pretraining$history, file.path(out_dir, "openeval_pretraining_history.csv"), row.names = FALSE)
write.csv(ref$joint_refinement$history, file.path(out_dir, "openeval_refinement_history.csv"), row.names = FALSE)
fit_label <- sprintf("openeval_H%d_G%s", H_fixed, G_label)
saveRDS(ref, file.path(out_dir, paste0(fit_label, "_fit.rds")))

factor_rows <- list()
example_rows <- list()
group_rows <- list()
feature_cols <- c(
  "is_instruction", "is_science", "is_mmlu_pro", "quantitative_symbolic",
  "negation_exception", "mechanism_causal", "long_context", "short_direct",
  "instruction_constraints"
)

for (h in seq_len(ref$H)) {
  loading_col <- paste0("loading_factor_", h)
  for (direction in c("positive", "negative")) {
    ord <- if (direction == "positive") {
      order(load_meta[[loading_col]], decreasing = TRUE)
    } else {
      order(load_meta[[loading_col]])
    }
    top <- load_meta[ord[seq_len(min(top_n, nrow(load_meta)))], ]
    bench_tab <- sort(table(top$benchmark), decreasing = TRUE)
    sem_tab <- sort(table(top$primary_semantic), decreasing = TRUE)
    factor_rows[[length(factor_rows) + 1L]] <- data.frame(
      factor = paste0("F", h),
      direction = direction,
      n_items = nrow(top),
      mean_signed_loading = mean(top[[loading_col]], na.rm = TRUE),
      mean_abs_loading = mean(abs(top[[loading_col]]), na.rm = TRUE),
      mean_item_accuracy = mean(top$empirical_accuracy, na.rm = TRUE),
      top_benchmarks = paste(paste0(names(bench_tab), "=", as.integer(bench_tab)), collapse = "; "),
      top_semantic_forms = paste(head(names(sem_tab), 4L), collapse = "; "),
      stringsAsFactors = FALSE
    )

    ex <- head(top, 12L)
    ex$factor <- paste0("F", h)
    ex$direction <- direction
    ex$signed_loading <- ex[[loading_col]]
    example_rows[[length(example_rows) + 1L]] <- ex[
      ,
      c("factor", "direction", "item_id", "benchmark", "category", "source",
        "empirical_accuracy", "alpha", "signed_loading", "primary_semantic",
        "question_snippet")
    ]

    top_items <- top$item_id
    tmp <- data.frame(
      model_id = rownames(X),
      factor = paste0("F", h),
      direction = direction,
      factor_group = scores[[paste0("group_factor_", h)]],
      factor_score = scores[[paste0("factor_", h)]],
      set_accuracy = rowMeans(X[, top_items, drop = FALSE]),
      overall_accuracy = scores$accuracy,
      stringsAsFactors = FALSE
    )
    gp <- aggregate(cbind(set_accuracy, overall_accuracy, factor_score) ~ factor_group, tmp, mean)
    gp$n_models <- as.integer(table(tmp$factor_group)[as.character(gp$factor_group)])
    gp$factor <- paste0("F", h)
    gp$direction <- direction
    group_rows[[length(group_rows) + 1L]] <- gp
  }
}

factor_summary <- do.call(rbind, factor_rows)
top_examples <- do.call(rbind, example_rows)
group_summary <- do.call(rbind, group_rows)
profile_summary <- aggregate(accuracy ~ profile_id, scores, function(z) {
  c(n = length(z), mean = mean(z), min = min(z), max = max(z))
})
profile_summary <- data.frame(
  profile_id = profile_summary$profile_id,
  n_models = profile_summary$accuracy[, "n"],
  mean_accuracy = profile_summary$accuracy[, "mean"],
  min_accuracy = profile_summary$accuracy[, "min"],
  max_accuracy = profile_summary$accuracy[, "max"]
)

write.csv(factor_summary, file.path(out_dir, "openeval_factor_interpretation_summary.csv"), row.names = FALSE)
write.csv(top_examples, file.path(out_dir, "openeval_top_loading_item_examples.csv"), row.names = FALSE)
write.csv(group_summary, file.path(out_dir, "openeval_factor_group_performance_on_top_items.csv"), row.names = FALSE)
write.csv(profile_summary, file.path(out_dir, "openeval_profile_summary.csv"), row.names = FALSE)

binary_ll <- binary_probit_loglik_alpha(ref$X, ref$F_hat, ref$Lambda_hat, ref$alpha_hat)
fit_summary <- data.frame(
  H = H_fixed,
  G_config = G_report,
  lambda_l1_penalty = refinement_lambda_l1_penalty,
  n_models = nrow(X),
  n_items = ncol(X),
  binary_loglik_per_response = binary_ll / length(X),
  pretraining_iterations = ref$pretraining$n_completed,
  refinement_iterations = ref$joint_refinement$n_completed,
  refinement_selected_iteration = ref$joint_refinement$selected_iteration,
  elapsed_sec = elapsed
)
write.csv(fit_summary, file.path(out_dir, paste0(fit_label, "_fit_summary.csv")), row.names = FALSE)

plot_factor_scatter(scores, out_dir)
plot_factor_scores_heatmap(scores, ref$H, out_dir)
plot_group_heatmap(scores, ref$H, out_dir)
plot_top_loading_heatmap(load_meta, ref$H, out_dir)

cat("\nOpenEval H=", H_fixed, ", G=[", G_report, "] fit summary:\n", sep = "")
print(fit_summary)
cat("\nFactor interpretation summary:\n")
print(factor_summary)
cat("\nMixture group summary:\n")
print(summarize_mixture_profiles_ordered(ref))
cat("\nOutputs saved in: ", normalizePath(out_dir), "\n", sep = "")
