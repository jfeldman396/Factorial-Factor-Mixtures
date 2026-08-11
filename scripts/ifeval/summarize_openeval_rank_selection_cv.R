#!/usr/bin/env Rscript

# Summarize OpenEval held-out cell-CV rank selection for:
#   1. independent-mixture probit factor model, G = 2
#   2. independent-mixture probit factor model, G = 3
#   3. ordinary Gaussian probit factor model

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ifeval_rank_selection_comparison")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_method_summary <- function(path, method, fold_path = NULL) {
  d <- read.csv(path, check.names = FALSE)
  d$method <- method
  if (!"mean_fit_seconds" %in% names(d) && !is.null(fold_path) && file.exists(fold_path)) {
    folds <- read.csv(fold_path, check.names = FALSE)
    fit_time <- aggregate(fit_seconds ~ H, data = folds, FUN = mean)
    names(fit_time)[names(fit_time) == "fit_seconds"] <- "mean_fit_seconds"
    d <- merge(d, fit_time, by = "H", all.x = TRUE)
  }
  keep <- c(
    "method",
    "H",
    "mean_heldout_loglik_per_response",
    "sd_heldout_loglik_per_response",
    "se_heldout_loglik_per_response",
    "best_H",
    "selected_by_one_se",
    "one_se_threshold",
    "mean_fit_seconds"
  )
  missing <- setdiff(keep, names(d))
  for (col in missing) d[[col]] <- NA
  d[, keep, drop = FALSE]
}

summaries <- rbind(
  read_method_summary(
    Sys.getenv(
      "G2_H_SUMMARY",
      unset = file.path(bundle_root, "results", "openeval_ifeval_cv_H1_10_G2", "mmlu_final_missing_aware_cv_H_summary.csv")
    ),
    "mixture_probit_G2",
    Sys.getenv(
      "G2_FOLD_SCORES",
      unset = file.path(bundle_root, "results", "openeval_ifeval_cv_H1_10_G2", "mmlu_final_missing_aware_cv_fold_scores.csv")
    )
  ),
  read_method_summary(
    Sys.getenv(
      "G3_H_SUMMARY",
      unset = file.path(bundle_root, "results", "openeval_ifeval_cv_H1_10_G3", "mmlu_final_missing_aware_cv_H_summary.csv")
    ),
    "mixture_probit_G3",
    Sys.getenv(
      "G3_FOLD_SCORES",
      unset = file.path(bundle_root, "results", "openeval_ifeval_cv_H1_10_G3", "mmlu_final_missing_aware_cv_fold_scores.csv")
    )
  ),
  read_method_summary(
    Sys.getenv(
      "ORDINARY_H_SUMMARY",
      unset = file.path(bundle_root, "results", "openeval_ifeval_cv_H1_10_ordinary", "ordinary_probit_H_summary.csv")
    ),
    "ordinary_probit_factor",
    Sys.getenv(
      "ORDINARY_FOLD_SCORES",
      unset = file.path(bundle_root, "results", "openeval_ifeval_cv_H1_10_ordinary", "ordinary_probit_factor_fold_scores.csv")
    )
  )
)

summaries <- summaries[order(summaries$method, summaries$H), ]
write.csv(summaries, file.path(out_dir, "openeval_ifeval_rank_selection_all_H.csv"), row.names = FALSE)

selected <- do.call(rbind, lapply(split(summaries, summaries$method), function(d) {
  best <- d[which.max(d$mean_heldout_loglik_per_response), ]
  one_se <- d[d$selected_by_one_se, ][1L, ]
  data.frame(
    method = best$method,
    best_H = best$H,
    best_mean_heldout_loglik_per_response = best$mean_heldout_loglik_per_response,
    best_se_heldout_loglik_per_response = best$se_heldout_loglik_per_response,
    one_se_H = one_se$H,
    one_se_threshold = best$one_se_threshold,
    mean_fit_seconds_at_best_H = best$mean_fit_seconds
  )
}))
selected <- selected[order(selected$method), ]
write.csv(selected, file.path(out_dir, "openeval_ifeval_selected_rank_by_method.csv"), row.names = FALSE)

cols <- c(
  mixture_probit_G2 = "#2B6CB0",
  mixture_probit_G3 = "#C05621",
  ordinary_probit_factor = "#4A5568"
)

png(file.path(out_dir, "openeval_ifeval_heldout_loglik_by_H_method.png"), width = 1500, height = 950, res = 160)
plot(
  range(summaries$H),
  range(
    summaries$mean_heldout_loglik_per_response - summaries$se_heldout_loglik_per_response,
    summaries$mean_heldout_loglik_per_response + summaries$se_heldout_loglik_per_response
  ),
  type = "n",
  xlab = "number of factors H",
  ylab = "held-out log likelihood per response",
  main = "OpenEval rank selection by held-out predictive log score"
)
for (m in names(cols)) {
  d <- summaries[summaries$method == m, ]
  lines(d$H, d$mean_heldout_loglik_per_response, type = "b", pch = 19, col = cols[m])
  arrows(
    d$H,
    d$mean_heldout_loglik_per_response - d$se_heldout_loglik_per_response,
    d$H,
    d$mean_heldout_loglik_per_response + d$se_heldout_loglik_per_response,
    angle = 90,
    code = 3,
    length = 0.035,
    col = cols[m]
  )
}
legend("bottomleft", legend = names(cols), col = cols, pch = 19, lty = 1, bty = "n")
dev.off()

png(file.path(out_dir, "openeval_ifeval_runtime_by_H_method.png"), width = 1400, height = 900, res = 160)
plot(
  range(summaries$H),
  range(summaries$mean_fit_seconds, na.rm = TRUE),
  type = "n",
  xlab = "number of factors H",
  ylab = "mean fit seconds per fold",
  main = "OpenEval CV runtime by H"
)
for (m in names(cols)) {
  d <- summaries[summaries$method == m, ]
  lines(d$H, d$mean_fit_seconds, type = "b", pch = 19, col = cols[m])
}
legend("topleft", legend = names(cols), col = cols, pch = 19, lty = 1, bty = "n")
dev.off()

accuracy_by_fold_path <- Sys.getenv(
  "H3_ACCURACY_BY_FOLD",
  unset = file.path(
    bundle_root,
    "results",
    "openeval_ifeval_rank_selection_comparison",
    "openeval_ifeval_H3_heldout_accuracy_by_fold.csv"
  )
)
if (file.exists(accuracy_by_fold_path)) {
  acc_fold <- read.csv(accuracy_by_fold_path, check.names = FALSE)
  write.csv(
    acc_fold,
    file.path(out_dir, "openeval_ifeval_H3_heldout_accuracy_by_fold.csv"),
    row.names = FALSE
  )

  acc_summary <- do.call(rbind, lapply(split(acc_fold, acc_fold$method), function(d) {
    data.frame(
      method = d$method[1L],
      mean_accuracy = mean(d$acc),
      se_accuracy = stats::sd(d$acc) / sqrt(nrow(d)),
      mean_brier = mean(d$brier),
      se_brier = stats::sd(d$brier) / sqrt(nrow(d)),
      mean_loglik = mean(d$loglik),
      se_loglik = stats::sd(d$loglik) / sqrt(nrow(d))
    )
  }))
  acc_summary <- acc_summary[order(acc_summary$method), ]
  write.csv(
    acc_summary,
    file.path(out_dir, "openeval_ifeval_H3_heldout_accuracy_summary.csv"),
    row.names = FALSE
  )

  png(file.path(out_dir, "openeval_ifeval_H3_predictive_metrics_by_method.png"), width = 1500, height = 650, res = 160)
  op <- par(mfrow = c(1, 3), mar = c(6, 5, 3, 1))
  methods <- acc_summary$method
  barplot(
    acc_summary$mean_accuracy,
    names.arg = methods,
    las = 2,
    ylim = range(c(acc_summary$mean_accuracy - acc_summary$se_accuracy, acc_summary$mean_accuracy + acc_summary$se_accuracy)),
    ylab = "held-out accuracy",
    main = "Accuracy"
  )
  arrows(
    seq_along(methods) * 1.2 - 0.5,
    acc_summary$mean_accuracy - acc_summary$se_accuracy,
    seq_along(methods) * 1.2 - 0.5,
    acc_summary$mean_accuracy + acc_summary$se_accuracy,
    angle = 90,
    code = 3,
    length = 0.035
  )
  barplot(
    acc_summary$mean_brier,
    names.arg = methods,
    las = 2,
    ylab = "held-out Brier score",
    main = "Brier"
  )
  barplot(
    acc_summary$mean_loglik,
    names.arg = methods,
    las = 2,
    ylim = range(c(acc_summary$mean_loglik - acc_summary$se_loglik, acc_summary$mean_loglik + acc_summary$se_loglik)),
    ylab = "held-out log likelihood per response",
    main = "Log likelihood"
  )
  arrows(
    seq_along(methods) * 1.2 - 0.5,
    acc_summary$mean_loglik - acc_summary$se_loglik,
    seq_along(methods) * 1.2 - 0.5,
    acc_summary$mean_loglik + acc_summary$se_loglik,
    angle = 90,
    code = 3,
    length = 0.035
  )
  par(op)
  dev.off()
}

cat("\nSelected ranks by held-out log score:\n")
print(selected)
cat("\nAll-H summaries saved in: ", normalizePath(out_dir), "\n", sep = "")
