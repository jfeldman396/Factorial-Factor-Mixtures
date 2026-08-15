#!/usr/bin/env Rscript

# Plot empirical factor-score marginals with the fitted one-dimensional
# Gaussian mixture components used for the independent mixture prior.

options(stringsAsFactors = FALSE)

fit_path <- Sys.getenv(
  "FIT_PATH",
  file.path(
    "results",
    "ifeval_columnwise_G_cv_20260812",
    "selected_mixture_H4_Gconfig3-3-3-2_lambda4_iter20_guarded_mix200",
    "openeval_H4_G3-3-3-2_fit.rds"
  )
)
out_dir <- Sys.getenv("OUT_DIR", dirname(fit_path))
out_file <- file.path(out_dir, "openeval_factor_marginal_mixtures.png")

fit <- readRDS(fit_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols <- c("#2b6cb0", "#c53030", "#2f855a", "#805ad5", "#b7791f")

png(out_file, width = 1800, height = 1300, res = 180)
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar), add = TRUE)

n_col <- min(2L, fit$H)
n_row <- ceiling(fit$H / n_col)
par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

for (h in seq_len(fit$H)) {
  y <- fit$F_hat[, h]
  mf <- fit$mixture_fits[[h]]
  cls <- fit$class_map[, h]
  ord <- order(mf$mu)
  pi_ord <- mf$pi[ord]
  mu_ord <- mf$mu[ord]
  var_ord <- mf$var[ord]
  score_sd <- stats::sd(y)
  if (!is.finite(score_sd) || score_sd <= 0) score_sd <- 1
  xg <- seq(min(y) - 0.5 * score_sd, max(y) + 0.5 * score_sd, length.out = 600)
  comp <- sapply(seq_along(pi_ord), function(g) {
    pi_ord[g] * stats::dnorm(xg, mean = mu_ord[g], sd = sqrt(var_ord[g]))
  })
  if (is.null(dim(comp))) comp <- matrix(comp, ncol = 1L)
  mix <- rowSums(comp)

  hist(
    y,
    breaks = 18,
    freq = FALSE,
    col = "#f3f4f6",
    border = "#d1d5db",
    main = paste0("F", h, " marginal mixture"),
    xlab = paste0("F", h, " score"),
    ylab = "density"
  )
  for (g in seq_along(pi_ord)) {
    lines(xg, comp[, g], col = cols[g], lwd = 2)
    abline(v = mu_ord[g], col = cols[g], lwd = 1.5, lty = 3)
    yy <- y[cls == g]
    if (length(yy) > 0L) {
      rug(yy, col = adjustcolor(cols[g], alpha.f = 0.7), ticksize = 0.045)
    }
  }

  legend(
    "topright",
    legend = paste0("cluster ", seq_along(pi_ord)),
    col = cols[seq_along(pi_ord)],
    lty = rep(1, length(pi_ord)),
    lwd = 2,
    bty = "n",
    cex = 0.75
  )
}

mtext("IFEval factor score marginals with fitted clusters", outer = TRUE, cex = 1.2, font = 2)

dev.off()
cat(normalizePath(out_file), "\n")
