#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE))
} else {
  getwd()
}
bundle_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

mixture_dir <- Sys.getenv(
  "MIXTURE_DIR",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ifeval_H3_G3_interpretation")
)
ordinary_dir <- Sys.getenv(
  "ORDINARY_DIR",
  unset = file.path(bundle_root, "results", "reproduced_openeval_ifeval_ordinary_probit_H3_visualization")
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = file.path(bundle_root, "results", "reproduced_ifeval_ordinary_vs_mixture_factor_visualization")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], all_permutations(x[-i]))
  }))
}

best_alignment <- function(ord_mat, mix_mat) {
  H <- ncol(mix_mat)
  C <- cor(ord_mat, mix_mat)
  perms <- all_permutations(seq_len(H))
  scores <- apply(perms, 1L, function(perm) sum(abs(C[cbind(perm, seq_len(H))])))
  perm <- perms[which.max(scores), ]
  signs <- sign(C[cbind(perm, seq_len(H))])
  signs[signs == 0] <- 1
  aligned <- sweep(ord_mat[, perm, drop = FALSE], 2L, signs, "*")
  colnames(aligned) <- paste0("ordinary_aligned_F", seq_len(H))
  list(perm = perm, signs = signs, aligned = aligned, cor = cor(aligned, mix_mat))
}

project3 <- function(mat) {
  data.frame(
    x = mat[, 1L] + 0.42 * mat[, 3L],
    y = mat[, 2L] + 0.24 * mat[, 3L]
  )
}

draw_panel <- function(mat, profiles, accuracy, title, xlab = "projected F1/F3", ylab = "projected F2/F3") {
  xy <- project3(mat)
  profile_levels <- sort(unique(as.character(profiles)))
  pal <- grDevices::hcl.colors(length(profile_levels), palette = "Dark 3")
  cols <- pal[match(as.character(profiles), profile_levels)]
  cex <- 0.7 + 1.7 * (accuracy - min(accuracy)) / pmax(diff(range(accuracy)), 1e-12)
  plot(
    xy$x, xy$y,
    pch = 19,
    col = cols,
    cex = cex,
    xlab = xlab,
    ylab = ylab,
    main = title
  )
  abline(h = 0, v = 0, col = "gray80", lty = 3)
}

mix <- read.csv(file.path(mixture_dir, "openeval_model_factor_scores_profiles.csv"), check.names = FALSE)
ordinary <- read.csv(file.path(ordinary_dir, "ordinary_probit_factor_scores.csv"), check.names = FALSE)
names(ordinary)[1L] <- "model_id"
names(ordinary)[3:5] <- paste0("ordinary_raw_F", 1:3)

dat <- merge(mix, ordinary, by = c("model_id", "accuracy"), all = FALSE)
mix_mat <- as.matrix(dat[, paste0("factor_", 1:3)])
ord_mat <- as.matrix(dat[, paste0("ordinary_raw_F", 1:3)])
aligned <- best_alignment(ord_mat, mix_mat)

for (h in seq_len(3L)) {
  dat[[paste0("ordinary_aligned_F", h)]] <- aligned$aligned[, h]
}

align_df <- data.frame(
  mixture_factor = paste0("F", seq_len(3L)),
  ordinary_factor_used = paste0("ordinary_raw_F", aligned$perm),
  sign_applied = aligned$signs,
  score_correlation = diag(aligned$cor)
)
write.csv(align_df, file.path(out_dir, "ifeval_ordinary_to_mixture_factor_alignment.csv"), row.names = FALSE)
write.csv(dat, file.path(out_dir, "ifeval_ordinary_aligned_factor_coordinates.csv"), row.names = FALSE)

png(file.path(out_dir, "ifeval_mixture_vs_ordinary_probit_3d_profiles.png"), width = 1800, height = 850, res = 160)
op <- par(mfrow = c(1, 2), mar = c(5, 5, 4, 1), oma = c(0, 0, 3, 0))
draw_panel(mix_mat, dat$profile_id, dat$accuracy, "Mixture factors")
draw_panel(aligned$aligned, dat$profile_id, dat$accuracy, "Ordinary probit factors, aligned")
mtext("IFEval factor spaces colored by mixture MAP profile", outer = TRUE, cex = 1.2, font = 2)
par(op)
dev.off()

cat("Saved base-R ordinary-vs-mixture comparison in: ", normalizePath(out_dir), "\n", sep = "")
print(align_df)
