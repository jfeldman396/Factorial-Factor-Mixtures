#!/usr/bin/env Rscript

# Paper-ready heatmaps for the DGP loading matrices used by the sample-size
# simulation. The simulation redraws Lambda in each replication; these figures
# use the same generator with fixed figure seeds to show representative
# realizations of the two loading designs.

options(stringsAsFactors = FALSE)

get_env <- function(name, default, FUN = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  FUN(value)
}

parse_int_csv <- function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])

out_dir <- get_env(
  "OUT_DIR",
  file.path(
    "..",
    "results",
    "moderate_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP_intercepts_centered"
  ),
  as.character
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p <- get_env("P", 500L, as.integer)
H_values <- get_env("H_VALUES", c(3L, 4L), parse_int_csv)
figure_seed <- get_env("FIGURE_SEED", 20260813L, as.integer)

loading_designs <- c(
  sparse = "balanced_moderate_few_positive_cross",
  cross = "balanced_moderate_dense_signed_cross"
)
loading_titles <- c(
  sparse = 'Loadings = "Sparse"',
  cross = 'Loadings = "Cross"'
)

balanced_block_sizes <- function(p, H) {
  block_sizes <- rep(floor(p / H), H)
  remainder <- p - sum(block_sizes)
  if (remainder > 0L) block_sizes[seq_len(remainder)] <- block_sizes[seq_len(remainder)] + 1L
  block_sizes
}

block_sign_matrix <- function(H) {
  signs <- matrix(sample(c(-1, 1), H * H, replace = TRUE), H, H)
  diag(signs) <- sample(c(-1, 1), H, replace = TRUE)
  signs
}

make_strong_same_sign_loadings <- function(design_name, p, H) {
  block_sizes <- balanced_block_sizes(p, H)
  block_id <- rep(seq_len(H), times = block_sizes)
  signs <- block_sign_matrix(H)
  Lambda <- matrix(0, p, H)

  if (design_name == "balanced_moderate_few_positive_cross") {
    for (j in seq_len(p)) {
      h <- block_id[j]
      Lambda[j, h] <- runif(1, 0.75, 1.25)
      for (k in setdiff(seq_len(H), h)) {
        if (runif(1) < 0.035) Lambda[j, k] <- runif(1, 0.12, 0.28)
      }
    }
  } else if (design_name == "balanced_moderate_dense_signed_cross") {
    for (j in seq_len(p)) {
      h <- block_id[j]
      Lambda[j, h] <- runif(1, 0.75, 1.25)
      for (k in setdiff(seq_len(H), h)) {
        if (runif(1) < 0.25) Lambda[j, k] <- sample(c(-1, 1), 1L) * runif(1, 0.20, 0.60)
      }
    }
  } else {
    stop("Unsupported loading design: ", design_name)
  }

  list(Lambda = Lambda, block_id = block_id, block_sizes = block_sizes)
}

make_design_lambda <- function(design_key, H, p, figure_seed) {
  set.seed(figure_seed + 1000L * H + match(design_key, names(loading_designs)))
  make_strong_same_sign_loadings(loading_designs[[design_key]], p = p, H = H)
}

plot_lambda_heatmap <- function(Lambda, block_id, title, out_file, zlim = c(-1.25, 1.25)) {
  p <- nrow(Lambda)
  H <- ncol(Lambda)
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
  breaks <- seq(zlim[1L], zlim[2L], length.out = length(pal) + 1L)
  z <- Lambda
  z[z < zlim[1L]] <- zlim[1L]
  z[z > zlim[2L]] <- zlim[2L]
  block_starts <- which(!duplicated(block_id))
  block_ends <- c(block_starts[-1L] - 1L, p)
  block_mids <- (block_starts + block_ends) / 2

  png(out_file, width = 920, height = 1240, res = 180)
  op <- par(mar = c(4.6, 5.2, 4.1, 7.1), xpd = FALSE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  image(
    x = seq_len(H),
    y = seq_len(p),
    z = t(z[p:1, , drop = FALSE]),
    col = pal,
    breaks = breaks,
    xaxt = "n",
    yaxt = "n",
    xlab = "factor",
    ylab = "item",
    main = title,
    useRaster = TRUE
  )
  axis(1, at = seq_len(H), labels = paste0("F", seq_len(H)), tick = FALSE)
  axis(2, at = p - block_mids + 1L, labels = paste0("block ", seq_len(H)), las = 1, tick = FALSE)
  abline(h = p - block_ends[-length(block_ends)] + 0.5, col = "#3A3A3A", lwd = 0.8)
  abline(v = seq_len(H) + 0.5, col = "#E8E8E8", lwd = 0.5)
  box()

  par(xpd = NA)
  legend_y <- seq(0.18, 0.82, length.out = length(pal))
  legend_x <- rep(H + 0.95, length(pal))
  points(legend_x, legend_y * p, pch = 15, col = pal, cex = 1.1)
  axis(
    side = 4,
    at = c(0.18, 0.34, 0.50, 0.66, 0.82) * p,
    labels = sprintf("%.2f", seq(zlim[1L], zlim[2L], length.out = 5L)),
    las = 1,
    line = 2.3,
    tick = FALSE
  )
  mtext("lambda", side = 4, line = 5.8)

  invisible(TRUE)
}

write_lambda_table <- function(Lambda, block_id, design_key, H, out_file) {
  tab <- data.frame(
    item = seq_len(nrow(Lambda)),
    block = block_id,
    loading_design = loading_titles[[design_key]],
    H = H,
    Lambda,
    check.names = FALSE
  )
  names(tab)[seq_len(H) + 4L] <- paste0("F", seq_len(H))
  write.csv(tab, out_file, row.names = FALSE)
}

out_files <- character(0)
for (H in H_values) {
  for (design_key in names(loading_designs)) {
    loading <- make_design_lambda(design_key, H = H, p = p, figure_seed = figure_seed)
    tag <- paste0(design_key, "_H", H)
    title <- sprintf("%s | H=%d DGP Lambda", loading_titles[[design_key]], H)
    png_file <- file.path(out_dir, paste0("dgp_lambda_heatmap_", tag, ".png"))
    csv_file <- file.path(out_dir, paste0("dgp_lambda_matrix_", tag, ".csv"))
    plot_lambda_heatmap(loading$Lambda, loading$block_id, title, png_file)
    write_lambda_table(loading$Lambda, loading$block_id, design_key, H, csv_file)
    out_files <- c(out_files, png_file, csv_file)
  }
}

cat("Wrote DGP loading heatmaps and matrices:\n")
cat(paste(normalizePath(out_files, mustWork = FALSE), collapse = "\n"), "\n")
