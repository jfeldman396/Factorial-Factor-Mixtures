#!/usr/bin/env Rscript

# Paper-ready heatmaps for the DGP loading matrices used by the sample-size
# simulation. The simulation redraws Lambda in each replication; these figures
# use the same generator with fixed figure seeds to show representative
# realizations of the two loading designs.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "../.."))
source(file.path(repo_root, "R", "sample_size_dgp.R"))

get_env <- function(name, default, FUN = identity) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  FUN(value)
}

parse_int_csv <- function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])

out_dir <- get_env(
  "OUT_DIR",
  file.path(
    "results",
    "full",
    "sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts"
  ),
  as.character
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p <- get_env("P", 500L, as.integer)
H_values <- get_env("H_VALUES", c(3L, 4L), parse_int_csv)
figure_seed <- get_env("FIGURE_SEED", 20260813L, as.integer)
block_size_mode <- get_env("BLOCK_SIZE_MODE", "balanced", as.character)

loading_designs <- c(
  sparse = "balanced_moderate_few_positive_cross",
  cross = "balanced_moderate_dense_signed_cross"
)
loading_designs <- normalize_sample_size_loading_design(loading_designs)
loading_titles <- c(
  sparse = 'Loadings = "Sparse"',
  cross = 'Loadings = "Cross"'
)
block_size_titles <- c(
  balanced = "balanced blocks",
  ifeval_like = "strongly unbalanced blocks",
  moderate_ifeval_like = "unbalanced blocks"
)

make_dgp_loadings <- function(design_name, p, H) {
  make_sample_size_loadings(
    design = design_name,
    p = p,
    H = H,
    block_size_mode = block_size_mode,
    loading_sign_mode = "block"
  )
}

make_design_lambda <- function(design_key, H, p, figure_seed) {
  set.seed(figure_seed + 1000L * H + match(design_key, names(loading_designs)))
  make_dgp_loadings(loading_designs[[design_key]], p = p, H = H)
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
    block_size_mode = block_size_mode,
    H = H,
    Lambda,
    check.names = FALSE
  )
  names(tab)[seq_len(H) + 5L] <- paste0("F", seq_len(H))
  write.csv(tab, out_file, row.names = FALSE)
}

out_files <- character(0)
for (H in H_values) {
  for (design_key in names(loading_designs)) {
    loading <- make_design_lambda(design_key, H = H, p = p, figure_seed = figure_seed)
    tag <- paste0(design_key, "_H", H)
    block_title <- if (block_size_mode %in% names(block_size_titles)) {
      block_size_titles[[block_size_mode]]
    } else {
      block_size_mode
    }
    title <- sprintf("%s | %s | H=%d DGP Lambda", loading_titles[[design_key]], block_title, H)
    png_file <- file.path(out_dir, paste0("dgp_lambda_heatmap_", tag, ".png"))
    csv_file <- file.path(out_dir, paste0("dgp_lambda_matrix_", tag, ".csv"))
    plot_lambda_heatmap(loading$Lambda, loading$block_id, title, png_file)
    write_lambda_table(loading$Lambda, loading$block_id, design_key, H, csv_file)
    out_files <- c(out_files, png_file, csv_file)
  }
}

cat("Wrote DGP loading heatmaps and matrices:\n")
cat(paste(normalizePath(out_files, mustWork = FALSE), collapse = "\n"), "\n")
