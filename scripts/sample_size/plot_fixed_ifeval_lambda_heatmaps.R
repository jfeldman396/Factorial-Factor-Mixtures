#!/usr/bin/env Rscript

# True loading heatmaps for the fixed-DGP IFEval-like simulation.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(repo_root, "R", "sample_size_dgp.R"))

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

split_csv <- function(x) {
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

parse_ints <- function(x) as.integer(split_csv(x))
parse_nums <- function(x) as.numeric(split_csv(x))

safe_token <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "value")
}

plot_lambda_heatmap <- function(Lambda, block_id, title, out_file, zlim) {
  p <- nrow(Lambda)
  H <- ncol(Lambda)
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
  z <- pmin(pmax(Lambda, zlim[1L]), zlim[2L])
  block_starts <- which(!duplicated(block_id))
  block_ends <- c(block_starts[-1L] - 1L, p)
  block_mids <- (block_starts + block_ends) / 2

  png(out_file, width = 1800, height = 1500, res = 190)
  op <- par(mar = c(5.0, 5.5, 5.0, 7.4), xpd = FALSE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  image(
    x = seq_len(H),
    y = seq_len(p),
    z = t(z[p:1, , drop = FALSE]),
    col = pal,
    breaks = seq(zlim[1L], zlim[2L], length.out = length(pal) + 1L),
    xlim = c(0.5, H + 0.5),
    xaxt = "n",
    yaxt = "n",
    xlab = "factor",
    ylab = "item block",
    main = title,
    cex.main = 0.85,
    useRaster = TRUE
  )
  axis(1, at = seq_len(H), labels = paste0("F", seq_len(H)), tick = FALSE)
  axis(2, at = p - block_mids + 1L, labels = paste0("block ", seq_len(H)), las = 1, tick = FALSE)
  abline(h = p - block_ends[-length(block_ends)] + 0.5, col = "#444444", lwd = 0.75, lty = 3)
  abline(v = seq_len(H) + 0.5, col = "#eeeeee", lwd = 0.5)
  box()

  par(xpd = NA)
  legend_y <- seq(0.16, 0.84, length.out = length(pal))
  points(rep(H + 0.95, length(pal)), legend_y * p, pch = 15, col = pal, cex = 1.05)
  axis(
    side = 4,
    at = c(0.16, 0.33, 0.50, 0.67, 0.84) * p,
    labels = sprintf("%.1f", seq(zlim[1L], zlim[2L], length.out = 5L)),
    las = 1,
    line = 2.3,
    tick = FALSE
  )
  mtext("loading", side = 4, line = 5.7)
  invisible(TRUE)
}

write_lambda_table <- function(loading, settings, out_file) {
  tab <- data.frame(
    item = seq_len(nrow(loading$Lambda)),
    block = loading$block_id,
    settings,
    loading$Lambda,
    check.names = FALSE
  )
  loading_cols <- seq_len(ncol(loading$Lambda)) + ncol(settings) + 2L
  names(tab)[loading_cols] <- paste0("F", seq_len(ncol(loading$Lambda)))
  write.csv(tab, out_file, row.names = FALSE)
}

seed_base <- as.integer(get_env("SEED", "20260731"))
p_values <- parse_ints(get_env("P_VALUES", "500,1000,1500,2000"))
p_master <- max(p_values)
H_values <- parse_ints(get_env("H_VALUES", "5,10"))
primary_range <- parse_nums(get_env("PRIMARY_LOADING_RANGE", "2,3"))
cross_range <- parse_nums(get_env("CROSS_LOADING_RANGE", "2,3"))
cross_prob <- as.numeric(get_env("CROSS_LOADING_PROB", "0.05"))
sep <- as.numeric(get_env("SEPARATION", "1"))
out_dir <- get_env(
  "OUT_DIR",
  file.path(repo_root, "results", "selected_plots", "sample_size", "fixed_ifeval_lambda_min30_u2_3_cp0_05_h5_h10", "true_lambda_heatmaps")
)
table_dir <- get_env(
  "TABLE_DIR",
  file.path(repo_root, "results", "selected_tables", "sample_size", "fixed_ifeval_lambda_min30_u2_3_cp0_05_h5_h10", "true_lambda")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

loading_design <- "balanced_moderate_dense_signed_cross"
loading_index <- 1L
block_index <- match("ifeval_min30", c("balanced", "ifeval_like", "moderate_ifeval_like", "ifeval_min30"))
strength_index <- match("strong", c("default", "weak", "strong"))
cross_prob_offset <- as.integer(round(10000 * cross_prob))

for (H in H_values) {
  loading_seed <- seed_base + 100000L * H +
    5000L * loading_index + 1000L * block_index +
    100L * strength_index + cross_prob_offset +
    as.integer(round(1000 * sep))
  set.seed(loading_seed)
  master <- make_sample_size_loadings(
    design = loading_design,
    p = p_master,
    H = H,
    block_size_mode = "ifeval_min30",
    loading_sign_mode = "block",
    primary_loading_range = primary_range,
    cross_loading_range = cross_range,
    cross_loading_prob = cross_prob,
    cross_sign_mode = "random"
  )
  for (p in p_values) {
    loading <- subset_sample_size_loading_output(
      master,
      p = p,
      H = H,
      block_size_mode = "ifeval_min30"
    )
    tag <- paste0("lambda_ifeval_min30_u2_3_cp", safe_token(cross_prob), "_H", H, "_p", p)
    settings <- data.frame(
      H = H,
      p = p,
      p_master = p_master,
      loading_seed = loading_seed,
      block_size_mode = "ifeval_min30",
      primary_loading_min = primary_range[1L],
      primary_loading_max = primary_range[2L],
      cross_loading_min = cross_range[1L],
      cross_loading_max = cross_range[2L],
      cross_loading_prob = cross_prob,
      stringsAsFactors = FALSE
    )
    title <- sprintf(
      "Fixed IFEval-like Lambda, min block size 30\nH=%d, p=%d, |nonzero loading| ~ Uniform(2,3), Pr(cross)=%.2f",
      H,
      p,
      cross_prob
    )
    plot_lambda_heatmap(loading$Lambda, loading$block_id, title, file.path(out_dir, paste0(tag, ".png")), c(-3, 3))
    write_lambda_table(loading, settings, file.path(table_dir, paste0(tag, ".csv")))
  }
}

cat("Wrote fixed-DGP Lambda heatmaps to:\n", normalizePath(out_dir, mustWork = FALSE), "\n")
cat("Wrote fixed-DGP Lambda matrices to:\n", normalizePath(table_dir, mustWork = FALSE), "\n")
