#!/usr/bin/env Rscript

# Representative DGP loading heatmaps for the final signal-support simulation.
#
# The simulation redraws Lambda in each replication.  These plots use fixed
# figure seeds and the same generator so the design can be inspected directly.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "../.."))
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

format_num <- function(x) format(x, scientific = FALSE, trim = TRUE)

loading_ranges <- function(strength) {
  strength <- match.arg(strength, c("weak", "strong"))
  switch(
    strength,
    weak = c(1.25, 1.75),
    strong = c(2.50, 3.00)
  )
}

plot_lambda_heatmap <- function(Lambda, block_id, title, out_file, zlim) {
  p <- nrow(Lambda)
  H <- ncol(Lambda)
  pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
  breaks <- seq(zlim[1L], zlim[2L], length.out = length(pal) + 1L)
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
    breaks = breaks,
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

write_lambda_table <- function(Lambda, block_id, settings, out_file) {
  tab <- data.frame(
    item = seq_len(nrow(Lambda)),
    block = block_id,
    settings,
    Lambda,
    check.names = FALSE
  )
  loading_cols <- seq_len(ncol(Lambda)) + ncol(settings) + 2L
  names(tab)[loading_cols] <- paste0("F", seq_len(ncol(Lambda)))
  write.csv(tab, out_file, row.names = FALSE)
}

p <- as.integer(get_env("P", "500"))
H_values <- parse_ints(get_env("H_VALUES", "5,10,15,20"))
loading_strengths <- split_csv(get_env("LOADING_STRENGTHS", "weak,strong"))
cross_loading_probs <- parse_nums(get_env("CROSS_LOADING_PROBS", "0.075,0.2"))
block_modes <- split_csv(get_env("BLOCK_SIZE_MODES", "balanced,ifeval_like"))
figure_seed <- as.integer(get_env("FIGURE_SEED", "20260813"))
out_dir <- get_env(
  "OUT_DIR",
  file.path(repo_root, "results", "selected_plots", "sample_size", "signal_support_grid", "dgp_heatmaps")
)
table_dir <- get_env(
  "TABLE_DIR",
  file.path(repo_root, "results", "selected_tables", "sample_size", "signal_support_grid_dgp")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

out_files <- character(0)
for (block_mode in block_modes) {
  for (strength in loading_strengths) {
    range <- loading_ranges(strength)
    for (cross_prob in cross_loading_probs) {
      for (H in H_values) {
        set.seed(figure_seed + 100000L * H + 1000L * match(block_mode, block_modes) +
                   100L * match(strength, loading_strengths) +
                   as.integer(round(1000 * cross_prob)))
        loading <- make_sample_size_loadings(
          design = "balanced_moderate_dense_signed_cross",
          p = p,
          H = H,
          block_size_mode = block_mode,
          loading_sign_mode = "block",
          primary_loading_range = range,
          cross_loading_range = range,
          cross_loading_prob = cross_prob,
          cross_sign_mode = "random"
        )
        block_title <- if (block_mode == "ifeval_like") "unbalanced blocks" else "balanced blocks"
        title <- sprintf(
          "Cross/IFEval-like Lambda\n%s | %s loadings | Pr(cross)=%.3f | H=%d",
          block_title,
          strength,
          cross_prob,
          H
        )
        tag <- paste(
          "lambda",
          safe_token(block_mode),
          safe_token(strength),
          paste0("cp", safe_token(format_num(cross_prob))),
          paste0("H", H),
          sep = "_"
        )
        zlim <- c(-max(abs(range)), max(abs(range)))
        png_file <- file.path(out_dir, paste0(tag, ".png"))
        csv_file <- file.path(table_dir, paste0(tag, ".csv"))
        settings <- data.frame(
          block_size_mode = block_mode,
          loading_strength = strength,
          primary_loading_min = range[1L],
          primary_loading_max = range[2L],
          cross_loading_min = range[1L],
          cross_loading_max = range[2L],
          cross_loading_prob = cross_prob,
          cross_sign_mode = "random",
          H = H,
          p = p,
          stringsAsFactors = FALSE
        )
        plot_lambda_heatmap(loading$Lambda, loading$block_id, title, png_file, zlim = zlim)
        write_lambda_table(loading$Lambda, loading$block_id, settings, csv_file)
        out_files <- c(out_files, png_file, csv_file)
      }
    }
  }
}

cat("Wrote DGP loading heatmaps and matrices:\n")
cat(paste(normalizePath(out_files, mustWork = FALSE), collapse = "\n"), "\n")
