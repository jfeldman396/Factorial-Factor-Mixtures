#!/usr/bin/env Rscript

# Paper-ready figures for the fixed IFEval-like simulation DGP.

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)])
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(file_arg[1L])) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(repo_root, "R", "sample_size_dgp.R"))

figure_dir <- Sys.getenv(
  "FIGURE_DIR",
  file.path(
    repo_root, "results", "selected_plots", "sample_size",
    "paper_simulation_section"
  )
)
table_dir <- Sys.getenv(
  "TABLE_DIR",
  file.path(
    repo_root, "results", "selected_tables", "sample_size",
    "paper_simulation_section"
  )
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

loading_seed <- 20772531L
p_master <- 2000L
p_display <- 500L
H_values <- c(5L, 10L)
loading_design <- "balanced_moderate_dense_signed_cross"
block_size_mode <- "ifeval_min30"
primary_loading_range <- c(1, 2)
cross_loading_range <- c(1, 2)
cross_loading_prob <- 0.05

loading_design_normalized <- normalize_sample_size_loading_design(loading_design)
design_levels <- sort(unique(normalize_sample_size_loading_design(
  names(sample_size_loading_design_aliases())
)))
design_id <- match(loading_design_normalized, design_levels)
block_id_code <- match(
  block_size_mode,
  c("balanced", "ifeval_like", "moderate_ifeval_like", "ifeval_min30")
)

make_current_loading <- function(H) {
  set.seed(
    loading_seed + 1000003L * H + 1009L * design_id + 101L * block_id_code
  )
  master <- make_sample_size_loadings(
    design = loading_design,
    p = p_master,
    H = H,
    block_size_mode = block_size_mode,
    loading_sign_mode = "block",
    primary_loading_range = primary_loading_range,
    cross_loading_range = cross_loading_range,
    cross_loading_prob = cross_loading_prob,
    cross_sign_mode = "random"
  )
  subset_sample_size_loading_output(
    master,
    p = p_display,
    H = H,
    block_size_mode = block_size_mode
  )
}

loading_outputs <- setNames(lapply(H_values, make_current_loading), paste0("H", H_values))

for (H in H_values) {
  loading <- loading_outputs[[paste0("H", H)]]
  loading_table <- data.frame(
    item = seq_len(nrow(loading$Lambda)),
    primary_block = loading$block_id,
    loading$Lambda,
    check.names = FALSE
  )
  names(loading_table)[-(1:2)] <- paste0("F", seq_len(H))
  write.csv(
    loading_table,
    file.path(table_dir, sprintf("dgp_loadings_p500_H%d.csv", H)),
    row.names = FALSE
  )
}

heat_colors <- colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(401)
heat_breaks <- seq(-2, 2, length.out = length(heat_colors) + 1L)

draw_loading_panel <- function(loading, H) {
  Lambda <- loading$Lambda
  p <- nrow(Lambda)
  block_id <- loading$block_id
  block_starts <- which(!duplicated(block_id))
  block_ends <- c(block_starts[-1L] - 1L, p)
  block_mids <- (block_starts + block_ends) / 2

  par(
    mar = c(4.3, 5.6, 3.8, 1.0),
    cex.axis = 1.14,
    cex.lab = 1.22
  )
  image(
    x = seq_len(H),
    y = seq_len(p),
    z = t(Lambda[p:1L, , drop = FALSE]),
    col = heat_colors,
    breaks = heat_breaks,
    xlim = c(0.5, H + 0.5),
    xaxt = "n",
    yaxt = "n",
    xlab = "Factor",
    ylab = "",
    main = sprintf("H = %d", H),
    cex.main = 1.30,
    font.main = 2,
    useRaster = TRUE
  )
  axis(
    1,
    at = seq_len(H),
    labels = paste0("F", seq_len(H)),
    tick = FALSE,
    cex.axis = 1.14
  )
  axis(
    2,
    at = p - block_mids + 1L,
    labels = paste0("B", seq_len(H), " (", loading$block_sizes, ")"),
    las = 1,
    tick = FALSE,
    cex.axis = 0.98
  )
  if (length(block_ends) > 1L) {
    abline(
      h = p - block_ends[-length(block_ends)] + 0.5,
      col = "#4D4D4D",
      lwd = 0.7,
      lty = 3
    )
  }
  abline(v = seq_len(H) + 0.5, col = "#E5E5E5", lwd = 0.45)
  box(col = "#4D4D4D")
}

draw_loading_figure <- function() {
  layout(matrix(c(1, 2, 3), nrow = 1L), widths = c(1, 1, 0.16))
  par(oma = c(1.8, 2.6, 2.5, 0.5), family = "sans")
  draw_loading_panel(loading_outputs$H5, 5L)
  draw_loading_panel(loading_outputs$H10, 10L)

  par(mar = c(4.5, 0.3, 3.8, 3.3))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(-2, 2), xaxs = "i", yaxs = "i")
  color_edges <- seq(-2, 2, length.out = length(heat_colors) + 1L)
  rect(
    0,
    color_edges[-length(color_edges)],
    1,
    color_edges[-1L],
    col = heat_colors,
    border = NA
  )
  axis(4, at = seq(-2, 2, by = 1), las = 1, cex.axis = 1.08)
  mtext("Loading", side = 4, line = 2.8, cex = 1.12)
  box(col = "#4D4D4D")

  mtext(
    "Data Generating Loadings (p = 500)",
    side = 3,
    outer = TRUE,
    line = 0.75,
    font = 2,
    cex = 1.42
  )
  mtext(
    "Primary item block",
    side = 2,
    outer = TRUE,
    line = 0.55,
    cex = 1.18
  )
}

write_figure <- function(stem, width, height, draw) {
  png(
    file.path(figure_dir, paste0(stem, ".png")),
    width = width,
    height = height,
    units = "in",
    res = 300,
    type = "quartz"
  )
  draw()
  dev.off()

  pdf(
    file.path(figure_dir, paste0(stem, ".pdf")),
    width = width,
    height = height,
    family = "Helvetica",
    useDingbats = FALSE
  )
  draw()
  dev.off()
}

write_figure("dgp_loadings_p500", 10.6, 6.3, draw_loading_figure)

mixture_specs <- list(
  separated = list(
    label = "Separated symmetric",
    g2 = list(pi = c(0.50, 0.50), multiplier = 1.35, sd = c(0.45, 0.45)),
    g3 = list(pi = c(0.30, 0.40, 0.30), multiplier = 1.35, sd = c(0.45, 0.65, 0.45))
  ),
  asymmetric_pi = list(
    label = "Asymmetric weights",
    g2 = list(pi = c(0.65, 0.35), multiplier = 1.35, sd = c(0.45, 0.45)),
    g3 = list(pi = c(0.20, 0.50, 0.30), multiplier = 1.35, sd = c(0.45, 0.65, 0.45))
  ),
  moderate_overlap = list(
    label = "Moderate overlap",
    g2 = list(pi = c(0.50, 0.50), multiplier = 1.10, sd = c(0.60, 0.60)),
    g3 = list(pi = c(0.30, 0.40, 0.30), multiplier = 1.10, sd = c(0.60, 0.75, 0.60))
  )
)
separation <- 2

canonicalize_mixture <- function(pi, mu, sd) {
  marginal_mean <- sum(pi * mu)
  marginal_var <- sum(pi * (sd^2 + mu^2)) - marginal_mean^2
  marginal_sd <- sqrt(marginal_var)
  list(
    pi = pi,
    mu = (mu - marginal_mean) / marginal_sd,
    sd = sd / marginal_sd,
    raw_mu = mu,
    raw_sd = sd,
    raw_marginal_mean = marginal_mean,
    raw_marginal_sd = marginal_sd
  )
}

canonical_mixtures <- list()
parameter_rows <- list()
row_index <- 0L
for (scenario in names(mixture_specs)) {
  for (G in c(2L, 3L)) {
    key <- paste0("g", G)
    spec <- mixture_specs[[scenario]][[key]]
    raw_mu <- separation * spec$multiplier * if (G == 2L) c(-1, 1) else c(-1, 0, 1)
    canonical <- canonicalize_mixture(spec$pi, raw_mu, spec$sd)
    canonical_mixtures[[paste(scenario, G, sep = "_")]] <- canonical
    for (g in seq_len(G)) {
      row_index <- row_index + 1L
      parameter_rows[[row_index]] <- data.frame(
        scenario = scenario,
        scenario_label = mixture_specs[[scenario]]$label,
        G = G,
        component = g,
        weight = canonical$pi[g],
        mean_multiplier = spec$multiplier,
        raw_mean = canonical$raw_mu[g],
        raw_sd = canonical$raw_sd[g],
        canonical_mean = canonical$mu[g],
        canonical_sd = canonical$sd[g],
        raw_marginal_mean = canonical$raw_marginal_mean,
        raw_marginal_sd = canonical$raw_marginal_sd,
        stringsAsFactors = FALSE
      )
    }
  }
}
mixture_parameter_table <- do.call(rbind, parameter_rows)
write.csv(
  mixture_parameter_table,
  file.path(table_dir, "dgp_mixture_parameters_raw_and_canonical.csv"),
  row.names = FALSE
)

component_colors <- c("#0072B2", "#D55E00", "#009E73")
x_grid <- seq(-2.25, 2.25, length.out = 1000L)

panel_density_max <- function(mix) {
  components <- vapply(
    seq_along(mix$pi),
    function(g) mix$pi[g] * dnorm(x_grid, mix$mu[g], mix$sd[g]),
    numeric(length(x_grid))
  )
  max(rowSums(components), components)
}
y_max <- 1.06 * max(vapply(canonical_mixtures, panel_density_max, numeric(1L)))

draw_mixture_panel <- function(
    scenario,
    G,
    show_x,
    show_y,
    show_title = (G == 2L),
    panel_y_max = y_max) {
  mix <- canonical_mixtures[[paste(scenario, G, sep = "_")]]
  components <- vapply(
    seq_len(G),
    function(g) mix$pi[g] * dnorm(x_grid, mix$mu[g], mix$sd[g]),
    numeric(length(x_grid))
  )
  total <- rowSums(components)

  par(
    mar = c(if (show_x) 4.5 else 2.3, if (show_y) 4.5 else 2.4, 3.5, 1.1),
    cex.axis = 1.13,
    cex.lab = 1.20
  )
  plot(
    x_grid,
    total,
    type = "n",
    xlim = range(x_grid),
    ylim = c(0, panel_y_max),
    xaxt = if (show_x) "s" else "n",
    yaxt = if (show_y) "s" else "n",
    xlab = if (show_x) "Canonical factor score" else "",
    ylab = "",
    bty = "n"
  )
  abline(h = pretty(c(0, panel_y_max)), col = "#EBEBEB", lwd = 0.65)
  abline(v = 0, col = "#B8B8B8", lwd = 0.75, lty = 3)
  for (g in seq_len(G)) {
    polygon(
      c(x_grid, rev(x_grid)),
      c(components[, g], rep(0, length(x_grid))),
      col = adjustcolor(component_colors[g], alpha.f = 0.17),
      border = NA
    )
    lines(x_grid, components[, g], col = component_colors[g], lwd = 1.8)
  }
  lines(x_grid, total, col = "#202020", lwd = 2.3)
  box(col = "#666666")

  title(
    main = if (show_title) mixture_specs[[scenario]]$label else "",
    font.main = 2,
    cex.main = 1.20,
    line = 1.0
  )
  text(
    -2.12,
    0.93 * panel_y_max,
    labels = paste0("pi = (", paste(sprintf("%.2f", mix$pi), collapse = ", "), ")"),
    adj = c(0, 1),
    cex = 0.94,
    col = "#4D4D4D"
  )
  if (scenario == "separated") {
    mtext(paste0("G = ", G), side = 2, line = 3.0, font = 2, cex = 1.08)
  }
}

draw_mixture_figure <- function() {
  layout(
    matrix(c(1, 2, 3, 4, 5, 6, 7, 7, 7), nrow = 3L, byrow = TRUE),
    heights = c(1, 1, 0.17)
  )
  par(oma = c(0.4, 2.0, 2.8, 0.4), family = "sans")
  scenario_names <- names(mixture_specs)
  for (scenario in scenario_names) {
    draw_mixture_panel(
      scenario,
      G = 2L,
      show_x = FALSE,
      show_y = identical(scenario, scenario_names[1L])
    )
  }
  for (scenario in scenario_names) {
    draw_mixture_panel(
      scenario,
      G = 3L,
      show_x = TRUE,
      show_y = identical(scenario, scenario_names[1L])
    )
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = c("Component 1", "Component 2", "Component 3", "Marginal mixture"),
    col = c(component_colors, "#202020"),
    lwd = c(2, 2, 2, 2.5),
    horiz = TRUE,
    bty = "n",
    cex = 1.04,
    x.intersp = 0.8,
    seg.len = 2.5
  )

  mtext(
    "Data-generating factor marginals",
    side = 3,
    outer = TRUE,
    line = 0.95,
    font = 2,
    cex = 1.42
  )
  mtext("Weighted density", side = 2, outer = TRUE, line = 0.70, cex = 1.12)
}

write_figure("dgp_factor_marginals_G2_G3", 10.6, 7.2, draw_mixture_figure)

g3_mixtures <- canonical_mixtures[grepl("_3$", names(canonical_mixtures))]
y_max_g3 <- 1.08 * max(vapply(g3_mixtures, panel_density_max, numeric(1L)))

draw_mixture_g3_figure <- function() {
  layout(
    matrix(c(1, 2, 3, 4, 4, 4), nrow = 2L, byrow = TRUE),
    heights = c(1, 0.18)
  )
  par(oma = c(0.4, 2.5, 2.9, 0.4), family = "sans")
  scenario_names <- names(mixture_specs)
  for (scenario in scenario_names) {
    draw_mixture_panel(
      scenario,
      G = 3L,
      show_x = TRUE,
      show_y = identical(scenario, scenario_names[1L]),
      show_title = TRUE,
      panel_y_max = y_max_g3
    )
  }

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = c("Component 1", "Component 2", "Component 3", "Marginal mixture"),
    col = c(component_colors, "#202020"),
    lwd = c(2, 2, 2, 2.5),
    horiz = TRUE,
    bty = "n",
    cex = 1.04,
    x.intersp = 0.8,
    seg.len = 2.5
  )

  mtext(
    "Data-generating factor marginals (G = 3)",
    side = 3,
    outer = TRUE,
    line = 0.95,
    font = 2,
    cex = 1.42
  )
  mtext("Weighted density", side = 2, outer = TRUE, line = 0.70, cex = 1.12)
}

write_figure("dgp_factor_marginals_G3", 10.6, 4.8, draw_mixture_g3_figure)

cat("Wrote paper-ready DGP figures to:\n", normalizePath(figure_dir), "\n")
cat("Wrote plotted DGP parameters to:\n", normalizePath(table_dir), "\n")
