#!/usr/bin/env Rscript

# Canonical data-generating utilities for the sample-size simulations.
#
# Keep simulation structure here, not copied into individual scripts.  This
# avoids a dangerous failure mode: two scripts using the same label, such as
# "Cross", while silently drawing different loading matrices.

sample_size_loading_design_aliases <- function() {
  c(
    sparse = "balanced_moderate_few_positive_cross",
    few_positive_cross = "balanced_moderate_few_positive_cross",
    balanced_moderate_few_positive_cross = "balanced_moderate_few_positive_cross",
    cross = "balanced_moderate_dense_signed_cross",
    dense_signed_cross = "balanced_moderate_dense_signed_cross",
    balanced_moderate_dense_signed_cross = "balanced_moderate_dense_signed_cross",
    block_sparse = "block_sparse",
    neighbor_cross = "block_sparse_multisigned",
    block_sparse_multisigned = "block_sparse_multisigned",
    # Historical command-line shorthand.  In the reproducible sample-size study,
    # "Cross" means the dense signed cross-loading design below.  Use
    # "neighbor_cross" for the older two-neighbor pattern.
    block_cross = "balanced_moderate_dense_signed_cross"
  )
}

normalize_sample_size_loading_design <- function(design) {
  aliases <- sample_size_loading_design_aliases()
  out <- unname(aliases[as.character(design)])
  if (any(is.na(out))) {
    bad <- unique(as.character(design)[is.na(out)])
    stop(
      "Unknown loading design(s): ", paste(bad, collapse = ", "),
      ". Valid names/aliases are: ", paste(names(aliases), collapse = ", "),
      call. = FALSE
    )
  }
  out
}

sample_size_loading_design_label <- function(design) {
  design <- normalize_sample_size_loading_design(design)
  labels <- c(
    block_sparse = "Block sparse",
    block_sparse_multisigned = "Neighbor cross",
    balanced_moderate_few_positive_cross = "Sparse",
    balanced_moderate_dense_signed_cross = "Cross"
  )
  unname(labels[design])
}

balanced_block_sizes <- function(p, H) {
  block_sizes <- rep(floor(p / H), H)
  remainder <- p - sum(block_sizes)
  if (remainder > 0L) block_sizes[seq_len(remainder)] <- block_sizes[seq_len(remainder)] + 1L
  block_sizes
}

ifeval_like_block_proportions <- function(H) {
  # The H=3 and H=4 proportions match the original IFEval-inspired plots.  For
  # larger H, use a smooth decreasing profile instead of silently falling back
  # to balanced blocks.
  switch(
    as.character(H),
    "3" = c(425, 31, 44) / 500,
    "4" = c(415, 35, 18, 32) / 500,
    {
      raw <- 1 / (seq_len(H)^1.25)
      raw / sum(raw)
    }
  )
}

integer_block_sizes_from_proportions <- function(p, proportions) {
  proportions <- proportions / sum(proportions)
  raw_sizes <- floor(p * proportions)
  remainder <- p - sum(raw_sizes)
  if (remainder > 0L) {
    order_idx <- order(p * proportions - raw_sizes, decreasing = TRUE)
    raw_sizes[order_idx[seq_len(remainder)]] <- raw_sizes[order_idx[seq_len(remainder)]] + 1L
  }
  pmax(raw_sizes, 1L)
}

ifeval_like_block_sizes <- function(p, H) {
  integer_block_sizes_from_proportions(p, ifeval_like_block_proportions(H))
}

moderate_ifeval_like_block_sizes <- function(p, H, blend_to_balanced = 0.50) {
  ifeval_proportions <- ifeval_like_block_proportions(H)
  balanced_proportions <- rep(1 / H, H)
  proportions <- (1 - blend_to_balanced) * ifeval_proportions +
    blend_to_balanced * balanced_proportions
  integer_block_sizes_from_proportions(p, proportions)
}

make_sample_size_block_sizes <- function(p, H, mode = "balanced") {
  mode <- match.arg(mode, c("balanced", "ifeval_like", "moderate_ifeval_like"))
  if (mode == "ifeval_like") return(ifeval_like_block_sizes(p, H))
  if (mode == "moderate_ifeval_like") return(moderate_ifeval_like_block_sizes(p, H))
  balanced_block_sizes(p, H)
}

sample_size_block_id <- function(block_sizes) {
  rep(seq_along(block_sizes), times = block_sizes)
}

sample_size_block_sign_matrix <- function(H) {
  signs <- matrix(sample(c(-1, 1), H * H, replace = TRUE), H, H)
  diag(signs) <- sample(c(-1, 1), H, replace = TRUE)
  signs
}

neighboring_factors_global <- function(h, H, n_cross) {
  if (H <= 1L || n_cross == 0L) return(integer(0))
  offsets <- as.vector(rbind(seq_len(H - 1L), -seq_len(H - 1L)))
  out <- integer(0)
  for (offset in offsets) {
    candidate <- ((h - 1L + offset) %% H) + 1L
    if (!(candidate %in% out) && candidate != h) out <- c(out, candidate)
    if (length(out) == n_cross) break
  }
  out
}

make_sample_size_loadings <- function(
    design,
    p,
    H,
    block_size_mode = "balanced",
    loading_sign_mode = c("block", "smoke"),
    block_sizes = NULL,
    primary_loading_range = NULL,
    cross_loading_range = NULL,
    cross_loading_prob = NULL,
    cross_sign_mode = NULL) {
  design <- normalize_sample_size_loading_design(design)
  loading_sign_mode <- match.arg(loading_sign_mode)
  if (is.null(block_sizes)) block_sizes <- make_sample_size_block_sizes(p, H, block_size_mode)
  if (sum(block_sizes) != p) stop("block_sizes must sum to p.", call. = FALSE)

  block_id <- sample_size_block_id(block_sizes)
  signs <- sample_size_block_sign_matrix(H)
  Lambda <- matrix(0, p, H)
  default <- switch(
    design,
    block_sparse = list(primary = c(0.75, 1.25), cross = c(0.12, 0.28), prob = 0.035, sign = "block"),
    block_sparse_multisigned = list(primary = c(0.75, 1.25), cross = c(0.55, 0.95), prob = NA_real_, sign = "block"),
    balanced_moderate_few_positive_cross = list(primary = c(0.75, 1.25), cross = c(0.12, 0.28), prob = 0.035, sign = "positive"),
    balanced_moderate_dense_signed_cross = list(primary = c(0.75, 1.25), cross = c(0.20, 0.60), prob = 0.25, sign = "random")
  )
  primary_loading_range <- if (is.null(primary_loading_range)) default$primary else sort(as.numeric(primary_loading_range))
  cross_loading_range <- if (is.null(cross_loading_range)) default$cross else sort(as.numeric(cross_loading_range))
  cross_loading_prob <- if (is.null(cross_loading_prob)) default$prob else as.numeric(cross_loading_prob)
  cross_sign_mode <- if (is.null(cross_sign_mode)) default$sign else as.character(cross_sign_mode)
  if (length(primary_loading_range) != 2L || any(!is.finite(primary_loading_range)) ||
      any(primary_loading_range < 0)) {
    stop("primary_loading_range must contain two non-negative finite values.", call. = FALSE)
  }
  if (length(cross_loading_range) != 2L || any(!is.finite(cross_loading_range)) ||
      any(cross_loading_range < 0)) {
    stop("cross_loading_range must contain two non-negative finite values.", call. = FALSE)
  }
  if (is.finite(cross_loading_prob) &&
      (!is.finite(cross_loading_prob) || cross_loading_prob < 0 || cross_loading_prob > 1)) {
    stop("cross_loading_prob must be in [0, 1].", call. = FALSE)
  }
  if (!cross_sign_mode %in% c("positive", "negative", "random", "block")) {
    stop("cross_sign_mode must be one of positive, negative, random, or block.", call. = FALSE)
  }

  draw_cross_sign <- function(h, k) {
    if (cross_sign_mode == "positive") return(1)
    if (cross_sign_mode == "negative") return(-1)
    if (cross_sign_mode == "random") return(sample(c(-1, 1), 1L))
    if (loading_sign_mode == "smoke") sample(c(-1, 1), 1L) else signs[h, k]
  }

  if (design == "block_sparse") {
    for (j in seq_len(p)) {
      h <- block_id[j]
      primary_sign <- if (loading_sign_mode == "smoke") 1 else signs[h, h]
      Lambda[j, h] <- primary_sign * runif(1, primary_loading_range[1L], primary_loading_range[2L])
      for (k in setdiff(seq_len(H), h)) {
        if (runif(1) < cross_loading_prob) {
          Lambda[j, k] <- draw_cross_sign(h, k) * runif(1, cross_loading_range[1L], cross_loading_range[2L])
        }
      }
    }
  } else if (design == "block_sparse_multisigned") {
    n_cross <- min(2L, max(0L, H - 1L))
    for (h in seq_len(H)) {
      block_rows <- which(block_id == h)
      cross_factors <- neighboring_factors_global(h, H, n_cross)
      for (j in block_rows) {
        primary_sign <- if (loading_sign_mode == "smoke") 1 else signs[h, h]
        Lambda[j, h] <- primary_sign * runif(1, primary_loading_range[1L], primary_loading_range[2L])
        for (k in cross_factors) {
          Lambda[j, k] <- draw_cross_sign(h, k) * runif(1, cross_loading_range[1L], cross_loading_range[2L])
        }
      }
    }
  } else if (design == "balanced_moderate_few_positive_cross") {
    for (j in seq_len(p)) {
      h <- block_id[j]
      Lambda[j, h] <- runif(1, primary_loading_range[1L], primary_loading_range[2L])
      for (k in setdiff(seq_len(H), h)) {
        if (runif(1) < cross_loading_prob) {
          Lambda[j, k] <- draw_cross_sign(h, k) * runif(1, cross_loading_range[1L], cross_loading_range[2L])
        }
      }
    }
  } else if (design == "balanced_moderate_dense_signed_cross") {
    for (j in seq_len(p)) {
      h <- block_id[j]
      Lambda[j, h] <- runif(1, primary_loading_range[1L], primary_loading_range[2L])
      for (k in setdiff(seq_len(H), h)) {
        if (runif(1) < cross_loading_prob) {
          Lambda[j, k] <- draw_cross_sign(h, k) * runif(1, cross_loading_range[1L], cross_loading_range[2L])
        }
      }
    }
  } else {
    stop("Unsupported loading design: ", design, call. = FALSE)
  }

  attr(Lambda, "loading_design") <- design
  list(Lambda = Lambda, block_id = block_id, block_sizes = block_sizes, loading_design = design)
}

make_sample_size_item_intercepts <- function(
    p,
    H,
    block_id,
    seed,
    mode = "none",
    intercept_sd = 0.45,
    intercept_block_span = 1.6,
    intercept_clip = 1.75) {
  set.seed(seed)
  mode <- match.arg(mode, c("none", "ifeval_like", "random", "block", "viroli_smoke"))
  if (mode == "none") return(rep(0, p))
  if (mode == "viroli_smoke") {
    block_shift <- seq(-0.65, 0.65, length.out = H)
    alpha <- block_shift[block_id] + rnorm(p, sd = 0.20)
    return(pmax(pmin(alpha, 1.50), -1.50))
  }

  block_means <- if (H == 1L) {
    0
  } else {
    seq(intercept_block_span / 2, -intercept_block_span / 2, length.out = H)
  }
  alpha <- block_means[block_id]
  if (mode %in% c("ifeval_like", "random")) {
    alpha <- alpha + rnorm(p, 0, intercept_sd)
  }
  if (mode == "random") {
    alpha <- rnorm(p, 0, intercept_sd)
  }
  pmin(pmax(alpha, -intercept_clip), intercept_clip)
}
