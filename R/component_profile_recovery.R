# Observation-level recovery metrics for independent factor-mixture components.

component_label_permutations <- function(x) {
  x <- as.integer(x)
  if (length(x) <= 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    rest <- component_label_permutations(x[-i])
    cbind(x[i], rest)
  }))
}

component_adjusted_rand_index <- function(a, b) {
  keep <- is.finite(a) & is.finite(b)
  a <- a[keep]
  b <- b[keep]
  if (length(a) < 2L) return(NA_real_)
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  total_pairs <- choose2(sum(tab))
  if (total_pairs <= 0) return(NA_real_)
  observed <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  col_pairs <- sum(choose2(colSums(tab)))
  expected <- row_pairs * col_pairs / total_pairs
  upper <- 0.5 * (row_pairs + col_pairs)
  if (abs(upper - expected) < 1e-12) {
    return(if (abs(observed - upper) < 1e-12) 1 else 0)
  }
  (observed - expected) / (upper - expected)
}

align_component_labels <- function(true_label, estimated_label, G = NULL) {
  true_label <- as.integer(true_label)
  estimated_label <- as.integer(estimated_label)
  keep <- is.finite(true_label) & is.finite(estimated_label)
  if (is.null(G)) {
    G <- max(c(true_label[keep], estimated_label[keep]), na.rm = TRUE)
  }
  G <- as.integer(G)
  if (!length(G) || !is.finite(G) || G < 1L || !any(keep)) {
    return(list(mapped = rep(NA_integer_, length(true_label)), mapping = integer(0), accuracy = NA_real_))
  }
  if (any(true_label[keep] < 1L | true_label[keep] > G) ||
      any(estimated_label[keep] < 1L | estimated_label[keep] > G)) {
    return(list(mapped = rep(NA_integer_, length(true_label)), mapping = integer(0), accuracy = NA_real_))
  }

  permutations <- component_label_permutations(seq_len(G))
  scores <- apply(permutations, 1L, function(mapping) {
    mean(mapping[estimated_label[keep]] == true_label[keep])
  })
  best <- which(scores == max(scores, na.rm = TRUE))[1L]
  mapping <- as.integer(permutations[best, ])
  mapped <- rep(NA_integer_, length(true_label))
  mapped[keep] <- mapping[estimated_label[keep]]
  list(mapped = mapped, mapping = mapping, accuracy = scores[best])
}

component_profile_recovery <- function(
    true_component,
    estimated_responsibilities,
    estimated_factor_index = NULL,
    estimated_component = NULL,
    eps = 1e-12) {
  true_component <- as.matrix(true_component)
  n <- nrow(true_component)
  H <- ncol(true_component)
  if (is.null(estimated_factor_index)) estimated_factor_index <- seq_len(H)
  estimated_factor_index <- as.integer(estimated_factor_index)

  empty_result <- function() {
    data.frame(
      component_profile_hamming_accuracy = NA_real_,
      component_profile_hamming_median = NA_real_,
      component_profile_hamming_p10 = NA_real_,
      component_profile_exact_accuracy = NA_real_,
      mean_component_ari = NA_real_,
      min_component_ari = NA_real_,
      mean_true_component_probability = NA_real_,
      component_brier_score = NA_real_,
      component_log_loss = NA_real_,
      mean_component_entropy = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  has_responsibilities <- is.list(estimated_responsibilities) &&
    length(estimated_factor_index) >= H
  estimated_component <- if (is.null(estimated_component)) NULL else as.matrix(estimated_component)
  has_hard_components <- !is.null(estimated_component) &&
    nrow(estimated_component) == n &&
    ncol(estimated_component) >= max(estimated_factor_index, na.rm = TRUE)
  if (!has_responsibilities && !has_hard_components) {
    return(empty_result())
  }

  mapped_class <- matrix(NA_integer_, n, H)
  correctness <- matrix(NA_real_, n, H)
  true_probability <- matrix(NA_real_, n, H)
  brier <- matrix(NA_real_, n, H)
  log_loss <- matrix(NA_real_, n, H)
  entropy <- matrix(NA_real_, n, H)
  ari <- rep(NA_real_, H)

  for (h in seq_len(H)) {
    est_h <- estimated_factor_index[h]
    if (!is.finite(est_h) || est_h < 1L) next
    responsibility <- NULL
    if (has_responsibilities && est_h <= length(estimated_responsibilities)) {
      responsibility <- as.matrix(estimated_responsibilities[[est_h]])
      if (nrow(responsibility) != n || ncol(responsibility) < 1L) {
        responsibility <- NULL
      }
    }
    if (has_hard_components && est_h <= ncol(estimated_component)) {
      # Gibbs comparisons may provide the posterior modal joint profile while
      # also providing marginal allocation probabilities for proper scores.
      estimated_class <- as.integer(estimated_component[, est_h])
    } else if (!is.null(responsibility)) {
      responsibility[!is.finite(responsibility) | responsibility < 0] <- 0
      row_total <- rowSums(responsibility)
      valid_rows <- is.finite(row_total) & row_total > 0
      responsibility[valid_rows, ] <- responsibility[valid_rows, , drop = FALSE] / row_total[valid_rows]
      responsibility[!valid_rows, ] <- NA_real_
      estimated_class <- max.col(responsibility, ties.method = "first")
    } else {
      next
    }

    G_true <- max(true_component[, h], na.rm = TRUE)
    if (!is.finite(G_true)) next
    if (!is.null(responsibility) && ncol(responsibility) != G_true) {
      ari[h] <- component_adjusted_rand_index(true_component[, h], estimated_class)
      next
    }

    aligned <- align_component_labels(true_component[, h], estimated_class, G = G_true)
    if (!length(aligned$mapping)) next
    mapped_class[, h] <- aligned$mapped
    correctness[, h] <- as.numeric(mapped_class[, h] == true_component[, h])
    ari[h] <- component_adjusted_rand_index(true_component[, h], estimated_class)

    if (is.null(responsibility)) next

    responsibility_true_order <- matrix(0, n, G_true)
    for (estimated_g in seq_len(G_true)) {
      true_g <- aligned$mapping[estimated_g]
      responsibility_true_order[, true_g] <- responsibility[, estimated_g]
    }
    true_index <- cbind(seq_len(n), true_component[, h])
    p_true <- responsibility_true_order[true_index]
    true_probability[, h] <- p_true
    log_loss[, h] <- -log(pmax(p_true, eps))

    target <- matrix(0, n, G_true)
    target[true_index] <- 1
    brier[, h] <- rowSums((responsibility_true_order - target)^2)
    entropy_h <- -rowSums(responsibility_true_order * log(pmax(responsibility_true_order, eps)))
    entropy[, h] <- if (G_true > 1L) entropy_h / log(G_true) else 0
  }

  observation_accuracy <- rowMeans(correctness, na.rm = TRUE)
  observation_accuracy[rowSums(is.finite(correctness)) != H] <- NA_real_
  exact_profile <- rowSums(correctness, na.rm = TRUE) == H
  exact_profile[rowSums(is.finite(correctness)) != H] <- NA

  safe_mean <- function(x) if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
  safe_min <- function(x) if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_
  safe_quantile <- function(x, probability) {
    if (any(is.finite(x))) as.numeric(quantile(x, probability, na.rm = TRUE, names = FALSE)) else NA_real_
  }

  data.frame(
    component_profile_hamming_accuracy = safe_mean(observation_accuracy),
    component_profile_hamming_median = safe_quantile(observation_accuracy, 0.5),
    component_profile_hamming_p10 = safe_quantile(observation_accuracy, 0.1),
    component_profile_exact_accuracy = safe_mean(as.numeric(exact_profile)),
    mean_component_ari = safe_mean(ari),
    min_component_ari = safe_min(ari),
    mean_true_component_probability = safe_mean(true_probability),
    component_brier_score = safe_mean(brier),
    component_log_loss = safe_mean(log_loss),
    mean_component_entropy = safe_mean(entropy),
    stringsAsFactors = FALSE
  )
}
