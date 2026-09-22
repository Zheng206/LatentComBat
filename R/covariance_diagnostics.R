#' Compute a Covariance or Correlation Matrix
#'
#' Computes either a covariance matrix or a correlation matrix from a numeric
#' data matrix. Observations are expected in rows and features in columns.
#'
#' @param x Numeric matrix or data frame with observations in rows and features
#'   in columns.
#' @param type Type of matrix to compute. Either `"correlation"` or
#'   `"covariance"`.
#'
#' @return A numeric square matrix containing the feature-wise correlation or
#'   covariance structure.
#'
#' @details
#' For `type = "correlation"`, correlations are computed using pairwise complete
#' observations. For `type = "covariance"`, covariances are computed using the
#' same missing-data rule.
#'
#' @keywords internal
compute_cov_matrix <- function(x, type = c("correlation", "covariance")) {

  type <- match.arg(type)
  x <- as.matrix(x)

  if (type == "correlation") {
    stats::cor(x, use = "pairwise.complete.obs")
  } else {
    stats::cov(x, use = "pairwise.complete.obs")
  }
}


#' Compute Pairwise Batch Covariance Distances
#'
#' Computes pairwise distances between batch-specific covariance or correlation
#' matrices.
#'
#' @param data Numeric matrix or data frame with observations in rows and
#'   features in columns.
#' @param bat Factor or vector specifying batch membership.
#' @param type Matrix type to compare, either `"correlation"` or `"covariance"`.
#' @param metric Distance metric. One of `"frobenius"`, `"mse"`,
#'   `"spectral"`, or `"eigen"`.
#' @param relative Logical indicating whether Frobenius or spectral distances
#'   are normalized by the average matrix norm of the two batches.
#'
#' @return A data frame containing one row per pair of batches with columns
#'   `batch1`, `batch2`, and `distance`.
#'
#' @export
batch_covariance_distance <- function(
    data,
    bat,
    type = c("correlation", "covariance"),
    metric = c("frobenius", "mse", "spectral", "eigen"),
    relative = TRUE
) {

  type <- match.arg(type)
  metric <- match.arg(metric)
  Y <- as.matrix(data)
  bat <- droplevels(as.factor(bat))

  if (nrow(Y) != length(bat)) {stop("`bat` must have one entry per observation.", call. = FALSE)}
  lev <- levels(bat)
  if (length(lev) < 2L) {stop("At least two batches are required.", call. = FALSE)}

  mats <- lapply(lev, function(b) compute_cov_matrix(Y[bat == b, , drop = FALSE], type = type))
  names(mats) <- lev
  pairs <- utils::combn(lev, 2, simplify = FALSE)

  out <- lapply(pairs, function(z) {
      A <- mats[[z[1]]]
      B <- mats[[z[2]]]
      d <- switch(metric, frobenius = frobenius_norm(A, B), mse = mse_cov(A, B), spectral = spectral_norm(A, B), eigen = eigen_error(A, B))

      if (isTRUE(relative) && metric %in% c("frobenius", "spectral")) {
        denom <- switch(metric,
          frobenius = 0.5 * (sqrt(sum(A^2)) + sqrt(sum(B^2))),
          spectral = 0.5 * (max(abs(eigen(A, symmetric = TRUE, only.values = TRUE)$values)) + max(abs(eigen(B, symmetric = TRUE, only.values = TRUE)$values)))
        )

        if (is.finite(denom) && denom > 0) {d <- d / denom}
      }

      data.frame(
        batch1 = z[1],
        batch2 = z[2],
        distance = d,
        stringsAsFactors = FALSE
      )
    }
  )

  dplyr::bind_rows(out)
}


#' Plot Batch Covariance Distance Before and After Harmonization
#'
#' Compares pairwise batch-specific covariance or correlation distances before
#' and after harmonization.
#'
#' @param raw_data Numeric matrix or data frame containing the original data.
#' @param harm_data Numeric matrix or data frame containing harmonized data.
#' @param bat Factor or vector specifying batch membership.
#' @param type Matrix type to compare, either `"correlation"` or `"covariance"`.
#' @param metric Distance metric. One of `"frobenius"`, `"mse"`,
#'   `"spectral"`, or `"eigen"`.
#' @param relative Logical indicating whether applicable distances are
#'   normalized by the average matrix norm of the corresponding batch pair.
#'
#' @return A ggplot object.
#'
#' @export
plot_batch_distance <- function(
    raw_data,
    harm_data,
    bat,
    type = c("correlation", "covariance"),
    metric = c("frobenius", "mse", "spectral", "eigen"),
    relative = TRUE
) {

  type <- match.arg(type)
  metric <- match.arg(metric)

  before <- batch_covariance_distance(
    data = raw_data,
    bat = bat,
    type = type,
    metric = metric,
    relative = relative
  )

  after <- batch_covariance_distance(
    data = harm_data,
    bat = bat,
    type = type,
    metric = metric,
    relative = relative
  )

  plot_df <- before |>
    dplyr::select(
      batch1 = .data$batch1,
      batch2 = .data$batch2,
      distance_before = .data$distance
    ) |>
    dplyr::left_join(
      after |>
        dplyr::select(
          batch1 = .data$batch1,
          batch2 = .data$batch2,
          distance_after = .data$distance
        ),
      by = c("batch1", "batch2")
    ) |>
    dplyr::mutate(
      pair = paste(.data$batch1, .data$batch2, sep = " vs "),
      change = .data$distance_after - .data$distance_before
    ) |>
    dplyr::arrange(.data$distance_before) |>
    dplyr::mutate(
      pair = factor(.data$pair, levels = .data$pair)
    )

  long_df <- plot_df |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(c(
        "distance_before",
        "distance_after"
      )),
      names_to = "stage",
      values_to = "distance"
    ) |>
    dplyr::mutate(
      stage = factor(
        .data$stage,
        levels = c(
          "distance_before",
          "distance_after"
        ),
        labels = c(
          "Before",
          "After"
        )
      )
    )

  point_cols <- c(
    "Before" = "#B8BDC7",
    "After" = "#4F86A6"
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = plot_df,
      ggplot2::aes(
        x = .data$distance_before,
        xend = .data$distance_after,
        y = .data$pair,
        yend = .data$pair
      ),
      linewidth = 1.1,
      color = "#D6DBE3"
    ) +
    ggplot2::geom_point(
      data = long_df,
      ggplot2::aes(
        x = .data$distance,
        y = .data$pair,
        fill = .data$stage
      ),
      shape = 21,
      size = 3.4,
      stroke = 0.8,
      color = "white"
    ) +
    ggplot2::scale_fill_manual(
      values = point_cols
    ) +
    ggplot2::labs(
      x = if (relative) {
        "Relative batch distance"
      } else {
        "Batch distance"
      },
      y = NULL,
      fill = NULL,
      subtitle = paste(
        tools::toTitleCase(type),
        "|",
        tools::toTitleCase(metric),
        "distance"
      )
    ) +
    ggplot2::theme_classic(
      base_size = 13
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 15
      ),
      plot.subtitle = ggplot2::element_text(
        size = 15,
        color = "grey35"
      ),
      plot.caption = ggplot2::element_text(
        size = 10,
        color = "grey40"
      ),
      axis.text.y = ggplot2::element_text(
        size = 10.5,
        color = "grey20",
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        size = 10.5,
        color = "grey20",
        face = "bold"
      ),
      axis.title.x = ggplot2::element_text(
        face = "bold",
        size = 12.5
      ),
      legend.position = "top",
      legend.text = ggplot2::element_text(
        face = "bold"
      ),
      panel.grid.major.x = ggplot2::element_line(
        color = "grey92",
        linewidth = 0.4
      ),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (relative) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = 1,
        linetype = "dashed",
        linewidth = 0.6,
        color = "grey55"
      )
  }

  p
}

#' Summarize Batch Covariance Homogeneity
#'
#' Summarizes pairwise batch-specific covariance or correlation distances.
#'
#' @inheritParams batch_covariance_distance
#'
#' @return A list containing pairwise distances and their mean and median.
#'
#' @export
summarize_batch_covariance <- function(
    data,
    bat,
    type = c("correlation", "covariance"),
    metric = c("frobenius", "mse", "spectral", "eigen"),
    relative = TRUE
) {

  d <- batch_covariance_distance(data = data, bat = bat, type = type, metric = metric, relative = relative)

  list(
    pairwise = d,
    mean_distance = mean(d$distance, na.rm = TRUE),
    median_distance = stats::median(d$distance, na.rm = TRUE)
  )
}

#' Evaluate Reduction in Batch Covariance Heterogeneity
#'
#' Compares average pairwise batch covariance or correlation distance before
#' and after harmonization.
#'
#' @inheritParams plot_batch_distance
#'
#' @return A named numeric vector containing the before and after distances and
#'   their relative reduction.
#'
#' @export
covariance_batch_reduction <- function(
    raw_data,
    harm_data,
    bat,
    type = c("correlation", "covariance"),
    metric = c("frobenius", "mse", "spectral", "eigen"),
    relative = TRUE
) {

  type <- match.arg(type)
  metric <- match.arg(metric)

  raw <- summarize_batch_covariance(raw_data, bat, type, metric, relative)
  harm <- summarize_batch_covariance(harm_data, bat, type, metric, relative)

  before <- raw$mean_distance
  after <- harm$mean_distance

  reduction <- if (is.finite(before) && before > 0) {
    1 - after / before
  } else {
    NA_real_
  }
  c(before = before, after = after, relative_reduction = reduction)
}


#' Evaluate Pooled Covariance Structure Perturbation
#'
#' Measures the difference between pooled covariance or correlation structure
#' before and after harmonization.
#'
#' @param raw_data Numeric matrix or data frame containing original data.
#' @param harm_data Numeric matrix or data frame containing harmonized data.
#' @param type Matrix type, either `"correlation"` or `"covariance"`.
#'
#' @return A list containing the raw and harmonized matrices and several
#'   discrepancy measures.
#'
#' @export
evaluate_structure_perturbation <- function(
    raw_data,
    harm_data,
    type = c("correlation", "covariance")
) {

  type <- match.arg(type)

  raw_data <- as.matrix(raw_data)
  harm_data <- as.matrix(harm_data)

  if (!all(dim(raw_data) == dim(harm_data))) {
    stop("`raw_data` and `harm_data` must have the same dimensions.", call. = FALSE)
  }

  raw_mat <- compute_cov_matrix(raw_data, type = type)
  harm_mat <- compute_cov_matrix(harm_data, type = type)

  frob <- frobenius_norm(raw_mat, harm_mat)
  denom <- sqrt(sum(raw_mat^2))

  relative_frobenius <- if (is.finite(denom) && denom > 0) {frob / denom} else {NA_real_}

  list(
    raw = raw_mat,
    harmonized = harm_mat,
    difference = harm_mat - raw_mat,
    Frobenius = frob,
    RelativeFrobenius = relative_frobenius,
    MSE = mse_cov(raw_mat, harm_mat),
    Spectral = spectral_norm(raw_mat, harm_mat),
    EigenError = eigen_error(raw_mat, harm_mat)
  )
}


#' Plot Pooled Structure Preservation
#'
#' Displays pooled covariance or correlation matrices before and after
#' harmonization together with their difference.
#'
#' @param raw_data Numeric matrix or data frame containing original data.
#' @param harm_data Numeric matrix or data frame containing harmonized data.
#' @param type Matrix type, either `"correlation"` or `"covariance"`.
#' @param labels Optional feature labels.
#'
#' @return A ggplot object.
#'
#' @export
plot_structure_preservation <- function(
    raw_data,
    harm_data,
    type = c("correlation", "covariance"),
    labels = NULL
) {

  type <- match.arg(type)
  result <- evaluate_structure_perturbation(raw_data, harm_data, type = type)
  p <- nrow(result$raw)

  if (is.null(labels)) {labels <- colnames(raw_data)}
  if (is.null(labels)) {labels <- paste0("V", seq_len(p))}

  make_long <- function(M, label) {
    data.frame(
      row = rep(seq_len(p), times = p),
      col = rep(seq_len(p), each = p),
      value = as.vector(M),
      matrix = label
    )
  }

  plot_df <- dplyr::bind_rows(
    make_long(result$raw, "Before"),
    make_long(result$harmonized, "After"),
    make_long(result$difference, "After - Before")
  )

  plot_df$matrix <- factor(plot_df$matrix, levels = c("Before", "After", "After - Before"))
  plot_df$row <- factor(plot_df$row, levels = rev(seq_len(p)), labels = rev(labels))
  plot_df$col <- factor(plot_df$col, levels = seq_len(p), labels = labels)

  if (type == "correlation") {fill_limit <- 1} else {fill_limit <- max(abs(c(result$raw, result$harmonized, result$difference)), na.rm = TRUE)}

  ggplot(plot_df, aes(x = .data[["col"]], y = .data[["row"]], fill = .data[["value"]])) +
    geom_tile() +
    facet_wrap(vars(.data[["matrix"]]), nrow = 1) +
    scale_fill_gradient2(midpoint = 0, limits = c(-fill_limit, fill_limit), oob = scales::squish) +
    coord_equal() +
    labs(
      x = NULL,
      y = NULL,
      fill = tools::toTitleCase(type),
      title = "Pooled Structure Preservation",
      subtitle = sprintf("Relative Frobenius perturbation: %.3f", result$RelativeFrobenius)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, margin = margin(b = 10)),
      legend.position = "right",
      plot.margin = margin(10, 15, 10, 10)
    )
}


#' Plot Batch-Specific Covariance Structure
#'
#' Displays covariance or correlation matrices separately by batch.
#'
#' @param data Numeric matrix or data frame with observations in rows and
#'   features in columns.
#' @param bat Factor or vector specifying batch membership.
#' @param type Matrix type, either `"correlation"` or `"covariance"`.
#' @param labels Optional feature labels.
#'
#' @return A ggplot object.
#'
#' @export
plot_batch_covariance <- function(
    data,
    bat,
    type = c("correlation", "covariance"),
    labels = NULL
) {

  type <- match.arg(type)
  Y <- as.matrix(data)
  bat <- droplevels(as.factor(bat))

  if (nrow(Y) != length(bat)) {
    stop("`bat` must have one entry per observation.", call. = FALSE)
  }

  p <- ncol(Y)

  if (is.null(labels)) {labels <- colnames(Y)}
  if (is.null(labels)) {labels <- paste0("V", seq_len(p))}

  plot_df <- lapply(levels(bat), function(b) {
      M <- compute_cov_matrix(Y[bat == b, , drop = FALSE], type = type)
      data.frame(
        row = rep(seq_len(p), times = p),
        col = rep(seq_len(p), each = p),
        value = as.vector(M),
        batch = b
      )
    }
  )

  plot_df <- dplyr::bind_rows(plot_df)

  plot_df$row <- factor(plot_df$row, levels = rev(seq_len(p)), labels = rev(labels))
  plot_df$col <- factor(plot_df$col, levels = seq_len(p), labels = labels)

  if (type == "correlation") {fill_limit <- 1} else {
    fill_limit <- max(abs(plot_df$value), na.rm = TRUE)
    if (!is.finite(fill_limit) || fill_limit == 0) {fill_limit <- 1}
  }

  ggplot(plot_df, aes(x = .data[["col"]], y = .data[["row"]], fill = .data[["value"]])) +
    geom_tile() +
    facet_wrap(vars(.data[["batch"]]), nrow = 1) +
    scale_fill_gradient2(midpoint = 0, limits = c(-fill_limit, fill_limit), oob = scales::squish) +
    coord_equal() +
    labs(x = NULL, y = NULL,
      fill = tools::toTitleCase(type),
      title = paste0("Batch-Specific ", tools::toTitleCase(type), " Structure"),
      subtitle = paste0("Feature-wise ", type, " matrices estimated separately within each batch")) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, margin = ggplot2::margin(b = 10)),
      legend.position = "right"
    )

}


#' Box's M test
#'
#' Tests equality of covariance matrices across groups using the
#' chi-square approximation with Box's small-sample correction.
#'
#' @param Y Numeric matrix n x p.
#' @param group Factor or coercible vector of group labels.
#' @param ridge Small diagonal ridge added to each group covariance matrix.
#'
#' @return A list with elements `statistic`, `parameter`, `p.value`,
#'   and `method`.
#'
#' @export
boxM_simple <- function(
    Y,
    group,
    ridge = 1e-8
) {

  Y <- as.matrix(Y)
  if (!is.factor(group)) {group <- factor(group)}
  ok <- stats::complete.cases(Y, group)
  Y <- Y[ok, , drop = FALSE]
  group <- droplevels(group[ok])
  p <- ncol(Y)
  k <- nlevels(group)
  ni <- as.integer(table(group))
  N <- sum(ni)

  if (k < 2L) {
    stop("At least two groups are required.", call. = FALSE)
  }

  if (any(ni <= 1L)) {
    stop("Each group must contain at least two observations.", call. = FALSE)
  }

  cov_i <- lapply(split(as.data.frame(Y), group), function(df) {stats::cov(df) + diag(ridge, p)})

  Sp <- Reduce(`+`, Map(function(S, n) {(n - 1) * S}, cov_i, ni)) / (N - k)
  logdet <- function(M) {as.numeric(determinant(M, logarithm = TRUE)$modulus)}
  Mstat <- (N - k) * logdet(Sp) - sum((ni - 1) * vapply(cov_i, logdet, numeric(1)))
  cfac <- ((2 * p^2 + 3 * p - 1) /(6 * (p + 1) * (k - 1))) * (sum(1 / (ni - 1)) - 1 / (N - k))

  X2 <- Mstat * (1 - cfac)

  df <- (k - 1) * p * (p + 1) / 2

  pval <- stats::pchisq(X2, df = df, lower.tail = FALSE)

  structure(
    list(
      statistic = unname(X2),
      parameter = unname(df),
      p.value = unname(pval),
      method = "Box's M test"
    )
  )
}

.prepare_batch_test_pca <- function(
    data,
    bat,
    n_pc = NULL,
    var_explained = 0.90
) {

  Y <- as.matrix(data)
  bat <- droplevels(as.factor(bat))

  if (nrow(Y) != length(bat)) {
    stop(
      "`bat` must have one entry per observation.",
      call. = FALSE
    )
  }

  if (nlevels(bat) < 2L) {
    stop(
      "At least two batches are required.",
      call. = FALSE
    )
  }

  if (anyNA(Y)) {
    stop(
      "PCA-based batch structure testing currently requires complete data.",
      call. = FALSE
    )
  }

  batch_n <- table(bat)
  min_batch_n <- min(batch_n)

  ## Keep dimension safely below the smallest batch size.
  max_pc <- min(
    ncol(Y),
    nrow(Y) - 1L,
    min_batch_n - 2L
  )

  if (max_pc < 1L) {
    stop(
      "Batch sizes are too small for multivariate structure testing.",
      call. = FALSE
    )
  }

  pca <- stats::prcomp(
    Y,
    center = TRUE,
    scale. = FALSE
  )

  prop_var <- pca$sdev^2 /
    sum(pca$sdev^2)

  cum_var <- cumsum(prop_var)

  if (is.null(n_pc)) {

    n_pc_target <- which(
      cum_var >= var_explained
    )[1L]

    if (is.na(n_pc_target)) {
      n_pc_target <- length(cum_var)
    }

    n_pc <- min(
      n_pc_target,
      max_pc
    )

  } else {

    n_pc <- as.integer(n_pc)

    if (
      length(n_pc) != 1L ||
      is.na(n_pc) ||
      n_pc < 1L
    ) {
      stop(
        "`n_pc` must be a positive integer.",
        call. = FALSE
      )
    }

    n_pc <- min(
      n_pc,
      max_pc
    )
  }

  loadings <- pca$rotation[
    ,
    seq_len(n_pc),
    drop = FALSE
  ]

  scores <- pca$x[
    ,
    seq_len(n_pc),
    drop = FALSE
  ]

  colnames(scores) <- paste0(
    "PC",
    seq_len(n_pc)
  )

  list(
    center = pca$center,
    loadings = loadings,
    scores = scores,
    n_pc = n_pc,
    variance_explained = cum_var[n_pc],
    max_pc = max_pc
  )
}


.manova_batch_test <- function(
    scores,
    bat
) {

  scores <- as.matrix(scores)
  bat <- droplevels(as.factor(bat))

  df <- data.frame(
    scores,
    bat = bat,
    check.names = FALSE
  )

  response <- paste0(
    "cbind(",
    paste(
      colnames(scores),
      collapse = ", "
    ),
    ")"
  )

  formula <- stats::as.formula(
    paste(
      response,
      "~ bat"
    )
  )

  fit <- stats::manova(
    formula,
    data = df
  )

  tab <- summary(
    fit,
    test = "Pillai"
  )$stats

  data.frame(
    statistic = unname(
      tab["bat", "Pillai"]
    ),
    approx_F = unname(
      tab["bat", "approx F"]
    ),
    df1 = unname(
      tab["bat", "num Df"]
    ),
    df2 = unname(
      tab["bat", "den Df"]
    ),
    p_value = unname(
      tab["bat", "Pr(>F)"]
    ),
    row.names = NULL
  )
}


#' Test Batch-Specific Mean and Covariance Structure
#'
#' Tests whether batches differ in multivariate mean structure using
#' MANOVA and in covariance structure using Box's M test.
#'
#' For high-dimensional data, tests are performed in a PCA representation
#' estimated from the input data.
#'
#' @param data Numeric matrix with observations in rows and features in columns.
#' @param bat Batch membership.
#' @param n_pc Number of principal components to use. If `NULL`, selected
#'   according to `var_explained`, subject to batch-size constraints.
#' @param var_explained Target cumulative variance explained when selecting PCs.
#' @param ridge Ridge passed to `boxM_simple()`.
#'
#' @return An object of class `"batch_structure_test"`.
#'
#' @export
test_batch_structure <- function(
    data,
    bat,
    n_pc = NULL,
    var_explained = 0.90,
    ridge = 1e-8
) {

  Y <- as.matrix(data)
  bat <- droplevels(as.factor(bat))

  pca_info <- .prepare_batch_test_pca(
    data = Y,
    bat = bat,
    n_pc = n_pc,
    var_explained = var_explained
  )

  Z <- pca_info$scores

  mean_test <- .manova_batch_test(
    scores = Z,
    bat = bat
  )

  covariance_fit <- boxM_simple(
    Y = Z,
    group = bat,
    ridge = ridge
  )

  covariance_test <- data.frame(
    statistic = covariance_fit$statistic,
    df = covariance_fit$parameter,
    p_value = covariance_fit$p.value,
    row.names = NULL
  )

  out <- list(
    n_pc = pca_info$n_pc,
    variance_explained = pca_info$variance_explained,
    mean = mean_test,
    covariance = covariance_test,
    pca = pca_info
  )

  class(out) <- "batch_structure_test"

  out
}


#' Compare Batch Structure Before and After Harmonization
#'
#' Compares residual batch-specific multivariate mean and covariance structure
#' before and after harmonization using a common PCA basis estimated from the
#' raw data.
#'
#' @param raw_data Original data matrix.
#' @param harm_data Harmonized data matrix.
#' @param bat Batch membership.
#' @param n_pc Number of PCs. If `NULL`, selected automatically.
#' @param var_explained Target cumulative variance explained.
#' @param ridge Ridge used in Box's M test.
#' @param alpha Significance level used for the convenience `significant`
#'   indicator.
#'
#' @return An object of class `"batch_structure_comparison"`.
#'
#' @export
compare_batch_structure <- function(
    raw_data,
    harm_data,
    bat,
    n_pc = NULL,
    var_explained = 0.90,
    ridge = 1e-8,
    alpha = 0.05
) {

  raw_data <- as.matrix(raw_data)
  harm_data <- as.matrix(harm_data)

  if (!all(dim(raw_data) == dim(harm_data))) {
    stop(
      "`raw_data` and `harm_data` must have the same dimensions.",
      call. = FALSE
    )
  }

  bat <- droplevels(as.factor(bat))

  if (nrow(raw_data) != length(bat)) {
    stop(
      "`bat` must have one entry per observation.",
      call. = FALSE
    )
  }

  if (anyNA(raw_data) || anyNA(harm_data)) {
    stop(
      "PCA-based batch structure comparison currently requires complete data.",
      call. = FALSE
    )
  }

  ## --------------------------------------------------
  ## Estimate PCA basis from RAW data only
  ## --------------------------------------------------

  pca_info <- .prepare_batch_test_pca(
    data = raw_data,
    bat = bat,
    n_pc = n_pc,
    var_explained = var_explained
  )

  V <- pca_info$loadings
  center <- pca_info$center

  raw_centered <- sweep(
    raw_data,
    2,
    center,
    FUN = "-"
  )

  harm_centered <- sweep(
    harm_data,
    2,
    center,
    FUN = "-"
  )

  Z_raw <- raw_centered %*% V
  Z_harm <- harm_centered %*% V

  pc_names <- paste0(
    "PC",
    seq_len(ncol(Z_raw))
  )

  colnames(Z_raw) <- pc_names
  colnames(Z_harm) <- pc_names

  ## --------------------------------------------------
  ## MANOVA
  ## --------------------------------------------------

  raw_mean <- .manova_batch_test(
    scores = Z_raw,
    bat = bat
  )

  harm_mean <- .manova_batch_test(
    scores = Z_harm,
    bat = bat
  )

  ## --------------------------------------------------
  ## Box's M
  ## --------------------------------------------------

  raw_cov <- boxM_simple(
    Y = Z_raw,
    group = bat,
    ridge = ridge
  )

  harm_cov <- boxM_simple(
    Y = Z_harm,
    group = bat,
    ridge = ridge
  )

  ## --------------------------------------------------
  ## Compact user-facing result
  ## --------------------------------------------------

  result <- dplyr::bind_rows(
    data.frame(
      structure = "Mean",
      test = "MANOVA (Pillai)",
      stage = "Before",
      statistic = raw_mean$statistic,
      p_value = raw_mean$p_value
    ),
    data.frame(
      structure = "Mean",
      test = "MANOVA (Pillai)",
      stage = "After",
      statistic = harm_mean$statistic,
      p_value = harm_mean$p_value
    ),
    data.frame(
      structure = "Covariance",
      test = "Box's M",
      stage = "Before",
      statistic = raw_cov$statistic,
      p_value = raw_cov$p.value
    ),
    data.frame(
      structure = "Covariance",
      test = "Box's M",
      stage = "After",
      statistic = harm_cov$statistic,
      p_value = harm_cov$p.value
    )
  )

  result$significant <- result$p_value < alpha

  out <- list(
    result = result,
    n_pc = pca_info$n_pc,
    variance_explained = pca_info$variance_explained,
    mean_before = raw_mean,
    mean_after = harm_mean,
    covariance_before = raw_cov,
    covariance_after = harm_cov,
    pca = pca_info,
    scores_raw = Z_raw,
    scores_harmonized = Z_harm
  )

  class(out) <- "batch_structure_comparison"

  out
}


#' @export
print.batch_structure_comparison <- function(
    x,
    ...
) {

  cat(
    "Batch Structure Comparison\n"
  )

  cat(
    "--------------------------\n"
  )

  cat(
    "Principal components:",
    x$n_pc,
    "\n"
  )

  cat(
    "Variance represented:",
    sprintf(
      "%.1f%%",
      100 * x$variance_explained
    ),
    "\n\n"
  )

  print(
    x$result,
    row.names = FALSE
  )

  invisible(x)
}


