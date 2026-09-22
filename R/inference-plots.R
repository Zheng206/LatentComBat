#' Plot LatentComBat Inference Diagnostics
#'
#' Visualizes how uncertainty propagation changes point estimates, confidence
#' intervals, standard errors, p-values, and discoveries. Also supports
#' forest-style plots for a selected inference method.
#'
#' @param x An object returned by [harm_inference()].
#' @param type Diagnostic to draw. One of `"estimate"`, `"estimate_ci"`,
#'   `"forest"`, `"se"`, `"se_ratio"`, `"pvalue"`, `"significance"`,
#'   `"harmonization"`, or `"bootstrap_bias"`.
#' @param term Optional character vector selecting downstream terms.
#' @param feature Optional character vector selecting features.
#' @param inference Inference method used by `type = "forest"`. Either
#'   `"propagated"` for uncertainty-aware inference or `"reference"` for
#'   plug-in/reference inference.
#' @param order_by Ordering used by `type = "forest"`. One of `"estimate"`,
#'   `"abs_estimate"`, `"pvalue"`, `"feature"`, or `"none"`.
#' @param decreasing Logical indicating whether forest-plot ordering should be
#'   decreasing. Defaults to `TRUE`.
#' @param annotate_significant Logical indicating whether significant features
#'   should be marked in the forest plot.
#' @param adjusted Logical indicating whether significance annotations in the
#'   forest plot should use multiplicity-adjusted p-values. Defaults to `TRUE`.
#' @param adjust_method Multiplicity adjustment passed to [stats::p.adjust()].
#' @param alpha Significance level.
#' @param bins Number of histogram bins.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#'
#' @method plot latentcombat_inference
#' @export
plot.latentcombat_inference <- function(
    x,
    type = c("estimate", "estimate_ci", "forest", "se", "se_ratio", "pvalue",
             "significance", "harmonization", "bootstrap_bias"),
    term = NULL,
    feature = NULL,
    inference = c("propagated", "reference"),
    order_by = c("estimate", "abs_estimate", "pvalue", "feature", "none"),
    decreasing = TRUE,
    annotate_significant = TRUE,
    adjusted = TRUE,
    adjust_method = "BY",
    alpha = 0.05,
    bins = 30L,
    ...
) {
  .require_ggplot2()
  type <- match.arg(type)
  inference <- match.arg(inference)
  order_by <- match.arg(order_by)

  tab <- inference_table(x, adjust_method = adjust_method, alpha = alpha)
  reference_label <- attr(tab, "reference_label")
  tab <- .filter_inference_terms(tab, term)
  tab <- .filter_inference_features(tab, feature)

  bins <- as.integer(bins)
  if (length(bins) != 1L || is.na(bins) || bins < 1L) {
    stop("`bins` must be a positive integer.", call. = FALSE)
  }

  switch(
    type,
    estimate = .plot_inference_estimate(tab, reference_label),
    estimate_ci = .plot_inference_estimate_ci(tab, reference_label, alpha),
    forest = .plot_inference_forest(
      tab, reference_label, inference, order_by, decreasing,
      annotate_significant, adjusted, alpha
    ),
    se = .plot_inference_se(tab, reference_label),
    se_ratio = .plot_inference_se_ratio(tab, reference_label),
    pvalue = .plot_inference_pvalue(tab, reference_label),
    significance = .plot_inference_significance(tab, reference_label, adjusted, adjust_method, alpha),
    harmonization = .plot_harmonization_fraction(tab, bins),
    bootstrap_bias = .plot_bootstrap_bias(tab, bins)
  )
}

.plot_inference_estimate <- function(x, reference_label) {
  finite <- is.finite(x$reference_estimate) & is.finite(x$estimate)
  x_plot <- x[finite, , drop = FALSE]

  if (nrow(x_plot) == 0L) {
    stop("No finite point estimates available for plotting.", call. = FALSE)
  }

  lim <- range(
    c(x_plot$reference_estimate, x_plot$estimate),
    finite = TRUE
  )

  pad <- diff(lim) * 0.05
  if (!is.finite(pad) || pad == 0) {
    pad <- max(abs(lim), na.rm = TRUE) * 0.05
  }

  lim <- c(lim[1] - pad, lim[2] + pad)

  n_terms <- length(unique(as.character(x_plot$term)))
  one_term <- n_terms == 1L

  estimate_diff <- x_plot$estimate - x_plot$reference_estimate
  med_abs_diff <- .finite_median(abs(estimate_diff))

  diag_shade <- data.frame(
    x = c(lim[1], lim[1], lim[2]),
    y = c(lim[1], lim[2], lim[2])
  )

  p <- ggplot2::ggplot(
    x_plot,
    ggplot2::aes(
      x = .data$reference_estimate,
      y = .data$estimate
    )
  ) +
    ggplot2::geom_polygon(
      data = diag_shade,
      ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      fill = "#0072B2",
      alpha = 0.035
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      linewidth = 0.6,
      color = "grey50"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      color = "grey85"
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      color = "grey85"
    )

  if (one_term) {
    p <- p +
      ggplot2::geom_point(
        color = "#0072B2",
        fill = "#0072B2",
        shape = 21,
        size = 2.8,
        stroke = 0.4,
        alpha = 0.75,
        na.rm = TRUE
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(
          color = term,
          shape = term
        ),
        size = 2.6,
        stroke = 0.7,
        alpha = 0.75,
        na.rm = TRUE
      )
  }

  subtitle <- if (is.finite(med_abs_diff)) {
    paste0(
      "Points near the diagonal have similar estimates; ",
      "median absolute difference = ",
      format(round(med_abs_diff, 3), nsmall = 3)
    )
  } else {
    "Points near the diagonal have similar estimates"
  }

  p +
    ggplot2::coord_equal(
      xlim = lim,
      ylim = lim
    ) +
    ggplot2::labs(
      x = paste0(.title_case(reference_label), " estimate"),
      y = "Uncertainty-aware estimate",
      color = if (one_term) NULL else "Term",
      shape = if (one_term) NULL else "Term",
      title = if (one_term) {
        paste0(
          "Point-estimate comparison for ",
          unique(as.character(x_plot$term))
        )
      } else {
        "Downstream point estimates after uncertainty propagation"
      },
      subtitle = subtitle
    ) +
    .theme_inference() +
    ggplot2::theme(
      legend.position = if (one_term) "none" else "top",
      panel.grid = ggplot2::element_blank()
    )
}

.plot_inference_estimate_ci <- function(x, reference_label, alpha) {
  zcrit <- stats::qnorm(1 - alpha / 2)
  reference_low <- x$reference_estimate - zcrit * x$reference_se
  reference_high <- x$reference_estimate + zcrit * x$reference_se

  propagated <- data.frame(
    feature = as.character(x$feature),
    term = as.character(x$term),
    method = "Uncertainty-aware",
    estimate = x$estimate,
    conf_low = x$conf_low,
    conf_high = x$conf_high,
    stringsAsFactors = FALSE
  )

  reference <- data.frame(
    feature = as.character(x$feature),
    term = as.character(x$term),
    method = reference_label,
    estimate = x$reference_estimate,
    conf_low = reference_low,
    conf_high = reference_high,
    stringsAsFactors = FALSE
  )

  plot_data <- rbind(reference, propagated)
  plot_data$method <- factor(
    plot_data$method,
    levels = c(reference_label, "Uncertainty-aware")
  )

  feature_levels <- unique(as.character(x$feature))
  plot_data$feature <- factor(
    plot_data$feature,
    levels = rev(feature_levels)
  )

  dodge <- ggplot2::position_dodge(width = 0.55)
  method_cols <- stats::setNames(
    c("#6C757D", "#0072B2"),
    c(reference_label, "Uncertainty-aware")
  )
  method_shapes <- stats::setNames(
    c(1, 16),
    c(reference_label, "Uncertainty-aware")
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = feature,
      y = estimate,
      color = method,
      shape = method
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0, linetype = 2,
      linewidth = 0.5, color = "grey55"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = conf_low, ymax = conf_high),
      width = 0.12, linewidth = 0.75, alpha = 0.65,
      position = dodge, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = 2.8, stroke = 0.85,
      position = dodge, na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~term, scales = "free") +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(values = method_cols) +
    ggplot2::scale_shape_manual(values = method_shapes) +
    ggplot2::labs(
      x = NULL,
      y = "Estimated effect",
      color = NULL,
      shape = NULL,
      title = "Downstream effect estimates and confidence intervals",
      subtitle = paste0(
        .title_case(reference_label),
        " compared with uncertainty-aware inference"
      )
    ) +
    .theme_inference() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      panel.spacing = grid::unit(1, "lines")
    )
}

.plot_inference_forest <- function(
    x,
    reference_label,
    inference,
    order_by,
    decreasing,
    annotate_significant,
    adjusted,
    alpha
) {
  zcrit <- stats::qnorm(1 - alpha / 2)

  if (identical(inference, "reference")) {
    plot_data <- data.frame(
      feature = as.character(x$feature),
      term = as.character(x$term),
      estimate = x$reference_estimate,
      se = x$reference_se,
      conf_low = x$reference_estimate - zcrit * x$reference_se,
      conf_high = x$reference_estimate + zcrit * x$reference_se,
      p_value = x$reference_p_value,
      p_adjust = x$reference_p_adjust,
      significant = if (adjusted) {
        x$reference_significant_adjusted
      } else {
        x$reference_significant
      },
      stringsAsFactors = FALSE
    )
    method_label <- .title_case(reference_label)
  } else {
    plot_data <- data.frame(
      feature = as.character(x$feature),
      term = as.character(x$term),
      estimate = x$estimate,
      se = x$se,
      conf_low = x$conf_low,
      conf_high = x$conf_high,
      p_value = x$p_value,
      p_adjust = x$p_adjust,
      significant = if (adjusted) {
        x$significant_adjusted
      } else {
        x$significant
      },
      stringsAsFactors = FALSE
    )
    method_label <- "Uncertainty-aware"
  }

  plot_data$significant[is.na(plot_data$significant)] <- FALSE
  plot_data$significance <- ifelse(
    plot_data$significant,
    "Significant",
    "Not significant"
  )

  one_feature <- length(unique(plot_data$feature)) == 1L

  if (one_feature) {
    if (identical(order_by, "estimate")) {
      ord <- order(plot_data$estimate, decreasing = decreasing, na.last = TRUE)
    } else if (identical(order_by, "abs_estimate")) {
      ord <- order(abs(plot_data$estimate), decreasing = decreasing, na.last = TRUE)
    } else if (identical(order_by, "pvalue")) {
      p_ord <- ifelse(
        is.finite(plot_data$p_adjust),
        plot_data$p_adjust,
        plot_data$p_value
      )
      ord <- order(p_ord, decreasing = FALSE, na.last = TRUE)
    } else {
      ord <- seq_len(nrow(plot_data))
    }

    plot_data <- plot_data[ord, , drop = FALSE]
    plot_data$term <- factor(
      plot_data$term,
      levels = rev(as.character(plot_data$term))
    )

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = estimate,
        y = term,
        color = significance,
        shape = significance
      )
    ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = 2,
        linewidth = 0.5,
        color = "grey55"
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          xmin = conf_low,
          xmax = conf_high
        ),
        width = 0.16,
        orientation = "y",
        linewidth = 0.75,
        alpha = 0.75,
        na.rm = TRUE
      ) +
      ggplot2::geom_point(
        size = 3,
        stroke = 0.9,
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "Not significant" = "#A0A0A0",
          "Significant" = "#0072B2"
        )
      ) +
      ggplot2::scale_shape_manual(
        values = c(
          "Not significant" = 1,
          "Significant" = 16
        )
      ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.05, 0.12))
      ) +
      ggplot2::labs(
        x = "Estimated effect",
        y = NULL,
        color = NULL,
        shape = NULL,
        title = paste0(
          method_label,
          " effects for ",
          unique(plot_data$feature)
        ),
        subtitle = paste0(
          if (adjusted) "Adjusted" else "Nominal",
          " significance at alpha = ",
          alpha
        )
      ) +
      .theme_inference()

    if (isTRUE(annotate_significant)) {
      sig <- plot_data[plot_data$significant, , drop = FALSE]

      if (nrow(sig) > 0L) {
        p <- p +
          ggplot2::geom_text(
            data = sig,
            ggplot2::aes(
              x = conf_high,
              y = term,
              label = "*"
            ),
            inherit.aes = FALSE,
            hjust = -0.5,
            size = 4,
            fontface = "bold",
            na.rm = TRUE
          )
      }
    }

    return(p)
  }

  # Multiple-feature forest plot
  plot_data <- .order_forest_data(
    plot_data,
    order_by = order_by,
    decreasing = decreasing
  )

  plot_data$plot_label <- paste(
    plot_data$term,
    plot_data$feature,
    sep = "___"
  )

  plot_data$plot_label <- factor(
    plot_data$plot_label,
    levels = rev(unique(plot_data$plot_label))
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = estimate,
      y = plot_label,
      color = significance,
      shape = significance
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.5,
      color = "grey55"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = conf_low,
        xmax = conf_high
      ),
      width = 0.16,
      orientation = "y",
      linewidth = 0.75,
      alpha = 0.75,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = 2.8,
      stroke = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(
      ~term,
      scales = "free_y"
    ) +
    ggplot2::scale_y_discrete(
      labels = function(z) sub("^.*___", "", z)
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Not significant" = "#A0A0A0",
        "Significant" = "#0072B2"
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Not significant" = 1,
        "Significant" = 16
      )
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.12))
    ) +
    ggplot2::labs(
      x = "Estimated effect",
      y = NULL,
      color = NULL,
      shape = NULL,
      title = paste0(
        method_label,
        " effect estimates and confidence intervals"
      ),
      subtitle = paste0(
        if (adjusted) "Adjusted" else "Nominal",
        " significance at alpha = ",
        alpha
      )
    ) +
    .theme_inference()

  if (isTRUE(annotate_significant)) {
    sig <- plot_data[plot_data$significant, , drop = FALSE]

    if (nrow(sig) > 0L) {
      p <- p +
        ggplot2::geom_text(
          data = sig,
          ggplot2::aes(
            x = conf_high,
            y = plot_label,
            label = "*"
          ),
          inherit.aes = FALSE,
          hjust = -0.5,
          size = 4,
          fontface = "bold",
          na.rm = TRUE
        )
    }
  }

  p
}

.order_forest_data <- function(
    x,
    order_by = c("estimate", "abs_estimate", "pvalue", "feature", "none"),
    decreasing = TRUE
) {
  order_by <- match.arg(order_by)

  pieces <- split(
    x,
    x$term,
    drop = TRUE
  )

  pieces <- lapply(pieces, function(z) {
    if (identical(order_by, "estimate")) {
      ord <- order(z$estimate, decreasing = decreasing, na.last = TRUE)
    } else if (identical(order_by, "abs_estimate")) {
      ord <- order(abs(z$estimate), decreasing = decreasing, na.last = TRUE)
    } else if (identical(order_by, "pvalue")) {
      p <- ifelse(is.finite(z$p_adjust), z$p_adjust, z$p_value)
      ord <- order(
        p,
        decreasing = !decreasing,
        na.last = TRUE
      )
    } else if (identical(order_by, "feature")) {
      ord <- order(
        z$feature,
        decreasing = decreasing,
        na.last = TRUE
      )
    } else {
      ord <- seq_len(nrow(z))
    }

    z[ord, , drop = FALSE]
  })

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

.forest_order_label <- function(order_by) {
  switch(
    order_by,
    estimate = "effect estimate",
    abs_estimate = "absolute effect size",
    pvalue = "statistical significance",
    feature = "feature name",
    none = "input order"
  )
}

.plot_inference_se <- function(x, reference_label) {
  finite <- is.finite(x$reference_se) & is.finite(x$se)
  x_plot <- x[finite, , drop = FALSE]

  if (nrow(x_plot) == 0L) {
    stop("No finite standard errors available for plotting.", call. = FALSE)
  }

  lim <- range(
    c(x_plot$reference_se, x_plot$se),
    finite = TRUE
  )
  pad <- diff(lim) * 0.05
  if (!is.finite(pad) || pad == 0) pad <- max(lim, na.rm = TRUE) * 0.05
  lim <- c(max(0, lim[1] - pad), lim[2] + pad)

  med_ratio <- .finite_median(x_plot$se_ratio)
  n_terms <- length(unique(as.character(x_plot$term)))
  one_term <- n_terms == 1L

  diag_shade <- data.frame(
    x = c(lim[1], lim[1], lim[2]),
    y = c(lim[1], lim[2], lim[2])
  )

  p <- ggplot2::ggplot(
    x_plot,
    ggplot2::aes(
      x = .data$reference_se,
      y = .data$se
    )
  ) +
    ggplot2::geom_polygon(
      data = diag_shade,
      ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      fill = "#0072B2",
      alpha = 0.04
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      linewidth = 0.6,
      color = "grey50"
    )

  if (one_term) {
    p <- p +
      ggplot2::geom_point(
        color = "#0072B2",
        fill = "#0072B2",
        shape = 21,
        size = 2.8,
        stroke = 0.4,
        alpha = 0.72,
        na.rm = TRUE
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(
          color = term,
          shape = term
        ),
        size = 2.6,
        stroke = 0.7,
        alpha = 0.75,
        na.rm = TRUE
      )
  }

  subtitle <- if (is.finite(med_ratio)) {
    paste0(
      "Points above the diagonal have larger uncertainty-aware SEs; ",
      "median SE ratio = ",
      format(round(med_ratio, 2), nsmall = 2)
    )
  } else {
    "Points above the diagonal have larger uncertainty-aware SEs"
  }

  p +
    ggplot2::coord_equal(
      xlim = lim,
      ylim = lim
    ) +
    ggplot2::labs(
      x = paste0(.title_case(reference_label), " SE"),
      y = "Uncertainty-aware SE",
      color = if (one_term) NULL else "Term",
      shape = if (one_term) NULL else "Term",
      title = if (one_term) {
        paste0(
          "Standard-error propagation for ",
          unique(as.character(x_plot$term))
        )
      } else {
        "Standard errors after harmonization uncertainty propagation"
      },
      subtitle = subtitle
    ) +
    .theme_inference() +
    ggplot2::theme(
      legend.position = if (one_term) "none" else "top",
      panel.grid = ggplot2::element_blank()
    )
}

.plot_inference_se_ratio <- function(x, reference_label) {
  keep <- is.finite(x$se_ratio) & !is.na(x$term)
  x_plot <- x[keep, , drop = FALSE]

  if (nrow(x_plot) == 0L) {
    stop("No finite SE ratios available for plotting.", call. = FALSE)
  }

  med_tab <- stats::aggregate(
    se_ratio ~ term,
    data = x_plot,
    FUN = function(z) stats::median(z, na.rm = TRUE)
  )
  ord <- order(med_tab$se_ratio, decreasing = TRUE, na.last = TRUE)
  term_levels <- as.character(med_tab$term[ord])
  x_plot$term <- factor(as.character(x_plot$term), levels = term_levels)

  one_term <- length(unique(as.character(x_plot$term))) == 1L
  overall_med <- stats::median(x_plot$se_ratio, na.rm = TRUE)

  p <- ggplot2::ggplot(
    x_plot,
    ggplot2::aes(x = term, y = se_ratio)
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = 2,
      linewidth = 0.6,
      color = "grey50"
    )

  if (one_term) {
    p <- p +
      ggplot2::geom_violin(
        fill = "#0072B2",
        color = NA,
        alpha = 0.18,
        width = 0.9,
        na.rm = TRUE
      ) +
      ggplot2::geom_boxplot(
        width = 0.18,
        outlier.shape = NA,
        fill = "#0072B2",
        color = "#0072B2",
        alpha = 0.45,
        linewidth = 0.6,
        na.rm = TRUE
      ) +
      ggplot2::geom_jitter(
        width = 0.08,
        height = 0,
        shape = 16,
        size = 1.8,
        alpha = 0.35,
        color = "#0072B2",
        na.rm = TRUE
      )
  } else {
    p <- p +
      ggplot2::geom_violin(
        ggplot2::aes(fill = term),
        color = NA,
        alpha = 0.18,
        width = 0.9,
        na.rm = TRUE
      ) +
      ggplot2::geom_boxplot(
        ggplot2::aes(fill = term),
        width = 0.18,
        outlier.shape = NA,
        alpha = 0.45,
        linewidth = 0.6,
        na.rm = TRUE
      ) +
      ggplot2::geom_jitter(
        ggplot2::aes(color = term),
        width = 0.08,
        height = 0,
        shape = 16,
        size = 1.7,
        alpha = 0.30,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      ggplot2::guides(fill = "none")
  }

  p +
    ggplot2::stat_summary(
      fun = stats::median,
      geom = "point",
      shape = 23,
      size = 3,
      stroke = 0.8,
      fill = "white",
      color = "black"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "SE ratio",
      title = if (one_term) {
        paste0("SE inflation for ", unique(as.character(x_plot$term)))
      } else {
        "Distribution of SE inflation after uncertainty propagation"
      },
      subtitle = paste0(
        "SE ratio = uncertainty-aware SE / ",
        reference_label,
        " SE; dashed line at 1, overall median = ",
        format(round(overall_med, 2), nsmall = 2)
      )
    ) +
    .theme_inference() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(face = "bold")
    )
}

.plot_inference_pvalue <- function(x, reference_label) {
  keep <- is.finite(x$reference_p_value) &
    is.finite(x$p_value) &
    x$reference_p_value > 0 &
    x$p_value > 0

  x_plot <- x[keep, , drop = FALSE]

  if (nrow(x_plot) == 0L) {
    stop("No finite p-values available for plotting.", call. = FALSE)
  }

  x_plot$reference_logp <- -log10(x_plot$reference_p_value)
  x_plot$propagated_logp <- -log10(x_plot$p_value)

  n_terms <- length(unique(as.character(x_plot$term)))
  one_term <- n_terms == 1L

  lim <- range(
    c(x_plot$reference_logp, x_plot$propagated_logp),
    finite = TRUE
  )

  pad <- diff(lim) * 0.05
  if (!is.finite(pad) || pad == 0) {
    pad <- max(lim, na.rm = TRUE) * 0.05
  }

  lim <- c(
    max(0, lim[1] - pad),
    lim[2] + pad
  )

  # Region below diagonal:
  # propagated -log10(p) < reference -log10(p)
  # => propagated p-value is larger / less significant
  less_sig_region <- data.frame(
    x = c(lim[1], lim[2], lim[2]),
    y = c(lim[1], lim[1], lim[2])
  )

  delta_logp <- x_plot$propagated_logp - x_plot$reference_logp

  n_less <- sum(delta_logp < 0, na.rm = TRUE)
  n_more <- sum(delta_logp > 0, na.rm = TRUE)

  p <- ggplot2::ggplot(
    x_plot,
    ggplot2::aes(
      x = .data$reference_logp,
      y = .data$propagated_logp
    )
  ) +
    ggplot2::geom_polygon(
      data = less_sig_region,
      ggplot2::aes(x = .data$x, y = .data$y),
      inherit.aes = FALSE,
      fill = "#0072B2",
      alpha = 0.035
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      linewidth = 0.6,
      color = "grey50"
    )

  if (one_term) {
    p <- p +
      ggplot2::geom_point(
        color = "#0072B2",
        fill = "#0072B2",
        shape = 21,
        size = 2.8,
        stroke = 0.4,
        alpha = 0.75,
        na.rm = TRUE
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(
          color = term,
          shape = term
        ),
        size = 2.6,
        stroke = 0.7,
        alpha = 0.75,
        na.rm = TRUE
      )
  }

  subtitle <- paste0(
    "Below diagonal = weaker evidence after uncertainty propagation; ",
    n_less, " weaker, ",
    n_more, " stronger"
  )

  p +
    ggplot2::coord_equal(
      xlim = lim,
      ylim = lim
    ) +
    ggplot2::labs(
      x = paste0(
        .title_case(reference_label),
        " -log10(p)"
      ),
      y = "Uncertainty-aware -log10(p)",
      color = if (one_term) NULL else "Term",
      shape = if (one_term) NULL else "Term",
      title = if (one_term) {
        paste0(
          "P-value propagation for ",
          unique(as.character(x_plot$term))
        )
      } else {
        "P-values after harmonization uncertainty propagation"
      },
      subtitle = subtitle
    ) +
    .theme_inference() +
    ggplot2::theme(
      legend.position = if (one_term) "none" else "top",
      panel.grid = ggplot2::element_blank()
    )
}

.plot_inference_significance <- function(
    x,
    reference_label,
    adjusted,
    adjust_method,
    alpha
) {
  if (adjusted) {
    ref_sig <- x$reference_significant_adjusted
    prop_sig <- x$significant_adjusted
  } else {
    ref_sig <- x$reference_significant
    prop_sig <- x$significant
  }

  transition <- ifelse(
    ref_sig & prop_sig,
    "Retained",
    ifelse(
      ref_sig & !prop_sig,
      "Lost",
      ifelse(
        !ref_sig & prop_sig,
        "Gained",
        "Not significant"
      )
    )
  )

  plot_data <- data.frame(
    term = x$term,
    transition = transition,
    stringsAsFactors = FALSE
  )

  counts <- as.data.frame(
    table(
      term = plot_data$term,
      transition = plot_data$transition
    ),
    stringsAsFactors = FALSE
  )

  names(counts)[names(counts) == "Freq"] <- "n"
  counts <- counts[counts$n > 0L, , drop = FALSE]

  transition_levels <- c(
    "Retained",
    "Lost",
    "Gained",
    "Not significant"
  )

  counts$transition <- factor(
    counts$transition,
    levels = transition_levels
  )

  transition_cols <- c(
    "Retained" = "#4C78A8",
    "Lost" = "#D98C5F",
    "Gained" = "#72A0C1",
    "Not significant" = "#D3D3D3"
  )

  ggplot2::ggplot(
    counts,
    ggplot2::aes(
      x = term,
      y = n,
      fill = transition
    )
  ) +
    ggplot2::geom_col(
      width = 0.7,
      color = "white",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = n
      ),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white",
      fontface = "bold",
      size = 3.5
    ) +
    ggplot2::scale_fill_manual(
      values = transition_cols,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of features",
      fill = "Transition",
      title = "Discovery changes after uncertainty propagation",
      subtitle = if (adjusted) {
        paste0(
          adjust_method,
          "-adjusted significance at alpha = ",
          alpha,
          ": ",
          reference_label,
          " vs uncertainty-aware inference"
        )
      } else {
        paste0(
          "Nominal significance at alpha = ",
          alpha,
          ": ",
          reference_label,
          " vs uncertainty-aware inference"
        )
      }
    ) +
    .theme_inference()
}

.plot_harmonization_fraction <- function(x, bins = 30L) {
  keep <- is.finite(x$harmonization_fraction)
  x_plot <- x[keep, , drop = FALSE]

  if (nrow(x_plot) == 0L) {
    stop(
      "No finite harmonization fractions are available for plotting.",
      call. = FALSE
    )
  }

  one_term <- length(unique(as.character(x_plot$term))) == 1L

  med_frac <- stats::median(
    x_plot$harmonization_fraction,
    na.rm = TRUE
  )

  p <- ggplot2::ggplot(
    x_plot,
    ggplot2::aes(x = .data$harmonization_fraction)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins,
      boundary = 0,
      fill = "#0072B2",
      color = "white",
      linewidth = 0.25,
      alpha = 0.65,
      na.rm = TRUE
    ) +
    ggplot2::geom_density(
      color = "#0072B2",
      linewidth = 0.8,
      alpha = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = med_frac,
      linetype = 2,
      linewidth = 0.7,
      color = "grey35"
    ) +
    ggplot2::annotate(
      "text",
      x = med_frac,
      y = Inf,
      label = paste0(
        "Median = ",
        format(round(med_frac, 2), nsmall = 2)
      ),
      vjust = 1.6,
      hjust = if (med_frac < 0.5) -0.1 else 1.1,
      size = 3.5,
      color = "grey30"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      x = "Fraction of inferential variance due to harmonization",
      y = "Density",
      title = if (one_term) {
        paste0(
          "Harmonization uncertainty for ",
          unique(as.character(x_plot$term))
        )
      } else {
        "Contribution of harmonization uncertainty"
      },
      subtitle = paste0(
        "Higher values indicate a larger share of total uncertainty ",
        "arising from harmonization"
      )
    ) +
    .theme_inference() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )

  if (!one_term) {
    p <- p +
      ggplot2::facet_wrap(
        ~term,
        scales = "free_y"
      )
  }

  p
}

.plot_bootstrap_bias <- function(x, bins = 30L) {
  keep <- is.finite(x$bootstrap_bias)
  x_plot <- x[keep, , drop = FALSE]

  if (nrow(x_plot) == 0L) {
    stop(
      "No finite bootstrap bias values are available for plotting.",
      call. = FALSE
    )
  }

  one_term <- length(unique(as.character(x_plot$term))) == 1L

  med_bias <- stats::median(
    x_plot$bootstrap_bias,
    na.rm = TRUE
  )

  p <- ggplot2::ggplot(
    x_plot,
    ggplot2::aes(x = .data$bootstrap_bias)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins,
      fill = "#0072B2",
      color = "white",
      linewidth = 0.25,
      alpha = 0.60,
      na.rm = TRUE
    ) +
    ggplot2::geom_density(
      color = "#0072B2",
      linewidth = 0.8,
      alpha = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.7,
      color = "grey45"
    ) +
    ggplot2::geom_vline(
      xintercept = med_bias,
      linetype = 2,
      linewidth = 0.7,
      color = "#D98C5F"
    ) +
    ggplot2::labs(
      x = "Bootstrap bias",
      y = "Density",
      title = if (one_term) {
        paste0(
          "Bootstrap bias for ",
          unique(as.character(x_plot$term))
        )
      } else {
        "Bootstrap bias after sequential harmonization"
      },
      subtitle = paste0(
        "Solid line indicates zero bias; dashed line indicates median bias = ",
        format(round(med_bias, 3), nsmall = 3)
      )
    ) +
    .theme_inference() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )

  if (!one_term) {
    p <- p +
      ggplot2::facet_wrap(
        ~term,
        scales = "free_y"
      )
  }

  p
}

.filter_inference_terms <- function(x, term) {
  if (is.null(term)) return(x)

  term <- as.character(term)
  unknown <- setdiff(
    term,
    unique(as.character(x$term))
  )

  if (length(unknown) > 0L) {
    stop(
      "Unknown term(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  out <- x[x$term %in% term, , drop = FALSE]

  if (nrow(out) == 0L) {
    stop(
      "No inference rows remain after term filtering.",
      call. = FALSE
    )
  }

  out
}

.filter_inference_features <- function(x, feature) {
  if (is.null(feature)) return(x)

  feature <- as.character(feature)
  unknown <- setdiff(
    feature,
    unique(as.character(x$feature))
  )

  if (length(unknown) > 0L) {
    stop(
      "Unknown feature(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  out <- x[x$feature %in% feature, , drop = FALSE]

  if (nrow(out) == 0L) {
    stop(
      "No inference rows remain after feature filtering.",
      call. = FALSE
    )
  }

  out
}

.theme_inference <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "left",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey35")
    )
}

.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package `ggplot2` is required for inference plots. Install it with ",
      "`install.packages(\"ggplot2\")`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.title_case <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) {
    return("")
  }

  paste0(
    toupper(substr(x, 1L, 1L)),
    substr(x, 2L, nchar(x))
  )
}

utils::globalVariables(c(
  "bootstrap_bias",
  "conf_high",
  "conf_low",
  "estimate",
  "feature",
  "harmonization_fraction",
  "method",
  "n",
  "plot_label",
  "propagated_logp",
  "reference_estimate",
  "reference_logp",
  "reference_se",
  "se",
  "se_ratio",
  "significance",
  "term",
  "transition"
))


