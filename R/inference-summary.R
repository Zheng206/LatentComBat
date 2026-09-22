#' Create an Analysis-Ready LatentComBat Inference Table
#'
#' Adds multiplicity-adjusted p-values, naive (plug-in or conditional)
#' inferential quantities, standard-error ratios, and significance transitions
#' to an object returned by [harm_inference()]. Adjustment is performed
#' separately within each downstream model term.
#'
#' @param x An object returned by [harm_inference()].
#' @param adjust_method Method passed to [stats::p.adjust()]. Defaults to
#'   `"BY"`.
#' @param alpha Significance level used to classify discoveries. Defaults to
#'   `0.05`.
#'
#' @return A data frame with the original inference results and derived
#'   diagnostic columns.
#'
#' @export
inference_table <- function(x, adjust_method = "BY", alpha = 0.05) {
  .check_inference_object(x)
  .check_inference_alpha(alpha)
  adjust_method <- .match_adjust_method(adjust_method)

  out <- as.data.frame(x)
  reference <- .inference_reference(out)

  out$reference_estimate <- reference$estimate
  out$reference_se <- reference$se
  out$reference_statistic <- ifelse(
    is.finite(out$reference_estimate) &
      is.finite(out$reference_se) &
      out$reference_se > 0,
    out$reference_estimate / out$reference_se,
    NA_real_
  )
  out$reference_p_value <- 2 * stats::pnorm(-abs(out$reference_statistic))
  out$se_ratio <- ifelse(
    is.finite(out$se) &
      is.finite(out$reference_se) &
      out$reference_se > 0,
    out$se / out$reference_se,
    NA_real_
  )

  out$p_adjust <- .adjust_within_term(
    p = out$p_value,
    term = out$term,
    method = adjust_method
  )
  out$reference_p_adjust <- .adjust_within_term(
    p = out$reference_p_value,
    term = out$term,
    method = adjust_method
  )

  out$significant <- !is.na(out$p_value) & out$p_value < alpha
  out$significant_adjusted <- !is.na(out$p_adjust) & out$p_adjust < alpha
  out$reference_significant <-
    !is.na(out$reference_p_value) & out$reference_p_value < alpha
  out$reference_significant_adjusted <-
    !is.na(out$reference_p_adjust) & out$reference_p_adjust < alpha

  out$significance_transition <- .significance_transition(
    reference = out$reference_significant_adjusted,
    propagated = out$significant_adjusted
  )

  attr(out, "reference_label") <- reference$label
  attr(out, "adjust_method") <- adjust_method
  attr(out, "alpha") <- alpha
  out
}


#' Summarize Discovery Changes After Uncertainty Propagation
#'
#' Compares multiplicity-adjusted discoveries from naive inference with those
#' from uncertainty-aware inference.
#'
#' @inheritParams inference_table
#' @param detailed Logical. If `TRUE`, return one row per feature and term. If
#'   `FALSE`, return transition counts by term.
#'
#' @return A data frame containing either feature-level discovery transitions
#'   or their counts by term.
#'
#' @export
inference_significance_changes <- function(
    x,
    adjust_method = "BY",
    alpha = 0.05,
    detailed = FALSE
) {
  tab <- inference_table(x, adjust_method = adjust_method, alpha = alpha)
  keep <- c(
    "feature", "term", "reference_p_value", "p_value",
    "reference_p_adjust", "p_adjust",
    "reference_significant_adjusted", "significant_adjusted",
    "significance_transition"
  )
  detail <- tab[, keep, drop = FALSE]
  if (isTRUE(detailed)) return(detail)

  counts <- as.data.frame(
    table(
      term = detail$term,
      transition = detail$significance_transition,
      useNA = "no"
    ),
    stringsAsFactors = FALSE
  )
  names(counts)[names(counts) == "Freq"] <- "n"
  counts[counts$n > 0L, , drop = FALSE]
}


#' Summarize LatentComBat Uncertainty-Aware Inference
#'
#' Reports standard-error inflation, discovery changes, and, when available,
#' the contribution of harmonization uncertainty or bootstrap bias.
#'
#' @param object An object returned by [harm_inference()].
#' @inheritParams inference_table
#' @param ... Unused.
#'
#' @return An object of class `summary_latentcombat_inference`.
#'
#' @method summary latentcombat_inference
#' @export
summary.latentcombat_inference <- function(
    object,
    adjust_method = "BY",
    alpha = 0.05,
    ...
) {
  tab <- inference_table(
    object,
    adjust_method = adjust_method,
    alpha = alpha
  )

  pieces <- split(tab, tab$term, drop = TRUE)
  term_summary <- lapply(pieces, function(z) {
    transitions <- table(z$significance_transition)
    transition_n <- function(label) {
      if (label %in% names(transitions)) unname(transitions[[label]]) else 0L
    }

    data.frame(
      term = as.character(z$term[1L]),
      n_features = nrow(z),
      median_reference_se = .finite_median(z$reference_se),
      median_propagated_se = .finite_median(z$se),
      median_se_ratio = .finite_median(z$se_ratio),
      reference_discoveries = sum(z$reference_significant_adjusted, na.rm = TRUE),
      propagated_discoveries = sum(z$significant_adjusted, na.rm = TRUE),
      retained = transition_n("Retained"),
      lost = transition_n("Lost"),
      gained = transition_n("Gained"),
      median_harmonization_fraction = if (
        "harmonization_fraction" %in% names(z)
      ) .finite_median(z$harmonization_fraction) else NA_real_,
      median_absolute_bootstrap_bias = if (
        "bootstrap_bias" %in% names(z)
      ) .finite_median(abs(z$bootstrap_bias)) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  term_summary <- do.call(rbind, term_summary)
  rownames(term_summary) <- NULL

  inference_method <- unique(as.character(tab$inference_method))
  harmonization_method <- unique(as.character(tab$harmonization_method))

  out <- list(
    inference_method = inference_method,
    harmonization_method = harmonization_method,
    reference_label = attr(tab, "reference_label"),
    adjust_method = adjust_method,
    alpha = alpha,
    n_rows = nrow(tab),
    n_features = length(unique(tab$feature)),
    terms = unique(as.character(tab$term)),
    term_summary = term_summary,
    results = tab
  )
  class(out) <- "summary_latentcombat_inference"
  out
}


#' Print a LatentComBat Inference Summary
#'
#' @param x An object returned by [summary.latentcombat_inference()].
#' @param digits Number of digits shown for numeric summaries.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#'
#' @method print summary_latentcombat_inference
#' @export
print.summary_latentcombat_inference <- function(x, digits = 3L, ...) {
  cat("LatentComBat uncertainty-aware inference summary\n")
  cat("  Harmonization: ", paste(x$harmonization_method, collapse = ", "), "\n", sep = "")
  cat("  Inference:     ", paste(x$inference_method, collapse = ", "), "\n", sep = "")
  cat("  Comparison:    ", x$reference_label, " vs uncertainty-aware\n", sep = "")
  cat("  Multiplicity:  ", x$adjust_method, " within term\n", sep = "")
  cat("  Alpha:         ", format(x$alpha), "\n", sep = "")
  cat("  Features:      ", x$n_features, "\n\n", sep = "")

  shown <- x$term_summary
  numeric_cols <- vapply(shown, is.numeric, logical(1L))
  shown[numeric_cols] <- lapply(shown[numeric_cols], function(z) round(z, digits))
  print(shown, row.names = FALSE)
  invisible(x)
}


.check_inference_object <- function(x) {
  if (!inherits(x, "latentcombat_inference")) {
    stop("`x` must be an object returned by `harm_inference()`.", call. = FALSE)
  }
  required <- c(
    "feature", "term", "estimate", "se", "p_value",
    "inference_method", "harmonization_method"
  )
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(
      "Inference object is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


.check_inference_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}


.match_adjust_method <- function(method) {
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("`adjust_method` must be one method accepted by `stats::p.adjust()`.", call. = FALSE)
  }
  choices <- stats::p.adjust.methods
  matched <- choices[tolower(choices) == tolower(method)]
  if (length(matched) != 1L) {
    stop(
      "Unknown `adjust_method`. Choose one of: ",
      paste(choices, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  matched
}


.adjust_within_term <- function(p, term, method) {
  out <- rep(NA_real_, length(p))
  groups <- split(seq_along(p), term, drop = TRUE)
  for (idx in groups) out[idx] <- stats::p.adjust(p[idx], method = method)
  out
}


.inference_reference <- function(x) {
  if (all(c("plugin_estimate", "plugin_se") %in% names(x))) {
    return(list(
      estimate = x$plugin_estimate,
      se = x$plugin_se,
      label = "plug-in"
    ))
  }
  if ("conditional_se" %in% names(x)) {
    return(list(
      estimate = x$estimate,
      se = x$conditional_se,
      label = "conditional"
    ))
  }
  stop(
    "Naive reference inference is unavailable. Expected `plugin_estimate` and ",
    "`plugin_se`, or `conditional_se`.",
    call. = FALSE
  )
}


.significance_transition <- function(reference, propagated) {
  out <- ifelse(
    reference & propagated,
    "Retained",
    ifelse(
      reference & !propagated,
      "Lost",
      ifelse(!reference & propagated, "Gained", "Never significant")
    )
  )
  factor(
    out,
    levels = c("Retained", "Lost", "Gained", "Never significant")
  )
}


.finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else stats::median(x)
}
