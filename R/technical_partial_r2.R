#' Compute Technical-Factor Partial R-Squared
#'
#' Computes the additional residual variation explained by each technical
#' factor, conditional on biological covariates and (by default) primary batch,
#' with optional adjustment for the other secondary technical factors.
#'
#' @param data Numeric matrix or data frame with observations in rows and
#'   imaging features in columns.
#' @param technical Data frame of secondary technical factors. Its columns are
#'   tested one at a time. Whether other technical columns enter the models
#'   is controlled by `test`.
#' @param covar Optional data frame of protected biological covariates.
#' @param bat Vector identifying the primary batch. Required when
#'   `adjust_batch = TRUE`.
#' @param model Model-fitting function, such as [stats::lm()] or
#'   [mgcv::gam()].
#' @param formula Formula for the biological covariates, with response `y`.
#'   For example, `y ~ s(age, k = 5) + sex`. Predictors must belong to `covar`;
#'   `.` expands to all biological covariates.
#' @param adjust_batch Whether to condition on primary batch. Defaults to `TRUE`.
#' @param test Either `"separate"` (default), testing each technical factor
#'   conditional on all other technical factors, or `"joint"`, comparing biology
#'   plus optional batch against biology plus optional batch plus that factor
#'   alone. In `"joint"` mode, no other secondary technical factors are adjusted
#'   for; this option does not test the entire technical block together.
#' @param ... Additional arguments passed to `model`. Fits must retain the
#'   supplied analysis sample; arguments that drop rows are not supported.
#'
#' @return A data frame with one row per feature and technical factor, containing
#'   `partial_R2`, `rss_reduced`, `rss_full`, `n_obs`, `rank_reduced`, `rank_full`,
#'   `rank_added`, `rank_deficient_reduced`, `rank_deficient_full`,
#'   `rank_deficient`, `estimable`, `negative`, `status`, and `error`.
#'   `partial_R2` is `NA` for failed fits, non-estimable factors, or a zero
#'   reduced RSS. Rank deficiency in unrelated nuisance terms does not by itself
#'   invalidate the comparison. `estimable` means the tested term adds at least
#'   one design dimension; it does not guarantee every coefficient is identified.
#'   Unknown diagnostics are `NA`. Original metadata names are retained.
#'   `status` is one of `"ok"`, `"negative"`, `"non_estimable"`, `"zero_rss"`,
#'   or `"fit_error"`. The `error` column preserves fitting error messages.
#'   Both modes return one row per feature and technical factor, preserving the
#'   original technical names.
#'   Attributes `test` and `technical_factors` record the test mode and original
#'   technical column names. Both modes work with [summarize_r2()] and
#'   [plot_technical_r2()]. The `formula` attribute contains the full-model
#'   formula in separate mode or a named list of per-factor full-model formulas
#'   in joint mode.
#'
#' @details
#' The measure is `(rss_reduced - rss_full) / rss_reduced`. Both models use the
#' same complete cases of the variables in their full formula. Unused covariates
#' do not affect this sample.
#'
#' In separate mode, the reduced model drops only the tested technical factor.
#' All technical factors enter the full model, which is fitted once per feature.
#' Complete cases are common to all factors within that feature.
#'
#' In joint mode, each factor is added individually to a baseline containing
#' only biology and optional primary batch. For factor `T_k`, the comparison is
#' `y ~ biology + batch` versus `y ~ biology + batch + T_k` (omitting batch
#' when `adjust_batch = FALSE`). This includes associations shared with other
#' secondary technical factors. Complete cases are selected per feature and
#' tested factor; missingness in other technical columns does not remove rows.
#' Contributions of correlated technical factors can overlap and should not be
#' summed. In either mode, a factor adding no design dimension gives `NA`.
#'
#' For nested ordinary least-squares models this is conventional partial
#' R-squared. For penalized GAMs, smoothing parameters may change between fits,
#' so this is an empirical RSS comparison rather than conventional OLS partial
#' R-squared. Negative values are retained and flagged, never clipped. Design
#' ranks describe unpenalized column spaces, not GAM effective degrees of freedom
#' or all penalty-related identifiability issues.
#'
#' Conditioning on batch measures secondary technical variation beyond batch;
#' it can be insensitive to removal of additive batch means. Setting
#' `adjust_batch = FALSE` also includes variation shared with batch. Neither
#' measure is causal. Before/after changes can reflect changes in both explained
#' and residual variation, so assess biological preservation separately.
#'
#' @seealso [summarize_r2()] for readable feature summaries and
#'   [plot_technical_r2()] for before/after plots.
#' @examples
#' set.seed(1)
#' technical <- data.frame(scanner = rnorm(30), field = rnorm(30))
#' data <- cbind(feature1 = technical$scanner + rnorm(30), feature2 = rnorm(30))
#' joint <- technical_partial_r2(data, technical, adjust_batch = FALSE,
#'                               test = "joint")
#' summarize_r2(joint, method = "Raw")
#' @export
technical_partial_r2 <- function(
    data, technical, covar = NULL, bat = NULL, model = stats::lm,
    formula = y ~ 1, adjust_batch = TRUE, test = c("separate", "joint"), ...
) {
  test <- match.arg(test)
  data <- as.data.frame(data)
  technical <- as.data.frame(technical)
  covar <- if (is.null(covar)) data.frame(row.names = seq_len(nrow(data))) else as.data.frame(covar)
  if (nrow(data) == 0L || ncol(data) == 0L) {
    stop("`data` must contain at least one observation and feature.", call. = FALSE)
  }
  if (!all(vapply(data, is.numeric, logical(1)))) {
    stop("All features in `data` must be numeric.", call. = FALSE)
  }
  if (!is.logical(adjust_batch) || length(adjust_batch) != 1L || is.na(adjust_batch)) {
    stop("`adjust_batch` must be TRUE or FALSE.", call. = FALSE)
  }
  if (adjust_batch && is.null(bat)) {
    stop("`bat` is required when `adjust_batch = TRUE`.", call. = FALSE)
  }
  if (nrow(technical) != nrow(data) || nrow(covar) != nrow(data) ||
      (adjust_batch && length(bat) != nrow(data))) {
    stop("`data`, `technical`, `covar`, and (when adjusted for) `bat` must have the same number of observations.", call. = FALSE)
  }
  if (ncol(technical) == 0L) {
    stop("`technical` must contain at least one factor.", call. = FALSE)
  }
  metadata_names <- c(names(covar), names(technical))
  if (anyDuplicated(names(data)) || anyDuplicated(metadata_names) ||
      anyNA(metadata_names) || any(!nzchar(metadata_names))) {
    stop("Feature names must be unique; covariate and technical names must be nonempty and unique across both data frames.", call. = FALSE)
  }
  if (any(metadata_names %in% c("y", ".primary_batch"))) {
    stop("Metadata names `y` and `.primary_batch` are reserved.", call. = FALSE)
  }
  if (!inherits(formula, "formula") || length(formula) != 3L || !identical(formula[[2L]], as.name("y"))) {
    stop("`formula` must use `y` as its response.", call. = FALSE)
  }
  # Expand dots using only biological covariates; build added terms as symbols
  # so non-syntactic metadata names never need parsing or renaming.
  biological_formula <- stats::formula(stats::terms(
    formula, data = cbind(covar, y = rep(0, nrow(data)))
  ))
  biological_vars <- all.vars(biological_formula[[3L]])
  if (!all(biological_vars %in% names(covar))) {
    stop("All biological formula predictors must be columns of `covar`.", call. = FALSE)
  }
  add_terms <- function(base, variables) {
    for (variable in variables) base[[3L]] <- call("+", base[[3L]], as.name(variable))
    base
  }
  technical_names <- names(technical)
  batch_terms <- if (adjust_batch) ".primary_batch" else character()
  base_data <- cbind(covar, technical)
  if (adjust_batch) base_data$.primary_batch <- as.factor(bat)
  full_formula <- add_terms(biological_formula, c(batch_terms, technical_names))
  used_vars <- all.vars(full_formula)

  fit_model <- function(fit_formula, frame) {
    tryCatch({
      fit <- model(fit_formula, data = frame, ...)
      design <- stats::model.matrix(fit)
      residual <- stats::residuals(fit)
      if (nrow(design) != nrow(frame) || length(residual) != nrow(frame) ||
          !identical(rownames(design), rownames(frame))) {
        stop("Fit did not retain the common analysis sample.")
      }
      if (any(!is.finite(design)) || any(!is.finite(residual))) {
        stop("Fit produced non-finite design entries or residuals.")
      }
      rank <- qr(design)$rank
      rss <- sum(residual^2)
      if (!is.finite(rss)) stop("Fit produced a non-finite residual sum of squares.")
      list(rank = rank, deficient = rank < ncol(design),
           rss = rss, error = NA_character_)
    }, error = function(e) list(rank = NA_integer_, deficient = NA,
                               rss = NA_real_, error = conditionMessage(e)))
  }
  full_formulas <- if (test == "joint") stats::setNames(lapply(technical_names, function(name) {
    add_terms(biological_formula, c(batch_terms, name))
  }), technical_names) else NULL
  rows <- vector("list", ncol(data) * length(technical_names))
  idx <- 1L
  for (feature_name in names(data)) {
    model_data <- base_data
    model_data$y <- data[[feature_name]]
    if (test == "separate") {
      frame <- model_data[stats::complete.cases(model_data[, used_vars, drop = FALSE]), , drop = FALSE]
      full <- if (nrow(frame)) fit_model(full_formula, frame) else
        list(rank = NA_integer_, deficient = NA, rss = NA_real_, error = "No complete observations for the full model.")
    }
    for (technical_name in technical_names) {
      if (test == "joint") {
        factor_formula <- full_formulas[[technical_name]]
        factor_vars <- all.vars(factor_formula)
        frame <- model_data[stats::complete.cases(model_data[, factor_vars, drop = FALSE]), , drop = FALSE]
        full <- if (nrow(frame)) fit_model(factor_formula, frame) else
          list(rank = NA_integer_, deficient = NA, rss = NA_real_, error = "No complete observations for the full model.")
      }
      remaining_technical <- if (test == "joint") character() else setdiff(technical_names, technical_name)
      reduced_formula <- add_terms(biological_formula, c(batch_terms, remaining_technical))
      reduced <- if (nrow(frame)) fit_model(reduced_formula, frame) else full
      errors <- c(if (!is.na(full$error)) paste0("Full model: ", full$error),
                  if (!is.na(reduced$error)) paste0("Reduced model: ", reduced$error))
      rank_added <- full$rank - reduced$rank
      estimable <- if (is.na(rank_added)) NA else rank_added > 0L
      partial_r2 <- NA_real_
      status <- if (length(errors)) "fit_error" else if (!estimable) "non_estimable" else if (reduced$rss <= 0) "zero_rss" else "ok"
      if (status == "ok") partial_r2 <- (reduced$rss - full$rss) / reduced$rss
      negative <- if (is.na(partial_r2)) NA else partial_r2 < 0
      if (isTRUE(negative)) status <- "negative"
      rows[[idx]] <- data.frame(
        Feature = feature_name, Technical = technical_name, partial_R2 = partial_r2,
        rss_reduced = reduced$rss, rss_full = full$rss, n_obs = nrow(frame),
        rank_reduced = reduced$rank, rank_full = full$rank, rank_added = rank_added,
        rank_deficient_reduced = reduced$deficient, rank_deficient_full = full$deficient,
        rank_deficient = reduced$deficient | full$deficient, estimable = estimable,
        negative = negative, status = status,
        error = if (length(errors)) paste(errors, collapse = "; ") else NA_character_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  out <- do.call(rbind, rows)
  attr(out, "formula") <- if (test == "joint") full_formulas else full_formula
  attr(out, "model") <- model
  attr(out, "adjust_batch") <- adjust_batch
  attr(out, "test") <- test
  attr(out, "technical_factors") <- technical_names
  out
}

#' Summarize Technical-Factor Partial R-Squared
#'
#' Summarizes feature-wise partial R-squared estimates for each technical factor.
#'
#' @param x A data frame returned by [technical_partial_r2()], or a data frame
#'   containing `Technical` and numeric `partial_R2` columns.
#' @param method Optional nonempty string identifying the harmonization method
#'   or data stage, such as `"Raw"` or `"CALM"`. If supplied, adds `Method`.
#'
#' @return A data frame with one row per observed technical factor, in order
#'   of first appearance, with columns `Technical`, `median_R2`, `q25_R2`,
#'   `q75_R2`, `n_valid`, and `n_total`, followed by `Method` when supplied.
#' @details Statistics use finite estimates, retaining negative values. Groups
#'   with no valid estimates remain in the output with `NA` statistics. Quartiles use
#'   [stats::quantile()] with its default interpolation (`type = 7`). These
#'   statistics describe variation across features, not uncertainty in a fitted
#'   effect. Counts refer to rows of `x`; each row should represent one feature
#'   and technical factor, as in the output of [technical_partial_r2()].
#' @examples
#' x <- data.frame(
#'   Technical = c("scanner", "scanner", "scanner", "field_strength"),
#'   partial_R2 = c(0.1, 0.3, NA, 0.2)
#' )
#' summarize_r2(x, method = "Raw")
#' @export
summarize_r2 <- function(x, method = NULL) {
  if (!is.data.frame(x) || !all(c("Technical", "partial_R2") %in% names(x)) ||
      !is.numeric(x$partial_R2) || !is.null(dim(x$partial_R2))) {
    stop("`x` must be a data frame with `Technical` and numeric `partial_R2` columns.", call. = FALSE)
  }
  if (!is.null(method) && (!is.character(method) || length(method) != 1L ||
                          is.na(method) || !nzchar(method))) {
    stop("`method` must be NULL or a single nonempty string.", call. = FALSE)
  }
  technical <- as.character(x$Technical)
  if (anyNA(technical) || any(!nzchar(technical))) {
    stop("`Technical` must contain nonmissing, nonempty names.", call. = FALSE)
  }
  out <- data.frame(
    Technical = character(), median_R2 = numeric(), q25_R2 = numeric(),
    q75_R2 = numeric(), n_valid = integer(), n_total = integer()
  )
  rows <- lapply(unique(technical), function(name) {
    values <- x$partial_R2[technical == name]
    valid <- values[is.finite(values)]
    quartiles <- if (length(valid)) stats::quantile(valid, c(0.25, 0.75), names = FALSE) else c(NA_real_, NA_real_)
    data.frame(
      Technical = name,
      median_R2 = if (length(valid)) stats::median(valid) else NA_real_,
      q25_R2 = quartiles[1L], q75_R2 = quartiles[2L],
      n_valid = length(valid), n_total = length(values)
    )
  })
  if (length(rows)) out <- do.call(rbind, rows)
  if (!is.null(method)) out$Method <- rep(method, nrow(out))
  out
}

#' Plot Technical-Factor Partial R-Squared
#'
#' Visualizes median changes in technical partial R-squared relative to raw
#' data, or feature distributions before and after harmonization.
#'
#' @param raw_result Result from [technical_partial_r2()] for raw data.
#' @param harmonized_result Result from [technical_partial_r2()] for harmonized
#'   data, or a nonempty named list of results for several harmonization methods.
#' @param type Either `"aggregate"` (default: median paired change across features)
#'   or `"feature"` (boxplots of absolute partial R-squared across features,
#'   including raw data and each harmonization method).
#' @param labels Distinct labels for raw data followed by each harmonized result.
#'   Defaults to `c("Before", "After")` for a single result. For a named list,
#'   the default is `c("Before", names(harmonized_result))`.
#' @param limit Optional positive number specifying the symmetric heatmap color
#'   limits in percentage points. Defaults to the largest absolute displayed
#'   change, or 1 if all changes are zero or unavailable. Values outside the
#'   color limits are squished; aggregate numeric labels retain actual changes.
#'   Ignored for `type = "feature"`.
#'
#' @return A `ggplot` object. Aggregate data contain `Technical`, `Method`,
#'   `median_delta_r2` (percentage points), `n_valid`, `n_total`, and `n_failed`.
#'   Feature boxplot data contain `Feature`, `Technical`, `Method`, and
#'   `partial_R2` on its original scale.
#' @details Changes are computed as
#'   `100 * (harmonized partial_R2 - raw partial_R2)` for each matched feature
#'   and technical factor. This is a percentage-point difference, not a percent
#'   reduction. Negative values indicate reduced residual technical variation.
#'
#'   The aggregate heatmap has methods in columns and technical factors in rows,
#'   with signed labels showing the median of finite paired changes. This is
#'   not a difference of unpaired medians. Missing or non-finite pairs remain in
#'   the counts but do not enter the median. Groups without valid pairs are grey
#'   and labeled `NA`. Counts refer to the union of raw and harmonized keys for
#'   each method.
#'
#'   The feature boxplots show absolute partial R-squared on its original scale,
#'   with raw data and harmonization methods in columns and one facet per
#'   technical factor. Boxes summarize finite estimates within each group;
#'   unavailable estimates are excluded, and negative values are retained.
#'   Each method uses all its finite feature estimates, without requiring a
#'   matching raw estimate. Whiskers extend to the most extreme estimates within
#'   1.5 times the interquartile range; more extreme estimates are plotted as
#'   outliers. No individual feature labels are drawn. Feature names and
#'   unavailable estimates remain in the plot's `data` component. Groups without
#'   finite estimates have no box. See [summarize_r2()] for valid/total counts.
#'
#'   Duplicate feature/technical-factor keys are rejected. Use the same model
#'   specification, test mode (`"separate"` or `"joint"`), and analysis sample
#'   before and after harmonization to make the comparison meaningful. Both test
#'   modes give one heatmap row or boxplot facet per technical factor.
#'   See [summarize_r2()] for tabular summaries.
#' @export
plot_technical_r2 <- function(raw_result, harmonized_result,
                              type = c("aggregate", "feature"), labels = c("Before", "After"),
                              limit = NULL) {
  type <- match.arg(type)
  required <- c("Feature", "Technical", "partial_R2")
  multiple <- is.list(harmonized_result) && !is.data.frame(harmonized_result)
  if (multiple) {
    methods <- names(harmonized_result)
    if (!length(harmonized_result) || is.null(methods) || anyNA(methods) ||
        any(!nzchar(methods)) || anyDuplicated(methods)) {
      stop("Multiple harmonized results must be supplied as a nonempty list with distinct method names.", call. = FALSE)
    }
    results <- harmonized_result
    if (missing(labels)) labels <- c("Before", methods)
  } else {
    results <- list(harmonized_result)
  }
  all_results <- c(list(raw_result), results)
  valid_result <- function(x) {
    is.data.frame(x) && all(required %in% names(x)) &&
      is.numeric(x$partial_R2) && is.null(dim(x$partial_R2))
  }
  if (!all(vapply(all_results, valid_result, logical(1)))) {
    stop("All results must be data frames returned by `technical_partial_r2()`.", call. = FALSE)
  }
  if (!is.character(labels) || length(labels) != length(all_results) || anyNA(labels) ||
      any(!nzchar(labels)) || anyDuplicated(labels)) {
    stop("`labels` must contain distinct nonempty strings for raw data and each harmonized result.", call. = FALSE)
  }
  if (type == "aggregate" && !is.null(limit) && (!is.numeric(limit) || length(limit) != 1L ||
                         !is.finite(limit) || limit <= 0)) {
    stop("`limit` must be NULL or a positive finite number in percentage points.", call. = FALSE)
  }
  keys <- c("Feature", "Technical")
  prepare <- function(x) {
    out <- data.frame(Feature = as.character(x$Feature), Technical = as.character(x$Technical),
                      partial_R2 = x$partial_R2)
    if (anyNA(out[, keys]) || any(vapply(out[, keys], function(v) any(!nzchar(v)), logical(1)))) {
      stop("Feature and technical-factor names must be nonmissing and nonempty.", call. = FALSE)
    }
    if (anyDuplicated(out[, keys])) {
      stop("Each feature and technical factor must identify a unique row in each result.", call. = FALSE)
    }
    out
  }
  prepared <- lapply(all_results, prepare)
  technical_levels <- unique(unlist(lapply(prepared, function(x) x$Technical), use.names = FALSE))
  if (type == "feature") {
    feature_data <- do.call(rbind, lapply(seq_along(prepared), function(i) {
      out <- prepared[[i]]
      out$Method <- rep(labels[i], nrow(out))
      out
    }))
    feature_data$Method <- factor(feature_data$Method, levels = labels)
    feature_data$Technical <- factor(feature_data$Technical, levels = technical_levels)
    # Keep every input row available in plot$data, while excluding non-finite
    # values from boxplot statistics.
    finite_data <- feature_data
    finite_data$partial_R2[!is.finite(finite_data$partial_R2)] <- NA_real_
    method_colors <- c("#4F86A6", "#66A69C", "#C77963", "#9B8FB5", "#C4A35A", "#7F9BA3")
    if (length(results) > length(method_colors)) {
      method_colors <- grDevices::hcl.colors(length(results), palette = "Dark 3")
    }
    fill_colors <- stats::setNames(c("#B8BDC5", method_colors[seq_along(results)]), labels)
    return(ggplot2::ggplot(feature_data, ggplot2::aes(x = .data[["Method"]], y = .data[["partial_R2"]], fill = .data[["Method"]])) +
      ggplot2::geom_boxplot(data = finite_data, na.rm = TRUE, width = 0.65,
                            colour = "#35404A", linewidth = 0.5, outlier.size = 1.4, outlier.alpha = 0.6) +
      ggplot2::facet_wrap(~ Technical, nrow = 1, drop = FALSE) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::scale_fill_manual(values = fill_colors, drop = FALSE) +
      ggplot2::labs(x = NULL, y = "Partial R2") +
      ggplot2::theme_classic(base_size = 13) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(face = "bold", angle = 45, hjust = 1),
        axis.text.y = ggplot2::element_text(face = "bold"), axis.title = ggplot2::element_text(face = "bold"),
        legend.position = "none"
      ))
  }
  raw <- prepared[[1L]]
  names(raw)[names(raw) == "partial_R2"] <- "raw_R2"
  paired <- lapply(seq_along(results), function(i) {
    harmonized <- prepared[[i + 1L]]
    names(harmonized)[names(harmonized) == "partial_R2"] <- "harmonized_R2"
    out <- merge(raw, harmonized, by = keys, all = TRUE, sort = FALSE)
    out$delta_r2 <- 100 * (out$harmonized_R2 - out$raw_R2)
    out$delta_r2[!is.finite(out$delta_r2)] <- NA_real_
    out$Method <- rep(labels[i + 1L], nrow(out))
    out
  })
  paired <- do.call(rbind, paired)
  groups <- unique(paired[, c("Technical", "Method"), drop = FALSE])
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    group <- groups[i, , drop = FALSE]
    values <- paired$delta_r2[paired$Technical == group$Technical & paired$Method == group$Method]
    valid <- is.finite(values)
    data.frame(group, median_delta_r2 = if (any(valid)) stats::median(values[valid]) else NA_real_,
               n_valid = sum(valid), n_total = length(values), n_failed = sum(!valid))
  })
  heatmap_data <- if (length(rows)) do.call(rbind, rows) else
    data.frame(Technical = character(), Method = character(), median_delta_r2 = numeric(),
               n_valid = integer(), n_total = integer(), n_failed = integer())
  heatmap_data$Technical <- factor(heatmap_data$Technical, levels = rev(technical_levels))
  heatmap_data$change_label <- ifelse(is.finite(heatmap_data$median_delta_r2),
                                    sprintf("%+.2f", heatmap_data$median_delta_r2), "NA")
  heatmap_data$Method <- factor(heatmap_data$Method, levels = labels[-1L])
  changes <- heatmap_data$median_delta_r2
  plot <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = .data[["Method"]], y = .data[["Technical"]], fill = .data[["median_delta_r2"]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data[["change_label"]]), size = 3.2, fontface = "bold")
  if (is.null(limit)) {
    limit <- max(c(0, abs(changes[is.finite(changes)])))
    if (limit == 0) limit <- 1
  }
  plot <- plot +
    ggplot2::scale_fill_gradient2(low = "#4F86A6", mid = "#F7F7F7", high = "#C77963",
                                  midpoint = 0, limits = c(-limit, limit), oob = scales::squish,
                                  na.value = "grey85", name = expression(bold(Delta * R^2~"(pp) vs Raw"))) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"),
      axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold"), axis.text.y = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"), legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"), legend.text = ggplot2::element_text(face = "bold")
    )
  plot
}
