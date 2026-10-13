#' Extract Feature-wise Fixed Effects
#'
#' Extracts parametric biological coefficients from feature-wise linear models
#' or GAMs for exploratory comparisons before and after harmonization.
#'
#' @param data Numeric observation-by-feature matrix or data frame.
#' @param covar Data frame of biological covariates, aligned with rows of `data`.
#' @param model Model-fitting function returning an `lm` or `gam` object.
#' @param formula Formula with response `y`; predictors must belong to `covar`.
#'   The default `y ~ .` uses all covariates.
#' @param method Nonempty method label, such as `"Raw"`, `"ComBat"`, or `"CALM"`.
#' @param terms Optional character vector of exact coefficient names to retain,
#'   such as `"age"` or `"sexM"`. The intercept is always excluded.
#' @param ... Additional arguments passed to `model`. Arguments that drop or
#'   reorder the supplied complete-case sample are not supported.
#'   To set a GAM fitting method such as REML, supply a wrapper as `model`,
#'   because `method` here labels the resulting data rather than the GAM fit.
#'
#' @return A data frame with `Feature`, `Term`, `Estimate`, `SE`, `Statistic`,
#'   `P_value`, `Q_value`, `Method`, `N_obs`, `Status`, and `Error`. Failed or
#'   non-estimable coefficients are retained with unavailable values set to `NA`.
#'   The `feature_diagnostics` attribute records fit status, observation count,
#'   and errors for every feature, including when there are no retained terms.
#' @details Complete cases are selected per feature using only variables in the
#'   formula. Factors retain their global contrasts: fits losing observed factor
#'   levels after filtering are marked as failures to avoid changing references.
#'   Actual R coefficient names identify the fitted contrasts. GAMs contribute
#'   only their parametric coefficient table, excluding spline coefficients and
#'   smooth-term tests. A smooth-only GAM returns no coefficient rows.
#'
#'   BH adjustment is performed separately for each coefficient across features.
#'   The adjustment family includes all input features, including failed tests;
#'   unavailable p-values stay `NA`. Unknown requested terms are rejected when
#'   coefficients are available. If every fit fails and no terms were specified,
#'   extraction stops with the fitting errors.
#'
#'   These are conditional, conventional downstream model estimates and p-values;
#'   they do not propagate harmonization uncertainty. Use [harm_inference()] for
#'   uncertainty-aware inference, subject to its assumptions. Use matched subjects,
#'   features, formulas, units, and contrasts for before/after comparisons. Raw
#'   estimates may themselves be confounded, so coefficient stability alone does
#'   not establish accuracy. GAM smooth effects require separate curve diagnostics.
#' @seealso [plot_fixed_effect_comparison()], [technical_partial_r2()]
#' @examples
#' set.seed(1)
#' covar <- data.frame(age = rnorm(30))
#' data <- cbind(feature1 = covar$age + rnorm(30), feature2 = rnorm(30))
#' extract_fixed_effects(data, covar, formula = y ~ age)
#' @export
extract_fixed_effects <- function(data, covar, model = stats::lm,
                                  formula = y ~ ., method = "Raw", terms = NULL, ...) {
  data <- as.data.frame(data)
  covar <- as.data.frame(covar)
  if (!nrow(data) || !ncol(data) || !all(vapply(data, is.numeric, logical(1)))) {
    stop("`data` must contain numeric observations and features.", call. = FALSE)
  }
  if (nrow(data) != nrow(covar)) stop("`data` and `covar` must have the same number of rows.", call. = FALSE)
  if (anyDuplicated(names(data)) || anyNA(names(data)) || any(!nzchar(names(data))) ||
      anyDuplicated(names(covar)) || anyNA(names(covar)) || any(!nzchar(names(covar)))) {
    stop("Feature and covariate names must be nonempty and unique within each data frame.", call. = FALSE)
  }
  if ("y" %in% names(covar)) stop("Covariate name `y` is reserved for the response.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 3L || !identical(formula[[2L]], as.name("y"))) {
    stop("`formula` must use `y` as its response.", call. = FALSE)
  }
  if (!is.function(model)) stop("`model` must be a fitting function.", call. = FALSE)
  if (!is.character(method) || length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("`method` must be a single nonempty string.", call. = FALSE)
  }
  if (!is.null(terms) && (!is.character(terms) || !length(terms) || anyNA(terms) ||
                        any(!nzchar(terms)) || anyDuplicated(terms) || "(Intercept)" %in% terms)) {
    stop("`terms` must contain distinct nonempty coefficient names, excluding the intercept.", call. = FALSE)
  }
  covar[] <- lapply(covar, function(x) {
    if (is.character(x)) as.factor(x) else if (is.factor(x)) droplevels(x) else x
  })
  fit_formula <- stats::formula(stats::terms(formula, data = cbind(covar, y = rep(0, nrow(data)))))
  predictors <- all.vars(fit_formula[[3L]])
  if (!all(predictors %in% names(covar))) {
    stop("All formula predictors must be columns of `covar`.", call. = FALSE)
  }
  empty_table <- data.frame(Term = character(), Estimate = numeric(), SE = numeric(),
                            Statistic = numeric(), P_value = numeric())
  model_args <- list(...)
  tables <- vector("list", ncol(data))
  diagnostics <- vector("list", ncol(data))
  for (g in seq_along(data)) {
    frame <- covar
    frame$y <- data[[g]]
    frame <- frame[stats::complete.cases(frame[, c("y", predictors), drop = FALSE]), , drop = FALSE]
    error <- NA_character_
    table <- tryCatch({
      if (!nrow(frame)) stop("No complete observations for the model.")
      for (name in predictors) {
        if (is.factor(frame[[name]]) && !all(levels(covar[[name]]) %in% as.character(frame[[name]]))) {
          stop("Factor `", name, "` lost levels after filtering; coefficient contrasts would change.")
        }
      }
      fit <- do.call(model, c(list(formula = fit_formula, data = frame), model_args))
      if (!inherits(fit, "lm") && !inherits(fit, "gam")) stop("Only lm and gam fits are supported.")
      design <- stats::model.matrix(fit)
      if (nrow(design) != nrow(frame) || !identical(rownames(design), rownames(frame)) ||
          length(stats::residuals(fit)) != nrow(frame)) stop("Fit did not retain the analysis sample.")
      if (inherits(fit, "gam")) {
        tab <- summary(fit)$p.table
        coefficient_names <- rownames(tab)
      } else {
        tab <- summary(fit)$coefficients
        coefficient_names <- names(stats::coef(fit))
      }
      coefficient_names <- setdiff(coefficient_names, "(Intercept)")
      if (!length(coefficient_names)) {
        empty_table
      } else {
        p_col <- grep("^Pr\\(", colnames(tab))
        if (length(p_col) != 1L) stop("P-value column is unavailable.")
        idx <- match(coefficient_names, rownames(tab))
        data.frame(Term = coefficient_names, Estimate = as.numeric(tab[idx, 1L]),
                   SE = as.numeric(tab[idx, 2L]), Statistic = as.numeric(tab[idx, 3L]),
                   P_value = as.numeric(tab[idx, p_col]))
      }
    }, error = function(e) {
      error <<- conditionMessage(e)
      empty_table
    })
    tables[[g]] <- table
    diagnostics[[g]] <- data.frame(Feature = names(data)[g], N_obs = nrow(frame),
                                   Status = if (is.na(error)) "ok" else "fit_error", Error = error)
  }
  diagnostics <- do.call(rbind, diagnostics)
  available <- unique(unlist(lapply(tables, function(x) x$Term), use.names = FALSE))
  if (is.null(terms)) {
    terms <- available
    if (!length(terms) && all(diagnostics$Status == "fit_error")) {
      stop("No coefficients could be extracted: ", paste(unique(diagnostics$Error), collapse = "; "), call. = FALSE)
    }
  } else if (any(diagnostics$Status == "ok") && !all(terms %in% available)) {
    stop("Unknown coefficient names in `terms`: ", paste(setdiff(terms, available), collapse = ", "), call. = FALSE)
  }
  result <- data.frame(Feature = character(), Term = character(), Estimate = numeric(),
                        SE = numeric(), Statistic = numeric(), P_value = numeric(), Q_value = numeric(),
                        Method = character(), N_obs = integer(), Status = character(), Error = character())
  if (length(terms)) {
    rows <- lapply(seq_along(tables), function(g) {
      table <- tables[[g]][match(terms, tables[[g]]$Term), , drop = FALSE]
      status <- if (diagnostics$Status[g] == "fit_error") rep("fit_error", length(terms)) else
        ifelse(is.finite(table$Estimate), "ok", "non_estimable")
      data.frame(Feature = names(data)[g], Term = terms, Estimate = table$Estimate, SE = table$SE,
                 Statistic = table$Statistic, P_value = table$P_value, Q_value = NA_real_, Method = method,
                 N_obs = diagnostics$N_obs[g], Status = status, Error = diagnostics$Error[g])
    })
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    for (term in terms) {
      idx <- which(result$Term == term & is.finite(result$P_value))
      result$Q_value[idx] <- stats::p.adjust(result$P_value[idx], method = "BH", n = ncol(data))
    }
  }
  attr(result, "feature_diagnostics") <- diagnostics
  attr(result, "formula") <- fit_formula
  result
}

#' Compare Fixed Effects Across Harmonization Methods
#'
#' Compares biological coefficient estimates, paired changes, or detection rates.
#'
#' @param raw_result Output of [extract_fixed_effects()] for raw data.
#' @param harmonized_result One extracted result or a nonempty named list of
#'   results for harmonization methods. List names identify the plotted methods.
#' @param type `"estimate"` for raw-versus-harmonized scatterplots with identity
#'   lines; `"change"` for boxplots of paired coefficient changes; or
#'   `"significance"` for proportions of features with `Q_value < alpha`.
#' @param alpha Significance cutoff, strictly between zero and one.
#'
#' @return A `ggplot` object. Estimate/change plot data retain all paired keys,
#'   including unavailable pairs, with `Raw_Estimate`, `Estimate`, and `Delta`.
#'   Significance plot data include `significant` (the valid-test proportion),
#'   `n_significant`, `n_valid`, `n_total`, and `n_failed` per term and method.
#' @details Raw data are labeled `"Raw"`; a single harmonized result uses its
#'   unique `Method` label when available, otherwise `"Harmonized"`. Significance
#'   comparisons use the union of feature/term keys across all methods so missing
#'   tests remain in the counts. A group with no valid tests has an unavailable
#'   proportion, not zero. Labels show significant/valid feature counts (`0/0`
#'   when no tests are valid). Total and failed counts remain in the plot data.
#'   Detection rates
#'   are descriptive: changes can reflect standard errors, sample sizes, and
#'   residual variance, as well as coefficient changes. Estimates may remain
#'   close to confounded raw estimates without being accurate. Ordinary p-values
#'   do not include harmonization uncertainty; see [harm_inference()].
#'
#'   Estimate/change plots match by feature and coefficient name. Only finite
#'   pairs enter the geometry. Compare matched subjects, feature families, units,
#'   formulas, and contrasts; in simulations also compare against known truth.
#'   The extractor excludes GAM smooth effects, which need curve comparisons.
#' @examples
#' raw <- data.frame(Feature = c("a", "b"), Term = "age",
#'                   Estimate = c(0.5, 0.8), Q_value = c(0.01, 0.2))
#' adjusted <- raw
#' adjusted$Estimate <- c(0.48, 0.78)
#' plot_fixed_effect_comparison(raw, list(CALM = adjusted), type = "estimate")
#' @export
plot_fixed_effect_comparison <- function(raw_result, harmonized_result,
                                         type = c("estimate", "change", "significance"), alpha = 0.05) {
  type <- match.arg(type)
  if (is.data.frame(harmonized_result)) {
    method <- if ("Method" %in% names(harmonized_result)) unique(as.character(harmonized_result$Method)) else character()
    if (length(method) != 1L || is.na(method) || !nzchar(method) || method == "Raw") method <- "Harmonized"
    harmonized_result <- stats::setNames(list(harmonized_result), method)
  }
  method_names <- names(harmonized_result)
  if (!is.list(harmonized_result) || !length(harmonized_result) || is.null(method_names) ||
      anyNA(method_names) || any(!nzchar(method_names)) || anyDuplicated(method_names) || "Raw" %in% method_names) {
    stop("Supply a nonempty named list with distinct harmonized method names other than `Raw`.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be between zero and one.", call. = FALSE)
  }
  all_results <- c(list(Raw = raw_result), harmonized_result)
  required <- c("Feature", "Term", "Estimate", "Q_value")
  prepared <- lapply(all_results, function(x) {
    if (!is.data.frame(x) || !all(required %in% names(x)) || !is.numeric(x$Estimate) || !is.numeric(x$Q_value)) {
      stop("All results need Feature, Term, numeric Estimate, and numeric Q_value.", call. = FALSE)
    }
    x <- x[, required, drop = FALSE]
    x$Feature <- as.character(x$Feature)
    x$Term <- as.character(x$Term)
    if (anyNA(x[, c("Feature", "Term")]) || any(!nzchar(x$Feature)) || any(!nzchar(x$Term))) {
      stop("Feature and term names must be nonmissing and nonempty.", call. = FALSE)
    }
    if (anyDuplicated(x[, c("Feature", "Term")])) stop("Duplicate Feature/Term keys.", call. = FALSE)
    if (any(is.finite(x$Q_value) & (x$Q_value < 0 | x$Q_value > 1))) stop("Finite Q_value entries must lie between zero and one.", call. = FALSE)
    x
  })
  keys <- unique(do.call(rbind, lapply(prepared, function(x) x[, c("Feature", "Term"), drop = FALSE])))
  if (!nrow(keys)) stop("No feature/term results are available to plot.", call. = FALSE)
  term_names <- unique(keys$Term)
  method_levels <- names(all_results)
  colors <- c("#4F86A6", "#66A69C", "#C77963", "#9B8FB5", "#C4A35A", "#7F9BA3")
  if (length(method_names) > length(colors)) colors <- grDevices::hcl.colors(length(method_names), "Dark 3")
  colors <- stats::setNames(c("#B8BDC5", colors[seq_along(method_names)]), method_levels)
  theme <- ggplot2::theme_classic(base_size = 13) + ggplot2::theme(
    strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")
  if (type == "significance") {
    long <- do.call(rbind, lapply(seq_along(prepared), function(i) {
      out <- merge(keys, prepared[[i]], by = c("Feature", "Term"), all.x = TRUE, sort = FALSE)
      out$Method <- rep(method_levels[i], nrow(out))
      out
    }))
    groups <- unique(long[, c("Term", "Method"), drop = FALSE])
    rows <- lapply(seq_len(nrow(groups)), function(i) {
      group <- groups[i, , drop = FALSE]
      q <- long$Q_value[long$Term == group$Term & long$Method == group$Method]
      valid <- is.finite(q)
      n_sig <- sum(q[valid] < alpha)
      data.frame(group, significant = if (any(valid)) n_sig / sum(valid) else NA_real_,
                 n_significant = n_sig, n_valid = sum(valid), n_total = length(q), n_failed = sum(!valid))
    })
    out <- if (length(rows)) do.call(rbind, rows) else data.frame(Term = character(), Method = character(),
      significant = numeric(), n_significant = integer(), n_valid = integer(), n_total = integer(), n_failed = integer())
    out$Method <- factor(out$Method, levels = method_levels)
    out$Term <- factor(out$Term, levels = term_names)
    out$count_label <- paste0(out$n_significant, "/", out$n_valid)
    out$label_y <- ifelse(is.finite(out$significant), out$significant, 0)
    return(ggplot2::ggplot(out, ggplot2::aes(x = .data[["Method"]], y = .data[["significant"]], fill = .data[["Method"]])) +
      ggplot2::geom_col(width = 0.65, na.rm = TRUE) +
      ggplot2::geom_text(ggplot2::aes(y = .data[["label_y"]], label = .data[["count_label"]]), vjust = -0.4, size = 3.5) +
      ggplot2::facet_wrap(~ Term, drop = FALSE) + ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.08),
                                 breaks = seq(0, 1, 0.25), expand = ggplot2::expansion(mult = 0)) +
      ggplot2::labs(x = NULL, y = "Proportion of significant features") + theme)
  }
  raw <- prepared[[1L]][, c("Feature", "Term", "Estimate"), drop = FALSE]
  names(raw)[3L] <- "Raw_Estimate"
  paired <- do.call(rbind, lapply(seq_along(harmonized_result), function(i) {
    out <- merge(raw, prepared[[i + 1L]], by = c("Feature", "Term"), all = TRUE, sort = FALSE)
    out$Delta <- out$Estimate - out$Raw_Estimate
    out$Method <- rep(method_names[i], nrow(out))
    out
  }))
  paired$Method <- factor(paired$Method, levels = method_names)
  paired$Term <- factor(paired$Term, levels = term_names)
  finite <- paired[is.finite(paired$Estimate) & is.finite(paired$Raw_Estimate) & is.finite(paired$Delta), , drop = FALSE]
  if (type == "estimate") {
    return(ggplot2::ggplot(paired, ggplot2::aes(x = .data[["Raw_Estimate"]], y = .data[["Estimate"]], color = .data[["Method"]])) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_point(data = finite, alpha = 0.5, size = 1.5) +
      ggplot2::facet_grid(Term ~ Method, scales = "free", drop = FALSE) +
      ggplot2::scale_color_manual(values = colors, drop = FALSE) +
      ggplot2::labs(x = "Raw coefficient estimate", y = "Harmonized coefficient estimate") + theme)
  }
  ggplot2::ggplot(paired, ggplot2::aes(x = .data[["Method"]], y = .data[["Delta"]], fill = .data[["Method"]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_boxplot(data = finite, width = 0.65, na.rm = TRUE, colour = "#35404A", outlier.alpha = 0.6) +
    ggplot2::facet_wrap(~ Term, scales = "free_y", drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Change in coefficient estimate") + theme
}
