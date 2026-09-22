#' Retrieve Posterior Harmonized-Data Draws
#'
#' Returns original-scale posterior harmonized-data draws retained by a
#' Bayesian LatentComBat fit.
#'
#' @param fit A fitted object returned by [com_harm()] with
#'   `method = "bayesian"` and `posterior_draws > 0`.
#'
#' @return A three-dimensional numeric array with dimensions
#'   posterior draw by observation by feature.
#'
#' @export
get_harm_draws <- function(fit) {
  if (!inherits(fit, "latentcombat_bayesian")) {
    stop("`get_harm_draws()` requires a Bayesian LatentComBat fit.", call. = FALSE)
  }

  draws <- fit$posterior_harmonization$harm_data_draws
  if (is.null(draws)) {
    stop(
      "No posterior harmonized-data draws are stored. Refit with `posterior_draws > 0`.",
      call. = FALSE
    )
  }

  draws
}


#' Combine Inference Across Bayesian Harmonization Draws
#'
#' Combines repeated downstream estimates obtained from posterior harmonized-
#' data draws using a law-of-total-variance decomposition. This is intended for
#' posterior Monte Carlo draws and therefore does not use the conventional
#' multiple-imputation `(1 + 1/S)` inflation factor.
#'
#' @param results Data frame containing one row per posterior draw and inferential
#'   target.
#' @param group_cols Character vector identifying the inferential target, such as
#'   `c("feature", "term")`.
#' @param estimate_col Name of the estimate column. Defaults to `"estimate"`.
#' @param se_col Name of the conditional standard-error column. Defaults to
#'   `"se"`.
#' @param conf_level Confidence level for normal-approximation intervals.
#'   Defaults to `0.95`.
#'
#' @return A data frame containing posterior-averaged estimates, within-
#'   harmonization variance, between-harmonization variance, total variance,
#'   propagated standard errors, confidence intervals, and two-sided p-values.
#'
#' @export
combine_harm_draw_inference <- function(
    results,
    group_cols = c("feature", "term"),
    estimate_col = "estimate",
    se_col = "se",
    conf_level = 0.95
) {
  if (!is.data.frame(results)) {
    stop("`results` must be a data frame.", call. = FALSE)
  }

  required <- c(group_cols, estimate_col, se_col)
  missing <- setdiff(required, names(results))
  if (length(missing) > 0L) {
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  split_key <- interaction(results[group_cols], drop = TRUE, lex.order = TRUE)
  pieces <- split(results, split_key)

  out <- lapply(pieces, function(x) {
    est <- as.numeric(x[[estimate_col]])
    se <- as.numeric(x[[se_col]])
    ok <- is.finite(est) & is.finite(se) & se >= 0
    est <- est[ok]
    se <- se[ok]

    base <- x[1L, group_cols, drop = FALSE]
    if (length(est) == 0L) {
      base$n_draws <- 0L
      base$estimate <- NA_real_
      base$within_var <- NA_real_
      base$between_var <- NA_real_
      base$total_var <- NA_real_
      base$se <- NA_real_
      base$statistic <- NA_real_
      base$p_value <- NA_real_
      base$conf_low <- NA_real_
      base$conf_high <- NA_real_
      return(base)
    }

    estimate <- mean(est)
    within_var <- mean(se^2)
    between_var <- if (length(est) > 1L) stats::var(est) else 0
    total_var <- within_var + between_var
    total_se <- sqrt(total_var)
    statistic <- if (is.finite(total_se) && total_se > 0) estimate / total_se else NA_real_

    base$n_draws <- length(est)
    base$estimate <- estimate
    base$within_var <- within_var
    base$between_var <- between_var
    base$total_var <- total_var
    base$se <- total_se
    base$statistic <- statistic
    base$p_value <- if (is.finite(statistic)) 2 * stats::pnorm(-abs(statistic)) else NA_real_
    base$conf_low <- estimate - zcrit * total_se
    base$conf_high <- estimate + zcrit * total_se
    base$harmonization_fraction <- if (is.finite(total_var) && total_var > 0) {
      between_var / total_var
    } else {
      NA_real_
    }
    base
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}


# Internal helper: fit the same downstream feature-wise model to every column.
.fit_feature_inference <- function(
    Y,
    covar,
    model = stats::lm,
    formula,
    terms = NULL,
    ...
) {
  Y <- as.matrix(Y)
  if (is.null(formula)) {stop("A downstream `formula` is required.", call. = FALSE)}
  if (is.null(covar)) {covar <- data.frame(row.names = seq_len(nrow(Y)))} else {covar <- as.data.frame(covar)}
  if (nrow(covar) != nrow(Y)) {stop("`covar` must have the same number of rows as the harmonized data.", call. = FALSE)}

  feature_names <- colnames(Y)
  if (is.null(feature_names)) feature_names <- paste0("feature_", seq_len(ncol(Y)))

  out <- vector("list", ncol(Y))

  for (g in seq_len(ncol(Y))) {
    dat <- data.frame(y = Y[, g], covar, check.names = FALSE)
    fit_g <- do.call(model, c(list(formula = formula, data = dat), list(...)))

    if (inherits(fit_g, "gam")) {
      tab <- summary(fit_g)$p.table
      if (is.null(tab) || nrow(tab) == 0L) {
        tmp <- data.frame(term = character(0), estimate = numeric(0), se = numeric(0))
      } else {
        tmp <- data.frame(
          term = rownames(tab),
          estimate = as.numeric(tab[, 1L]),
          se = as.numeric(tab[, 2L]),
          stringsAsFactors = FALSE
        )
      }
    } else if (inherits(fit_g, "merMod")) {
      est <- lme4::fixef(fit_g)
      vc <- stats::vcov(fit_g)
      tmp <- data.frame(
        term = names(est),
        estimate = as.numeric(est),
        se = sqrt(diag(as.matrix(vc))),
        stringsAsFactors = FALSE
      )
    } else if (inherits(fit_g, "lm")) {
      tab <- summary(fit_g)$coefficients
      tmp <- data.frame(
        term = rownames(tab),
        estimate = as.numeric(tab[, "Estimate"]),
        se = as.numeric(tab[, "Std. Error"]),
        stringsAsFactors = FALSE
      )
    } else {
      stop(
        "Unsupported downstream model class. Built-in inference currently supports `lm`, `gam`/`bam`, and `lmer` fixed effects.",
        call. = FALSE
      )
    }

    if (is.null(terms)) {
      tmp <- tmp[tmp$term != "(Intercept)", , drop = FALSE]
    } else {
      tmp <- tmp[tmp$term %in% terms, , drop = FALSE]
    }

    tmp$feature <- feature_names[g]
    tmp <- tmp[, c("feature", "term", "estimate", "se"), drop = FALSE]
    out[[g]] <- tmp
  }

  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}


.bootstrap_indices_within_batch <- function(bat) {
  bat <- droplevels(as.factor(bat))
  idx_by_batch <- split(seq_along(bat), bat)
  unlist(
    lapply(idx_by_batch, function(idx) sample(idx, length(idx), replace = TRUE)),
    use.names = FALSE
  )
}


.refit_sequential_bootstrap <- function(fit, idx) {
  info <- fit$inference_data
  if (is.null(info)) {
    stop("The sequential fit does not contain inference inputs. Refit with `keep_inference_data = TRUE`.", call. = FALSE)
  }

  data_b <- info$data[idx, , drop = FALSE]
  bat_b <- droplevels(as.factor(info$bat[idx]))
  covar_b <- if (is.null(info$covar)) NULL else info$covar[idx, , drop = FALSE]
  protected_b <- if (is.null(info$protected_covar)) {
    NULL
  } else {
    as.data.frame(info$protected_covar)[idx, , drop = FALSE]
  }

  args <- list(
    bat = bat_b,
    data = data_b,
    covar = covar_b,
    model = info$model,
    formula = info$formula,
    ref.batch = info$ref.batch,
    eb = info$eb,
    stan = info$stan,
    robust.LS = info$robust.LS,
    cov = info$cov,
    var_thresh = info$var_thresh,
    min_rblock = info$min_rblock,
    max_rblock = info$max_rblock,
    sva = info$sva,
    var_stable = info$var_stable,
    protected_covar = protected_b,
    keep_inference_data = FALSE
  )

  extra <- info$extra_args
  if (length(extra) > 0L) {
    # Preserve observation-specific weights under bootstrap resampling.
    if (!is.null(extra$weights)) {
      if (is.matrix(extra$weights) && nrow(extra$weights) == nrow(info$data)) {
        extra$weights <- extra$weights[idx, , drop = FALSE]
      } else if (!is.matrix(extra$weights) && length(extra$weights) == nrow(info$data)) {
        extra$weights <- extra$weights[idx]
      }
    }

    duplicate <- intersect(names(extra), names(args))
    if (length(duplicate) > 0L) extra[duplicate] <- NULL
    args <- c(args, extra)
  }

  do.call(com_harm.sequential, args)
}


#' Uncertainty-Aware Post-Harmonization Inference
#'
#' Fits feature-wise downstream models while propagating uncertainty from the
#' LatentComBat harmonization stage. Bayesian fits use retained posterior
#' harmonized-data draws. Sequential fits use a full-pipeline nonparametric
#' bootstrap that resamples observations within batch, reruns the complete
#' sequential harmonization workflow, and refits the downstream model.
#'
#' @param fit A fitted object returned by [com_harm()].
#' @param covar Downstream-analysis covariate data frame. If `NULL`, the
#'   covariates stored with the harmonization fit are used.
#' @param formula Downstream model formula written with response `y`, for
#'   example `y ~ age + sex + diagnosis`. If `NULL`, the formula stored with
#'   the harmonization fit is used.
#' @param model Downstream model-fitting function. If `NULL`, the model stored
#'   with the harmonization fit is used. Built-in extraction supports
#'   [stats::lm()], [mgcv::gam()]/[mgcv::bam()], and fixed effects from
#'   [lme4::lmer()].
#' @param terms Optional character vector of coefficient names to return. If
#'   `NULL`, all parametric/fixed-effect coefficients except the intercept are
#'   returned.
#' @param B Number of full-pipeline bootstrap replicates for sequential fits.
#'   Defaults to `200`. Ignored for Bayesian fits.
#' @param seed Random seed for the sequential bootstrap. Defaults to `1`.
#' @param conf_level Confidence level. Defaults to `0.95`.
#' @param progress Logical indicating whether bootstrap progress is printed.
#' @param ... Additional arguments passed to the downstream model-fitting
#'   function.
#'
#' @return A data frame with one row per feature and coefficient. For Bayesian
#'   fits, uncertainty is combined as mean conditional variance plus variance
#'   of the coefficient across posterior harmonization draws. For sequential
#'   fits, standard errors are the empirical standard deviations of estimates
#'   from the full-pipeline bootstrap; percentile confidence intervals are also
#'   returned.
#'
#' @details
#' The Bayesian calculation is posterior propagation of harmonization
#' uncertainty, not a fully Bayesian downstream model. The sequential
#' calculation is a subject/observation-level bootstrap and currently assumes
#' independent observational units. For longitudinal repeated-measures
#' analyses, users should use a cluster-aware resampling strategy rather than
#' the built-in sequential bootstrap.
#'
#' Smooth terms from GAMs are not reduced to a single coefficient by this
#' helper; only the parametric coefficient table is returned. Fixed effects are
#' returned for mixed models.
#'
#' @export
harm_inference <- function(
    fit,
    covar = NULL,
    formula = NULL,
    model = NULL,
    terms = NULL,
    B = 200L,
    seed = 1L,
    conf_level = 0.95,
    progress = TRUE,
    ...
) {
  if (!inherits(fit, "latentcombat_fit")) {
    stop("`fit` must inherit from `latentcombat_fit`.", call. = FALSE)
  }

  stored <- fit$inference_data
  if (is.null(covar)) {
    if (is.null(stored)) stop("Supply downstream `covar`.", call. = FALSE)
    covar <- stored$covar
  }
  if (is.null(formula)) {
    if (is.null(stored) || is.null(stored$formula)) {
      stop("Supply a downstream `formula`.", call. = FALSE)
    }
    formula <- stored$formula
  }
  if (is.null(model)) {
    model <- if (!is.null(stored) && !is.null(stored$model)) stored$model else stats::lm
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  if (inherits(fit, "latentcombat_bayesian")) {
    Y_draws <- get_harm_draws(fit)
    S <- dim(Y_draws)[1L]

    draw_results <- vector("list", S)
    for (s in seq_len(S)) {
      draw_results[[s]] <- .fit_feature_inference(
        Y = Y_draws[s, , ],
        covar = covar,
        model = model,
        formula = formula,
        terms = terms,
        ...
      )
      draw_results[[s]]$draw <- s
    }

    draw_results <- do.call(rbind, draw_results)
    combined <- combine_harm_draw_inference(
      draw_results,
      group_cols = c("feature", "term"),
      estimate_col = "estimate",
      se_col = "se",
      conf_level = conf_level
    )

    plugin <- .fit_feature_inference(
      Y = fit$harm_data,
      covar = covar,
      model = model,
      formula = formula,
      terms = terms,
      ...
    )
    names(plugin)[names(plugin) == "estimate"] <- "plugin_estimate"
    names(plugin)[names(plugin) == "se"] <- "plugin_se"

    out <- merge(combined, plugin, by = c("feature", "term"), all.x = TRUE, sort = FALSE)
    out$inference_method <- "posterior_propagation"
    out$harmonization_method <- "bayesian"
    class(out) <- c("latentcombat_inference", "data.frame")
    return(out)
  }

  if (!inherits(fit, "latentcombat_sequential")) {
    stop("Unsupported LatentComBat fit class.", call. = FALSE)
  }
  if (is.null(stored)) {
    stop(
      "Sequential uncertainty propagation requires stored original inputs. Refit with `keep_inference_data = TRUE`.",
      call. = FALSE
    )
  }

  formula_text <- paste(deparse(formula), collapse = " ")
  if (identical(model, lme4::lmer) || grepl("\\|", formula_text)) {
    stop(
      "The built-in sequential bootstrap currently assumes independent observational units and does not bootstrap repeated-measures/mixed-model data. Use a cluster-aware bootstrap for longitudinal inference.",
      call. = FALSE
    )
  }

  B <- as.integer(B)
  if (length(B) != 1L || is.na(B) || B < 2L) {
    stop("`B` must be an integer of at least 2 for sequential inference.", call. = FALSE)
  }

  point <- .fit_feature_inference(
    Y = fit$harm_data,
    covar = covar,
    model = model,
    formula = formula,
    terms = terms,
    ...
  )

  set.seed(as.integer(seed))
  boot_results <- vector("list", B)

  for (b in seq_len(B)) {
    idx <- .bootstrap_indices_within_batch(stored$bat)
    boot_fit <- .refit_sequential_bootstrap(fit, idx)
    covar_b <- if (is.null(covar)) NULL else as.data.frame(covar)[idx, , drop = FALSE]

    boot_results[[b]] <- .fit_feature_inference(
      Y = boot_fit$harm_data,
      covar = covar_b,
      model = model,
      formula = formula,
      terms = terms,
      ...
    )
    boot_results[[b]]$bootstrap <- b

    if (isTRUE(progress) && (b == 1L || b == B || b %% max(1L, floor(B / 10L)) == 0L)) {
      message("Sequential inference bootstrap: ", b, "/", B)
    }
  }

  boots <- do.call(rbind, boot_results)
  key <- interaction(boots[c("feature", "term")], drop = TRUE, lex.order = TRUE)
  pieces <- split(boots, key)
  alpha <- 1 - conf_level
  zcrit <- stats::qnorm(1 - alpha / 2)

  boot_summary <- lapply(pieces, function(x) {
    vals <- x$estimate[is.finite(x$estimate)]
    data.frame(
      feature = x$feature[1L],
      term = x$term[1L],
      n_boot = length(vals),
      bootstrap_mean = if (length(vals)) mean(vals) else NA_real_,
      bootstrap_se = if (length(vals) > 1L) stats::sd(vals) else NA_real_,
      percentile_low = if (length(vals)) as.numeric(stats::quantile(vals, alpha / 2, na.rm = TRUE, names = FALSE)) else NA_real_,
      percentile_high = if (length(vals)) as.numeric(stats::quantile(vals, 1 - alpha / 2, na.rm = TRUE, names = FALSE)) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  boot_summary <- do.call(rbind, boot_summary)

  out <- merge(point, boot_summary, by = c("feature", "term"), all.x = TRUE, sort = FALSE)
  out$conditional_se <- out$se
  out$se <- out$bootstrap_se
  out$bootstrap_bias <- out$bootstrap_mean - out$estimate
  out$statistic <- ifelse(is.finite(out$se) & out$se > 0, out$estimate / out$se, NA_real_)
  out$p_value <- 2 * stats::pnorm(-abs(out$statistic))
  out$conf_low <- out$percentile_low
  out$conf_high <- out$percentile_high
  out$normal_conf_low <- out$estimate - zcrit * out$se
  out$normal_conf_high <- out$estimate + zcrit * out$se
  out$inference_method <- "full_pipeline_bootstrap"
  out$harmonization_method <- "sequential"

  class(out) <- c("latentcombat_inference", "data.frame")
  out
}


#' Print LatentComBat Inference Results
#'
#' @param x An object returned by [harm_inference()].
#' @param ... Additional arguments passed to `print.data.frame`.
#'
#' @return `x`, invisibly.
#'
#' @method print latentcombat_inference
#' @export
print.latentcombat_inference <- function(x, ...) {
  method <- unique(x$inference_method)
  if (length(method) == 1L) {
    cat("LatentComBat uncertainty-aware inference\n")
    cat("  Method: ", method, "\n", sep = "")
  }
  NextMethod("print", x, ...)
  invisible(x)
}
