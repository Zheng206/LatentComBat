# Lightweight Bayesian convergence diagnostics for LatentComBat
#
# Suggested package location: R/bayesian_diagnostics.R
#
# Design goals:
#   - evaluate whether Stan sampled the fitted posterior reliably;
#   - support both Joint Bayesian LatentComBat fits and univariate
#     Bayesian ComBat fits used in the sequential harmonization workflow;
#   - assess sampler geometry using divergences, maximum treedepth hits,
#     and E-BFMI;
#   - assess convergence of identifiable / scientifically relevant parameter
#     families using R-hat and effective sample size;
#   - automatically detect whether the supplied fit is a joint latent model
#     or a univariate Stan model and select appropriate diagnostic parameters;
#   - accept a high-level LatentComBat fit object, including
#     `latentcombat_bayesian` and `latentcombat_sequential` objects, as well
#     as a raw CmdStanMCMC object;
#   - for latent-factor models, deliberately exclude element-wise latent
#     coordinates such as H, Lambda, kappa, and factor-specific scale
#     parameters from the overall convergence decision because latent-factor
#     sign / permutation / rotation ambiguity can make their element-wise
#     R-hat and ESS misleading even when the fitted latent surface is stable;
#   - for univariate Bayesian ComBat models, directly assess identifiable
#     hierarchical parameters such as gamma_i, sigma_psi, mu_ig, and
#     sigma_delta.
#
# Main API:
#   diag <- bayes_diagnostics(fit)
#   diag
#   summary(diag)
#   plot(diag, type = "rhat")
#   plot(diag, type = "ess")
#   plot(diag, type = "trace")
#
# Examples:
#
#   # Joint Bayesian LatentComBat
#   diag <- bayes_diagnostics(fit_bayes)
#
#   # Sequential LatentComBat using the univariate Stan engine
#   diag <- bayes_diagnostics(harm_result)
#
#   # Raw CmdStanR fit
#   diag <- bayes_diagnostics(cmdstan_fit)
#
# This file provides lightweight MCMC convergence diagnostics only.
# More computationally intensive checks of derived quantities such as
# latent_hat, latent_remove, R_adj, posterior predictive quantities, and
# harmonization stability should be implemented separately.

.get_bayes_stan_fit <- function(x) {
  if (inherits(x, "CmdStanMCMC")) {
    return(x)
  }

  if (
    inherits(x, "latentcombat_bayesian") &&
    !is.null(x$stan_result$fit) &&
    inherits(x$stan_result$fit, "CmdStanMCMC")
  ) {
    return(x$stan_result$fit)
  }

  if (
    inherits(x, "latentcombat_sequential") &&
    !is.null(x$stan_result$fit) &&
    inherits(x$stan_result$fit, "CmdStanMCMC")
  ) {
    return(x$stan_result$fit)
  }

  if (
    is.list(x) &&
    !is.null(x$stan_result) &&
    is.list(x$stan_result) &&
    !is.null(x$stan_result$fit) &&
    inherits(x$stan_result$fit, "CmdStanMCMC")
  ) {
    return(x$stan_result$fit)
  }

  if (
    is.list(x) &&
    !is.null(x$fit) &&
    inherits(x$fit, "CmdStanMCMC")
  ) {
    return(x$fit)
  }

  stop(
    paste0(
      "Could not locate a CmdStanMCMC object in `fit`. ",
      "Expected a raw CmdStanMCMC object, a `latentcombat_bayesian` fit, ",
      "a `latentcombat_sequential` fit with `stan_result$fit`, ",
      "or a list containing `$fit` or `$stan_result$fit`."
    ),
    call. = FALSE
  )
}


.stan_has_variable <- function(stan_fit, variable) {
  ok <- tryCatch({
    z <- stan_fit$summary(variables = variable)
    is.data.frame(z) && nrow(z) > 0L
  }, error = function(e) FALSE)

  isTRUE(ok)
}


.detect_bayes_model <- function(stan_fit) {
  if (
    .stan_has_variable(stan_fit, "mu_g") ||
    .stan_has_variable(stan_fit, "alpha")
  ) {
    return("latent")
  }

  if (
    .stan_has_variable(stan_fit, "gamma_i") &&
    .stan_has_variable(stan_fit, "mu_ig") &&
    .stan_has_variable(stan_fit, "sigma_delta")
  ) {
    return("univariate")
  }

  "unknown"
}


.default_core_variables <- function(stan_fit, model_type = NULL) {
  if (is.null(model_type)) {
    model_type <- .detect_bayes_model(stan_fit)
  }

  if (identical(model_type, "latent")) {
    candidates <- c(
      "mu_g",
      "tau_b",
      "sigma0",
      "sigma_logdelta_b",
      # Longitudinal analogues, when present
      "alpha",
      "tau_alpha",
      "sigma_u"
    )
  } else if (identical(model_type, "univariate")) {
    candidates <- c(
      "gamma_i",
      "sigma_psi",
      "mu_ig",
      "sigma_delta"
    )
  } else {
    stop("Unable to determine the Bayesian Stan model type.", call. = FALSE)
  }

  keep <- vapply(
    candidates,
    function(v) .stan_has_variable(stan_fit, v),
    logical(1)
  )

  candidates[keep]
}

.classify_core_parameter <- function(x) {
  out <- rep("Other", length(x))

  # Joint Bayesian LatentComBat
  out[grepl("^mu_g\\[|^alpha\\[", x)] <- "Mean structure"
  out[grepl("^tau_b\\[|^tau_g\\[|^tau_alpha\\[", x)] <- "Batch / subject scale"
  out[grepl("^sigma0\\[", x)] <- "Residual scale"
  out[grepl("^sigma_logdelta_b\\[|^sigma_logdelta$", x)] <- "Batch-scale heterogeneity"
  out[grepl("^sigma_u\\[|^sigma_u$", x)] <- "Longitudinal latent scale"

  # Univariate Bayesian ComBat
  out[grepl("^gamma_i\\[", x)] <- "Batch mean"
  out[grepl("^sigma_psi\\[", x)] <- "Batch-effect heterogeneity"
  out[grepl("^mu_ig\\[", x)] <- "Feature-by-batch mean"
  out[grepl("^sigma_delta\\[", x)] <- "Feature-by-batch residual scale"

  factor(
    out,
    levels = c(
      "Mean structure",
      "Batch / subject scale",
      "Residual scale",
      "Batch-scale heterogeneity",
      "Longitudinal latent scale",
      "Batch mean",
      "Batch-effect heterogeneity",
      "Feature-by-batch mean",
      "Feature-by-batch residual scale",
      "Other"
    )
  )
}


.safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}


.safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}


.safe_prop <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}


.status_rank <- function(x) {
  match(x, c("PASS", "WARNING", "FAIL"))
}


.worst_status <- function(...) {
  x <- c(...)
  x <- x[!is.na(x)]
  if (length(x) == 0L) return("WARNING")
  c("PASS", "WARNING", "FAIL")[max(.status_rank(x))]
}


.extract_sampler_diagnostics <- function(stan_fit, max_treedepth = 12L) {
  d <- stan_fit$sampler_diagnostics(format = "draws_array")

  if (length(dim(d)) != 3L) {
    stop("Unexpected sampler-diagnostics format returned by CmdStanR.", call. = FALSE)
  }

  n_iter <- dim(d)[1L]
  n_chain <- dim(d)[2L]
  diag_names <- dimnames(d)[[3L]]

  get_diag <- function(name) {
    j <- match(name, diag_names)
    if (is.na(j)) return(NULL)
    d[, , j, drop = TRUE]
  }

  divergent <- get_diag("divergent__")
  treedepth <- get_diag("treedepth__")
  energy <- get_diag("energy__")

  chain <- lapply(seq_len(n_chain), function(ch) {
    div_ch <- if (is.null(divergent)) {
      NA_real_
    } else {
      sum(divergent[, ch] > 0, na.rm = TRUE)
    }

    tree_ch <- if (is.null(treedepth)) {
      NA_real_
    } else {
      sum(treedepth[, ch] >= max_treedepth, na.rm = TRUE)
    }

    ebfmi_ch <- NA_real_
    if (!is.null(energy)) {
      e <- as.numeric(energy[, ch])
      e <- e[is.finite(e)]
      if (length(e) >= 3L && stats::var(e) > 0) {
        ebfmi_ch <- mean(diff(e)^2) / stats::var(e)
      }
    }

    data.frame(
      chain = ch,
      n_draws = n_iter,
      divergences = div_ch,
      treedepth_hits = tree_ch,
      ebfmi = ebfmi_ch,
      stringsAsFactors = FALSE
    )
  })

  chain <- do.call(rbind, chain)

  list(
    chain = chain,
    n_iter = n_iter,
    n_chains = n_chain,
    n_draws_total = n_iter * n_chain,
    divergences = if (all(is.na(chain$divergences))) {
      NA_real_
    } else {
      sum(chain$divergences, na.rm = TRUE)
    },
    treedepth_hits = if (all(is.na(chain$treedepth_hits))) {
      NA_real_
    } else {
      sum(chain$treedepth_hits, na.rm = TRUE)
    },
    ebfmi = chain$ebfmi,
    min_ebfmi = .safe_min(chain$ebfmi),
    divergence_prop = if (n_iter * n_chain > 0) {
      sum(chain$divergences, na.rm = TRUE) / (n_iter * n_chain)
    } else {
      NA_real_
    },
    treedepth_prop = if (n_iter * n_chain > 0) {
      sum(chain$treedepth_hits, na.rm = TRUE) / (n_iter * n_chain)
    } else {
      NA_real_
    },
    max_treedepth = as.integer(max_treedepth)
  )
}


.sampler_status <- function(
    sampler,
    ebfmi_threshold = 0.30,
    ebfmi_fail_threshold = 0.20,
    divergence_fail_prop = 0.01,
    treedepth_warning_prop = 0.01,
    treedepth_fail_prop = 0.05
) {
  n_draws <- sampler$n_draws_total

  divergence_prop <- if (is.finite(n_draws) && n_draws > 0) {
    sampler$divergences / n_draws
  } else {
    NA_real_
  }

  treedepth_prop <- if (is.finite(n_draws) && n_draws > 0) {
    sampler$treedepth_hits / n_draws
  } else {
    NA_real_
  }

  min_ebfmi <- sampler$min_ebfmi

  if (
    (is.finite(divergence_prop) && divergence_prop > divergence_fail_prop) ||
    (is.finite(treedepth_prop) && treedepth_prop > treedepth_fail_prop) ||
    (is.finite(min_ebfmi) && min_ebfmi < ebfmi_fail_threshold)
  ) {
    return("FAIL")
  }

  if (
    (!is.na(sampler$divergences) && sampler$divergences > 0L) ||
    (is.finite(treedepth_prop) && treedepth_prop > treedepth_warning_prop) ||
    (is.finite(min_ebfmi) && min_ebfmi < ebfmi_threshold)
  ) {
    return("WARNING")
  }

  "PASS"
}


.parameter_status <- function(
    max_rhat,
    min_bulk_ess,
    min_tail_ess,
    rhat_threshold = 1.01,
    ess_threshold = 100
) {
  if ((!is.na(max_rhat) && max_rhat > 1.05) ||
      (!is.na(min_bulk_ess) && min_bulk_ess < 50) ||
      (!is.na(min_tail_ess) && min_tail_ess < 50)) {
    return("FAIL")
  }

  if ((!is.na(max_rhat) && max_rhat > rhat_threshold) ||
      (!is.na(min_bulk_ess) && min_bulk_ess < ess_threshold) ||
      (!is.na(min_tail_ess) && min_tail_ess < ess_threshold)) {
    return("WARNING")
  }

  "PASS"
}


# -----------------------------------------------------------------------------
# Main constructor
# -----------------------------------------------------------------------------

#' Lightweight Diagnostics for Bayesian LatentComBat / Univariate Stan Fits
#'
#' Evaluates the most important MCMC diagnostics for a Joint Bayesian
#' LatentComBat fit while deliberately excluding element-wise latent-factor
#' coordinates such as `H` and `Lambda` from the overall convergence decision.
#' Those coordinates can exhibit sign, permutation, or rotation ambiguity even
#' when the fitted latent surface is stable.
#'
#' The diagnostic has two components:
#' \itemize{
#'   \item sampler geometry: divergences, treedepth hits, and E-BFMI;
#'   \item convergence of a small set of identifiable / structural parameter
#'   families using R-hat and effective sample size.
#' }
#'
#' @param fit A fitted `latentcombat_bayesian` object, a univariate Stan result
#'   containing a CmdStanR `CmdStanMCMC` object in `$fit`, or a raw
#'   `CmdStanMCMC` object.
#' @param variables Optional character vector of Stan base variable names to
#'   monitor. If `NULL`, a lightweight model-aware default is used.
#' @param rhat_threshold R-hat warning threshold. Defaults to `1.01`.
#' @param ess_threshold Minimum bulk/tail ESS used as a warning threshold.
#'   Defaults to `100`.
#' @param ebfmi_threshold E-BFMI warning threshold. Defaults to `0.30`.
#' @param ebfmi_fail_threshold E-BFMI failure threshold. Defaults to `0.20`.
#' @param divergence_fail_prop Proportion of post-warmup draws with divergences
#'   above which sampler geometry is classified as `FAIL`. Defaults to `0.01`.
#' @param treedepth_warning_prop Proportion of post-warmup draws reaching the
#'   maximum tree depth above which sampler geometry is classified as `WARNING`.
#'   Defaults to `0.01`.
#' @param treedepth_fail_prop Proportion of post-warmup draws reaching the
#'   maximum tree depth above which sampler geometry is classified as `FAIL`.
#'   Defaults to `0.05`.
#' @param max_treedepth Maximum NUTS tree depth used when counting treedepth
#'   hits. Defaults to `12`.
#'
#' @return An object of class `latentcombat_bayes_diagnostics`.
#'
#' @export
bayes_diagnostics <- function(
    fit,
    variables = NULL,
    rhat_threshold = 1.01,
    ess_threshold = 100,
    ebfmi_threshold = 0.30,
    ebfmi_fail_threshold = 0.20,
    divergence_fail_prop = 0.01,
    treedepth_warning_prop = 0.01,
    treedepth_fail_prop = 0.05,
    max_treedepth = 12L
) {
  stan_fit <- .get_bayes_stan_fit(fit)
  model_type <- .detect_bayes_model(stan_fit)

  if (identical(model_type, "unknown")) {
    stop(
      "Unable to determine whether the Stan fit is latent or univariate.",
      call. = FALSE
    )
  }

  if (is.null(variables)) {
    variables <- .default_core_variables(
      stan_fit,
      model_type = model_type
    )
  }

  if (length(variables) == 0L) {
    stop("No core diagnostic variables were found in the Stan fit.", call. = FALSE)
  }

  smry <- stan_fit$summary(variables = variables)

  needed <- c("variable", "rhat", "ess_bulk", "ess_tail")
  missing_cols <- setdiff(needed, names(smry))
  if (length(missing_cols) > 0L) {
    stop(
      "CmdStanR summary is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  smry$group <- .classify_core_parameter(smry$variable)
  smry$rhat_flag <- is.finite(smry$rhat) & smry$rhat > rhat_threshold
  smry$ess_flag <-
    (is.finite(smry$ess_bulk) & smry$ess_bulk < ess_threshold) |
    (is.finite(smry$ess_tail) & smry$ess_tail < ess_threshold)

  sampler <- .extract_sampler_diagnostics(
    stan_fit,
    max_treedepth = as.integer(max_treedepth)
  )

  max_rhat <- .safe_max(smry$rhat)
  min_bulk_ess <- .safe_min(smry$ess_bulk)
  min_tail_ess <- .safe_min(smry$ess_tail)

  sampling_status <- .sampler_status(
    sampler,
    ebfmi_threshold = ebfmi_threshold,
    ebfmi_fail_threshold = ebfmi_fail_threshold,
    divergence_fail_prop = divergence_fail_prop,
    treedepth_warning_prop = treedepth_warning_prop,
    treedepth_fail_prop = treedepth_fail_prop
  )

  parameter_status <- .parameter_status(
    max_rhat = max_rhat,
    min_bulk_ess = min_bulk_ess,
    min_tail_ess = min_tail_ess,
    rhat_threshold = rhat_threshold,
    ess_threshold = ess_threshold
  )

  overall_status <- .worst_status(sampling_status, parameter_status)

  latent <- list(K_target = NULL, K_bayes = NULL)
  if (
    identical(model_type, "latent") &&
    inherits(fit, "latentcombat_bayesian")
  ) {
    latent$K_target <- fit$latent$K_target
    latent$K_bayes <- fit$latent$K_bayes
  }

  out <- list(
    model_type = model_type,
    status = overall_status,
    sampling_status = sampling_status,
    parameter_status = parameter_status,
    parameter_summary = smry,
    sampler = sampler,
    convergence = list(
      max_rhat = max_rhat,
      n_rhat_bad = sum(smry$rhat_flag, na.rm = TRUE),
      prop_rhat_bad = .safe_prop(smry$rhat_flag),
      min_bulk_ess = min_bulk_ess,
      min_tail_ess = min_tail_ess,
      n_ess_bad = sum(smry$ess_flag, na.rm = TRUE),
      prop_ess_bad = .safe_prop(smry$ess_flag)
    ),
    thresholds = list(
      rhat = rhat_threshold,
      ess = ess_threshold,
      ebfmi = ebfmi_threshold,
      ebfmi_fail = ebfmi_fail_threshold,
      divergence_fail_prop = divergence_fail_prop,
      treedepth_warning_prop = treedepth_warning_prop,
      treedepth_fail_prop = treedepth_fail_prop,
      max_treedepth = as.integer(max_treedepth)
    ),
    latent = latent,
    variables = variables,
    stan_fit = stan_fit
  )

  class(out) <- "latentcombat_bayes_diagnostics"
  out
}


# -----------------------------------------------------------------------------
# Print / summary
# -----------------------------------------------------------------------------

#' @param x An object returned by [bayes_diagnostics()].
#' @param ... Additional arguments passed to or from methods.
#'
#' @rdname bayes_diagnostics
#' @method print latentcombat_bayes_diagnostics
#' @export
print.latentcombat_bayes_diagnostics <- function(x, ...) {
  title <- switch(
    x$model_type,
    latent = "Joint Bayesian LatentComBat diagnostics",
    univariate = "Univariate Bayesian ComBat diagnostics",
    "Bayesian Stan diagnostics"
  )

  cat(title, "
")
  cat(strrep("-", nchar(title)), "
", sep = "")
  cat("Overall status:              ", x$status, "
", sep = "")
  cat("  Sampler geometry:          ", x$sampling_status, "
", sep = "")
  cat("  Core parameter convergence:", x$parameter_status, "

", sep = "")

  cat("Sampling
")
  cat("  Chains:                    ", x$sampler$n_chains, "
", sep = "")
  cat("  Post-warmup draws:         ", x$sampler$n_draws_total, "
", sep = "")
  cat(
    "  Divergences:               ", x$sampler$divergences,
    " (", sprintf("%.2f%%", 100 * x$sampler$divergence_prop), ")
",
    sep = ""
  )
  cat(
    "  Max treedepth hits:        ", x$sampler$treedepth_hits,
    " (", sprintf("%.2f%%", 100 * x$sampler$treedepth_prop), ")
",
    sep = ""
  )
  cat(
    "  E-BFMI by chain:           ",
    paste(sprintf("%.2f", x$sampler$ebfmi), collapse = ", "),
    "

",
    sep = ""
  )

  cat("Core parameter convergence
")
  cat("  Maximum R-hat:             ", sprintf("%.3f", x$convergence$max_rhat), "
", sep = "")
  cat(
    "  Above ", x$thresholds$rhat, ":                ",
    x$convergence$n_rhat_bad, " / ", nrow(x$parameter_summary), "
",
    sep = ""
  )
  cat("  Minimum bulk ESS:          ", round(x$convergence$min_bulk_ess), "
", sep = "")
  cat("  Minimum tail ESS:          ", round(x$convergence$min_tail_ess), "
", sep = "")

  if (!is.null(x$latent$K_bayes)) {
    cat("
Latent dimension
")
    cat("  K target:                  ", x$latent$K_target, "
", sep = "")
    cat("  K fitted:                  ", x$latent$K_bayes, "
", sep = "")
  }

  if (identical(x$model_type, "latent")) {
    cat("
Note: element-wise H, Lambda, kappa, and latent-factor-specific
")
    cat("scale parameters are excluded from the overall status because latent
")
    cat("factor sign/permutation/rotation ambiguity can inflate their R-hat.
")
  }

  invisible(x)
}


#' @param object An object returned by [bayes_diagnostics()].
#'
#' @rdname bayes_diagnostics
#' @method summary latentcombat_bayes_diagnostics
#' @export
summary.latentcombat_bayes_diagnostics <- function(object, ...) {
  group_summary <- lapply(
    split(object$parameter_summary, object$parameter_summary$group, drop = TRUE),
    function(z) {
      data.frame(
        group = as.character(z$group[1L]),
        n_parameters = nrow(z),
        max_rhat = .safe_max(z$rhat),
        prop_rhat_bad = .safe_prop(z$rhat_flag),
        min_bulk_ess = .safe_min(z$ess_bulk),
        min_tail_ess = .safe_min(z$ess_tail),
        prop_ess_bad = .safe_prop(z$ess_flag),
        stringsAsFactors = FALSE
      )
    }
  )

  group_summary <- do.call(rbind, group_summary)
  rownames(group_summary) <- NULL

  out <- list(
    model_type = object$model_type,
    status = object$status,
    sampling_status = object$sampling_status,
    parameter_status = object$parameter_status,
    group_summary = group_summary,
    sampler = object$sampler,
    thresholds = object$thresholds
  )

  class(out) <- "summary_latentcombat_bayes_diagnostics"
  out
}


#' @rdname bayes_diagnostics
#' @method print summary_latentcombat_bayes_diagnostics
#' @export
print.summary_latentcombat_bayes_diagnostics <- function(x, ...) {
  title <- switch(
    x$model_type,
    latent = "Joint Bayesian diagnostic summary",
    univariate = "Univariate Bayesian diagnostic summary",
    "Bayesian diagnostic summary"
  )

  cat(title, "
")
  cat("Overall status:               ", x$status, "
", sep = "")
  cat("Sampler geometry:             ", x$sampling_status, "
", sep = "")
  cat("Core parameter convergence:   ", x$parameter_status, "

", sep = "")
  print(x$group_summary, row.names = FALSE)
  invisible(x)
}

# Plot helpers
.plot_bayes_rhat <- function(x) {
  df <- x$parameter_summary
  df <- df[is.finite(df$rhat), , drop = FALSE]

  subtitle <- if (identical(x$model_type, "latent")) {
    "Core identifiable parameters; latent-factor coordinates excluded"
  } else {
    "Core parameters from the univariate Bayesian model"
  }

  df$status <- ifelse(df$rhat > x$thresholds$rhat, "Above threshold", "Within threshold")

  rhat_upper <- max(c(df$rhat, x$thresholds$rhat), na.rm = TRUE)

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["rhat"]],
      y = .data[["group"]]
    )
  ) +
    ggplot2::annotate(
      "rect",
      xmin = x$thresholds$rhat,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.05
    ) +
    ggplot2::geom_vline(
      xintercept = x$thresholds$rhat,
      linetype = 2,
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = .data[["status"]],
        alpha = .data[["status"]]
      ),
      size = 2.4,
      position = ggplot2::position_jitter(
        height = 0.09,
        width = 0
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Within threshold" = 16,
        "Above threshold" = 17
      )
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        "Within threshold" = 0.45,
        "Above threshold" = 1
      )
    ) +
    ggplot2::labs(
      title = "R-hat convergence diagnostic",
      subtitle = subtitle,
      x = expression(hat(R)),
      y = NULL,
      shape = NULL,
      alpha = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(color = "grey25"),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = ggplot2::element_text(
        color = "grey35",
        size = 10.5
      ),
      legend.position = "top"
    )
}


.plot_bayes_ess <- function(x) {
  df <- x$parameter_summary

  subtitle <- if (identical(x$model_type, "latent")) {
    "Core identifiable parameters; latent-factor coordinates excluded"
  } else {
    "Core parameters from the univariate Bayesian model"
  }

  df_long <- rbind(
    data.frame(
      variable = df$variable,
      group = df$group,
      ess = df$ess_bulk,
      type = "Bulk ESS",
      stringsAsFactors = FALSE
    ),
    data.frame(
      variable = df$variable,
      group = df$group,
      ess = df$ess_tail,
      type = "Tail ESS",
      stringsAsFactors = FALSE
    )
  )

  df_long <- df_long[
    is.finite(df_long$ess),
    ,
    drop = FALSE
  ]

  df_long$status <- ifelse(
    df_long$ess < x$thresholds$ess,
    "Below threshold",
    "Within threshold"
  )

  ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = .data[["ess"]],
      y = .data[["group"]]
    )
  ) +
    ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = x$thresholds$ess,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.05
    ) +
    ggplot2::geom_vline(
      xintercept = x$thresholds$ess,
      linetype = 2,
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = .data[["status"]],
        alpha = .data[["status"]]
      ),
      size = 2.4,
      position = ggplot2::position_jitter(
        height = 0.09,
        width = 0
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "Within threshold" = 16,
        "Below threshold" = 17
      )
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        "Within threshold" = 0.45,
        "Below threshold" = 1
      )
    ) +
    ggplot2::facet_wrap(
      ~type,
      ncol = 1
    ) +
    ggplot2::labs(
      title = "Effective sample size diagnostic",
      subtitle = subtitle,
      x = "Effective sample size",
      y = NULL,
      shape = NULL,
      alpha = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(color = "grey25"),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = ggplot2::element_text(
        color = "grey35",
        size = 10.5
      ),
      legend.position = "top"
    )
}


.select_trace_variables <- function(
    x,
    variables = NULL,
    n = 4L
) {
  if (!is.null(variables)) {
    return(as.character(variables))
  }

  smry <- x$parameter_summary
  n <- max(1L, as.integer(n))

  worst_ess <- pmin(
    smry$ess_bulk,
    smry$ess_tail,
    na.rm = TRUE
  )

  ord <- order(
    -ifelse(
      is.finite(smry$rhat),
      smry$rhat,
      -Inf
    ),
    ifelse(
      is.finite(worst_ess),
      worst_ess,
      Inf
    )
  )

  unique(smry$variable[ord])[
    seq_len(min(n, nrow(smry)))
  ]
}


.plot_bayes_trace <- function(
    x,
    variables = NULL,
    n = 4L
) {
  auto_selected <- is.null(variables)

  vars <- .select_trace_variables(
    x,
    variables = variables,
    n = n
  )

  d <- x$stan_fit$draws(
    variables = vars,
    format = "draws_array"
  )

  iter <- dim(d)[1L]
  chains <- dim(d)[2L]
  vnames <- dimnames(d)[[3L]]

  pieces <- vector(
    "list",
    length(vnames) * chains
  )

  z <- 1L

  for (j in seq_along(vnames)) {
    for (ch in seq_len(chains)) {
      pieces[[z]] <- data.frame(
        iteration = seq_len(iter),
        chain = factor(ch),
        variable = vnames[j],
        value = as.numeric(d[, ch, j]),
        stringsAsFactors = FALSE
      )

      z <- z + 1L
    }
  }

  df <- do.call(rbind, pieces)

  subtitle <- if (auto_selected) {
    "Parameters prioritized by high R-hat and low effective sample size"
  } else {
    "User-selected parameters"
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[["iteration"]],
      y = .data[["value"]],
      color = .data[["chain"]]
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.4,
      alpha = 0.75
    ) +
    ggplot2::facet_wrap(
      ~variable,
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::labs(
      title = "MCMC trace diagnostic",
      subtitle = subtitle,
      x = "Post-warmup iteration",
      y = NULL,
      color = "Chain"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = ggplot2::element_text(
        color = "grey35",
        size = 10.5
      ),
      legend.position = "top"
    )
}

.plot_bayes_worst <- function(x, n = 10L) {
  df <- x$parameter_summary
  n <- max(1L, as.integer(n))

  # Keep parameters for which at least one diagnostic is available.
  keep <-
    is.finite(df$rhat) |
    is.finite(df$ess_bulk) |
    is.finite(df$ess_tail)

  df <- df[keep, , drop = FALSE]

  if (nrow(df) == 0L) {
    stop(
      "No finite convergence diagnostics are available.",
      call. = FALSE
    )
  }

  # Worst ESS across bulk and tail.
  df$min_ess <- mapply(
    function(bulk, tail) {
      z <- c(bulk, tail)
      z <- z[is.finite(z)]

      if (length(z) == 0L) {
        return(NA_real_)
      }

      min(z)
    },
    df$ess_bulk,
    df$ess_tail
  )

  # ---------------------------------------------------------
  # Rank parameters
  #
  # Priority:
  #   1. parameters violating either R-hat or ESS threshold
  #   2. larger R-hat
  #   3. smaller ESS
  # ---------------------------------------------------------

  df$flagged <-
    (is.finite(df$rhat) &
       df$rhat > x$thresholds$rhat) |
    (is.finite(df$min_ess) &
       df$min_ess < x$thresholds$ess)

  ord <- order(
    -as.integer(df$flagged),
    -ifelse(is.finite(df$rhat), df$rhat, -Inf),
    ifelse(is.finite(df$min_ess), df$min_ess, Inf)
  )

  df <- df[ord, , drop = FALSE]
  df <- utils::head(df, n)

  # ---------------------------------------------------------
  # Convert R-hat and ESS to threshold-relative severity.
  #
  # A score of:
  #   1 = exactly at warning threshold
  #   >1 = worse than threshold
  #   <1 = within threshold
  #
  # R-hat:
  #   (Rhat - 1) / (threshold - 1)
  #
  # ESS:
  #   ESS threshold / observed ESS
  # ---------------------------------------------------------

  df$rhat_score <- ifelse(
    is.finite(df$rhat),
    pmax(
      0,
      (df$rhat - 1) /
        (x$thresholds$rhat - 1)
    ),
    NA_real_
  )

  df$ess_score <- ifelse(
    is.finite(df$min_ess) & df$min_ess > 0,
    x$thresholds$ess / df$min_ess,
    NA_real_
  )

  # Long format for plotting.
  plot_df <- rbind(
    data.frame(
      variable = df$variable,
      group = df$group,
      diagnostic = "R-hat",
      severity = df$rhat_score,
      value = df$rhat,
      stringsAsFactors = FALSE
    ),
    data.frame(
      variable = df$variable,
      group = df$group,
      diagnostic = "Minimum ESS",
      severity = df$ess_score,
      value = df$min_ess,
      stringsAsFactors = FALSE
    )
  )

  plot_df <- plot_df[
    is.finite(plot_df$severity),
    ,
    drop = FALSE
  ]

  # Preserve ranking from worst to best.
  variable_levels <- rev(unique(df$variable))

  plot_df$variable <- factor(
    plot_df$variable,
    levels = variable_levels
  )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[["severity"]],
      y = .data[["variable"]],
      shape = .data[["diagnostic"]]
    )
  ) +
    ggplot2::annotate(
      "rect",
      xmin = 1,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.05
    ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = 2,
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      size = 2.8,
      alpha = 0.8
    ) +
    ggplot2::labs(
      title = "Parameters with the weakest MCMC diagnostics",
      subtitle = paste0(
        "Top ", nrow(df),
        " parameters ranked by R-hat and effective sample size"
      ),
      x = "Diagnostic severity relative to warning threshold",
      y = NULL,
      shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        color = "grey25"
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = ggplot2::element_text(
        color = "grey35",
        size = 10.5
      ),
      legend.position = "top"
    )
}


# -----------------------------------------------------------------------------
# Plot method
# -----------------------------------------------------------------------------

#' Plot Bayesian LatentComBat / Univariate Stan Diagnostics
#'
#' @param x An object returned by [bayes_diagnostics()].
#' @param type Diagnostic plot type: `"rhat"`, `"ess"`, `"trace"` or `"worst`.
#' @param variables Optional scalar Stan variable names for `type = "trace"`.
#' @param n Number of automatically selected trace plots when `variables` is
#'   `NULL`. Defaults to `4`.
#' @param ... Additional arguments currently ignored.
#'
#' @return A ggplot object.
#'
#' @method plot latentcombat_bayes_diagnostics
#' @export
plot.latentcombat_bayes_diagnostics <- function(
    x,
    type = c("rhat", "ess", "trace", "worst"),
    variables = NULL,
    n = 4L,
    ...
) {
  type <- match.arg(type)

  switch(
    type,
    rhat = .plot_bayes_rhat(x),
    ess = .plot_bayes_ess(x),
    trace = .plot_bayes_trace(
      x,
      variables = variables,
      n = n
    ),
    worst = .plot_bayes_worst(
      x,
      n = n
    )
  )
}
