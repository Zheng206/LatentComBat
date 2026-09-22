#' Sequential Latent Harmonization Followed by ComBat
#'
#' Applies the sequential harmonization workflow, with optional surrogate-
#' variable adjustment for latent mean effects and optional latent variance
#' stabilization, followed by primary-batch ComBat and optional covariance
#' harmonization.
#'
#' @inheritParams com_harm
#'
#' @return An object of class `"latentcombat_sequential"` containing:
#' \describe{
#'   \item{harm_data}{
#'     Harmonized data matrix.
#'   }
#'   \item{resid}{
#'     Residual matrix from the final ComBat fit.
#'   }
#'   \item{data_stand_result}{
#'     Standardization results from the final ComBat step.
#'   }
#'   \item{latent}{
#'     Information from surrogate-variable estimation when `sva = TRUE`;
#'     otherwise `NULL`.
#'   }
#'   \item{protection}{
#'     Information on protected covariate adjustment when applicable;
#'     otherwise `NULL`.
#'   }
#'   \item{variance}{
#'     Results from latent variance stabilization when `var_stable = TRUE`;
#'     otherwise `NULL`.
#'   }
#'   \item{inference_data}{
#'     Original inputs and harmonization settings used by the optional
#'     full-pipeline bootstrap in [harm_inference()] when
#'     `keep_inference_data = TRUE`.
#'   }
#'   \item{method}{
#'     Character string `"sequential"`.
#'   }
#'   \item{eb_result}{
#'     Empirical-Bayes results when the final ComBat step uses the EB engine.
#'   }
#'   \item{stan_result}{
#'     Stan results when the final ComBat step uses the Stan engine.
#'   }
#' }
#'
#' @details
#' When `sva = TRUE`, latent mean variation is estimated and removed before
#' batch harmonization. When `var_stable = TRUE`, latent variance
#' stabilization is subsequently applied. These options are independent; if
#' both are enabled, latent mean adjustment is performed first, followed by
#' variance stabilization.
#'
#' The resulting data are then harmonized for the primary batch effect using
#' ComBat. If `cov = TRUE`, covariance harmonization is additionally applied
#' using the selected principal-component variance threshold.
#'
#' When `stan = TRUE`, the final ComBat step is fit with the Stan/NUTS engine.
#' The arguments `adapt_delta`, `iter_warmup`, `iter_sampling`, and
#' `max_treedepth` control the NUTS sampler and are ignored when the
#' empirical-Bayes engine is used.
#'
#' @export
com_harm.sequential <- function(
    bat,
    data,
    covar,
    model = lm,
    formula = NULL,
    ref.batch = NULL,
    eb = TRUE,
    stan = FALSE,
    robust.LS = FALSE,
    cov = FALSE,
    var_thresh = 0.95,
    min_rblock = 1,
    max_rblock = Inf,
    sva = FALSE,
    var_stable = FALSE,
    protected_covar = NULL,
    adapt_delta = 0.98,
    iter_sampling = 1000,
    iter_warmup = 1000,
    max_treedepth = 12,
    keep_inference_data = TRUE,
    ...
) {
  original_data <- as.matrix(data)
  extra_args <- list(...)

  batch_result <- batch_matrix(bat, ref.batch = ref.batch)

  latent_info <- NULL
  protection_info <- NULL
  variance_info <- NULL

  if (isTRUE(sva)) {
    fitted_model <- model_fitting(
      data = data,
      batch = batch_result$batch_matrix,
      covar = covar,
      model = model,
      formula = formula,
      batch_include = TRUE,
      ...
    )

    sva_fit <- .com_harm_sva_mean(
      data = data,
      bat = bat,
      covar = covar,
      fitted_model = fitted_model,
      model = model,
      formula = formula,
      protected_covar = protected_covar,
      B = 50L,
      alpha = 0.05
    )

    data <- sva_fit$data
    latent_info <- sva_fit$latent
    protection_info <- sva_fit$protection
  }

  if (isTRUE(var_stable)) {
    variance_fit <- .stabilize_latent_variance(
      data = data,
      bat = bat,
      covar = covar,
      batch_result = batch_result,
      model = model,
      formula = formula,
      ...
    )

    data <- variance_fit$data
    variance_info <- variance_fit
  }

  combat_fit <- .com_harm_combat(
    data = data,
    bat = bat,
    covar = covar,
    batch_result = batch_result,
    model = model,
    formula = formula,
    eb = eb,
    stan = stan,
    robust.LS = robust.LS,
    ref.batch = ref.batch,
    cov = cov,
    var_thresh = var_thresh,
    min_rblock = min_rblock,
    max_rblock = max_rblock,
    adapt_delta = adapt_delta,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    max_treedepth = max_treedepth,
    ...
  )

  inference_data <- NULL
  if (isTRUE(keep_inference_data)) {
    inference_data <- list(
      data = original_data,
      bat = bat,
      covar = covar,
      model = model,
      formula = formula,
      ref.batch = ref.batch,
      eb = eb,
      stan = stan,
      robust.LS = robust.LS,
      cov = cov,
      var_thresh = var_thresh,
      min_rblock = min_rblock,
      max_rblock = max_rblock,
      sva = sva,
      var_stable = var_stable,
      protected_covar = protected_covar,
      adapt_delta = adapt_delta,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      max_treedepth = max_treedepth,
      extra_args = extra_args
    )
  }

  out <- list(
    harm_data = combat_fit$harm_data,
    resid = combat_fit$resid,
    data_stand_result = combat_fit$data_stand_result,
    latent = latent_info,
    protection = protection_info,
    variance = variance_info,
    inference_data = inference_data,
    method = "sequential"
  )

  if (combat_fit$engine == "eb") out$eb_result <- combat_fit$fit
  if (combat_fit$engine == "stan") out$stan_result <- combat_fit$fit

  class(out) <- c("latentcombat_sequential", "latentcombat_fit", "list")
  out
}
