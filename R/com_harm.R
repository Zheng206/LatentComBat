#' Harmonize Data Across Batches
#'
#' Main user-facing function for batch harmonization. Depending on `method`,
#' the function applies either a sequential SVA/ComBat workflow or a joint
#' Bayesian latent ComBat model.
#'
#' @param bat Batch or site factor with one value per observation.
#' @param data Numeric matrix or data frame with observations in rows and
#'   features in columns.
#' @param covar Optional data frame of observed covariates to preserve.
#' @param model Mean-model fitting function, such as [stats::lm()] or
#'   [mgcv::gam()].
#' @param formula Optional model formula.
#' @param ref.batch Optional reference batch.
#' @param eb Logical indicating whether empirical-Bayes shrinkage is used in
#'   the sequential ComBat step. Defaults to `TRUE`.
#' @param stan Logical indicating whether Stan-based estimation is used in the
#'   sequential workflow. Defaults to `FALSE`.
#' @param robust.LS Logical indicating whether robust location and scale
#'   estimates are used. Defaults to `FALSE`.
#' @param cov Logical indicating whether covariance harmonization is applied.
#'   Defaults to `FALSE`.
#' @param var_thresh Cumulative variance threshold used to select principal
#'   components for covariance harmonization. Defaults to `0.95`.
#' @param min_rblock Minimum number of principal components to harmonize.
#'   Defaults to `1`.
#' @param max_rblock Maximum number of principal components to harmonize.
#'   Defaults to `Inf`.
#' @param sva Logical indicating whether surrogate-variable adjustment is
#'   performed before ComBat in the sequential workflow. Defaults to `FALSE`.
#' @param bayes_sva Deprecated compatibility flag. If `TRUE` and `method` is
#'   not specified, the Bayesian workflow is selected.
#' @param method Harmonization strategy, either `"sequential"` or `"bayesian"`.
#'   If `NULL`, the sequential workflow is used unless `bayes_sva = TRUE`.
#' @param K_max Maximum number of latent factors considered in latent-variable
#'   estimation. Defaults to `10`.
#' @param K_buffer Additional latent dimensions considered during latent-factor
#'   selection. Defaults to `2`.
#' @param subj Optional subject identifier for longitudinal or repeated-measures
#'   data.
#' @param var_stable Logical indicating whether latent variance stabilization
#'   is applied before the final ComBat step. This option is independent of
#'   `sva`. Defaults to `FALSE`.
#' @param protected_covar Optional additional covariates or protected subspace
#'   used to reduce removal of preserved biological variation.
#' @param posterior_draws Number of posterior harmonized-data draws to retain
#'   when `method = "bayesian"`. Defaults to `0`, which preserves the usual
#'   point-harmonization workflow without storing draws.
#' @param posterior_seed Random seed used only to select retained posterior
#'   draws from the fitted Bayesian model. Defaults to `1`.
#' @param keep_inference_data Logical indicating whether the original inputs and
#'   harmonization settings needed by [harm_inference()] should be stored in the
#'   fitted object. Defaults to `TRUE`.
#' @param adapt_delta Target acceptance probability for the NUTS sampler.
#'   Larger values use smaller leapfrog step sizes and can reduce divergent
#'   transitions at the cost of additional computation. Defaults to `0.98`.
#' @param iter_sampling Number of post-warmup MCMC iterations retained per
#'   chain when `method = "bayesian"`. Defaults to `1000`.
#' @param iter_warmup Number of warmup iterations per chain used for NUTS
#'   adaptation when `method = "bayesian"`. Warmup draws are not retained for
#'   posterior inference. Defaults to `1000`.
#' @param max_treedepth Maximum tree depth allowed for each NUTS trajectory.
#'   Larger values permit longer trajectories when the posterior geometry is
#'   difficult to explore, at the cost of additional computation. Defaults to
#'   `12`.
#' @param ... Additional arguments passed to lower-level model-fitting and
#'   harmonization functions.
#'
#' @return A harmonization result returned by either [com_harm.sequential()]
#'   or [com_harm.bayesian()], depending on the selected method.
#'
#' @details
#' With `method = "sequential"`, optional surrogate-variable adjustment and
#' latent variance stabilization are performed before the final ComBat step.
#' When both `sva = TRUE` and `var_stable = TRUE`, latent mean adjustment is
#' performed first, followed by variance stabilization.
#'
#' With `method = "bayesian"`, latent variation and batch effects are estimated
#' jointly using the Bayesian latent ComBat model.
#'
#' @export
com_harm <- function(
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
    bayes_sva = FALSE,
    method = NULL,
    K_max = 10,
    K_buffer = 2,
    subj = NULL,
    var_stable = FALSE,
    protected_covar = NULL,
    posterior_draws = 0L,
    posterior_seed = 1L,
    adapt_delta = 0.98,
    iter_sampling = 1000,
    iter_warmup = 1000,
    max_treedepth = 12,
    keep_inference_data = TRUE,
    ...
) {
  if (is.null(method)) {
    method <- if (isTRUE(bayes_sva)) "bayesian" else "sequential"
  } else {
    method <- match.arg(method, c("sequential", "bayesian"))
    if (isTRUE(bayes_sva) && method != "bayesian") {
      warning("`method` overrides legacy `bayes_sva = TRUE`.", call. = FALSE)
    }
  }

  if (method == "bayesian") {
    return(com_harm.bayesian(
      bat = bat, data = data, covar = covar, model = model,
      formula = formula, ref.batch = ref.batch,
      robust.LS = robust.LS, K_max = K_max, K_buffer = K_buffer,
      subj = subj, protected_covar = protected_covar,
      posterior_draws = posterior_draws, posterior_seed = posterior_seed,
      keep_inference_data = keep_inference_data, adapt_delta = adapt_delta, iter_sampling = iter_sampling,
      iter_warmup = iter_warmup, max_treedepth = max_treedepth, ...))
  }

  com_harm.sequential(
    bat = bat, data = data, covar = covar, model = model,
    formula = formula, ref.batch = ref.batch, eb = eb, stan = stan,
    robust.LS = robust.LS, cov = cov, var_thresh = var_thresh,
    min_rblock = min_rblock, max_rblock = max_rblock, sva = sva,
    var_stable = var_stable, protected_covar = protected_covar,
    keep_inference_data = keep_inference_data, adapt_delta = adapt_delta, iter_sampling = iter_sampling,
    iter_warmup = iter_warmup, max_treedepth = max_treedepth, ...)
}
