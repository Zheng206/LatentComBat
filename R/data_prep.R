#' Construct Batch Design Information
#'
#' Creates a batch design matrix and associated indexing information used by
#' the harmonization procedures.
#'
#' @param bat Factor specifying batch or site membership.
#' @param ref.batch Optional reference batch.
#'
#' @return A list containing the batch factor, batch design matrix, observation
#'   indices for each batch, batch sizes, and a reference-batch indicator.
#'
#' @export
batch_matrix <- function(bat, ref.batch = NULL){
  if (!is.factor(bat)) {
    bat <- factor(bat)
  }
  bat <- droplevels(bat)
  batch <- model.matrix(~ -1 + bat)
  batch_index <- lapply(levels(bat), function(x) which(bat == x))
  names(batch_index) <- levels(bat)
  n_batches <- sapply(batch_index, length)

  if (!is.null(ref.batch)) {
    if (!(ref.batch %in% levels(bat))) {
      stop("Reference batch must be in the batch levels")
    }
    ref <- bat == ref.batch
  }else{ref <- NULL}
  return(list("batch_vector" = bat, "batch_matrix" = batch, "batch_index" = batch_index, "n_batches" = n_batches, "ref" = ref))
}


#' Fit Feature-Wise Mean Models
#'
#' Fits a separate regression model to each feature using observed covariates
#' and, optionally, batch indicators.
#'
#' @param data Numeric matrix or data frame with observations in rows and
#'   features in columns.
#' @param batch Batch design matrix.
#' @param covar Optional data frame of observed covariates.
#' @param model Model-fitting function, such as [stats::lm()], [mgcv::gam()],
#'   or [lme4::lmer()].
#' @param formula Optional model formula.
#' @param batch_include Logical indicating whether batch terms are included in
#'   the fitted model. Defaults to `TRUE`.
#' @param weights Optional observation-level weight vector or feature-specific
#'   weight matrix.
#' @param ... Additional arguments passed to the model-fitting function.
#'
#' @return A list containing the original data, model data frame, formula, and
#'   fitted feature-wise models.
#'
#' @export
model_fitting <- function(data, batch, covar, model, formula = NULL,
                          batch_include = TRUE, weights = NULL, ...) {

  if (is.null(covar)) {
    mod <- data.frame(batch = I(batch))
    formula <- y ~ 1
    warning("No covariates were provided. A null model has been fitted.")
  } else {
    mod <- data.frame(covar, batch = I(batch))
    if (is.null(formula)) {
      stop("Please provide a formula to fit the model!")
    }
  }

  if (batch_include) {
    f_use <- update(formula, ~ . + batch + -1)
  } else {
    f_use <- update(formula, ~ . + -1)
  }

  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      if (nrow(weights) != nrow(data)) {
        stop("Weight matrix must have same number of rows as data")
      }
      if (ncol(weights) != ncol(data)) {
        stop("Weight matrix must have same number of columns as data")
      }
      use_feature_weights <- TRUE
    } else {
      if (length(weights) != nrow(data)) {
        stop("Weight vector must have length equal to number of observations")
      }
      use_feature_weights <- FALSE
    }
  } else {
    use_feature_weights <- FALSE
  }

  fits <- lapply(seq_len(ncol(data)), function(g) {
    y <- data[, g]
    names(y) <- "y"
    dat <- data.frame(y = y, mod)
    args <- list(formula = f_use, data = dat)

    if (!is.null(weights)) {
      if (use_feature_weights) {
        args$weights <- weights[, g]
      } else {
        args$weights <- weights
      }
    }

    extra_args <- list(...)
    extra_args$weights <- NULL
    args <- c(args, extra_args)
    do.call(model, args)
  })

  names(fits) <- colnames(data)

  return(list(
    data = data,
    mod = mod,
    formula = formula,
    fits = fits
  ))
}

#' Standardize Data for ComBat Harmonization
#'
#' Constructs the expected mean structure under a pooled or reference-batch
#' configuration and standardizes each feature using pooled residual variance.
#'
#' @param model A fitted model object returned by [model_fitting()].
#' @param batch_result Batch information returned by [batch_matrix()].
#' @param robust.LS Logical indicating whether a robust variance estimator is
#'   used instead of the sample variance.
#' @param bayes_sva Logical indicating whether predictions should exclude
#'   random effects for mixed-effects models. Defaults to `FALSE`.
#' @param keep_cols Optional character vector of covariate columns whose values
#'   should be retained when constructing the standardization mean.
#'
#' @return A list containing the standardized data (`data_stand`), fitted
#'   standardization mean (`stand_mean`), and feature-wise scaling matrix
#'   (`sd_mat`).
#'
#' @export
standardize_data <- function(model, batch_result, robust.LS = FALSE, bayes_sva = FALSE, keep_cols = NULL){
  pmod <- model$mod
  n = length(batch_result$batch_vector)
  pmod$batch[] <- matrix(batch_result$n_batches/n, n, nlevels(batch_result$batch_vector), byrow = TRUE)
  if (!is.null(batch_result$ref)) {
    pmod$batch[] <- 0
    pmod$batch[,batch_result$batch_vector[which(batch_result$ref == 1)] |> unique()] <- 1
  }

  if (!is.null(keep_cols)) {
    all_vars <- colnames(pmod)
    nuisance_vars <- setdiff(all_vars, c("batch", keep_cols))
    for (v in nuisance_vars) {
      if (v %in% names(pmod)) {
        if (is.numeric(pmod[[v]])) {
          pmod[[v]] <- 0
        } else if (is.factor(pmod[[v]]) || is.character(pmod[[v]])) {
          ref_level <- levels(as.factor(pmod[[v]]))[1]
          pmod[[v]] <- factor(rep(ref_level, n), levels = levels(pmod[[v]]))
        }
      }
    }
  }

  if(bayes_sva){
    get_pred <- function(fit, newdata) {
      if (inherits(fit, "merMod")) {
        return(predict(fit, newdata = newdata, type = "response", re.form = NA))
      } else {
        return(predict(fit, newdata = newdata, type = "response"))
      }
    }
    stand_mean <- sapply(model$fits, get_pred, newdata = pmod)
    resid_mean <- sapply(model$fits, get_pred, newdata = model$mod)
  }else{
    stand_mean <- sapply(model$fits, predict, newdata = pmod, type = "response")
    resid_mean <- sapply(model$fits, predict, newdata = model$mod, type = "response")
  }

  if(robust.LS){
    scale = biweight_midvar
  }else{
    scale = var
  }

  if (!is.null(batch_result$ref)) {
    nref <- sum(batch_result$ref)
    var_pooled <- apply((model$data - resid_mean)[batch_result$ref, , drop = FALSE], 2, scale) *
      (nref - 1)/nref
  } else {
    var_pooled <- apply(model$data - resid_mean, 2, scale) * (n - 1)/n
  }
  sd_mat <- sapply(sqrt(var_pooled), rep, n)
  data_stand <- (model$data-stand_mean)/sd_mat
  return(list("data_stand" = data_stand, "stand_mean" = stand_mean, "sd_mat" = sd_mat))
}


#' Predict the Fixed-Effect Mean
#'
#' Generates fitted mean values from a regression model, excluding random
#' effects for mixed-effects models.
#'
#' @param fit A fitted regression or mixed-effects model.
#' @param newdata Data used to generate predictions.
#'
#' @return A numeric vector of predicted mean values.
#'
#' @keywords internal
predict_mean_safe <- function(fit, newdata) {
  if (inherits(fit, "merMod")) {
    predict(fit, newdata = newdata, re.form = NA, type = "response")
  } else {
    predict(fit, newdata = newdata, type = "response")
  }
}

#' Extract Feature-Wise Fitted Means
#'
#' Computes the fitted mean for each feature from a feature-wise model object.
#'
#' @param fitted_model A fitted model object returned by [model_fitting()].
#'
#' @return A numeric matrix of fitted mean values with observations in rows and
#'   features in columns.
#'
#' @keywords internal
get_mu_hat <- function(fitted_model) {
  mu <- sapply(fitted_model$fits, predict_mean_safe, newdata = fitted_model$mod)
  mu <- as.matrix(mu)
  colnames(mu) <- names(fitted_model$fits)
  mu
}


#' Estimate Batch-Specific Residual Variance Effects
#'
#' Estimates feature- and batch-specific residual variance effects using
#' log-squared residuals.
#'
#' @param e Numeric residual matrix with observations in rows and features in
#'   columns.
#' @param bat Factor specifying batch membership.
#' @param eps Optional positive constant added before taking the logarithm of
#'   squared residuals.
#'
#' @return A list containing the fitted log-variance surface (`ell`), numerical
#'   offset (`eps`), pooled feature-specific log variance (`c_g`), batch-specific
#'   log-variance deviations (`theta_bg`), and corresponding log-scale effects
#'   (`log_delta_bg`).
#'
#' @keywords internal
variance_step_quasi <- function(e, bat, eps = NULL) {
  e <- as.matrix(e)
  N <- nrow(e); P <- ncol(e)
  bat <- droplevels(bat)
  if (is.null(eps)) {
    med <- stats::median(as.numeric(e^2), na.rm = TRUE)
    if (!is.finite(med) || med <= 0) med <- 1
    eps <- 1e-8 * med
  }
  Z <- log(e^2 + eps)
  Xb <- model.matrix(~ bat - 1)
  coef_bg <- qr.solve(Xb, Z)
  coef_bg <- t(coef_bg)
  if (nrow(coef_bg) != ncol(Xb)) coef_bg <- t(coef_bg)
  w_b <- as.numeric(table(bat)) / length(bat)
  mu_g <- as.numeric(t(w_b) %*% coef_bg)
  theta_bg <- sweep(coef_bg, 2, mu_g, "-")
  log_delta_bg <- 0.5 * theta_bg
  b_int <- as.integer(bat)
  ell <- matrix(rep(mu_g, each = N), N, P) + theta_bg[b_int, , drop = FALSE]
  list(ell = ell, eps = eps, c_g = mu_g, theta_bg = theta_bg, log_delta_bg = log_delta_bg)
}

#' Stabilize Batch-Specific Residual Variance
#'
#' Iteratively estimates and smooths batch-specific residual log-variance
#' effects while keeping the fitted mean structure fixed.
#'
#' @param data Numeric matrix or data frame with observations in rows and
#'   features in columns.
#' @param bat Factor specifying batch membership.
#' @param fitted_model A fitted mean-model object returned by [model_fitting()].
#' @param n_iter Number of variance-update iterations. Defaults to `4`.
#' @param damp Damping parameter controlling the update of the fitted
#'   log-variance surface. Defaults to `0.5`.
#' @param Kv Deprecated or unused latent-variance dimension parameter.
#' @param eps Optional positive constant added before taking the logarithm of
#'   squared residuals.
#' @param verbose Logical indicating whether iteration diagnostics are printed.
#'
#' @return A list containing the final fitted log-variance surface (`ell`),
#'   the most recent variance-fit result (`vfit`), and the fixed fitted mean
#'   matrix (`mu_hat`).
#'
#' @keywords internal
em_var_stabilize_only <- function(data, bat, fitted_model, n_iter = 4L,
                                  damp = 0.5, Kv = 0L, eps = NULL, verbose = FALSE) {
  Y <- as.matrix(data)
  N <- nrow(Y); P <- ncol(Y)
  mu_hat <- sapply(fitted_model$fits, predict_mean_safe, newdata = fitted_model$mod)
  mu_hat <- as.matrix(mu_hat)
  ell <- matrix(0, N, P)
  vfit <- NULL

  for (it in seq_len(n_iter)) {
    e <- Y - mu_hat
    vfit_new <- variance_step_quasi(e, bat, eps = eps)
    if (is.null(eps)) eps <- vfit_new$eps
    ell <- (1 - damp) * ell + damp * vfit_new$ell
    vfit <- vfit_new
    if (verbose) {
      obj <- 0.5 * sum(ell + (e^2) * exp(-ell))
      message(sprintf("  var-EM iter %02d | NLL~ %.4f", it, obj))
    }
  }
  list(ell = ell, vfit = vfit, mu_hat = mu_hat)
}


#' Tukey biweight mid-variance (univariate robust variance)
#'
#' Computes a robust variance estimate using Tukey's biweight with
#' median/MAD centering and scaling. Observations with \eqn{|u_i| \ge 1} receive
#' zero weight, where \eqn{u_i = (x_i - \mathrm{center})/(c \cdot \mathrm{MAD})}.
#'
#' @param data Numeric vector.
#' @param center Optional numeric center. Defaults to \code{median(data)}.
#' @param norm.unbiased Logical; if \code{TRUE} (default), uses
#'   \code{c = 9 / qnorm(0.75)} under
#'   normality. If \code{FALSE}, uses \code{c = 9}.
#'
#' @details
#' Let \eqn{d_i = x_i - \mathrm{center}}, \eqn{\mathrm{MAD} = \mathrm{median}(|d_i|)},
#' and \eqn{u_i = d_i/(c\,\mathrm{MAD})}. The estimator is
#' \deqn{\hat{\sigma}^2_{\mathrm{BI}} =
#' n \cdot \frac{\sum_i \mathbb{1}(|u_i|<1)\, d_i^2 (1-u_i^2)^4}
#' {\left\{\sum_i \mathbb{1}(|u_i|<1)\, (1-u_i^2)(1-5u_i^2)\right\}^2}.}
#'
#' @return A single numeric: robust variance estimate.
#'
#' @note
#' If \code{MAD == 0} or the denominator is zero (e.g., all points trimmed),
#' the result may be \code{NaN}/\code{Inf}. Consider falling back to
#' \code{var(data)} or adding a small jitter in such edge cases.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, 0, 1), 10)  # one outlier
#' var(x)
#' biweight_midvar(x)           # robust to the outlier
#'
#' @seealso \code{\link[stats]{var}}
#' @export

biweight_midvar <- function(data, center=NULL, norm.unbiased = TRUE) {
  if (is.null(center)) {
    center <- median(data)
  }

  mad <- median(abs(data - center))
  d <- data - center
  c <- ifelse(norm.unbiased, 9/qnorm(0.75), 9)
  u <- d/(c*mad)

  n <- length(data)
  indic <- abs(u) < 1

  num <- sum(indic * d^2 * (1 - u^2)^4)
  dem <- sum(indic * (1 - u^2) * (1 - 5*u^2))^2

  n * num/dem
}




