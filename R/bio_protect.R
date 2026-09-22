#' Create an empty matrix
#'
#' Creates an empty numeric matrix with a specified number of rows and
#' zero columns.
#'
#' @param n Integer. Number of rows in the matrix.
#'
#' @return A numeric matrix with \code{n} rows and 0 columns.
#'
#' @keywords internal
empty_matrix <- function(n) {
  matrix(numeric(0), nrow = n, ncol = 0L)
}

#' Reduce a design matrix to full column rank using QR decomposition
#'
#' Removes numerically zero columns from a design matrix and uses a
#' pivoted QR decomposition to retain a linearly independent subset of
#' columns. The function also returns an orthonormal basis for the column
#' space of the reduced design matrix.
#'
#' @param Z Numeric matrix. Design matrix to be reduced.
#' @param tol Numeric. Numerical tolerance used to identify near-zero
#'   columns and determine the rank in the QR decomposition.
#'   Defaults to \code{1e-8}.
#'
#' @return A list with three components:
#' \describe{
#'   \item{\code{Z}}{The reduced design matrix containing a linearly
#'     independent subset of the original columns.}
#'   \item{\code{Q}}{An orthonormal basis for the column space of
#'     \code{Z}.}
#'   \item{\code{rank}}{The numerical rank of the design matrix.}
#' }
#'
#' If the input contains no columns, all columns are numerically zero, or
#' the resulting matrix has rank zero, \code{Z} and \code{Q} are returned
#' as matrices with \code{nrow(Z)} rows and zero columns, and
#' \code{rank = 0}.
#'
#' @keywords internal
reduce_design_qr <- function(Z, tol = 1e-8) {
  Z <- as.matrix(Z)
  n <- nrow(Z)

  if (ncol(Z) == 0L) {
    Z0 <- empty_matrix(nrow(Z))
    return(
      list(Z = Z0, Q = Z0, rank = 0L)
    )
  }

  keep <- colSums(Z^2) > tol^2
  Z <- Z[ , keep, drop = FALSE]

  if (ncol(Z) == 0L) {
    Z0 <- empty_matrix(nrow(Z))
    return(
      list(Z = Z0, Q = Z0, rank = 0L)
    )
  }

  qrZ <- qr(Z,tol = tol)
  r <- qrZ$rank

  if (r == 0L) {
    Z0 <- empty_matrix(nrow(Z))
    return(
      list(Z = Z0, Q = Z0, rank = 0L)
    )
  }

  keep_idx <- qrZ$pivot[seq_len(r)]
  Z_reduced <- Z[ , keep_idx, drop = FALSE]
  Q <- qr.Q(qrZ, complete = FALSE)[ , seq_len(r), drop = FALSE]

  return(
    list(Z = Z_reduced, Q = Q, rank = r)
  )
}

#' Construct a protected design matrix for latent-factor estimation
#'
#' Constructs a design matrix representing observed variation that should
#' be protected when estimating additional latent factors during
#' harmonization. The primary biological structure to preserve is defined by
#' the fixed effects specified in \code{formula}, corresponding to the
#' biological signals protected in the main ComBat model. The protected
#' design can additionally include the primary ComBat batch variable and
#' user-specified covariates whose variation should not be absorbed by the
#' latent factors.
#'
#' For mixed-effects models fitted with \code{lme4::lmer} or
#' \code{lme4::glmer}, only the fixed-effect component of the supplied
#' formula is used. For \code{mgcv::gam} and \code{mgcv::bam}, the
#' corresponding GAM design matrix is constructed so that smooth biological
#' effects are represented in the protected space. For other model types,
#' the design matrix is constructed using \code{stats::model.matrix}.
#'
#' If \code{preserve_batch = TRUE}, the primary batch variable is also
#' included in the protected design. This prevents the latent-factor
#' estimation step from absorbing variation attributable to the known batch
#' effect, which is handled separately by ComBat. Additional covariates may
#' be supplied through \code{protected_covar} to protect observed variation
#' not already represented by the main ComBat model.
#'
#' The combined protected design is passed to \code{reduce_design_qr} to
#' determine its numerical rank, remove redundant columns, and construct an
#' orthonormal basis for its column space.
#'
#' @param model Model-fitting function used to construct the protected design
#'   matrix. Special handling is provided for \code{lme4::lmer},
#'   \code{lme4::glmer}, \code{mgcv::gam}, and \code{mgcv::bam}.
#'
#' @param formula Model formula specifying the fixed effects to preserve
#'   during the ComBat harmonization step. These effects represent the main
#'   observed biological signals that should be protected from removal when
#'   estimating additional latent variation. For mixed-effects models, only
#'   the fixed-effect component of the formula is used.
#'
#' @param covar Data frame or object coercible to a data frame containing the
#'   covariates referenced in \code{formula}.
#'
#' @param bat Primary batch variable specified for the ComBat harmonization
#'   step. Its length must equal the number of observations in \code{covar}.
#'
#' @param preserve_batch Logical. If \code{TRUE}, the primary batch variable
#'   \code{bat} is also included in the protected design during latent-factor
#'   estimation. This prevents latent factors from capturing variation
#'   attributable to the known batch effect that will be handled separately
#'   by ComBat. Defaults to \code{TRUE}.
#'
#' @param protected_covar Optional data frame containing additional observed
#'   covariates whose variation should be protected during latent-factor
#'   estimation beyond the biological effects specified in \code{formula}.
#'   Constant columns are removed automatically. Missing values are not
#'   allowed.
#'
#' @param qr_tol Numeric. Numerical tolerance used by
#'   \code{reduce_design_qr} to identify near-zero columns and determine the
#'   numerical rank of the protected design matrix. Defaults to
#'   \code{1e-8}.
#'
#' @param return_Q Logical. If \code{FALSE}, returns the combined protected
#'   design matrix before QR reduction. If \code{TRUE}, returns the raw
#'   design matrix together with its full-rank representation, orthonormal
#'   basis, and numerical rank. Defaults to \code{FALSE}.
#'
#' @return If \code{return_Q = FALSE}, a numeric matrix containing the
#'   combined protected design before QR reduction.
#'
#'   If \code{return_Q = TRUE}, a list with the following components:
#' \describe{
#'   \item{\code{Z_raw}}{The combined protected design matrix before QR
#'     reduction.}
#'   \item{\code{Z}}{A full-column-rank subset of \code{Z_raw}, obtained
#'     using pivoted QR decomposition.}
#'   \item{\code{Q}}{An orthonormal basis for the column space of the
#'     protected design.}
#'   \item{\code{rank}}{The numerical rank of the protected design matrix.}
#' }
#'
#' @keywords internal
build_Z_preserve <- function(
    model,
    formula,
    covar,
    bat,
    preserve_batch = TRUE,
    protected_covar = NULL,
    qr_tol = 1e-8,
    return_Q = FALSE
) {
  covar <- as.data.frame(covar)
  n <- nrow(covar)

  if (length(bat) != n) {
    stop("`bat` and `covar` must have the same number of observations.")
  }

  if (is.null(formula)) {
    stop("`formula` is required to construct the protected design.")
  }

  f_rhs <- formula

  if (identical(model, lme4::lmer) || identical(model, lme4::glmer)) {f_rhs <- lme4::nobars(f_rhs)}

  if (identical(model, mgcv::gam) || identical(model, mgcv::bam)) {
    tmp <- data.frame(y = rep(0, n), covar, bat = bat)
    gam_obj <- mgcv::gam(f_rhs, data = tmp, fit = FALSE)
    Z_cov <- gam_obj$X
  } else {
    tmp <- data.frame(covar, bat = bat)
    Z_cov <- stats::model.matrix(stats::delete.response(stats::terms(f_rhs)), data = tmp)
  }

  Z_parts <- list(Z_cov)

  if (isTRUE(preserve_batch)) {
    Z_bat <- stats::model.matrix( ~ bat - 1)
    Z_parts[[length(Z_parts) + 1L]] <- Z_bat
  }

  if (!is.null(protected_covar)) {
    protected_covar <- as.data.frame(protected_covar)
    if (nrow(protected_covar) != n) {
      stop("`protected_covar` must have the same number of rows as `covar`.")
    }

    if (anyNA(protected_covar)) {
      stop("`protected_covar` contains missing values. ", "Please impute or otherwise handle missing values first.")
    }

    keep_extra <- vapply(protected_covar, function(x) {length(unique(x)) > 1L}, logical(1))
    protected_covar <- protected_covar[ , keep_extra, drop = FALSE]

    if (ncol(protected_covar) > 0L) {
      Z_extra <- stats::model.matrix( ~ ., data = protected_covar)
      if ("(Intercept)" %in% colnames(Z_extra)) {Z_extra <- Z_extra[ , colnames(Z_extra) != "(Intercept)", drop = FALSE]}
      if (ncol(Z_extra) > 0L) {Z_parts[[length(Z_parts) + 1L]] <- Z_extra}
    }
  }

  Z_raw <- as.matrix(do.call(cbind, Z_parts))

  if (any(!is.finite(Z_raw))) {
    stop("Protected design contains non-finite values.")
  }

  qr_out <- reduce_design_qr(Z_raw, tol = qr_tol)

  if (isTRUE(return_Q)) {
    return(
      list(
        Z_raw = Z_raw,
        Z = qr_out$Z,
        Q = qr_out$Q,
        rank = qr_out$rank
      )
    )
  }

  Z_raw
}

#' Orthogonalize latent factors against a protected subspace
#'
#' Removes from the latent factor matrix \code{H} any variation that lies in
#' the protected subspace spanned by the columns of \code{Q}. This ensures
#' that the estimated latent factors are orthogonal to the observed variation
#' that has been designated for preservation, such as the biological effects
#' specified in the main ComBat model and, when requested, the primary batch
#' variable or additional protected covariates.
#'
#' The orthogonalization is performed by subtracting the projection of
#' \code{H} onto the column space of \code{Q}. The columns of \code{Q} are
#' assumed to form an orthonormal basis, as returned by
#' \code{reduce_design_qr} or \code{build_Z_preserve}.
#'
#' If \code{Q} has zero columns, no protected subspace is present and
#' \code{H} is returned unchanged.
#'
#' @param H Numeric matrix containing estimated latent factors or surrogate
#'   variables, with observations in rows and latent factors in columns.
#'
#' @param Q Numeric matrix whose columns form an orthonormal basis for the
#'   protected subspace. The number of rows must match the number of rows in
#'   \code{H}. This matrix is typically obtained from the QR decomposition of
#'   the protected design matrix.
#'
#' @return A numeric matrix with the same dimensions as \code{H}, containing
#'   latent factors after removing their projections onto the protected
#'   subspace spanned by \code{Q}.
#'
#' @keywords internal
orthogonalize_latent <- function(H, Q) {

  H <- as.matrix(H)
  Q <- as.matrix(Q)

  if (nrow(H) != nrow(Q)) {
    stop("`H` and `Q` must have the same number of rows.")
  }

  if (ncol(Q) == 0L) {
    return(H)
  }

  H - Q %*% crossprod(Q, H)
}
