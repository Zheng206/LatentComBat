#' Print a LatentComBat Fit
#'
#' Prints a concise summary of a fitted LatentComBat object, including the
#' harmonization method, data dimensions, and latent-dimension information
#' when available.
#'
#' @param x An object inheriting from class `"latentcombat_fit"`.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return The input object `x`, returned invisibly.
#'
#' @method print latentcombat_fit
#' @export
print.latentcombat_fit <- function(x, ...) {
  cat("LatentComBat fit\n")
  cat("  Method: ", x$method, "\n", sep = "")
  if (!is.null(x$harm_data)) {
    cat("  Data:   ", nrow(x$harm_data), " x ", ncol(x$harm_data), "\n", sep = "")
  }
  if (!is.null(x$latent$K_target)) {
    cat("  K target: ", x$latent$K_target, "\n", sep = "")
    cat("  K fitted: ", x$latent$K_bayes, "\n", sep = "")
  } else if (!is.null(x$latent$H)) {
    k <- if (is.null(x$latent$H)) 0L else ncol(x$latent$H)
    cat("  SVA factors: ", k, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a LatentComBat Fit
#'
#' Constructs a compact summary of a fitted LatentComBat object, including the
#' harmonization method, data dimensions, latent-model information, protected
#' subspace rank, and whether latent variance stabilization was applied.
#'
#' @param object An object inheriting from class `"latentcombat_fit"`.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return An object of class `"summary_latentcombat_fit"` containing:
#' \describe{
#'   \item{method}{
#'     Harmonization method used in the fitted model.
#'   }
#'   \item{n}{
#'     Number of observations in the harmonized data.
#'   }
#'   \item{p}{
#'     Number of features in the harmonized data.
#'   }
#'   \item{latent}{
#'     Latent-model information stored in the fitted object.
#'   }
#'   \item{protection_rank}{
#'     Rank of the protected covariate subspace, when available.
#'   }
#'   \item{has_variance_stabilization}{
#'     Logical indicating whether the sequential latent variance-stabilization
#'     step was applied.
#'   }
#' }
#'
#' @method summary latentcombat_fit
#' @export
summary.latentcombat_fit <- function(object, ...) {
  out <- list(
    method = object$method,
    n = if (!is.null(object$harm_data)) nrow(object$harm_data) else NA_integer_,
    p = if (!is.null(object$harm_data)) ncol(object$harm_data) else NA_integer_,
    latent = object$latent,
    protection_rank = if (!is.null(object$protection)) object$protection$rank else NULL,
    has_variance_stabilization = !is.null(object$variance)
  )
  class(out) <- "summary_latentcombat_fit"
  out
}

#' Print a LatentComBat Summary
#'
#' Prints a compact summary of a `"summary_latentcombat_fit"` object.
#'
#' @param x An object of class `"summary_latentcombat_fit"`.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return The input object `x`, returned invisibly.
#'
#' @method print summary_latentcombat_fit
#' @export
print.summary_latentcombat_fit <- function(x, ...) {
  cat("LatentComBat summary\n")
  cat("  Method: ", x$method, "\n", sep = "")
  cat("  Dimensions: ", x$n, " x ", x$p, "\n", sep = "")
  if (!is.null(x$protection_rank)) {
    cat("  Protected-space rank: ", x$protection_rank, "\n", sep = "")
  }
  cat("  Latent variance stabilization: ",
      if (isTRUE(x$has_variance_stabilization)) "yes" else "no", "\n", sep = "")
  invisible(x)
}
