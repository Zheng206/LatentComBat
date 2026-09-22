#' Select the Number of Principal Components by Explained Variance
#'
#' Selects the smallest number of principal components (PCs) required to
#' explain at least a specified proportion of the total variance, subject to
#' user-defined lower and upper bounds on the number of retained components.
#'
#' @param pca_fit A fitted PCA object containing a numeric `sdev` component,
#'   such as an object returned by [stats::prcomp()]. The squared values of
#'   `sdev` are used to calculate the proportion of variance explained by
#'   each principal component.
#'
#' @param var_thresh Numeric scalar between 0 and 1 specifying the minimum
#'   cumulative proportion of variance to be explained by the retained
#'   principal components. Defaults to `0.95`.
#'
#' @param min_r Integer specifying the minimum number of principal components
#'   to retain. Defaults to `5`.
#'
#' @param max_r Integer or `Inf` specifying the maximum number of principal
#'   components to retain. Defaults to `Inf`.
#'
#' @return An integer giving the selected number of principal components.
#'   The selected value is the smallest number of PCs whose cumulative
#'   explained variance reaches `var_thresh`, constrained to be at least
#'   `min_r`, no greater than `max_r`, and no greater than the total number
#'   of available principal components.
#'
#' @details
#' Let \eqn{d_k} denote the standard deviation of the \eqn{k}-th principal
#' component. The cumulative proportion of variance explained by the first
#' \eqn{r} components is
#'
#' \deqn{
#' \frac{\sum_{k=1}^{r} d_k^2}
#'      {\sum_{k=1}^{K} d_k^2},
#' }
#'
#' where \eqn{K} is the total number of available principal components.
#' The function first selects the smallest \eqn{r} for which this cumulative
#' proportion is at least `var_thresh`, and then applies the `min_r` and
#' `max_r` constraints.
#'
#' @examples
#' x <- scale(USArrests)
#' pca_fit <- prcomp(x)
#'
#' pick_r_from_pc(
#'   pca_fit,
#'   var_thresh = 0.95,
#'   min_r = 2
#' )
#'
#' pick_r_from_pc(
#'   pca_fit,
#'   var_thresh = 0.95,
#'   min_r = 2,
#'   max_r = 3
#' )
#'
#' @export
pick_r_from_pc <- function(pca_fit, var_thresh = 0.95, min_r = 5, max_r = Inf) {
  total_variance <- sum(pca_fit$sdev^2)
  if (!is.finite(total_variance) || total_variance <= 0) {
    stop("PCA fit contains no positive variance.", call. = FALSE)
  }
  cumve <- cumsum(pca_fit$sdev^2) / total_variance
  r <- which(cumve >= var_thresh)[1]
  r <- max(min_r, r)
  r <- min(r, length(pca_fit$sdev))
  r <- min(r, max_r)
  return(r)
}


#' Covariance Batch Harmonization Using Principal Components
#'
#' Applies a CovBat-style harmonization procedure to reduce batch- or
#' site-related differences in covariance structure. The input data are
#' decomposed using principal component analysis (PCA), a selected subset of
#' leading principal component scores is harmonized across sites using
#' [com_harm()], and the harmonized data are reconstructed in the original
#' feature space.
#'
#' @param R A numeric matrix or data frame with observations in rows and
#'   features in columns. In typical use, `R` contains residualized data after
#'   removal of biological or other covariate effects that should be preserved.
#'
#' @param site A vector specifying the batch or site membership of each
#'   observation. Its length must equal the number of rows of `R`.
#'
#' @param center Logical indicating whether the variables should be centered
#'   before PCA. Passed to [stats::prcomp()]. Defaults to `TRUE`.
#'
#' @param scale_scores Logical indicating whether the variables should be
#'   scaled to unit variance before PCA. Passed to the `scale.` argument of
#'   [stats::prcomp()]. Defaults to `TRUE`.
#'
#' @param var_thresh Numeric scalar between 0 and 1 specifying the cumulative
#'   proportion of variance used to determine the number of leading principal
#'   components to harmonize. Defaults to `0.95`.
#'
#' @param min_rblock Integer specifying the minimum number of leading principal
#'   components to harmonize. Defaults to `1`.
#'
#' @param max_rblock Integer or `Inf` specifying the maximum number of leading
#'   principal components to harmonize. Defaults to `Inf`.
#'
#' @param ref.batch Optional reference batch passed to [com_harm()]. If
#'   specified, harmonization is performed relative to this batch. Defaults
#'   to `NULL`.
#'
#' @return A numeric matrix with the same dimensions as `R`, containing the
#'   reconstructed data after harmonization of the selected leading principal
#'   component scores.
#'
#' @details
#' Let the PCA decomposition of the centered and optionally scaled data be
#'
#' \deqn{
#' R_s = U D V^\top,
#' }
#'
#' where the rows of \eqn{UD} correspond to subject-level principal component
#' scores and the columns of \eqn{V} are the principal component loading
#' vectors.
#'
#' The number of components to harmonize, \eqn{r}, is selected using
#' [pick_r_from_pc()] as the smallest number of leading components whose
#' cumulative explained variance reaches `var_thresh`, subject to the
#' `min_rblock` and `max_rblock` constraints.
#'
#' The first \eqn{r} principal component score vectors are then harmonized
#' across sites using [com_harm()] with empirical-Bayes shrinkage disabled.
#' Principal components beyond \eqn{r} are left unchanged. The complete score
#' matrix is then projected back to the original feature space using the PCA
#' loadings, after which the original centering and scaling transformations are
#' reversed.
#'
#' Thus, only variation represented by the selected leading PCA subspace is
#' directly adjusted for site effects, while variation in the remaining
#' principal components is retained unchanged.
#'
#' @seealso
#' [pick_r_from_pc()], [com_harm()], [stats::prcomp()]
#'
#' @examples
#' set.seed(123)
#'
#' R <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' site <- rep(c("Site1", "Site2"), each = 50)
#'
#' R_harm <- covbat(
#'   R = R,
#'   site = site,
#'   var_thresh = 0.95,
#'   min_rblock = 1
#' )
#'
#' @export
covbat <- function(R, site, center = TRUE,
                   scale_scores = TRUE, var_thresh = 0.95, min_rblock = 1, max_rblock = Inf,
                   ref.batch = NULL){
  R <- as.matrix(R)
  n <- nrow(R)
  p <- ncol(R)
  site <- as.factor(site)
  pf <- prcomp(R, center = center, scale. = scale_scores)
  r <- pick_r_from_pc(pf, var_thresh, min_r = min_rblock, max_r = max_rblock)
  scores <- pf$x[,1:r]
  scores_com <- com_harm(bat = site, data = scores, covar = NULL, eb = FALSE, ref.batch = ref.batch, cov = FALSE)
  full_scores <- pf$x
  full_scores[,1:r] <- scores_com$harm_data
  if (scale_scores) {
    data_covbat <- full_scores %*% t(pf$rotation) *
      matrix(pf$scale, n, p, byrow = TRUE) +
      matrix(pf$center, n, p, byrow = TRUE)
  } else {
    data_covbat <- full_scores %*% t(pf$rotation) +
      matrix(pf$center, n, p, byrow = TRUE)
  }
  return(data_covbat)
}

