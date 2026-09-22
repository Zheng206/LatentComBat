make_lc_fixture <- function(n_per_batch = 12L, p = 6L, seed = 123) {
  set.seed(seed)
  n <- 2L * n_per_batch
  bat <- factor(rep(c("A", "B"), each = n_per_batch))
  x <- seq(-1, 1, length.out = n)
  sex <- factor(rep(c("F", "M"), length.out = n))
  covar <- data.frame(x = x, sex = sex)

  beta_x <- seq(0.3, 0.8, length.out = p)
  batch_shift <- seq(0.5, 1.0, length.out = p)
  Y <- matrix(rnorm(n * p, sd = 0.6), nrow = n, ncol = p)
  Y <- Y + x %o% beta_x
  Y[bat == "B", ] <- sweep(Y[bat == "B", , drop = FALSE], 2, batch_shift, "+")
  colnames(Y) <- paste0("feature", seq_len(p))

  list(
    data = Y,
    bat = bat,
    covar = covar,
    formula = y ~ x + sex
  )
}

expect_finite_matrix <- function(x, nrow = NULL, ncol = NULL) {
  expect_true(is.matrix(x))
  if (!is.null(nrow)) expect_equal(base::nrow(x), nrow)
  if (!is.null(ncol)) expect_equal(base::ncol(x), ncol)
  expect_true(all(is.finite(x)))
}
