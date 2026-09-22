test_that("covbat preserves dimensions and finite values", {
  set.seed(11)
  n <- 30L
  p <- 5L
  site <- factor(rep(c("A", "B"), each = n / 2))
  R <- matrix(rnorm(n * p), n, p)
  R[site == "B", 1:2] <- R[site == "B", 1:2] + 0.8

  out <- suppressWarnings(covbat(
    R = R,
    site = site,
    var_thresh = 0.8,
    min_rblock = 1,
    max_rblock = 2
  ))

  expect_equal(dim(out), dim(R))
  expect_true(all(is.finite(out)))
})

test_that("pick_r_from_pc rejects zero-variance PCA fits", {
  expect_error(
    pick_r_from_pc(list(sdev = c(0, 0))),
    "no positive variance"
  )
})
