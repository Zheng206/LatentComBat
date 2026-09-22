test_that("residualize_one removes design-space components", {
  x <- seq(-1, 1, length.out = 30)
  Z <- cbind(1, x)
  Y <- cbind(2 + 3 * x, -1 + 0.5 * x) + matrix(rnorm(60, sd = 0.05), 30, 2)
  R <- LatentComBat:::residualize_one(Y, Z)

  expect_lt(max(abs(crossprod(Z, R))), 1e-10)
})

test_that("residualize_one rejects missing values", {
  Y <- matrix(rnorm(20), 10, 2)
  Z <- cbind(1, rnorm(10))
  Y[1, 1] <- NA
  expect_error(LatentComBat:::residualize_one(Y, Z), "NA values")
})

test_that("within-subject permutation preserves each subject-feature multiset", {
  set.seed(5)
  R <- matrix(seq_len(24), 8, 3)
  subj <- rep(c("s1", "s2", "s3", "s4"), each = 2)
  Rp <- LatentComBat:::permute_within_subject_featurewise(R, subj)

  expect_equal(dim(Rp), dim(R))
  for (s in unique(subj)) {
    idx <- which(subj == s)
    for (g in seq_len(ncol(R))) {
      expect_setequal(Rp[idx, g], R[idx, g])
    }
  }
})

test_that("detect_K_one returns zero for a zero residual matrix", {
  R <- matrix(0, 12, 4)
  Z <- matrix(1, 12, 1)
  out <- LatentComBat:::detect_K_one(R, Z, B = 3L, alpha = 0.05)

  expect_identical(out$K, 0L)
  expect_true(all(is.na(out$pvals)))
})

test_that("detect_K_one returns a valid sequential p-value sequence", {
  set.seed(6)
  h <- rnorm(24)
  load <- c(2, 1.5, 1, 0.5, 0.2)
  R <- h %o% load + matrix(rnorm(24 * 5, sd = 0.25), 24, 5)
  Z <- matrix(1, 24, 1)
  out <- LatentComBat:::detect_K_one(R, Z, B = 8L, alpha = 0.2)

  expect_true(out$K >= 0L)
  expect_true(out$K <= 5L)
  expect_true(all(diff(out$pvals) >= -1e-12))
  expect_true(all(out$pvals >= 0 & out$pvals <= 1))
})

test_that("construct_svs_one returns requested number of factors", {
  set.seed(7)
  R <- matrix(rnorm(30 * 12), 30, 12)
  s <- svd(R)
  H <- LatentComBat:::construct_svs_one(R, K = 2L, svd_R = s)

  expect_equal(dim(H), c(30L, 2L))
  expect_identical(colnames(H), c("SV1", "SV2"))
  expect_true(all(is.finite(H)))
})
