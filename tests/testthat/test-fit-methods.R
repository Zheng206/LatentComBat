test_that("summary and print methods expose fit metadata", {
  x <- structure(
    list(
      harm_data = matrix(1:12, 4, 3),
      method = "sequential",
      latent = NULL,
      protection = list(rank = 2L),
      variance = list(K_var = 1L)
    ),
    class = c("latentcombat_sequential", "latentcombat_fit", "list")
  )

  expect_output(print(x), "LatentComBat fit")
  s <- summary(x)
  expect_s3_class(s, "summary_latentcombat_fit")
  expect_identical(s$method, "sequential")
  expect_equal(c(s$n, s$p), c(4L, 3L))
  expect_equal(s$protection_rank, 2L)
  expect_true(s$has_variance_stabilization)
  expect_output(print(s), "LatentComBat summary")
})
