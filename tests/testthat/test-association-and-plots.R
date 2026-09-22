test_that("get_r2 handles perfect, constant, and incomplete associations", {
  x <- 1:10
  expect_equal(get_r2(2 * x + 1, x), 1)
  expect_true(is.na(get_r2(rnorm(10), rep(1, 10))))
  expect_true(is.na(get_r2(c(1, NA), c(1, 2))))
})

test_that("make_model_matrix_block expands mixed metadata", {
  df <- data.frame(
    age = 1:6,
    sex = factor(rep(c("F", "M"), 3)),
    site = c("a", "a", "b", "b", "c", "c")
  )
  X <- make_model_matrix_block(df)
  expect_equal(nrow(X), 6L)
  expect_true(all(vapply(as.data.frame(X), is.numeric, logical(1))))
})

test_that("pc_assoc_table stores PCA metadata and groups", {
  fx <- make_lc_fixture(n_per_batch = 10L, p = 6L)
  meta <- data.frame(batch = fx$bat, x = fx$covar$x)
  groups <- c(batch = "Technical", x = "Biological")
  out <- pc_assoc_table(fx$data, meta, n_pc = 3, variable_groups = groups)

  expect_equal(nrow(out), 6L)
  expect_setequal(unique(out$Panel), c("Technical", "Biological"))
  expect_length(attr(out, "pc_variance"), 3L)
  expect_false(is.null(attr(out, "pca")))
})

test_that("pc_assoc_table validates dimensions and group names", {
  fx <- make_lc_fixture(n_per_batch = 6L, p = 4L)
  expect_error(
    pc_assoc_table(fx$data, data.frame(x = 1:5)),
    "same number of rows"
  )
  expect_error(
    pc_assoc_table(fx$data, data.frame(x = fx$covar$x), variable_groups = "Metadata"),
    "named vector"
  )
})

test_that("CCA diagnostic returns bounded canonical correlations", {
  set.seed(9)
  H <- matrix(rnorm(60), 20, 3)
  X <- data.frame(x = rnorm(20), grp = factor(rep(c("A", "B"), 10)))
  out <- cca_diagnostic(H, X)

  expect_true(all(out$cor >= 0 & out$cor <= 1 + 1e-12))
  expect_equal(nrow(out$U), 20L)
  expect_equal(nrow(out$V), 20L)
})

test_that("plot helpers return ggplot objects", {
  fx <- make_lc_fixture(n_per_batch = 8L, p = 5L)
  assoc <- pc_assoc_table(
    fx$data,
    data.frame(batch = fx$bat, x = fx$covar$x),
    n_pc = 2,
    variable_groups = c(batch = "Technical", x = "Biological")
  )
  expect_s3_class(plot_pc_assoc_heatmap(assoc), "ggplot")

  set.seed(10)
  H <- matrix(rnorm(48), 16, 3)
  cca1 <- cca_diagnostic(H, data.frame(x = rnorm(16), z = rnorm(16)))
  cca2 <- cca_diagnostic(H, data.frame(a = rnorm(16), b = rnorm(16)))
  expect_s3_class(plot_cca_cor(cca1, cca2), "ggplot")
  expect_s3_class(plot_cca_scores(cca1, factor(rep(c("A", "B"), 8))), "ggplot")
})
