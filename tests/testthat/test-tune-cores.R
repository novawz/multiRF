test_that("sanitize_mc_cores always returns a valid positive core count", {
  skip_if_not_installed("pkgload")
  pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)

  sanitize_mc_cores <- get("sanitize_mc_cores", envir = asNamespace("multiRF"))

  expect_gte(sanitize_mc_cores(NULL), 1L)
  expect_gte(sanitize_mc_cores(NA_integer_), 1L)
  expect_gte(sanitize_mc_cores(0L), 1L)
  expect_gte(sanitize_mc_cores(-5L), 1L)
  expect_identical(sanitize_mc_cores(1L), 1L)
})
