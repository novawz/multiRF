test_that("native supervised ytry default uses response dimension", {
  fit_mv_forest_cpp <- get("fit_mv_forest_cpp", envir = asNamespace("multiRF"))
  fit_forest <- get("fit_forest", envir = asNamespace("multiRF"))

  set.seed(1)
  X <- matrix(rnorm(40 * 30), nrow = 40, ncol = 30)
  Y <- matrix(rnorm(40 * 5), nrow = 40, ncol = 5)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  colnames(Y) <- paste0("y", seq_len(ncol(Y)))

  expected_ytry <- ceiling(ncol(Y) / 3)

  cpp_fit <- fit_mv_forest_cpp(
    X = X,
    Y = Y,
    ntree = 5L,
    mtry = 0L,
    ytry = 0L,
    seed = 11L,
    nthread = 1L,
    prox_mode = -1L
  )

  expect_equal(cpp_fit$ytry, expected_ytry)

  r_fit <- fit_forest(
    X = as.data.frame(X),
    Y = as.data.frame(Y),
    type = "regression",
    ntree = 5L,
    ytry = NULL,
    seed = 11L,
    engine = "native",
    proximity = "none"
  )

  expect_equal(r_fit$ytry, expected_ytry)
})
