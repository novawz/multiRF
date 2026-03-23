test_that("native classification fit returns class-compatible fields", {
  fit_forest <- get("fit_forest", envir = asNamespace("multiRF"))

  set.seed(1)
  X <- as.data.frame(matrix(rnorm(30 * 5), nrow = 30, ncol = 5))
  Y <- factor(rep(c("A", "B", "C"), each = 10))

  fit <- fit_forest(
    X = X,
    Y = Y,
    type = "classification",
    ntree = 10L,
    seed = 11L,
    engine = "native"
  )

  expect_equal(class(fit)[3], "class+")
  expect_equal(dim(fit$membership), c(nrow(X), fit$ntree))
  expect_equal(length(fit$predicted), nrow(X))
  expect_equal(nrow(fit$class.prob), nrow(X))
  expect_equal(colnames(fit$class.prob), levels(Y))
  expect_true(is.numeric(fit$err.rate))
  expect_true(is.finite(fit$err.rate))
})
