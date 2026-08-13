test_that("mrf3_tsne validates similarity matrices", {
  expect_error(
    mrf3_tsne(list(dat = matrix(1:6, nrow = 2))),
    "square matrix"
  )
  expect_error(
    mrf3_tsne(list(dat = matrix(0, 3, 3))),
    "positive off-diagonal"
  )

  non_finite <- diag(3)
  non_finite[1, 2] <- NA_real_
  expect_error(mrf3_tsne(list(dat = non_finite)), "finite values")

  negative <- matrix(1, 3, 3)
  negative[1, 2] <- -0.1
  expect_error(mrf3_tsne(list(dat = negative)), "non-negative")
})

test_that("zero-similarity rows are allowed when total edge mass is positive", {
  similarity <- matrix(0, 4, 4)
  similarity[2, 3] <- similarity[3, 2] <- 1

  embedding <- mrf3_tsne(
    list(dat = similarity), max_iter = 10, learning_rate = 10, seed = 31
  )

  expect_true(is.matrix(embedding))
  expect_equal(dim(embedding), c(4L, 2L))
  expect_true(all(is.finite(embedding)))
  expect_equal(as.numeric(colMeans(embedding)), c(0, 0), tolerance = 1e-12)
})

test_that("seed makes mrf3_tsne reproducible without changing caller RNG", {
  similarity <- matrix(c(
    0, 1, 0.5, 0.2,
    1, 0, 0.3, 0.4,
    0.5, 0.3, 0, 0.8,
    0.2, 0.4, 0.8, 0
  ), 4, 4, byrow = TRUE)

  set.seed(90210)
  rng_before <- .Random.seed
  first <- mrf3_tsne(
    list(dat = similarity), max_iter = 25, learning_rate = 10, seed = 17
  )
  expect_identical(.Random.seed, rng_before)
  second <- mrf3_tsne(
    list(dat = similarity), max_iter = 25, learning_rate = 10, seed = 17
  )

  expect_identical(first, second)
  expect_true(all(is.finite(first)))
  expect_true(is.logical(attr(first, "converged")))
  expect_true(is.numeric(attr(first, "cost_history")))
})

test_that("patience counts consecutive stable cost changes", {
  similarity <- matrix(c(
    0, 1, 0.4,
    1, 0, 0.6,
    0.4, 0.6, 0
  ), 3, 3, byrow = TRUE)

  embedding <- mrf3_tsne(
    list(dat = similarity), max_iter = 20, learning_rate = 10,
    tol = 1e9, patience = 2, seed = 4
  )

  expect_true(attr(embedding, "converged"))
  expect_equal(attr(embedding, "iterations"), 3L)
  expect_length(attr(embedding, "cost_history"), 3L)
})

test_that("mrf3_tsne validates optimization parameters", {
  similarity <- matrix(c(0, 1, 1, 0), 2, 2)

  expect_error(mrf3_tsne(list(dat = similarity), dims = 0), "dims")
  expect_error(mrf3_tsne(list(dat = similarity), max_iter = 1.5), "max_iter")
  expect_error(mrf3_tsne(list(dat = similarity), learning_rate = 0), "learning_rate")
  expect_error(mrf3_tsne(list(dat = similarity), patience = 0), "patience")
  expect_error(mrf3_tsne(list(dat = similarity), tol = -1), "tol")
  expect_error(mrf3_tsne(list(dat = similarity), seed = -1), "seed")
  expect_error(mrf3_tsne(list(dat = similarity), verbose = NA), "verbose")
})

test_that("C++ tSNE gradient guards dimensions and invalid probabilities", {
  gradient <- multiRF:::tsne_cost_gradient_cpp
  Y <- matrix(c(-1, 0, 1, 0), 2, 2, byrow = TRUE)
  P <- matrix(c(0, 0.5, 0.5, 0), 2, 2)

  expect_error(gradient(Y, matrix(1, 2, 3)), "square matrix")
  expect_error(gradient(Y, matrix(1, 3, 3)), "square matrix")

  P_na <- P
  P_na[1, 2] <- NA_real_
  expect_error(gradient(Y, P_na), "finite values")

  P_negative <- P
  P_negative[1, 2] <- -0.5
  expect_error(gradient(Y, P_negative), "non-negative")

  expect_error(gradient(Y, matrix(0, 2, 2)), "positive off-diagonal")
  expect_true(all(is.finite(gradient(Y, P))))
})
