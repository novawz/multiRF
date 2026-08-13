test_that("AO embedding corresponds to the final fusion weights", {
  empty <- matrix(0, nrow = 6L, ncol = 6L)
  W1 <- empty
  W1[1:3, 1L] <- 1
  W1[4:6, 2L] <- 1
  W2 <- empty
  W2[c(1L, 2L, 4L), 1L] <- 1
  W2[c(3L, 5L, 6L), 2L] <- 1
  W_models <- list(W1, W2)

  fit <- multiRF:::ao_fuse_similarity(
    W_models = W_models,
    k = 2L,
    gamma = 0.01,
    alpha_init = c(0.6, 0.4),
    max_iter = 1L
  )

  # One iteration deliberately moves alpha, exposing a stale-U return.
  expect_gt(sum(abs(fit$alpha - c(0.6, 0.4))), 0.5)

  laplacians <- lapply(
    W_models,
    function(W) multiRF:::build_similarity_laplacian(W)$L
  )
  L_final <- Reduce("+", Map(`*`, fit$alpha, laplacians))
  eig <- eigen(L_final, symmetric = TRUE)
  expected_U <- eig$vectors[
    , order(eig$values)[seq_len(2L)], drop = FALSE
  ]

  # Compare eigenspaces rather than eigenvectors because signs are arbitrary.
  expect_equal(
    fit$U %*% t(fit$U),
    expected_U %*% t(expected_U),
    tolerance = 1e-12
  )
})
