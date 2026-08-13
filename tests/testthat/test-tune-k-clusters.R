make_tune_similarity <- function(cluster_sizes = c(3L, 3L)) {
  membership <- rep(seq_along(cluster_sizes), cluster_sizes)
  out <- outer(
    membership,
    membership,
    function(a, b) ifelse(a == b, 0.9, 0.1)
  )
  diag(out) <- 1
  out
}

test_that("tune_k_clusters dispatches matrices to the default method", {
  similarity <- make_tune_similarity()

  tuned_k <- tune_k_clusters(similarity, method = "Spectral")
  tuned <- tune_k_clusters(
    similarity,
    method = "Spectral",
    return_cluster = TRUE
  )

  expect_true(is.function(getS3method("tune_k_clusters", "default")))
  expect_true(tuned_k %in% 2:5)
  expect_identical(tuned_k, tuned$best_k)
  expect_length(tuned$cl, nrow(similarity))
  expect_identical(tuned$k_tune, 2:5)
})

test_that("mrf3 dispatch reuses the default tuning path", {
  similarity <- make_tune_similarity()
  object <- structure(list(dat = similarity), class = c("mrf3", "cl"))

  from_matrix <- tune_k_clusters(
    similarity,
    method = "Spectral",
    return_cluster = TRUE
  )
  from_mrf3 <- tune_k_clusters(
    object,
    method = "Spectral",
    return_cluster = TRUE
  )

  expect_true(is.function(getS3method("tune_k_clusters", "mrf3")))
  expect_identical(from_mrf3$best_k, from_matrix$best_k)
  expect_identical(from_mrf3$cl, from_matrix$cl)
  expect_error(
    multiRF:::tune_k_clusters(structure(list(), class = "mrf3")),
    "containing a `dat` matrix",
    fixed = TRUE
  )
})

test_that("Spectral and PAM tuning use the observed sample size", {
  similarity <- make_tune_similarity()

  spectral <- multiRF:::tune_k_clusters(
    similarity,
    method = "Spectral",
    return_cluster = TRUE
  )
  pam <- multiRF:::tune_k_clusters(
    similarity,
    method = "PAM",
    prox = TRUE,
    return_cluster = TRUE
  )

  expect_identical(spectral$k_tune, 2:5)
  expect_identical(pam$k_tune, 2:5)
  expect_true(spectral$best_k %in% spectral$k_tune)
  expect_true(pam$best_k %in% pam$k_tune)
  expect_length(spectral$cl, 6L)
  expect_length(pam$cl, 6L)
})

test_that("PAM silhouette tuning uses average silhouette width", {
  set.seed(1)
  x <- rbind(
    cbind(stats::rnorm(20, -4, 0.35), stats::rnorm(20, 0, 0.35)),
    cbind(stats::rnorm(20, 0, 0.35), stats::rnorm(20, 4, 0.35)),
    cbind(stats::rnorm(20, 4, 0.35), stats::rnorm(20, 0, 0.35))
  )
  tuned <- pam_cl(x, k_tune = 2:6, diss = FALSE,
                  tune_method = "silhouette")
  d <- stats::dist(x)
  expected <- vapply(2:6, function(k) {
    cluster::pam(d, k = k, diss = TRUE)$silinfo$avg.width
  }, numeric(1))

  expect_equal(tuned$sil, expected)
  expect_identical(tuned$best_k, (2:6)[which.max(expected)])
})

test_that("three samples provide the single valid k candidate", {
  similarity <- matrix(0.2, nrow = 3L, ncol = 3L)
  diag(similarity) <- 1

  spectral <- multiRF:::tune_k_clusters(
    similarity,
    method = "Spectral",
    return_cluster = TRUE
  )
  pam <- multiRF:::tune_k_clusters(
    similarity,
    method = "PAM",
    prox = TRUE,
    return_cluster = TRUE
  )

  expect_identical(spectral$best_k, 2L)
  expect_identical(spectral$k_tune, 2L)
  expect_identical(pam$best_k, 2L)
  expect_identical(pam$k_tune, 2L)
  expect_length(pam$sil, 1L)

  expect_error(
    multiRF:::tune_k_clusters(diag(2L)),
    "At least three samples"
  )
})

test_that("plot_k uses the evaluated k grid and restores graphics state", {
  similarity <- make_tune_similarity(c(4L, 4L))
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)

  graphics::par(mfrow = c(2L, 2L), mar = c(3, 3, 2, 1), oma = c(1, 0, 0, 0))
  before <- graphics::par(c("mfrow", "mar", "oma", "xpd"))

  expect_error(
    multiRF:::tune_k_clusters(
      similarity,
      method = "Spectral",
      k_tune = c(2L, 4L, 6L),
      plot_k = TRUE
    ),
    NA
  )
  expect_equal(graphics::par(c("mfrow", "mar", "oma", "xpd")), before)

  expect_error(
    multiRF:::tune_k_clusters(
      similarity,
      method = "PAM",
      tune_method = "ratio",
      prox = TRUE,
      k_tune = c(2L, 4L, 6L),
      plot_k = TRUE
    ),
    NA
  )
  expect_equal(graphics::par(c("mfrow", "mar", "oma", "xpd")), before)
})

test_that("tuning arguments fail with targeted messages", {
  similarity <- make_tune_similarity()

  expect_error(
    multiRF:::tune_k_clusters(similarity, method = "ward"),
    "`method` must be one of"
  )
  expect_error(
    multiRF:::tune_k_clusters(similarity, tune_method = "gap"),
    "`tune_method` must be one of"
  )
  expect_error(
    multiRF:::tune_k_clusters(similarity, gap_w = "sqrt"),
    "`gap_w` must be one of"
  )
  expect_error(
    multiRF:::tune_k_clusters(similarity, plot_k = NA),
    "`plot_k` must be TRUE or FALSE"
  )
  expect_error(
    multiRF:::tune_k_clusters(similarity, prox = 1),
    "`prox` must be TRUE or FALSE"
  )
  expect_error(
    multiRF:::tune_k_clusters(similarity, k_tune = c(2, 2.5)),
    "must be an integer"
  )
})
