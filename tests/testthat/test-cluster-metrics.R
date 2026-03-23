test_that("cluster metrics return expected values for identical partitions", {
  pred <- c("A", "A", "B", "B")
  ref <- c("A", "A", "B", "B")

  expect_equal(cluster_ari(pred, ref), 1)
  expect_equal(cluster_jaccard(pred, ref), 1)
  expect_equal(cluster_nmi(pred, ref), 1)
  expect_equal(cluster_purity(pred, ref), 1)

  metrics <- cluster_metrics(pred, ref)
  expect_equal(metrics$n, 4)
  expect_equal(metrics$ari, 1)
  expect_equal(metrics$jaccard, 1)
  expect_equal(metrics$nmi, 1)
  expect_equal(metrics$purity, 1)
})

test_that("cluster metrics handle missing labels and aligned named vectors", {
  pred <- c(s1 = "A", s2 = NA, s3 = "B", s4 = "B")
  ref <- c(s1 = "A", s2 = "A", s3 = "B", s4 = "B")

  metrics <- cluster_metrics(pred, ref, na.rm = TRUE)
  expect_equal(metrics$n, 3)
  expect_equal(metrics$ari, 1)

  mats <- cluster_metric_matrix(
    list(
      fit1 = c(s1 = "1", s2 = "1", s3 = "2"),
      fit2 = c(s3 = "2", s2 = "1", s1 = "1")
    ),
    metric = "ari"
  )

  expect_equal(dim(mats), c(2, 2))
  expect_equal(mats["fit1", "fit2"], 1)
  expect_equal(mats["fit2", "fit1"], 1)
})

test_that("cluster metric matrix validates inputs", {
  expect_error(cluster_metric_matrix(list(a = c("1", "2"))), "at least 2")
  expect_error(
    cluster_metric_matrix(list(a = c(x = "1", y = "2"), b = c(z = "1", w = "2"))),
    "common sample ids"
  )
})
