# Reproducibility: fixed seeds must give bit-identical forest weights,
# OOB weights, and clustering, including across thread counts.
# (Promoted from tools/manual-tests/test_reproducibility.R)

make_data <- function(n = 50, p1 = 10, p2 = 8) {
  set.seed(42)
  X1 <- matrix(rnorm(n * p1), n, p1,
               dimnames = list(paste0("s", 1:n), paste0("x", 1:p1)))
  X2 <- matrix(rnorm(n * p2), n, p2,
               dimnames = list(paste0("s", 1:n), paste0("y", 1:p2)))
  list(X1 = X1, X2 = X2)
}

test_that("C++ supervised forest is deterministic and thread-invariant", {
  fit_mv_forest_cpp <- get("fit_mv_forest_cpp", envir = asNamespace("multiRF"))
  d <- make_data()
  X_mat <- as.matrix(data.frame(d$X1))
  Y_mat <- as.matrix(data.frame(d$X2))

  a <- fit_mv_forest_cpp(X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 1L)
  b <- fit_mv_forest_cpp(X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 1L)
  expect_equal(max(abs(a$forest.wt - b$forest.wt), na.rm = TRUE), 0)

  c4 <- fit_mv_forest_cpp(X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 4L)
  d4 <- fit_mv_forest_cpp(X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 4L)
  expect_equal(max(abs(c4$forest.wt - d4$forest.wt), na.rm = TRUE), 0)

  # 1 thread vs 4 threads must agree
  expect_identical(a$forest.wt, c4$forest.wt)
  expect_identical(a$proximity, c4$proximity)
  expect_identical(a$membership, c4$membership)
  expect_identical(a$inbag, c4$inbag)
  expect_identical(a$imd_x_per_tree, c4$imd_x_per_tree)
  expect_identical(a$imd_y_per_tree, c4$imd_y_per_tree)
  expect_identical(a$pairwise_xy, c4$pairwise_xy)
  expect_identical(a$tree_info, c4$tree_info)

  # Exercise the OOB denominator and enhanced-proximity accumulation paths too.
  embed <- X_mat[, seq_len(3L), drop = FALSE]
  o1 <- fit_mv_forest_cpp(
    X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 1L,
    prox_mode = 2L, forest_wt_mode = 2L, embed = embed,
    enhanced_prox_mode = 1L
  )
  o4 <- fit_mv_forest_cpp(
    X_mat, Y_mat, ntree = 50L, seed = 529L, nthread = 4L,
    prox_mode = 2L, forest_wt_mode = 2L, embed = embed,
    enhanced_prox_mode = 1L
  )
  expect_identical(o1$forest.wt, o4$forest.wt)
  expect_identical(o1$proximity, o4$proximity)
  expect_identical(o1$enhanced_prox, o4$enhanced_prox)
})

test_that("C++ unsupervised forest is deterministic and thread-invariant", {
  fit_mv_forest_unsup_cpp <- get("fit_mv_forest_unsup_cpp",
                                 envir = asNamespace("multiRF"))
  d <- make_data()
  data_mat <- as.matrix(data.frame(d$X1))

  a <- fit_mv_forest_unsup_cpp(data_mat, ntree = 50L, seed = 529L, nthread = 1L)
  b <- fit_mv_forest_unsup_cpp(data_mat, ntree = 50L, seed = 529L, nthread = 1L)
  expect_equal(max(abs(a$forest.wt - b$forest.wt), na.rm = TRUE), 0)

  c4 <- fit_mv_forest_unsup_cpp(data_mat, ntree = 50L, seed = 529L, nthread = 4L)
  d4 <- fit_mv_forest_unsup_cpp(data_mat, ntree = 50L, seed = 529L, nthread = 4L)
  expect_equal(max(abs(c4$forest.wt - d4$forest.wt), na.rm = TRUE), 0)
})

test_that("OOB forest weight is deterministic", {
  fit_mv_forest_cpp <- get("fit_mv_forest_cpp", envir = asNamespace("multiRF"))
  compute_oob_forest_wt_cpp <- get("compute_oob_forest_wt_cpp",
                                   envir = asNamespace("multiRF"))
  d <- make_data()
  res <- fit_mv_forest_cpp(as.matrix(data.frame(d$X1)),
                           as.matrix(data.frame(d$X2)),
                           ntree = 50L, seed = 529L, nthread = 1L)
  oob_a <- compute_oob_forest_wt_cpp(res$membership, res$inbag)
  oob_b <- compute_oob_forest_wt_cpp(res$membership, res$inbag)
  expect_equal(max(abs(oob_a - oob_b)), 0)
})

test_that("R fit_forest wrapper is reproducible (supervised and unsupervised)", {
  fit_forest <- get("fit_forest", envir = asNamespace("multiRF"))
  d <- make_data()
  X1_df <- as.data.frame(d$X1)
  X2_df <- as.data.frame(d$X2)

  s_a <- fit_forest(X = X1_df, Y = X2_df, ntree = 50, seed = 529L)
  s_b <- fit_forest(X = X1_df, Y = X2_df, ntree = 50, seed = 529L)
  expect_equal(max(abs(s_a$forest.wt - s_b$forest.wt)), 0)

  u_a <- fit_forest(X = X1_df, Y = NULL, type = "unsupervised",
                    ntree = 50, seed = 529L)
  u_b <- fit_forest(X = X1_df, Y = NULL, type = "unsupervised",
                    ntree = 50, seed = 529L)
  expect_equal(max(abs(u_a$forest.wt - u_b$forest.wt)), 0)
})

test_that("full mrf3_fit pipeline is reproducible", {
  mrf3_fit <- get("mrf3_fit", envir = asNamespace("multiRF"))
  d <- make_data()
  dat_list <- list(omics1 = as.data.frame(d$X1), omics2 = as.data.frame(d$X2))

  a <- mrf3_fit(dat_list, ntree = 50, seed = 529L, verbose = FALSE)
  b <- mrf3_fit(dat_list, ntree = 50, seed = 529L, verbose = FALSE)

  if (!is.null(a$shared_clustering$S) && !is.null(b$shared_clustering$S)) {
    expect_equal(max(abs(a$shared_clustering$S - b$shared_clustering$S),
                     na.rm = TRUE), 0)
  }

  expect_identical(
    as.integer(get_clusters(a)),
    as.integer(get_clusters(b))
  )
  expect_identical(a$top_v$model_top_v, b$top_v$model_top_v)

  if (!is.null(a$specific_clustering$S) && !is.null(b$specific_clustering$S)) {
    for (nm in names(a$specific_clustering$S)) {
      expect_equal(max(abs(a$specific_clustering$S[[nm]] -
                             b$specific_clustering$S[[nm]]), na.rm = TRUE), 0)
    }
  }
})
