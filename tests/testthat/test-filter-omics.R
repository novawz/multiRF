test_that("filter_omics preserves names and produces summaries", {
  dat_list <- list(
    mRNA = data.frame(g1 = c(1, 2, 3), g2 = c(1, 5, 9), g3 = c(2, 2, 2)),
    miRNA = data.frame(m1 = c(0, 0, 1), m2 = c(1, 3, 6), m3 = c(2, 2, 2))
  )

  out <- filter_omics(
    dat.list = dat_list,
    filter_mode = "manual",
    top_n_manual = c(mRNA = 2, miRNA = 1),
    return_summary = TRUE,
    verbose = FALSE
  )

  expect_named(out$dat_filtered, c("mRNA", "miRNA"))
  expect_equal(ncol(out$dat_filtered$mRNA), 2)
  expect_equal(ncol(out$dat_filtered$miRNA), 1)
  expect_equal(out$filter_summary$n_features_in, c(3, 3))
  expect_equal(out$filter_summary$n_features_out, c(2, 1))
})

test_that("filter_omics none mode keeps all features", {
  dat_list <- list(
    expr = data.frame(a = c(1, 2), b = c(3, 4))
  )

  out <- filter_omics(
    dat.list = dat_list,
    filter_mode = "none",
    return_summary = FALSE,
    verbose = FALSE
  )

  expect_equal(out$expr, dat_list$expr)
})

test_that("filter_omics rejects non-numeric blocks", {
  dat_list <- list(
    bad = data.frame(x = c("a", "b"), y = c("c", "d"))
  )

  expect_error(
    filter_omics(dat.list = dat_list, filter_mode = "none", verbose = FALSE),
    "must be numeric"
  )
})
