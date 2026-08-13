test_that("omics preprocessing enforces sample and numeric invariants", {
  a <- data.frame(`gene A` = 1:5, constant = rep(2, 5), check.names = FALSE)
  b <- data.frame(`protein B` = 6:10, check.names = FALSE)
  rownames(a) <- rownames(b) <- paste0("S", 1:5)

  out <- filter_omics(
    list(rna = a, protein = b),
    filter_mode = "none", return_summary = TRUE, verbose = FALSE
  )
  expect_identical(colnames(out$dat_filtered$rna), "gene A")
  expect_equal(out$filter_summary$n_zero_variance_removed, c(1, 0))

  b_bad_order <- b[5:1, , drop = FALSE]
  expect_error(
    filter_omics(list(rna = a, protein = b_bad_order),
                 filter_mode = "none", verbose = FALSE),
    "Sample order differs"
  )

  b_missing <- b
  rownames(b_missing)[1] <- "other"
  expect_error(
    filter_omics(list(rna = a, protein = b_missing),
                 filter_mode = "none", verbose = FALSE),
    "Sample names differ"
  )

  a_nonfinite <- a
  a_nonfinite[1, 1] <- NA_real_
  expect_error(
    filter_omics(list(rna = a_nonfinite),
                 filter_mode = "none", verbose = FALSE),
    "NA, NaN, or infinite"
  )
})

test_that("native wrapper forwards RF structure and records resolved settings", {
  set.seed(10)
  ids <- paste0("S", 1:36)
  X <- as.data.frame(matrix(rnorm(36 * 6), 36, 6), check.names = FALSE)
  Y <- as.data.frame(matrix(rnorm(36 * 4), 36, 4), check.names = FALSE)
  rownames(X) <- rownames(Y) <- ids
  colnames(X) <- paste0("x", 1:6)
  colnames(Y) <- paste0("y", 1:4)

  fit <- fit_forest(
    X, Y, ntree = 12, mtry = 3, ytry = 2, nsplit = 4,
    nodesize = 4, max_depth = 2, nthread = 1, seed = 91,
    forest.wt = "all", proximity = "none",
    xvar.wt = setNames(c(6, 5, 4, 3, 2, 1), colnames(X)),
    yvar.wt = setNames(c(1, 2, 3, 4), colnames(Y))
  )

  expect_equal(fit$structural_params$ntree, 12L)
  expect_equal(fit$structural_params$nsplit, 4L)
  expect_equal(fit$structural_params$nodesize, 4L)
  expect_equal(fit$structural_params$max_depth, 2L)
  expect_equal(fit$structural_params$nthread, 1L)
  expect_identical(fit$forest.wt.mode, "all")
  expect_equal(sum(fit$xvar.wt), 1)
  expect_equal(sum(fit$yvar.wt), 1)
  expect_true(all(vapply(fit$tree_info, function(ti) max(ti$depth) <= 2L,
                         logical(1))))

  expect_error(fit_forest(X, Y, split.wt = rep(1, ncol(X))),
               "does not yet implement")
  expect_error(fit_forest(X, Y, case.wt = rep(1, nrow(X))),
               "does not yet implement")
})

test_that("native wrappers reject invalid scalar controls before fitting", {
  ids <- paste0("S", seq_len(8L))
  X <- data.frame(x1 = seq_len(8L), x2 = rev(seq_len(8L)), row.names = ids)
  Y <- data.frame(y1 = seq_len(8L)^2, row.names = ids)

  expect_error(fit_forest(X, Y, ntree = 0), "ntree.*positive")
  expect_error(fit_forest(X, Y, ntree = NA_real_), "ntree.*positive")
  expect_error(fit_forest(X, Y, ntree = c(2, 3)), "ntree.*positive")
  expect_error(fit_forest(X, Y, nthread = -1), "nthread.*non-negative")
  expect_error(fit_forest(X, Y, nsplit = 1.5), "nsplit.*integer")
  expect_error(fit_forest(X, Y, nodesize = 0), "nodesize.*positive")
  expect_error(fit_forest(X, Y, seed = 1.5), "seed.*integer")
  expect_error(fit_forest(X, Y, enhanced_prox = NA), "enhanced_prox.*TRUE or FALSE")
  expect_error(fit_forest(X, Y, sibling_gamma = Inf), "sibling_gamma.*finite")
})

test_that("all, inbag, and OOB forest-weight modes have native semantics", {
  compute_oob_forest_wt_cpp <- get(
    "compute_oob_forest_wt_cpp", envir = asNamespace("multiRF")
  )
  set.seed(11)
  ids <- paste0("S", 1:42)
  X <- as.data.frame(matrix(rnorm(42 * 5), 42, 5))
  Y <- as.data.frame(matrix(rnorm(42 * 3), 42, 3))
  rownames(X) <- rownames(Y) <- ids

  args <- list(X = X, Y = Y, ntree = 40, seed = 17, nthread = 1,
               proximity = "none")
  fit_all <- do.call(fit_forest, c(args, list(forest.wt = "all")))
  fit_inbag <- do.call(fit_forest, c(args, list(forest.wt = "inbag")))
  fit_oob <- do.call(fit_forest, c(args, list(forest.wt = "oob")))

  expect_equal(unname(rowSums(fit_all$forest.wt)), rep(1, nrow(X)), tolerance = 1e-12)
  expect_equal(unname(rowSums(fit_inbag$forest.wt)), rep(1, nrow(X)), tolerance = 1e-12)
  expect_equal(
    fit_oob$forest.wt,
    compute_oob_forest_wt_cpp(fit_oob$membership, fit_oob$inbag),
    tolerance = 1e-12
  )
  expect_gt(max(abs(fit_all$forest.wt - fit_inbag$forest.wt)), 0)
})

test_that("OOB NMSE averages all predictor and response coordinates", {
  get_oob_nmse <- get("get_oob_nmse", envir = asNamespace("multiRF"))
  compute_oob_forest_wt <- get(
    "compute_oob_forest_wt", envir = asNamespace("multiRF")
  )
  set.seed(112)
  X <- as.data.frame(matrix(rnorm(48 * 5), 48, 5))
  Y <- as.data.frame(matrix(rnorm(48 * 2), 48, 2))
  fit <- fit_forest(
    X, Y, ntree = 35, seed = 117, nthread = 1,
    proximity = "none", forest.wt = "oob"
  )
  W <- compute_oob_forest_wt(fit)
  nmse_x <- colMeans((W %*% as.matrix(X) - as.matrix(X))^2, na.rm = TRUE) /
    apply(X, 2, stats::var)
  nmse_y <- colMeans((W %*% as.matrix(Y) - as.matrix(Y))^2, na.rm = TRUE) /
    apply(Y, 2, stats::var)

  expect_equal(get_oob_nmse(fit), mean(c(nmse_x, nmse_y)), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(
    get_oob_nmse(fit),
    (mean(nmse_x) + mean(nmse_y)) / 2
  )))
})

test_that("OOB NMSE excludes samples with no OOB prediction", {
  compute_oob_forest_wt <- get(
    "compute_oob_forest_wt", envir = asNamespace("multiRF")
  )
  get_oob_nmse <- get("get_oob_nmse", envir = asNamespace("multiRF"))

  no_oob <- list(
    membership = matrix(c(1L, 1L, 1L), ncol = 1L),
    inbag = matrix(c(1L, 1L, 1L), ncol = 1L),
    xvar = data.frame(x = c(1, 2, 4)),
    yvar = NULL
  )
  expect_true(all(is.na(compute_oob_forest_wt(no_oob))))
  expect_true(is.na(get_oob_nmse(no_oob)))

  some_oob <- list(
    membership = matrix(1L, nrow = 3L, ncol = 2L),
    inbag = matrix(
      c(1L, 1L, 0L,
        1L, 0L, 1L),
      nrow = 3L,
      ncol = 2L
    ),
    xvar = data.frame(x = c(1, 2, 4)),
    yvar = data.frame(y1 = c(2, 5, 7), y2 = c(3, 3, 9))
  )
  W <- compute_oob_forest_wt(some_oob)
  expect_true(all(is.na(W[1L, ])))
  expect_true(all(is.finite(W[2:3, ])))

  X <- as.matrix(some_oob$xvar)
  Y <- as.matrix(some_oob$yvar)
  expected <- mean(c(
    colMeans((W %*% X - X)^2, na.rm = TRUE) / apply(X, 2, stats::var),
    colMeans((W %*% Y - Y)^2, na.rm = TRUE) / apply(Y, 2, stats::var)
  ))
  expect_equal(get_oob_nmse(some_oob), expected, tolerance = 1e-12)
})

test_that("native unsupervised path receives structural settings without fallback", {
  set.seed(12)
  X <- as.data.frame(matrix(rnorm(30 * 7), 30, 7))
  rownames(X) <- paste0("S", 1:30)

  fit <- fit_forest(
    X, Y = NULL, type = "unsupervised", ntree = 9,
    ytry = 3, nsplit = 5, nodesize = 4, max_depth = 2,
    nthread = 1, forest.wt = "all", proximity = "none", seed = 23
  )

  expect_identical(fit$engine, "multiRF")
  expect_equal(dim(fit$inbag), c(nrow(X), 9L))
  expect_equal(fit$structural_params$nsplit, 5L)
  expect_equal(fit$structural_params$nodesize, 4L)
  expect_equal(fit$structural_params$max_depth, 2L)
  expect_identical(fit$structural_params$forest.wt, "all")
})

test_that("multi-forest models retain unambiguous connection metadata", {
  set.seed(13)
  ids <- paste0("S", 1:24)
  dat <- list(
    `rna_block` = data.frame(a = rnorm(24), b = rnorm(24), row.names = ids),
    `protein_block` = data.frame(c = rnorm(24), d = rnorm(24), row.names = ids)
  )
  fit <- fit_multi_forest(
    dat, connect_list = list(c("rna_block", "protein_block")),
    ntree = 5, nthread = 1, proximity = "none", seed = 31
  )
  expect_identical(fit[[1]]$connection$response, "rna_block")
  expect_identical(fit[[1]]$connection$predictor, "protein_block")

  one <- fit_multi_forest(
    dat["rna_block"], connect_list = list("rna_block"),
    ntree = 5, nthread = 1, proximity = "none", seed = 32
  )
  expect_identical(one[[1]]$engine, "multiRF")
  expect_identical(one[[1]]$connection$response, "rna_block")
  expect_null(one[[1]]$connection$predictor)
})

test_that("multi-forest display names are unique when underscores collide", {
  set.seed(1301)
  ids <- paste0("S", seq_len(12L))
  dat <- list(
    a_b = data.frame(v = rnorm(12L), row.names = ids),
    c = data.frame(v = rnorm(12L), row.names = ids),
    a = data.frame(v = rnorm(12L), row.names = ids),
    b_c = data.frame(v = rnorm(12L), row.names = ids)
  )
  connections <- list(c("a_b", "c"), c("a", "b_c"))
  fit <- fit_multi_forest(
    dat, connect_list = connections, ntree = 2L, nthread = 1L,
    max_depth = 1L, proximity = "none", seed = 43L
  )

  expect_identical(names(fit), c("a_b_c", "a_b_c.1"))
  expect_identical(fit[[1L]]$connection,
                   list(response = "a_b", predictor = "c"))
  expect_identical(fit[[2L]]$connection,
                   list(response = "a", predictor = "b_c"))
})

test_that("residual unsupervised forests inherit requested RF structure", {
  get_shared_specific_weights <- get(
    "get_shared_specific_weights", envir = asNamespace("multiRF")
  )
  set.seed(14)
  ids <- paste0("S", 1:24)
  X <- matrix(rnorm(24 * 4), 24, 4,
              dimnames = list(ids, paste0("g", 1:4)))
  W <- diag(0.2, 24) + matrix(0.8 / 24, 24, 24)
  dimnames(W) <- list(ids, ids)
  recon <- list(
    W = list(W_all = W, W_per_response = list(rna = W)),
    fused_mat = list(rna = W %*% X)
  )

  out <- get_shared_specific_weights(
    list(rna = X), recon,
    specific_seed = 41, specific_ntree = 7,
    specific_samptype = "swor", specific_ytry = 2,
    specific_proximity = "none", specific_nsplit = 4,
    specific_nthread = 1, specific_nodesize = 4,
    specific_max_depth = 2, specific_forest_wt = "all"
  )
  params <- out$specific$residual_mod$rna$structural_params
  expect_equal(params$ntree, 7L)
  expect_equal(params$ytry, 2L)
  expect_equal(params$nsplit, 4L)
  expect_equal(params$nodesize, 4L)
  expect_equal(params$max_depth, 2L)
  expect_identical(params$forest.wt, "all")
})

test_that("workflow maps main node settings into residual forests", {
  inherit_args <- get(
    ".inherit_specific_forest_args", envir = asNamespace("multiRF")
  )
  inherited <- inherit_args(
    shared_specific_args = list(specific_nodesize = 9L),
    ntree = 101L, samptype = "swr", ytry = 7L, proximity = "all",
    dots = list(nsplit = 6L, nthread = 2L, nodesize = 5L, nodedepth = 4L)
  )

  expect_equal(inherited$specific_ntree, 101L)
  expect_identical(inherited$specific_samptype, "swr")
  expect_equal(inherited$specific_ytry, 7L)
  expect_equal(inherited$specific_nsplit, 6L)
  expect_equal(inherited$specific_nthread, 2L)
  expect_equal(inherited$specific_max_depth, 4L)
  # Explicit residual overrides take precedence over inherited main settings.
  expect_equal(inherited$specific_nodesize, 9L)
  expect_identical(inherited$specific_forest_wt, "all")
})
