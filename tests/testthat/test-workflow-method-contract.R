test_that("fusion normalizes modularity within response blocks", {
  ids <- paste0("S", 1:3)
  mk <- function(response, predictor, value) {
    W <- diag(value, 3) + matrix((1 - value) / 3, 3, 3)
    dimnames(W) <- list(ids, ids)
    list(
      forest.wt = W,
      xvar = data.frame(x = 1:3, row.names = ids),
      yvar = data.frame(y = 3:1, row.names = ids),
      connection = list(response = response, predictor = predictor)
    )
  }
  models <- list(
    A_B = mk("A", "B", .2),
    A_C = mk("A", "C", .3),
    B_A = mk("B", "A", .4)
  )
  scores <- matrix(NA_real_, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  scores[cbind(c("A", "A", "B"), c("B", "C", "A"))] <- c(9, 1, 1)
  cache <- multiRF:::build_model_weight_cache_list(models, vmax = 3)
  out <- multiRF:::build_response_fused_weight_from_cache(
    cache, top_v = 3, connection_score = scores,
    model_list = models, response_blocks = c("A", "B")
  )
  mats <- lapply(cache, multiRF:::materialize_weight_from_cache, v = 3)
  expected <- 0.5 * (0.9 * mats$A_B + 0.1 * mats$A_C) + 0.5 * mats$B_A
  expect_equal(out, expected)
})

test_that("top-v grid reaches the no-truncation region", {
  expect_equal(multiRF:::infer_fused_tune_vmax(100), 100L)
  expect_gte(multiRF:::infer_fused_tune_vmax(7), ceiling(.8 * 7))
})

test_that("missing response weights fail instead of leaking global signal", {
  ids <- paste0("S", 1:4)
  W <- diag(4)
  dimnames(W) <- list(ids, ids)
  model <- list(
    forest.wt = W,
    xvar = data.frame(x = 1:4, row.names = ids),
    yvar = data.frame(y = 4:1, row.names = ids),
    connection = list(response = "A", predictor = "B")
  )
  expect_error(
    get_reconstr_matrix(
      rfit = list(A_B = model), model_top_v = 4
    ),
    "No fitted connection has block `B` as response"
  )
})

test_that("specific labels are exposed and cluster IMD receives shared labels", {
  ss <- list(
    shared = diag(3), shared_frac = data.frame(),
    specific = list(imd = NULL, imd_per_tree = NULL, A = diag(3))
  )
  sc <- list(
    shared = list(cl = setNames(c(1, 1, 2), paste0("S", 1:3))),
    specific = list(
      method = "Spectral", similarity_type = "second",
      by_omics = list(A = list(cl = setNames(c(2, 1, 2), paste0("S", 1:3))))
    )
  )
  out <- multiRF:::merge_cluster_outputs(ss, sc)
  expect_equal(out$clusters$A, sc$specific$by_omics$A$cl)
})

test_that("connection shorthand detects underscore ambiguity", {
  expect_error(
    multiRF:::normalize_connect_list(
      "a_b_c", valid_names = c("a_b", "c", "a", "b_c"), n_blocks = 4
    ),
    "ambiguous"
  )
  expect_equal(
    multiRF:::normalize_connect_list(
      list(c("a_b", "c")), valid_names = c("a_b", "c"), n_blocks = 2
    ),
    list(c("a_b", "c"))
  )
})

test_that("small-n clustering caps candidate k", {
  S <- diag(5)
  expect_error(
    cluster_similarity_matrix(S, k = 5, method = "Spectral"),
    "smaller than the number of samples"
  )
})

test_that("an entirely absent response block is never dropped from Eq. 8", {
  ids <- paste0("S", 1:4)
  W <- matrix(c(
    .2, .4, .3, .1,
    .3, .2, .4, .1,
    .1, .4, .2, .3,
    .4, .1, .3, .2
  ), 4, byrow = TRUE, dimnames = list(ids, ids))
  model <- list(
    forest.wt = W,
    xvar = data.frame(x = 1:4, row.names = ids),
    yvar = data.frame(y = 4:1, row.names = ids),
    connection = list(response = "A", predictor = "B")
  )
  model_b <- model
  model_b$connection <- list(response = "B", predictor = "A")
  expect_error(
    get_reconstr_matrix(
      rfit = list(A_B = model, B_A = model_b), model_top_v = 4,
      response_blocks = c("A", "B", "C")
    ),
    "block `C` as response"
  )
})

test_that("tuning cache reproduces final Eq. 6-8 reconstruction", {
  ids <- paste0("S", 1:4)
  mk <- function(response, predictor, shift) {
    W <- matrix(c(
      .25, .45, .20, .10,
      .15, .30, .45, .10,
      .10, .35, .25, .30,
      .40, .10, .30, .20
    ), 4, byrow = TRUE)
    W <- W[, ((seq_len(4) + shift - 1L) %% 4L) + 1L, drop = FALSE]
    dimnames(W) <- list(ids, ids)
    list(
      forest.wt = W,
      xvar = data.frame(x = 1:4, row.names = ids),
      yvar = data.frame(y = 4:1, row.names = ids),
      connection = list(response = response, predictor = predictor)
    )
  }
  models <- list(
    A_B = mk("A", "B", 0L),
    A_C = mk("A", "C", 1L),
    B_A = mk("B", "A", 2L),
    B_C = mk("B", "C", 3L)
  )
  scores <- matrix(NA_real_, 3, 3,
                   dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  scores[cbind(c("A", "A", "B", "B"), c("B", "C", "A", "C"))] <- c(4, 1, 3, .5)

  cache <- multiRF:::build_model_weight_cache_list(models, vmax = 2)
  tuned_W <- multiRF:::build_response_fused_weight_from_cache(
    cache_list = cache,
    top_v = 2,
    connection_score = scores,
    model_list = models,
    response_blocks = c("A", "B"),
    recon_fusion = "weighted",
    score_power = 2,
    score_floor = .25
  )
  final <- get_reconstr_matrix(
    rfit = models,
    model_top_v = 2,
    connection_score = scores,
    response_blocks = c("A", "B"),
    recon_fusion = "weighted",
    score_power = 2,
    score_floor = .25
  )

  expect_equal(tuned_W, final$W$W_all, tolerance = 1e-12)
  expect_equal(unname(rowSums(tuned_W)), rep(1, 4), tolerance = 1e-12)
  expect_equal(
    unname(rowSums(multiRF:::materialize_weight_from_cache(cache[[1]], 2))),
    rep(1, 4), tolerance = 1e-12
  )
})

test_that("80-percent rule maps fixed top-v to no truncation", {
  dat <- list(
    A = matrix(0, 5, 1, dimnames = list(paste0("S", 1:5), "a")),
    B = matrix(0, 5, 1, dimnames = list(paste0("S", 1:5), "b"))
  )
  out <- multiRF:::resolve_top_v_values(
    dat_input = dat,
    mod_input = list(),
    connection_input = list(c("A", "B"), c("B", "A")),
    connection_score = NULL,
    model_top_v_input = 4,
    fused_top_v_input = 4,
    verbose = FALSE
  )
  expect_equal(out$model_top_v, 5L)
  expect_null(out$fused_top_v)
})

test_that("sub-MRF inherits the main forest controls unless overridden", {
  inherited <- multiRF:::.inherit_sub_mrf_args(
    sub_mrf_args = list(n_sub = 8L, nodesize = 9L),
    ntree = 101L,
    samptype = "swr",
    ytry = 4L,
    dots = list(nsplit = 17L, nodesize = 3L, mtry = 6L, engine = "native")
  )
  expect_equal(inherited$ntree_full, 101L)
  expect_equal(inherited$ntree_per_sub, 13L)
  expect_equal(inherited$samptype, "swr")
  expect_equal(inherited$ytry, 4L)
  expect_equal(inherited$nsplit, 17L)
  expect_equal(inherited$mtry, 6L)
  expect_equal(inherited$engine, "native")
  expect_equal(inherited$nodesize, 9L)

  # Features absent from a replicate contribute zero to the ensemble mean.
  expect_equal(
    multiRF:::.average_sub_mrf_imd(c(a = .8, b = .2), n_sub = 4L),
    c(a = .2, b = .05)
  )
})

test_that("clustering defaults are spectral and similarities are hollow", {
  expect_identical(eval(formals(cluster_similarity_matrix)$method), c("Spectral", "PAM"))
  expect_identical(eval(formals(mrf3)$filter_mode), c("auto", "none", "manual"))

  W <- matrix(c(
    .2, .5, .2, .1,
    .4, .2, .3, .1,
    .1, .4, .2, .3,
    .3, .1, .4, .2
  ), 4, byrow = TRUE)
  S <- multiRF:::build_similarity_from_weights(W)
  expect_equal(diag(S), rep(0, 4))
  expect_equal(S, t(S))
  expect_equal(multiRF:::calc_modularity(diag(4)), 0)
})

test_that("direct clustering helpers filter invalid small-n k candidates", {
  S <- matrix(c(
    0, .9, .2, .1, .1,
    .9, 0, .2, .1, .1,
    .2, .2, 0, .8, .7,
    .1, .1, .8, 0, .8,
    .1, .1, .7, .8, 0
  ), 5, byrow = TRUE)
  sp <- spectral_cl(S)
  pa <- pam_cl(1 - S)
  expect_true(sp$best_k >= 2L && sp$best_k < 5L)
  expect_true(pa$best_k >= 2L && pa$best_k < 5L)
})

test_that("one-block reconstruction retains explicit unsupervised metadata", {
  ids <- paste0("S", 1:4)
  W <- matrix(c(
    .2, .4, .3, .1,
    .3, .2, .4, .1,
    .1, .4, .2, .3,
    .4, .1, .3, .2
  ), 4, byrow = TRUE, dimnames = list(ids, ids))
  model <- list(
    forest.wt = W,
    xvar = data.frame(x = 1:4, row.names = ids),
    yvar = NULL,
    connection = list(response = "rna_seq", predictor = NULL)
  )
  out <- get_reconstr_matrix(
    rfit = list(legacy_display_name = model),
    model_top_v = 4,
    response_blocks = "rna_seq"
  )
  expect_identical(names(out$fused_mat), "rna_seq")
  expect_identical(names(out$W$W_per_response), "rna_seq")
})
