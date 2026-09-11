brute_entropy <- function(W, v) {
  multiRF:::calc_fused_weight_entropy(
    multiRF:::postprocess_fused_weight(
      W, top_v = v, row_normalize = TRUE, keep_ties = TRUE
    )
  )
}

test_that("closed-form entropy curve matches truncate-and-recompute", {
  set.seed(11)
  n <- 40
  W <- matrix(runif(n * n), n, n)
  diag(W) <- 0
  W[W < 0.6] <- 0
  W[c(3, 9), ] <- 0
  W[5, ] <- sample(c(0, 0.1, 0.2, 0.5), n, TRUE)

  curve <- multiRF:::fused_entropy_curve(W)
  brute <- vapply(seq_len(n), function(v) brute_entropy(W, v), numeric(1))
  expect_equal(curve$v, seq_len(n))
  expect_equal(curve$entropy, brute, tolerance = 1e-12)
  expect_equal(curve$entropy_inf, brute_entropy(W, NULL), tolerance = 1e-12)
  expect_null(names(curve$entropy))
  expect_true(all(diff(curve$entropy) >= -1e-12))
})

test_that("closed-form entropy curve keeps tied neighbours at the cutoff", {
  W <- rbind(
    c(0, 0.5, 0.5, 0.2),
    c(0.3, 0, 0.3, 0.3),
    c(0.1, 0.1, 0, 0.1),
    c(0.9, 0, 0, 0)
  )
  curve <- multiRF:::fused_entropy_curve(W)
  brute <- vapply(1:4, function(v) brute_entropy(W, v), numeric(1))
  expect_equal(curve$entropy, brute, tolerance = 1e-12)
  # Row 1 truncated to top-1 keeps both tied 0.5 entries, so it is not 0.
  expect_gt(curve$entropy[1], 0)
})

test_that("saturation rule selects the first v reaching tau and interpolates", {
  df <- data.frame(
    obj = NA_real_, sil = NA_real_, diffe = NA_real_,
    entropy = c(0.20, 0.45, 0.62, 0.74, 0.80, 0.83),
    fused_top_v = c(10, 20, 30, 40, 50, Inf),
    k = NA_integer_,
    is_no_trunc = c(rep(FALSE, 5), TRUE)
  )
  out <- multiRF:::select_saturation_summary(df, "fused_top_v", tau = 0.9)
  expect_equal(nrow(out), 1L)
  # target 0.747 lies between v = 40 (0.74) and v = 50 (0.80): 41.2 -> 42
  expect_equal(out$fused_top_v, 42)
  expect_equal(out$saturation_rule, "interpolated")
  expect_true(out$saturation_interpolated)
  expect_equal(out$saturation_target, 0.9 * 0.83)
  expect_equal(out$entropy, 0.752)
  expect_equal(out$entropy_frac, 0.752 / 0.83)

  exact <- multiRF:::select_saturation_summary(
    transform(df, fused_top_v = c(1:5, Inf)), "fused_top_v", tau = 0.9
  )
  expect_equal(exact$fused_top_v, 5)
  expect_equal(exact$saturation_rule, "exact")

  first <- multiRF:::select_saturation_summary(df, "fused_top_v", tau = 0.2)
  expect_equal(first$fused_top_v, 10)
  expect_equal(first$saturation_rule, "first_candidate")

  unreached <- multiRF:::select_saturation_summary(df, "fused_top_v", tau = 0.99)
  expect_equal(unreached$fused_top_v, Inf)
  expect_equal(unreached$saturation_rule, "not_reached")

  flat <- multiRF:::select_saturation_summary(
    transform(df, entropy = 0), "fused_top_v", tau = 0.9
  )
  expect_equal(flat$fused_top_v, Inf)
  expect_equal(flat$saturation_rule, "flat_entropy")

  expect_error(
    multiRF:::select_saturation_summary(df, "fused_top_v", tau = 1),
    "tau"
  )
})

test_that("saturation selection is invariant to the candidate grid", {
  v <- 1:200
  entropy <- 0.8 * (1 - exp(-v / 30))
  make_df <- function(keep) {
    data.frame(
      obj = NA_real_, sil = NA_real_, diffe = NA_real_,
      entropy = c(entropy[keep], 0.8),
      model_top_v = c(v[keep], 200),
      k = NA_integer_
    )
  }
  fine <- multiRF:::select_saturation_summary(make_df(v), "model_top_v", tau = 0.9)
  coarse <- multiRF:::select_saturation_summary(
    make_df(seq(1, 200, by = 25)), "model_top_v", tau = 0.9
  )
  expect_equal(fine$saturation_rule, "exact")
  expect_equal(coarse$saturation_rule, "interpolated")
  expect_lte(abs(fine$model_top_v - coarse$model_top_v), 3)
})

test_that("tune_fused_top_v uses the closed-form curve at integer resolution", {
  ids <- paste0("S", 1:30)
  set.seed(3)
  mk <- function(response, predictor) {
    W <- matrix(runif(900), 30, 30)
    diag(W) <- 0
    dimnames(W) <- list(ids, ids)
    list(
      forest.wt = W,
      xvar = data.frame(x = rnorm(30), row.names = ids),
      yvar = data.frame(y = rnorm(30), row.names = ids),
      connection = list(response = response, predictor = predictor)
    )
  }
  models <- list(A_B = mk("A", "B"), B_A = mk("B", "A"))
  dat <- list(
    A = matrix(rnorm(60), 30, 2, dimnames = list(ids, c("a1", "a2"))),
    B = matrix(rnorm(60), 30, 2, dimnames = list(ids, c("b1", "b2")))
  )
  mod <- list(
    mod = models,
    connection = list(c("A", "B"), c("B", "A")),
    connection_score = NULL,
    recon_fusion = "uniform",
    response_blocks = c("A", "B")
  )

  out <- tune_fused_top_v(dat, mod, vmin = 5, model_top_v = 30, parallel = FALSE)
  expect_equal(nrow(out$object), (30 - 5 + 1) + 1)
  expect_true(all(diff(out$object$entropy[is.finite(out$object$fused_top_v)]) >= -1e-12))
  expect_true(all(c("entropy_frac", "is_no_trunc") %in% names(out$object)))
  expect_true(out$vtop_tb$fused_top_v >= 5)
  expect_true(!is.null(out$vtop_tb$saturation_rule))

  elbow <- tune_fused_top_v(
    dat, mod, vmin = 5, model_top_v = 30, parallel = FALSE,
    object = "entropy_elbow"
  )
  expect_true(!is.null(elbow$vtop_tb$elbow_rule))
})

test_that("top-v tuning uses the same predictor fallback as reconstruction", {
  ids <- paste0("S", 1:12)
  set.seed(17)
  W <- matrix(runif(12 * 12), 12, 12, dimnames = list(ids, ids))
  diag(W) <- 0
  model <- list(
    forest.wt = W,
    xvar = data.frame(x = rnorm(12), row.names = ids),
    yvar = data.frame(y = rnorm(12), row.names = ids),
    connection = list(response = "A", predictor = "B")
  )
  models <- list(A_B = model)
  dat <- list(
    A = matrix(rnorm(24), 12, 2, dimnames = list(ids, c("a1", "a2"))),
    B = matrix(rnorm(24), 12, 2, dimnames = list(ids, c("b1", "b2")))
  )
  mod <- list(
    mod = models,
    connection = list(c("A", "B")),
    connection_score = NULL,
    recon_fusion = "uniform",
    response_blocks = names(dat)
  )

  model_tune <- tune_model_top_v(
    dat, mod, tmin = 2, by = 2, max_candidates = 10,
    parallel = FALSE, object = "saturation"
  )
  model_entropy <- vapply(model_tune$object$model_top_v, function(v) {
    recon <- get_reconstr_matrix(
      models, model_top_v = v, recon_fusion = "uniform",
      response_blocks = names(dat)
    )
    multiRF:::calc_fused_weight_entropy(recon$W$W_all)
  }, numeric(1))
  expect_equal(model_tune$object$entropy, model_entropy, tolerance = 1e-12)

  fused_tune <- tune_fused_top_v(
    dat, mod, vmin = 2, model_top_v = 12,
    parallel = FALSE, object = "saturation"
  )
  fused_entropy <- vapply(fused_tune$object$fused_top_v, function(v) {
    recon <- get_reconstr_matrix(
      models, model_top_v = 12, recon_fusion = "uniform",
      response_blocks = names(dat),
      fused_top_v = if (is.finite(v)) v else NULL
    )
    multiRF:::calc_fused_weight_entropy(recon$W$W_all)
  }, numeric(1))
  expect_equal(fused_tune$object$entropy, fused_entropy, tolerance = 1e-12)
  expect_identical(mod$connection, list(c("A", "B")))
})

test_that("model top-v tuning honors pmin global fusion", {
  ids <- paste0("S", 1:10)
  set.seed(23)
  mk <- function(response, predictor) {
    W <- matrix(runif(100), 10, 10, dimnames = list(ids, ids))
    diag(W) <- 0
    list(
      forest.wt = W,
      xvar = data.frame(x = rnorm(10), row.names = ids),
      yvar = data.frame(y = rnorm(10), row.names = ids),
      connection = list(response = response, predictor = predictor)
    )
  }
  models <- list(A_B = mk("A", "B"), B_A = mk("B", "A"))
  dat <- list(
    A = matrix(rnorm(20), 10, 2, dimnames = list(ids, c("a1", "a2"))),
    B = matrix(rnorm(20), 10, 2, dimnames = list(ids, c("b1", "b2")))
  )
  mod <- list(
    mod = models,
    connection = list(c("A", "B"), c("B", "A")),
    recon_fusion = "uniform",
    global_fusion = "pmin",
    response_blocks = names(dat)
  )

  tuned <- tune_model_top_v(
    dat, mod, tmin = 2, by = 2, parallel = FALSE, object = "saturation"
  )
  expected <- vapply(tuned$object$model_top_v, function(v) {
    recon <- get_reconstr_matrix(
      models, model_top_v = v, recon_fusion = "uniform",
      global_fusion = "pmin", response_blocks = names(dat)
    )
    multiRF:::calc_fused_weight_entropy(recon$W$W_all)
  }, numeric(1))
  expect_equal(tuned$object$entropy, expected, tolerance = 1e-12)
  expect_error(
    tune_fused_top_v(
      dat, mod, vmin = 2, model_top_v = 10, parallel = FALSE
    ),
    "not used"
  )
})

test_that("workflow exposes the top-v method and tau in config", {
  expect_equal(
    eval(formals(mrf3_fit)$top_v_method)[1L],
    "saturation"
  )
  expect_equal(formals(mrf3_fit)$top_v_tau, 0.9)
  expect_equal(
    eval(formals(multiRF:::resolve_top_v_values)$top_v_method)[1L],
    "saturation"
  )
})
