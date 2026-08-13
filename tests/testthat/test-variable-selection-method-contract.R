test_that("fixed filter uses tau times raw-IMD SD", {
  w <- list(A = c(a = 0.05, b = 0.10, c = 0.40, d = 0.80))
  dat <- list(A = data.frame(a = 1:5, b = 6:10, c = 11:15, d = 16:20))
  mod <- structure(
    list(imd = w, connection = list(c("A", "A")), ntree = 20L,
         type = "regression"),
    class = "mrf3"
  )

  out <- mrf3_vs(mod, dat, method = "thres", se = 1, re_fit = FALSE)

  expect_equal(unname(out$thres[["A"]]), stats::sd(w$A))
  expect_equal(out$selected_vars$A, c("c", "d"))
  expect_equal(out$imd$A[c("a", "b")], c(a = 0, b = 0))
  # Selected variables retain their raw Eq. 8 IMD values.
  expect_equal(out$imd$A[c("c", "d")], w$A[c("c", "d")])
  expect_false(isTRUE(all.equal(sqrt(sum(out$imd$A^2)), 1)))
})

test_that("adaptive filter uses adjacent candidate plateaus", {
  dat <- list(
    A = data.frame(a1 = 1:6, a2 = 6:1),
    B = data.frame(b1 = c(1, 1, 2, 2, 3, 3), b2 = c(3, 2, 1, 3, 2, 1))
  )
  calls <- list()
  candidate_scores <- rep(c(0.505, 0.600, 0.601), each = 2L)
  fake_fit <- function(...) {
    args <- list(...)
    calls[[length(calls) + 1L]] <<- args
    list(list(score = candidate_scores[[length(calls)]]))
  }
  fake_score <- function(models) models[[1L]]$score
  baseline <- list(list(
    membership = matrix(1L, nrow = 1L, ncol = 1L),
    inbag = matrix(1L, nrow = 1L, ncol = 1L),
    score = 0.5
  ))

  out <- multiRF:::choose_thres2(
    weights = list(A = c(a1 = 0.2, a2 = 0.8),
                   B = c(b1 = 0.3, b2 = 0.9)),
    connection = list(c("B", "A")),
    new_dat = dat,
    ytry = 1L,
    ntree = 10L,
    type = "regression",
    models_init = baseline,
    k = 2L,
    tol = 0.01,
    tau_grid = seq(0, 0.2, by = 0.1),
    seed = 101L,
    .fit_fun = fake_fit,
    .score_fun = fake_score
  )

  expect_length(calls, 6L)
  expect_true(all(vapply(calls, function(x) identical(x$forest.wt, "oob"),
                         logical(1))))
  expect_equal(vapply(calls, `[[`, integer(1), "seed"),
               rep(c(101L, 102L), 3L))
  # The one-fit baseline (.500) is close to the first candidate (.505), but it
  # is diagnostic only. The candidate plateau starts at .1: .600 vs .601.
  expect_equal(attr(out, "tau"), 0.1)
  expect_equal(attr(out, "baseline_oob_nmse"), 0.5)
  expect_equal(attr(out, "oob_trace")$mean_oob_nmse,
               c(0.505, 0.600, 0.601))
  expect_equal(attr(out, "oob_trace")$tau, seq(0, 0.2, by = 0.1))
})

test_that("get_multi_weights defaults to raw forest IMD", {
  expect_identical(formals(get_multi_weights)$normalized, FALSE)
  expect_identical(formals(mrf3_vs)$normalized, FALSE)
  expect_identical(formals(cluster_imd)$imd_normalized_weights, FALSE)
})

test_that("native IMD maps sides from connection metadata and tolerates NULL Y", {
  block <- "rna_block_with_underscore"
  dat_unsup <- setNames(
    list(data.frame(g1 = 1:4, g2 = 5:8)),
    block
  )
  x_pt <- rbind(g1 = c(0.4, 0.6), g2 = c(0.1, 0.2))
  unsup_mod <- list(
    imd_weights = list(X = c(g1 = 0.5, g2 = 0.15)),
    imd_weights_per_tree = list(X = x_pt, Y = NULL),
    connection = list(response = block, predictor = NULL)
  )
  out_unsup <- get_multi_weights(
    setNames(list(unsup_mod), "legacy_name_that_must_not_be_parsed"),
    dat_unsup
  )

  expect_equal(out_unsup$weight_list[[block]], unsup_mod$imd_weights$X)
  expect_length(out_unsup$weight_list_init[[1L]], ncol(x_pt))
  expect_named(out_unsup$weight_list_init[[1L]][[1L]], block)
  expect_equal(out_unsup$weight_list_init[[1L]][[1L]][[block]], x_pt[, 1L])

  response <- "protein_response_block"
  predictor <- "rna_predictor_block"
  dat_directed <- list(
    protein_response_block = data.frame(p1 = 1:4, p2 = 4:1),
    rna_predictor_block = data.frame(g1 = 5:8, g2 = 8:5)
  )
  directed_mod <- list(
    imd_weights = list(
      X = c(g1 = 0.6, g2 = 0.2),
      Y = c(p1 = 0.7, p2 = 0.1)
    ),
    imd_weights_per_tree = list(
      X = rbind(g1 = c(0.5, 0.7), g2 = c(0.1, 0.3)),
      Y = rbind(p1 = c(0.6, 0.8), p2 = c(0.05, 0.15))
    ),
    connection = list(response = response, predictor = predictor)
  )
  out_directed <- get_multi_weights(
    setNames(list(directed_mod), "also_not_a_parseable_connection"),
    dat_directed
  )

  expect_equal(out_directed$weight_list[[predictor]], directed_mod$imd_weights$X)
  expect_equal(out_directed$weight_list[[response]], directed_mod$imd_weights$Y)
  expect_named(
    out_directed$weight_list_init[[1L]][[1L]],
    c(predictor, response)
  )
})

test_that("raw IMD aggregation aligns names and rejects incomplete feature sets", {
  dat <- list(
    A = data.frame(a = 1:4, b = 5:8),
    B = data.frame(y = 1:4),
    C = data.frame(z = 4:1)
  )
  models <- list(
    B_A = list(
      connection = list(response = "B", predictor = "A"),
      imd_weights = list(X = c(a = 0.2, b = 0.8), Y = c(y = 0.1))
    ),
    C_A = list(
      connection = list(response = "C", predictor = "A"),
      # Deliberately reversed: aggregation must be keyed by feature name.
      imd_weights = list(X = c(b = 0.6, a = 0.4), Y = c(z = 0.3))
    )
  )

  out <- get_multi_weights(models, dat, normalized = FALSE)
  expect_equal(out$weight_list$A, c(a = 0.3, b = 0.7))

  missing <- models
  missing$C_A$imd_weights$X <- c(a = 0.4)
  expect_error(
    get_multi_weights(missing, dat, normalized = FALSE),
    "missing: b"
  )

  extra <- models
  extra$C_A$imd_weights$X <- c(b = 0.6, a = 0.4, rogue = 0.2)
  expect_error(
    get_multi_weights(extra, dat, normalized = FALSE),
    "unexpected: rogue"
  )
})

test_that("Eq. 16 uses per-tree SE and Student t with ntree minus one df", {
  mat <- rbind(
    a = c(0.80, 0.90, 0.85, 0.90, 0.88),
    b = c(0.10, 0.20, 0.10, 0.20, 0.10),
    c = c(0.20, 0.10, 0.20, 0.10, 0.20)
  )
  out <- multiRF:::.imd_transformation_one(mat, alpha = 0.05)
  expected_t <- (rowMeans(mat) - mean(mat)) /
    (apply(mat, 1, stats::sd) / sqrt(ncol(mat)))

  expect_equal(out$ts, expected_t)
  expect_equal(
    out$pval,
    stats::pt(expected_t, df = ncol(mat) - 1L, lower.tail = FALSE)
  )
  expect_equal(out$keep_idx, c(a = 1L, b = 0L, c = 0L))
})

test_that("transformation aggregates connections and warns for unevaluated blocks", {
  A <- rbind(
    a = c(0.80, 0.90, 0.85, 0.90, 0.88),
    b = c(0.10, 0.20, 0.10, 0.20, 0.10),
    c = c(0.20, 0.10, 0.20, 0.10, 0.20)
  )
  B <- rbind(
    x = c(0.10, 0.20, 0.10, 0.20, 0.10),
    y = c(0.70, 0.80, 0.70, 0.80, 0.75)
  )
  make_connection <- function(A, B) {
    lapply(seq_len(ncol(A)), function(tree) list(A = A[, tree], B = B[, tree]))
  }
  expect_warning(
    out <- multiRF:::test_fn(
      wl = list(make_connection(A, B), make_connection(A * 0.9, B * 0.9)),
      connection = list(c("B", "A"), c("B", "A")),
      dat_names = c("A", "B", "C"),
      sig.thres = c(A = 0.05, B = 0.10, C = 0.20),
      feature_names = list(A = rownames(A), B = rownames(B), C = c("u", "v"))
    ),
    "No connected forest supplies per-tree IMD for block `C`"
  )

  expect_named(out$keep_idx, c("A", "B", "C"))
  expect_equal(out$n_models, c(A = 2L, B = 2L, C = 0L))
  expect_equal(names(out$ts$A), rownames(A))
  expect_equal(out$keep_idx$A, c(a = 1L, b = 0L, c = 0L))
  expect_equal(out$keep_idx$C, c(u = 1L, v = 1L))
})

test_that("transformation requires a strict majority across connections", {
  high <- rbind(a = c(0.9, 0.8, 0.9, 0.8, 0.9),
                b = c(0.1, 0.2, 0.1, 0.2, 0.1))
  low <- rbind(a = c(0.1, 0.2, 0.1, 0.2, 0.1),
               b = c(0.9, 0.8, 0.9, 0.8, 0.9))
  as_trees <- function(mat) {
    lapply(seq_len(ncol(mat)), function(tree) list(A = mat[, tree]))
  }
  out <- multiRF:::test_fn(
    wl = list(as_trees(high), as_trees(low)),
    connection = list(c("A", "A"), c("A", "A")),
    dat_names = "A",
    sig.thres = 0.05,
    feature_names = list(A = c("a", "b"))
  )

  # Each feature wins exactly one of two forests, which is not a majority.
  expect_equal(out$keep_idx$A, c(a = 0L, b = 0L))
})

test_that("IMD mixture has an explicit zero mass and ordered noise component", {
  x <- c(z1 = 0, z2 = 0, n1 = 0.01, n2 = 0.02, s1 = 0.70, s2 = 0.80)
  post <- multiRF:::.fit_imd_mixture(x, c1 = "normal", c2 = "normal")

  expect_equal(colnames(post), c("noise", "signal", "zero"))
  expect_equal(post[c("z1", "z2"), "zero"], c(z1 = 1, z2 = 1))
  expect_equal(unname(rowSums(post)), rep(1, length(x)))
  expect_gt(post["n1", "noise"], post["n1", "signal"])
  expect_lt(post["s1", "noise"], post["s1", "signal"])
})

test_that("unidentifiable positive IMD mixture fails conservatively", {
  expect_warning(
    post <- multiRF:::.fit_imd_mixture(c(z = 0, a = 0.01, b = 0.01, c = 0.01)),
    "do not contain enough variation"
  )
  expect_equal(post["z", "zero"], 1)
  expect_equal(unname(post[c("a", "b", "c"), "noise"]), rep(1, 3))
  expect_equal(unname(post[c("a", "b", "c"), "signal"]), rep(0, 3))
  expect_true(isTRUE(attr(post, "degenerate")))
})

test_that("optional truncated-normal and gamma mixture families use stable MLE updates", {
  x <- c(z = 0, seq(0.005, 0.08, length.out = 20),
         seq(0.45, 0.9, length.out = 20))

  for (families in list(c("truncn", "normal"), c("normal", "gamma"),
                        c("truncn", "gamma"))) {
    post <- multiRF:::.fit_imd_mixture(
      x, c1 = families[[1L]], c2 = families[[2L]], iter = 200L
    )
    expect_true(all(is.finite(post)))
    expect_equal(unname(rowSums(post)), rep(1, length(x)), tolerance = 1e-8)
    expect_true(all(diff(attr(post, "component_means")) >= 0))
    trace <- attr(post, "loglik_trace")
    if (length(trace) > 1L) {
      expect_true(all(diff(trace) >= -1e-6))
    }
  }

  # Near-degenerate low gamma means must not be moved upward by a scale floor.
  edge <- c(0, seq(0.0001, 0.0002, length.out = 20), rep(0.005, 30))
  edge_post <- multiRF:::.fit_imd_mixture(
    edge, c1 = "normal", c2 = "gamma", iter = 100L
  )
  edge_trace <- attr(edge_post, "loglik_trace")
  expect_true(all(is.finite(edge_post)))
  if (length(edge_trace) > 1L) {
    expect_true(all(diff(edge_trace) >= -1e-6))
  }

  separated <- c(rep(0, 10), seq(0.01, 0.05, length.out = 30),
                 0.5 + seq(-1e-4, 1e-4, length.out = 30))
  separated_post <- multiRF:::.fit_imd_mixture(
    separated, c1 = "normal", c2 = "gamma", iter = 100L
  )
  separated_trace <- attr(separated_post, "loglik_trace")
  expect_true(all(diff(separated_trace) >= -1e-6))
  expect_gt(mean(tail(separated_post[, "signal"], 30)), 0.99)

  # Extreme positive boundary masses previously caused a gamma M-step to drop
  # the observed-data log likelihood by more than 80 units.
  boundary <- c(
    rep(0, 17), rep(.Machine$double.xmin, 30), rep(1e-300, 26),
    rep(1e-12, 25), rep(0.5, 15), rep(1 - 1e-12, 20), rep(1, 28)
  )
  boundary_post <- multiRF:::.fit_imd_mixture(
    boundary, c1 = "truncn", c2 = "gamma", iter = 300L, eps = 1e-8
  )
  boundary_trace <- attr(boundary_post, "loglik_trace")
  expect_true(all(is.finite(boundary_post)))
  expect_equal(unname(rowSums(boundary_post)), rep(1, length(boundary)),
               tolerance = 1e-8)
  expect_true(all(diff(boundary_trace) >= 0))
})

test_that("signal all returns a union object and specific refit is explicit", {
  dat <- list(
    A = data.frame(a = rnorm(10), b = rnorm(10), c = rnorm(10)),
    B = data.frame(x = rnorm(10), y = rnorm(10), z = rnorm(10))
  )
  A_pt <- rbind(a = rep(0.8, 3), b = rep(0.05, 3), c = rep(0.1, 3))
  B_pt <- rbind(x = rep(0.05, 3), y = rep(0.8, 3), z = rep(0.1, 3))
  wf <- structure(
    list(
      data = dat,
      imd = list(A = c(a = 0.1, b = 0.7, c = 0.1),
                 B = c(x = 0.7, y = 0.1, z = 0.1)),
      imd_init = NULL,
      connection = list(c("B", "A")),
      config = list(ntree = 3L, ytry = 1L),
      type = "regression",
      oob_err = 1,
      models = list(),
      specific = list(
        imd = list(A = rowMeans(A_pt), B = rowMeans(B_pt)),
        imd_per_tree = list(A = A_pt, B = B_pt)
      )
    ),
    class = c("mrf3_fit", "list")
  )

  specific <- mrf3_vs(wf, method = "thres", signal = "specific",
                      se = 1, re_fit = FALSE)
  expect_identical(specific$mod, wf$models)
  expect_false(isTRUE(specific$refit_performed))

  all <- mrf3_vs(wf, method = "thres", signal = "all", se = 1, re_fit = FALSE)
  expect_s3_class(all, "mrf3_vs_all")
  expect_equal(all$selected_vars$A, c("a", "b"))
  expect_equal(all$selected_vars$B, c("x", "y"))
  expect_named(all$signal_results, c("shared", "specific"))
  expect_false(isTRUE(all$refit_performed))
})

test_that("signal all keeps adaptive filtering for the shared component", {
  dat <- list(
    A = data.frame(a = 1:8, b = 8:1),
    B = data.frame(x = rep(1:4, each = 2), y = rep(4:1, each = 2))
  )
  wf <- structure(
    list(
      data = dat,
      imd = list(A = c(a = 0.2, b = 0.8), B = c(x = 0.3, y = 0.9)),
      imd_init = NULL,
      connection = list(c("B", "A")),
      config = list(ntree = 3L, ytry = 1L), type = "regression",
      oob_err = NULL, models = list(),
      specific = list(
        imd = list(A = c(a = 0.3, b = 0.7), B = c(x = 0.4, y = 0.6)),
        imd_per_tree = NULL
      )
    ),
    class = c("mrf3_fit", "list")
  )

  calls <- new.env(parent = emptyenv())
  calls$adaptive <- 0L
  calls$fixed <- 0L
  testthat::local_mocked_bindings(
    choose_thres2 = function(...) {
      calls$adaptive <- calls$adaptive + 1L
      c(A = -1, B = -1)
    },
    chooss_thres3 = function(weights, se) {
      calls$fixed <- calls$fixed + 1L
      stats::setNames(rep(-1, length(weights)), names(weights))
    },
    .package = "multiRF"
  )

  expect_warning(
    out <- mrf3_vs(wf, method = "filter", signal = "all", re_fit = FALSE),
    "shared component will use adaptive OOB filtering"
  )
  expect_equal(c(adaptive = calls$adaptive, fixed = calls$fixed),
               c(adaptive = 1L, fixed = 1L))
  expect_identical(out$selection_method, "filter")
})

test_that("specific selection preserves consensus mean instead of last-run mean", {
  dat <- list(A = data.frame(a = 1:8, b = 8:1))
  wf <- structure(
    list(
      data = dat,
      imd = list(A = c(a = 0.5, b = 0.5)), imd_init = NULL,
      connection = list("A"), config = list(ntree = 3L, ytry = 1L),
      type = "unsupervised", models = list(),
      specific = list(
        imd = list(A = c(a = 0.8, b = 0.2)),
        # Represents a final consensus run disagreeing with the mean.
        imd_per_tree = list(A = rbind(a = rep(0.1, 3), b = rep(0.9, 3)))
      )
    ),
    class = c("mrf3_fit", "list")
  )
  out <- mrf3_vs(
    wf, method = "thres", signal = "specific", se = 1, re_fit = FALSE
  )
  expect_equal(out$imd$A, c(a = 0.8, b = 0))
  expect_equal(out$selected_vars$A, "a")
})

test_that("specific selection refits the original cross-omics connection", {
  dat <- list(
    A = data.frame(a = 1:8, b = 8:1),
    B = data.frame(x = rep(1:4, each = 2), y = rep(4:1, each = 2))
  )
  original_connection <- list(c("B", "A"))
  wf <- structure(
    list(
      data = dat, imd = list(A = c(a = 0.5, b = 0.5),
                             B = c(x = 0.5, y = 0.5)),
      imd_init = NULL, connection = original_connection,
      config = list(ntree = 3L, ytry = 1L), type = "regression",
      models = list(),
      specific = list(
        imd = list(A = c(a = 0.8, b = 0.2), B = c(x = 0.7, y = 0.3)),
        imd_per_tree = NULL
      )
    ),
    class = c("mrf3_fit", "list")
  )
  captured <- new.env(parent = emptyenv())
  captured$connection <- NULL
  testthat::local_mocked_bindings(
    fit_multi_forest = function(dat.list, connect_list, ...) {
      captured$connection <- connect_list
      list(list(refitted = TRUE))
    },
    .vs_mean_model_oob_nmse = function(models) 0.1,
    .package = "multiRF"
  )

  out <- mrf3_vs(
    wf, method = "thres", signal = "specific", se = 0,
    re_fit = TRUE, re_weights = TRUE
  )
  expect_true(out$refit_performed)
  expect_equal(captured$connection, original_connection)
  expect_true(out$mod[[1]]$refitted)
})

test_that("workflow honors explicit final refit outside robust mode", {
  workflow_body <- paste(deparse(body(mrf3_fit)), collapse = "\n")
  expect_match(
    workflow_body,
    "if \\(force_refit\\)[[:space:]]+final_vs_args\\$re_fit <- TRUE"
  )
  expect_false(grepl(
    "final_vs_args\\$re_fit <- force_refit",
    workflow_body,
    fixed = FALSE
  ))
})
