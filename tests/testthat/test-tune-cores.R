test_that("sanitize_mc_cores always returns a valid positive core count", {
  sanitize_mc_cores <- get("sanitize_mc_cores", envir = asNamespace("multiRF"))

  expect_gte(sanitize_mc_cores(NULL), 1L)
  expect_gte(sanitize_mc_cores(NA_integer_), 1L)
  expect_gte(sanitize_mc_cores(0L), 1L)
  expect_gte(sanitize_mc_cores(-5L), 1L)
  expect_identical(sanitize_mc_cores(1L), 1L)
})

test_that("PSOCK lapply keeps sequential semantics with one core", {
  parallel_lapply_psock <- get(
    "parallel_lapply_psock", envir = asNamespace("multiRF")
  )
  parent_pid <- Sys.getpid()

  result <- parallel_lapply_psock(
    4:1,
    function(i) list(value = i^2, pid = Sys.getpid()),
    cores = 1L
  )

  expect_identical(vapply(result, `[[`, numeric(1), "value"), c(16, 9, 4, 1))
  expect_true(all(vapply(result, `[[`, integer(1), "pid") == parent_pid))
})

test_that("PSOCK lapply uses fresh processes without forking", {
  skip_on_cran()
  available_cores <- parallel::detectCores()
  skip_if(is.na(available_cores) || available_cores < 2L)
  parallel_lapply_psock <- get(
    "parallel_lapply_psock", envir = asNamespace("multiRF")
  )
  parent_pid <- Sys.getpid()

  result <- parallel_lapply_psock(
    c(4L, 1L),
    function(i) {
      list(
        value = i^2,
        pid = Sys.getpid(),
        fork_child = parallel:::isChild()
      )
    },
    cores = 8L
  )

  expect_identical(vapply(result, `[[`, numeric(1), "value"), c(16, 1))
  expect_true(all(vapply(result, `[[`, integer(1), "pid") != parent_pid))
  expect_identical(length(unique(vapply(result, `[[`, integer(1), "pid"))), 2L)
  expect_false(any(vapply(result, `[[`, logical(1), "fork_child")))
})

test_that("PSOCK workers inherit library paths and multiRF options", {
  skip_on_cran()
  available_cores <- parallel::detectCores()
  skip_if(is.na(available_cores) || available_cores < 2L)
  parallel_lapply_psock <- get(
    "parallel_lapply_psock", envir = asNamespace("multiRF")
  )
  old_options <- options(multiRF.engine = "native", multiRF.nthread = 3L)
  on.exit(options(old_options), add = TRUE)

  result <- parallel_lapply_psock(
    1:2,
    function(i) {
      list(
        package_path = system.file(package = "multiRF"),
        engine = getOption("multiRF.engine"),
        nthread = getOption("multiRF.nthread")
      )
    },
    cores = 2L
  )

  expect_true(all(nzchar(vapply(result, `[[`, character(1), "package_path"))))
  expect_identical(vapply(result, `[[`, character(1), "engine"), rep("native", 2L))
  expect_identical(vapply(result, `[[`, integer(1), "nthread"), rep(3L, 2L))
})

test_that("parallel grid evaluation preserves candidate order", {
  skip_on_cran()
  available_cores <- parallel::detectCores()
  skip_if(is.na(available_cores) || available_cores < 2L)
  eval_grid <- get("eval_grid", envir = asNamespace("multiRF"))

  result <- eval_grid(
    list(7L, 2L, 9L),
    function(i) c(candidate = i, score = i / 10),
    parallel = TRUE,
    cores = 2L
  )

  expect_identical(vapply(result, `[[`, numeric(1), "candidate"), c(7, 2, 9))
  expect_equal(vapply(result, `[[`, numeric(1), "score"), c(0.7, 0.2, 0.9))

  expect_error(
    eval_grid(
      list(4L, 6L),
      function(i) {
        if (i == 6L) stop("deliberate failure")
        i
      },
      parallel = TRUE,
      cores = 2L
    ),
    "candidate value\\(s\\) 6: deliberate failure"
  )
})

test_that("forest post-processing agrees between sequential and PSOCK runs", {
  skip_on_cran()
  available_cores <- parallel::detectCores()
  skip_if(is.na(available_cores) || available_cores < 2L)

  set.seed(529)
  sample_ids <- paste0("s", seq_len(16L))
  X <- matrix(
    stats::rnorm(16L * 6L),
    nrow = 16L,
    dimnames = list(sample_ids, paste0("x", seq_len(6L)))
  )
  Y <- matrix(
    stats::rnorm(16L * 5L),
    nrow = 16L,
    dimnames = list(sample_ids, paste0("y", seq_len(5L)))
  )
  mod <- multiRF::fit_forest(
    X,
    Y,
    ntree = 6L,
    nthread = 1L,
    forest.wt = "all",
    proximity = "all",
    seed = 529L
  )
  get_imp_forest <- get("get_imp_forest", envir = asNamespace("multiRF"))
  cl_forest <- get("cl_forest", envir = asNamespace("multiRF"))

  imp_sequential <- get_imp_forest(mod, parallel = FALSE, calc = "Both")
  expect_warning(
    imp_psock <- get_imp_forest(
      mod,
      parallel = TRUE,
      cores = 2L,
      calc = "Both"
    ),
    NA
  )
  expect_equal(imp_psock$imp_ls, imp_sequential$imp_ls, tolerance = 0)

  # Exercise the PSOCK path first. A preceding sequential call would load the
  # Matrix S4 classes in the main process and could mask cold-session failures.
  prox_psock <- cl_forest(mod, parallel = TRUE, cores = 2L, symm = TRUE)
  prox_sequential <- cl_forest(mod, parallel = FALSE, symm = TRUE)
  expect_equal(prox_psock$prox, prox_sequential$prox, tolerance = 1e-14)
})

test_that("ordinary proximity can be reconstructed from stored membership", {
  set.seed(529)
  sample_ids <- paste0("s", seq_len(18L))
  X <- matrix(
    stats::rnorm(18L * 6L),
    nrow = 18L,
    dimnames = list(sample_ids, paste0("x", seq_len(6L)))
  )
  Y <- matrix(
    stats::rnorm(18L * 5L),
    nrow = 18L,
    dimnames = list(sample_ids, paste0("y", seq_len(5L)))
  )

  fit_none <- multiRF::fit_forest(
    X, Y,
    ntree = 8L,
    nthread = 1L,
    proximity = "none",
    seed = 529L
  )
  fit_all <- multiRF::fit_forest(
    X, Y,
    ntree = 8L,
    nthread = 1L,
    proximity = "all",
    seed = 529L
  )
  expect_identical(fit_none$membership, fit_all$membership)

  reconstruct <- get(
    "reconstruct_plain_proximity", envir = asNamespace("multiRF")
  )
  prox_sequential <- reconstruct(fit_none, parallel = FALSE)
  expect_equal(as.matrix(prox_sequential), fit_all$proximity, tolerance = 0)

  available_cores <- parallel::detectCores()
  if (!is.na(available_cores) && available_cores >= 2L && !identical(Sys.getenv("NOT_CRAN"), "false")) {
    prox_psock <- reconstruct(fit_none, parallel = TRUE, cores = 2L)
    expect_equal(as.matrix(prox_psock), fit_all$proximity, tolerance = 0)
  }

  clustering <- multiRF::mrf3_cl_prox(
    list(fit_none),
    k = 2L,
    enhanced = FALSE,
    parallel = FALSE
  )
  expect_equal(clustering$dat, fit_all$proximity, tolerance = 0)

  mixed <- multiRF::mrf3_cl_prox(
    list(fit_none, fit_all),
    k = 2L,
    enhanced = FALSE,
    parallel = FALSE
  )
  expect_equal(mixed$dat, fit_all$proximity, tolerance = 0)

  combine <- get(
    "combine_proximity_matrices", envir = asNamespace("multiRF")
  )
  reversed_ids <- rev(rownames(fit_all$proximity))
  reordered <- fit_all$proximity[reversed_ids, reversed_ids, drop = FALSE]
  expect_equal(
    combine(list(fit_all$proximity, reordered)) / 2,
    fit_all$proximity,
    tolerance = 0
  )
  colnames(reordered)[1L] <- "different-sample"
  expect_error(
    combine(list(fit_all$proximity, reordered)),
    "model 2.*sample names"
  )

  legacy <- fit_all
  legacy$proximity.mode <- NULL
  legacy_result <- suppressMessages(
    multiRF::mrf3_cl_prox(
      list(legacy), k = 2L, enhanced = FALSE, parallel = FALSE
    )
  )
  expect_equal(legacy_result$dat, fit_all$proximity, tolerance = 0)

  invalid_stored <- fit_all
  invalid_stored$proximity[,] <- 0
  expect_error(
    multiRF::mrf3_cl_prox(
      list(fit_all, invalid_stored),
      k = 2L,
      enhanced = FALSE,
      parallel = FALSE
    ),
    "model 2.*diagonal"
  )

  invalid <- fit_none
  invalid$membership[1L, 1L] <- NA_real_
  expect_error(
    multiRF::mrf3_cl_prox(
      list(invalid),
      k = 2L,
      enhanced = FALSE,
      parallel = FALSE
    ),
    "model 1.*membership.*missing"
  )

  invalid$membership <- fit_none$membership
  invalid$membership[1L, 1L] <- 1.5
  expect_error(
    multiRF::mrf3_cl_prox(
      list(invalid),
      k = 2L,
      enhanced = FALSE,
      parallel = FALSE
    ),
    "model 1.*integer terminal-node IDs"
  )
})

test_that("randomForestSRC fallback records and reconstructs omitted proximity", {
  skip_if_not_installed("randomForestSRC")

  set.seed(41)
  ids <- paste0("s", seq_len(12L))
  X <- matrix(
    rnorm(12L * 4L), nrow = 12L,
    dimnames = list(ids, paste0("x", seq_len(4L)))
  )
  Y <- matrix(
    rnorm(12L * 3L), nrow = 12L,
    dimnames = list(ids, paste0("y", seq_len(3L)))
  )
  fit <- multiRF::fit_forest(
    X, Y,
    ntree = 5L,
    proximity = "none",
    engine = "randomForestSRC",
    seed = 41L
  )

  expect_identical(fit$proximity.mode, "none")
  expect_true(is.null(fit$proximity) || identical(fit$proximity, FALSE))
  clustered <- suppressMessages(
    multiRF::mrf3_cl_prox(
      list(fit), k = 2L, enhanced = FALSE, parallel = FALSE
    )
  )
  expect_equal(dim(clustered$dat), c(12L, 12L))
  expect_true(all(is.finite(clustered$dat)))
})
