#' Fit random forest model and create model lists
#'
#' @param X A data frame that consider as predictor set
#' @param Y A data frame that consider as response set. If Y = NULL, an unsupervised RF is conducted
#' @param type Select the type of RF model. The default is regression. Can select from "regression", "classification", and "unsupervised"
#' @param nodedepth Backward-compatible alias for `max_depth`.
#' @param max_depth Maximum tree depth; `0` means unlimited.
#' @param nodesize Minimum terminal-node size. Defaults to 5 for regression,
#' 1 for classification, and 3 for unsupervised forests.
#' @param nthread Number of threads used by the native engine.
#' @param ntree Number of trees.
#' @param forest.wt Forest-weight output mode: `"all"` gives all-sample
#' terminal-node weights, `"inbag"` restricts donors to bootstrap
#' samples, and `"oob"` evaluates each target only on trees where it is OOB.
#' @param proximity Proximity output mode: `"all"`, `"inbag"`, `"oob"`, or
#' `"none"`.
#' @param mtry Number of candidate X variables per split. Default `NULL` =
#'   `ceiling(px/3)` for multivariate regression.
#' @param ytry Number of candidate Y variables per split. Default `NULL` =
#'   `ceiling(qy/3)` for supervised, `15` for unsupervised.
#' @param nsplit Number of candidate numeric cutpoints evaluated per variable.
#'   Native and `randomForestSRC` both default to `10`; set `0` to scan all.
#' @param samptype Sampling scheme: `"swor"` (without replacement) or `"swr"`
#'   (with replacement).
#' @param xvar.wt,yvar.wt Optional non-negative predictor/response sampling
#' weights used by the native multivariate forest.
#' @param split.wt,case.wt Optional split/case weights. These are forwarded to
#' `randomForestSRC`; the native engine fails explicitly because it does not
#' yet implement them.
#' @param seed Random seed passed to the selected engine.
#' @param engine Forest backend. Default is `getOption("multiRF.engine", "native")`.
#' Native is the default and recommended engine. `randomForestSRC` is used only
#' as a non-native fallback when explicitly requested.
#' @param enhanced_prox Logical; whether to compute enhanced proximity in the
#'   native engine.
#' @param sibling_gamma Strength of the sibling-leaf correction used by
#'   enhanced proximity.
#' @param leaf_embed_dim Embedding dimension used by the native enhanced
#'   proximity path.
#' @param ... Additional arguments passed to `randomForestSRC::rfsrc()`
#' when `engine != "native"`.
#' @return A model list
#' @details `fit_forest()` now defaults to the package-native engine for
#' classification, multivariate regression, and unsupervised fitting.
#' `randomForestSRC` is optional and is only used when `engine != "native"`.
#' If `type` is omitted and `Y = NULL`, unsupervised fitting is selected.
#' @export
fit_forest <- function(X, Y = NULL,
                       type = c("regression", "classification", "unsupervised"),
                       nodedepth = NULL, max_depth = NULL, nodesize = NULL,
                       ntree = 200, forest.wt = "all", proximity = "all",
                       mtry = NULL, ytry = NULL, nsplit = 10,
                       samptype = c("swor", "swr"),
                       nthread = getOption("multiRF.nthread", 0L),
                       xvar.wt = NULL, yvar.wt = NULL,
                       split.wt = NULL, case.wt = NULL,
                       seed = 529L, engine = getOption("multiRF.engine", "native"),
                       enhanced_prox = FALSE, sibling_gamma = 0.5,
                       leaf_embed_dim = 10L, ...){

  X <- as.data.frame(X, check.names = FALSE)
  if (missing(type) && is.null(Y)) type <- "unsupervised"
  type <- match.arg(type)
  samptype <- match.arg(samptype)
  forest.wt <- match.arg(forest.wt, c("all", "inbag", "oob"))
  proximity <- match.arg(proximity, c("all", "inbag", "oob", "none"))

  ntree <- .native_integer_scalar(ntree, "ntree", 1L)
  nsplit <- .native_integer_scalar(nsplit, "nsplit", 0L)
  nthread <- .native_integer_scalar(nthread, "nthread", 0L)
  seed <- .native_seed_scalar(seed)
  enhanced_prox <- .native_logical_scalar(enhanced_prox, "enhanced_prox")
  leaf_embed_dim <- .native_integer_scalar(leaf_embed_dim, "leaf_embed_dim", 1L)
  sibling_gamma <- .native_nonnegative_scalar(sibling_gamma, "sibling_gamma")

  if (!is.null(max_depth) && !is.null(nodedepth)) {
    max_depth_checked <- .native_integer_scalar(max_depth, "max_depth", 0L)
    nodedepth_checked <- .native_integer_scalar(nodedepth, "nodedepth", 0L)
    if (max_depth_checked != nodedepth_checked) {
      stop("Specify only one of `max_depth` and `nodedepth`, or give them the same value.")
    }
  }
  if (is.null(max_depth)) max_depth <- if (is.null(nodedepth)) 0L else nodedepth
  max_depth <- .native_integer_scalar(max_depth, "max_depth/nodedepth", 0L)
  if (is.null(nodesize)) {
    nodesize <- switch(type, regression = 5L, classification = 1L, unsupervised = 3L)
  }
  nodesize <- .native_integer_scalar(nodesize, "nodesize", 1L)

  if (identical(engine, "native")) {
    if (!is.null(split.wt) || !is.null(case.wt)) {
      stop("The native engine does not yet implement `split.wt` or `case.wt`.")
    }
    dot_args <- list(...)
    if (length(dot_args) > 0L) {
      dot_names <- names(dot_args)
      if (is.null(dot_names)) dot_names <- rep("<unnamed>", length(dot_args))
      dot_names[dot_names == ""] <- "<unnamed>"
      stop(
        "Unsupported native `fit_forest()` argument(s): ",
        paste(dot_names, collapse = ", "), "."
      )
    }
  }

  if (identical(engine, "native") && identical(type, "classification")) {
    return(fit_class_forest(
      X = X,
      Y = Y,
      ntree = as.integer(ntree),
      mtry = mtry,
      forest.wt = forest.wt,
      nsplit = as.integer(nsplit),
      nodesize = nodesize,
      max_depth = max_depth,
      seed = as.integer(seed),
      proximity = proximity,
      samptype = samptype,
      nthread = as.integer(nthread),
      xvar.wt = xvar.wt
    ))
  }

  # Use native C++ engine for multivariate regression (default)
  if (identical(engine, "native") && identical(type, "regression") && !is.null(Y)) {
    return(fit_mv_forest(
      X = X, Y = Y,
      ntree = as.integer(ntree),
      mtry = mtry,
      ytry = ytry,
      nsplit = as.integer(nsplit),
      forest.wt = forest.wt,
      nodesize = nodesize,
      max_depth = max_depth,
      seed = as.integer(seed),
      proximity = proximity,
      samptype = samptype,
      nthread = as.integer(nthread),
      xvar.wt = xvar.wt,
      yvar.wt = yvar.wt,
      enhanced_prox = enhanced_prox,
      sibling_gamma = sibling_gamma,
      leaf_embed_dim = leaf_embed_dim
    ))
  }

  # Keep the all-native unsupervised path independent of the optional
  # randomForestSRC package.
  if (identical(engine, "native") && identical(type, "unsupervised")) {
    if (!is.null(xvar.wt) || !is.null(yvar.wt)) {
      stop("Native unsupervised forests do not currently support variable weights.")
    }
    return(fit_mv_forest_unsup(
      X = X,
      ntree = as.integer(ntree),
      ytry = ytry,
      nsplit = as.integer(nsplit),
      forest.wt = forest.wt,
      nodesize = nodesize,
      max_depth = max_depth,
      seed = as.integer(seed),
      proximity = proximity,
      samptype = samptype,
      nthread = as.integer(nthread),
      enhanced_prox = enhanced_prox,
      sibling_gamma = sibling_gamma,
      leaf_embed_dim = leaf_embed_dim
    ))
  }

  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    stop(
      "`randomForestSRC` is only needed for non-native fallback paths. ",
      "Install it or use `engine = 'native'`.",
      call. = FALSE
    )
  }

  # rfsrc doesn't accept "none"; map to FALSE for the fallback path
  if (identical(proximity, "none")) proximity <- FALSE
  fallback_nodedepth <- if (identical(max_depth, 0L)) NULL else max_depth

  # rfsrc only accepts integer mtry/ytry; resolve formula strings here
  n.xvar <- ncol(data.frame(X))
  if (is.character(mtry)) mtry <- resolve_param(mtry, p = n.xvar, default = NULL, name = "mtry")
  if (is.character(ytry)) ytry <- resolve_param(ytry, p = ncol(data.frame(Y)), default = NULL, name = "ytry")

  if(type == "classification"){

    mrf <- randomForestSRC::rfsrc(
      Y ~ .,
      data = data.frame(Y = Y, X, check.names = FALSE),
      membership = TRUE,
      nodedepth = fallback_nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      mtry = mtry,
      na.action = "na.impute",
      statistics = TRUE,
      proximity = proximity,
      samptype = samptype,
      nodesize = nodesize,
      nsplit = nsplit,
      xvar.wt = xvar.wt,
      split.wt = split.wt,
      case.wt = case.wt,
      seed = seed,
      ...
    )
  }

  if(type == "regression") {

    Y <- as.data.frame(Y, check.names = FALSE)
    # Keep the alternative backend on the same response-subsampling default
    # as the native implementation (qtry = ceil(q/3)).
    if (is.null(ytry)) ytry <- ceiling(ncol(Y) / 3)

    mrf <- randomForestSRC::rfsrc(
      randomForestSRC::get.mv.formula(colnames(Y)),
      data = data.frame(X, Y, check.names = FALSE),
      membership = TRUE,
      nodedepth = fallback_nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      mtry = mtry,
      na.action = "na.impute",
      statistics = TRUE,
      proximity = proximity,
      samptype = samptype,
      ytry = ytry,
      nodesize = nodesize,
      nsplit = nsplit,
      xvar.wt = xvar.wt,
      yvar.wt = yvar.wt,
      split.wt = split.wt,
      case.wt = case.wt,
      seed = seed,
      ...
    )

  }
  if(type == "unsupervised"){
    if(is.null(ytry)) ytry <- 15L  # legacy rfsrc fallback default

    mrf <- randomForestSRC::rfsrc(
      data = X,
      membership = TRUE,
      nodedepth = fallback_nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      mtry = mtry,
      statistics = TRUE,
      proximity = proximity,
      samptype = samptype,
      ytry = ytry,
      nodesize = nodesize,
      nsplit = nsplit,
      xvar.wt = xvar.wt,
      split.wt = split.wt,
      case.wt = case.wt,
      seed = seed,
      ...
    )
  }

  return(mrf)
}


#' @param dat.list A list that contains multi-omics datasets with samples in rows and features in columns. Samples should be matched in each dataset.
#' @param connect_list A pre-defined connection list between datasets. If `NULL`,
#' all directed pairwise connections are enumerated.
#' @param var.wt Optional variable-weight list aligned with `dat.list`.
#' @param yprob Deprecated. Use `ytry` directly instead.
#' @param ytry Number of response variables randomly selected per split.
#'   Default `NULL` means the native engine uses `ceiling(qy/3)`.
#'   Set to a specific integer to override (e.g., `ytry = ncol(Y) / 2`).
#' @rdname fit_forest
#' @export
fit_multi_forest <- function(dat.list, connect_list = NULL, var.wt = NULL,
                             yprob = 1, ytry = NULL, seed = 529L,
                             parallel_connections = FALSE,
                             cores_connections = NULL, ...){

  connection_used <- if (is.null(connect_list)) {
    enumerate_connections(names(dat.list))
  } else {
    connect_list
  }

  fit_one <- function(d){
    dat_fit <- dat.list[d]
    if(!is.null(var.wt)) {
      varwt <- var.wt[d]
    } else {
      varwt <- NULL
    }

    if(length(dat_fit) == 1){
      dots_fit <- list(...)
      type_fit <- if (is.null(dots_fit$type)) "unsupervised" else dots_fit$type
      dots_fit$type <- NULL
      fit_args <- c(
        list(X = dat_fit[[1]], Y = NULL, type = type_fit,
             xvar.wt = varwt[[1]], ytry = ytry, seed = seed),
        dots_fit
      )
      mod <- do.call(fit_forest, fit_args)
    } else {
      mod <- fit_forest(dat_fit[[2]], dat_fit[[1]], xvar.wt = varwt[[2]], yvar.wt = varwt[[1]], ytry = ytry, seed = seed, ...)
    }
    mod
  }

  n_conn <- length(connection_used)
  if (parallel_connections && n_conn > 1L) {
    n_par <- sanitize_mc_cores(
      cores = if (!is.null(cores_connections)) cores_connections else parallel::detectCores(),
      fallback = 1L
    )
    n_par <- min(n_conn, n_par)
    mod_l <- parallel::mclapply(connection_used, fit_one, mc.cores = n_par)
  } else {
    mod_l <- lapply(connection_used, fit_one)
  }
  names(mod_l) <- make.unique(vapply(
    connection_used, paste, collapse = "_", FUN.VALUE = character(1)
  ))

  # Store unambiguous connection metadata so downstream code need not parse
  # block names containing underscores.
  if (length(mod_l) == length(connection_used)) {
    for (i in seq_along(mod_l)) {
      d <- as.character(connection_used[[i]])
      mod_l[[i]]$connection <- list(
        response = d[[1L]],
        predictor = if (length(d) >= 2L) d[[2L]] else NULL
      )
      mod_l[[i]]$response_block <- d[[1L]]
      mod_l[[i]]$predictor_block <- if (length(d) >= 2L) d[[2L]] else NULL
    }
  }

  return(mod_l)

}

# Backward compatibility aliases
fit_rfsrc <- function(X, Y = NULL, type = "regression", nodedepth = NULL,
                      max_depth = NULL, nodesize = NULL,
                      ntree = 200, forest.wt = "all", proximity = "all",
                      mtry = NULL, ytry = NULL, nsplit = 10,
                      samptype = c("swor", "swr"),
                      nthread = getOption("multiRF.nthread", 0L),
                      xvar.wt = NULL, yvar.wt = NULL,
                      split.wt = NULL, case.wt = NULL,
                      seed = 529L, engine = getOption("multiRF.engine", "native"),
                      enhanced_prox = FALSE, sibling_gamma = 0.5,
                      leaf_embed_dim = 10L, ...) {
  .Deprecated("fit_forest")
  fit_forest(
    X = X,
    Y = Y,
    type = type,
    nodedepth = nodedepth,
    max_depth = max_depth,
    nodesize = nodesize,
    ntree = ntree,
    forest.wt = forest.wt,
    proximity = proximity,
    mtry = mtry,
    ytry = ytry,
    nsplit = nsplit,
    samptype = samptype,
    nthread = nthread,
    xvar.wt = xvar.wt,
    yvar.wt = yvar.wt,
    split.wt = split.wt,
    case.wt = case.wt,
    seed = seed,
    engine = engine,
    enhanced_prox = enhanced_prox,
    sibling_gamma = sibling_gamma,
    leaf_embed_dim = leaf_embed_dim,
    ...
  )
}

fit_multi_rfsrc <- function(dat.list, connect_list = NULL, var.wt = NULL,
                            yprob = 1, ytry = NULL, seed = 529L, ...) {
  .Deprecated("fit_multi_forest")
  fit_multi_forest(
    dat.list = dat.list,
    connect_list = connect_list,
    var.wt = var.wt,
    yprob = yprob,
    ytry = ytry,
    seed = seed,
    ...
  )
}
