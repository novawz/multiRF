#' Fit random forest model and create model lists
#'
#' @param X A data frame that consider as predictor set
#' @param Y A data frame that consider as response set. If Y = NULL, an unsupervised RF is conducted
#' @param type Select the type of RF model. The default is regression. Can select from "regression", "classification", and "unsupervised"
#' @param nodedepth Maximum depth for the legacy `randomForestSRC` fallback.
#' Ignored by the native engine.
#' @param ntree Number of trees.
#' @param forest.wt Forest-weight output mode for the legacy
#' `randomForestSRC` fallback. Ignored by the native engine, which always
#' returns the full `n x n` forest-weight matrix.
#' @param proximity Proximity output mode for the legacy `randomForestSRC`
#' fallback. Ignored by the native engine, which always returns the full
#' `n x n` proximity matrix.
#' @param mtry Number of candidate X variables per split. Default `NULL` =
#'   `floor(sqrt(px))`. Passed to native engine; ignored by rfsrc fallback.
#' @param ytry Number of candidate Y variables per split. Default `NULL` =
#'   `min(qy, ceiling(px/3))` for supervised, `15` for unsupervised.
#' @param seed Random seed passed to the selected engine.
#' @param engine Forest backend. Default is `getOption("multiRF.engine", "native")`.
#' Native is the default and recommended engine. `randomForestSRC` is used only
#' as a non-native fallback when explicitly requested.
#' @param ... Additional arguments passed to `randomForestSRC::rfsrc()`
#' when `engine != "native"`.
#' @return A model list
#' @details `fit_forest()` now defaults to the package-native engine for
#' classification, multivariate regression, and unsupervised fitting.
#' `randomForestSRC` is optional and is only used when `engine != "native"`.
fit_forest <-  function(X, Y = NULL, type = "regression", nodedepth = NULL,
                       ntree = 200, forest.wt = "all", proximity = "all",
                       mtry = NULL, ytry = NULL,
                       samptype = c("swor", "swr"),
                       seed = -10, engine = getOption("multiRF.engine", "native"), ...){

  X <- data.frame(X)
  samptype <- match.arg(samptype)

  if (identical(engine, "native") && identical(type, "classification")) {
    return(fit_class_forest(
      X = X,
      Y = Y,
      ntree = as.integer(ntree),
      mtry = mtry,
      nodesize = 1L,
      seed = as.integer(seed),
      proximity = proximity,
      samptype = samptype
    ))
  }

  # Use native C++ engine for multivariate regression (default)
  if (identical(engine, "native") && identical(type, "regression") && !is.null(Y)) {
    return(fit_mv_forest(
      X = X, Y = Y,
      ntree = as.integer(ntree),
      mtry = mtry,
      ytry = ytry,
      nodesize = 5L,
      seed = as.integer(seed),
      proximity = proximity,
      samptype = samptype
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

  if(type == "classification"){

    mrf <- randomForestSRC::rfsrc(
      Y ~ .,
      data = data.frame(Y = Y, X),
      membership = TRUE,
      nodedepth = nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      na.action = "na.impute",
      statistics = TRUE,
      proximity = proximity,
      samptype = samptype,
      seed = seed,
      ...
    )
  }

  if(type == "regression") {

    Y <- data.frame(Y)

    mrf <- randomForestSRC::rfsrc(
      randomForestSRC::get.mv.formula(colnames(Y)),
      data = data.frame(X,Y),
      membership = TRUE,
      nodedepth = nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      na.action = "na.impute",
      statistics = TRUE,
      proximity = proximity,
      samptype = samptype,
      ytry = ytry,
      seed = seed,
      ...
    )

  }
  if(type == "unsupervised"){

    if (identical(engine, "native")) {
      return(fit_mv_forest_unsup(
        X = X,
        ntree = as.integer(ntree),
        ytry = ytry,
        nodesize = 3L,
        seed = as.integer(seed),
        proximity = proximity,
        samptype = samptype
      ))
    }

    if(is.null(ytry)) ytry <- 15L  # legacy rfsrc fallback default

    mrf <- randomForestSRC::rfsrc(
      data = X,
      membership = TRUE,
      nodedepth = nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      statistics = TRUE,
      proximity = proximity,
      samptype = samptype,
      ytry = ytry,
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
#'   Default `NULL` means the native engine uses `min(qy, ceiling(px/3))`.
#'   Set to a specific integer to override (e.g., `ytry = ncol(Y) / 2`).
#' @rdname fit_forest
#'
fit_multi_forest <- function(dat.list, connect_list = NULL, var.wt = NULL, yprob = 1, ytry = NULL, seed = -10, ...){

  if(is.null(connect_list)){
    ytry_use <- if (!is.null(ytry)) as.integer(ytry) else NULL
    mod_l <- full_connect(dat.list, ytry = ytry_use, seed = seed, ...)

  } else {

    mod_l <- plyr::llply(
      connect_list,
      .fun = function(d){

        dat_fit <- dat.list[d]
        if(!is.null(var.wt)) {
          varwt <- var.wt[d]
        } else {
          varwt <- NULL
        }

        # Resolve ytry: explicit > default (0 = let C++ decide)
        ytry_use <- if (!is.null(ytry)) as.integer(ytry) else 0L

        if(length(dat_fit) == 1){

          mod <- fit_forest(dat_fit[[1]], xvar.wt = varwt[[1]], ytry = ytry_use, seed = seed, ...)
        } else {

          mod <- fit_forest(dat_fit[[2]], dat_fit[[1]], xvar.wt = varwt[[2]], yvar.wt = varwt[[1]], ytry = ytry_use, seed = seed, ...)
        }


        return(mod)
      }
    )

    names(mod_l) <- purrr::map(connect_list, ~paste0(., collapse = "_"))
  }

  return(mod_l)

}

# Backward compatibility aliases
fit_rfsrc <- function(...) {
  .Deprecated("fit_forest")
  fit_forest(...)
}

fit_multi_rfsrc <- function(...) {
  .Deprecated("fit_multi_forest")
  fit_multi_forest(...)
}
