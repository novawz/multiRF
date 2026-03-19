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
#' @param ytry Optional `ytry` value passed to `randomForestSRC::rfsrc()`
#' when `engine != "native"`.
#' @param seed Random seed passed to the selected engine.
#' @param engine Forest backend. Default is `getOption("multiRF.engine", "native")`.
#' Native is the default and recommended engine. `randomForestSRC` is used only
#' as a non-native fallback when explicitly requested.
#' @param ... Additional arguments passed to `randomForestSRC::rfsrc()`
#' when `engine != "native"`.
#' @return A model list
#' @details `fit_rfsrc()` now defaults to the package-native engine for
#' classification, multivariate regression, and unsupervised fitting.
#' `randomForestSRC` is optional and is only used when `engine != "native"`.
fit_rfsrc <-  function(X, Y = NULL, type = "regression", nodedepth = NULL,
                       ntree = 200, forest.wt = "all", proximity = "all", ytry = NULL,
                       seed = -10, engine = getOption("multiRF.engine", "native"), ...){

  X <- data.frame(X)

  if (identical(engine, "native") && identical(type, "classification")) {
    return(fit_class_forest(
      X = X,
      Y = Y,
      ntree = as.integer(ntree),
      nodesize = 5L,
      seed = as.integer(seed)
    ))
  }

  # Use native C++ engine for multivariate regression (default)
  if (identical(engine, "native") && identical(type, "regression") && !is.null(Y)) {
    ytry_val <- if (is.null(ytry)) 0L else as.integer(ytry)
    return(fit_mv_forest(
      X = X, Y = Y,
      ntree = as.integer(ntree),
      ytry = ytry_val,
      nodesize = 5L,
      seed = as.integer(seed)
    ))
  }

  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    stop(
      "`randomForestSRC` is only needed for non-native fallback paths. ",
      "Install it or use `engine = 'native'`.",
      call. = FALSE
    )
  }

  if(type == "classification"){

    mrf <- randomForestSRC::rfsrc(
      Y ~ .,
      data = data.frame(Y = Y, X),
      membership = T,
      nodedepth = nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      na.action = "na.impute",
      statistics = T,
      proximity = proximity,
      seed = seed,
      ...
    )
  }

  if(type == "regression") {

    Y <- data.frame(Y)

    mrf <- randomForestSRC::rfsrc(
      randomForestSRC::get.mv.formula(colnames(Y)),
      data = data.frame(X,Y),
      membership = T,
      nodedepth = nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      na.action = "na.impute",
      statistics = T,
      proximity = proximity,
      ytry = ytry,
      seed = seed,
      ...
    )

  }
  if(type == "unsupervised"){

    if(is.null(ytry)) ytry <- 10

    if (identical(engine, "native")) {
      return(fit_mv_forest_unsup(
        X = X,
        ntree = as.integer(ntree),
        ytry = as.integer(ytry),
        nodesize = 5L,
        seed = as.integer(seed)
      ))
    }

    mrf <- randomForestSRC::rfsrc(
      data = X,
      membership = T,
      nodedepth = nodedepth,
      var.used = "by.tree",
      forest.wt = forest.wt,
      ntree = ntree,
      statistics = T,
      proximity = proximity,
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
#' @param yprob Proportion of response variables used for tuning `ytry`.
#' @rdname fit_rfsrc
#'
fit_multi_rfsrc <- function(dat.list, connect_list = NULL, var.wt = NULL, yprob = 1, seed = -10, ...){

  if(is.null(connect_list)){

    mod_l <- full_connect(dat.list, seed = seed, ...)

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

        if(length(dat_fit) == 1){

          ytry <- min(ceiling(ncol(dat_fit[[1]]) * yprob), 20)
          mod <- fit_rfsrc(dat_fit[[1]], xvar.wt = varwt[[1]], ytry = ytry, seed = seed, ...)
        } else {

          if(is.null(varwt)) {
            ytry <- min(ceiling(ncol(dat_fit[[1]]) * yprob), 5000)
          } else {
            ytry <- min(ceiling(ncol(dat_fit[[1]]) * yprob), length(varwt[[1]][varwt[[1]] != 0]))
          }

          if(ytry == ncol(dat_fit[[1]])) ytry <- NULL

          mod <- fit_rfsrc(dat_fit[[2]], dat_fit[[1]], xvar.wt = varwt[[2]], yvar.wt = varwt[[1]], ytry = ytry, seed = seed, ...)
        }


        return(mod)
      }
    )

    names(mod_l) <- purrr::map(connect_list, ~paste0(., collapse = "_"))
  }

  return(mod_l)

}
