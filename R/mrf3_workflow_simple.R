#' Simplified mrf3 Entry Point
#'
#' A compact wrapper around `mrf3_fit()` for common use cases.
#' It keeps only the most frequently used parameters and forwards advanced
#' options through `...`.
#' By default, this wrapper runs shared/specific weighting and specific
#' clustering branches.
#'
#' @param dat.list A named list of omics matrices (samples in rows, features in columns).
#' @param k Optional cluster number used for both shared and specific clustering.
#' If `NULL`, `k` is tuned in downstream clustering routines.
#' @param ntree Number of trees for RF fitting.
#' @param top_v Optional unified top-v cutoff.
#' If set, applies to both `model_top_v` and `fused_top_v`.
#' @param samptype Sampling scheme passed through to `mrf3_fit()`.
#' @param main_clustering Global clustering strategy applied consistently to
#' shared, specific, and robust branches: `"similarity"` (default),
#' `"proximity"`, or `"enhanced_proximity"`.
#' @param filter_mode Filtering mode passed to `mrf3_fit()`.
#' @param seed Random seed.
#' @param ... Advanced options passed to `mrf3_fit()`, e.g.
#' `model_top_v_tune_args`, `fused_top_v_tune_args`,
#' `clustering_args`, `run_imd`, `imd_args`,
#' `run_variable_selection`, `variable_selection_args`,
#' `cluster_imd_args`, `run_robust_clustering`, `compact_output`.
#'
#' @return An object of class `"mrf3_fit"`.
#' @export
mrf3 <- function(dat.list,
                 k = NULL,
                 ntree = 500,
                 top_v = NULL,
                 samptype = c("swor", "swr"),
                 main_clustering = c("similarity", "proximity", "enhanced_proximity"),
                 filter_mode = c("none", "auto", "manual"),
                 seed = 529,
                 ...) {

  if (is.character(main_clustering) && length(main_clustering) == 1L &&
      identical(main_clustering, "off")) {
    stop(
      "`main_clustering = 'off'` has been removed. ",
      "Use one of: 'similarity', 'proximity', 'enhanced_proximity'."
    )
  }
  main_clustering <- match.arg(main_clustering)
  filter_mode <- match.arg(filter_mode)

  dots <- list(...)
  if (!is.null(dots$specific_method)) dots$specific_method <- NULL
  clustering_args <- list(
    shared_mode = "average",
    shared_method = "Spectral",
    shared_similarity_type = "second",
    specific_method = "Spectral",
    specific_similarity_type = "second",
    specific_prox_method_cl = "PAM"
  )

  if (!is.null(k)) {
    if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2) {
      stop("`k` must be NULL or a single integer >= 2.")
    }
    k <- as.integer(k)
    clustering_args$shared_k <- k
    clustering_args$specific_k <- k
  }

  if (!is.null(dots$clustering_args)) {
    clustering_args <- utils::modifyList(clustering_args, dots$clustering_args)
  }
  dots$clustering_args <- clustering_args

  if (!is.null(top_v)) {
    if (!is.numeric(top_v) || length(top_v) != 1L || !is.finite(top_v) || top_v <= 0) {
      stop("`top_v` must be NULL or a single positive numeric value.")
    }
    top_v <- as.integer(top_v)
  }

  base_args <- list(
    dat.list = dat.list,
    ntree = ntree,
    scale = TRUE,
    ytry = NULL,
    samptype = match.arg(samptype),
    connect_list = NULL,
    filter_mode = filter_mode,
    main_clustering = main_clustering,
    top_v = top_v,
    seed = seed
  )

  final_args <- utils::modifyList(base_args, dots)
  do.call(mrf3_fit, final_args)
}

#' Full mrf3 entry point
#'
#' @inheritParams mrf3_fit
#' @return An object of class `"mrf3_fit"`.
#' @export
mrf3_full <- function(...) {
  mrf3_fit(...)
}
