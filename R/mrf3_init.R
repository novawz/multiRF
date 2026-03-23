#' Fit an initial MRF model
#'
#' @param dat.list A list containing multi-omics datasets with samples in rows and features in columns.
#' @param ntree Number of trees for fitting MRF model. Default is 300.
#' @param scale Whether to z-standardize each feature. Default is TRUE.
#' @param ytry Number of response variables sampled at each split. `NULL`
#'   delegates to the engine default.
#' @param samptype Sampling scheme passed to forest fitting: `"swor"` or `"swr"`.
#' @param connect_list Optional pre-defined connection list. If `NULL`, the
#' function fits all directed pairwise models first, then uses `find_connection()`
#' to select one direction per omics pair.
#' @param filter_mode Feature filtering mode passed to `filter_omics()`:
#' `"auto"`, `"none"`, or `"manual"`.
#' @param filter_method Feature dispersion metric passed to `filter_omics()`:
#' `"mad"` or `"variance"`.
#' @param top_n_by_type Optional auto-filter overrides passed to `filter_omics()`.
#' @param top_n_manual Optional manual top-n configuration passed to `filter_omics()`.
#' @param filter_verbose Logical; whether to print filtering diagnostics.
#' @param return_data Whether to return the data list. Default is FALSE.
#' @param sub_mrf  Logical; if `TRUE`, uses `fit_sub_multi_rfsrc()` instead of
#'   `fit_multi_forest()` for the final model fit.  This sub-samples response and
#'   predictor features per connection, fits smaller MRFs, and averages the
#'   resulting n x n matrices.  Useful when p is large.  Default is `FALSE`.
#' @param sub_mrf_args  Named list of arguments forwarded to
#'   `fit_sub_multi_rfsrc()` when `sub_mrf = TRUE`.  Common options:
#'   `n_sub` (number of replicates, default 15), `frac_response` (default 0.2),
#'   `frac_predictor` (default 0.2), `ntree_per_sub` (default 25),
#'   `min_response_for_sub` (default 500; blocks smaller than this use all features),
#'   `min_predictor_for_sub` (default 500), `ytry` (default 0.5).
#' @param verbose Logical; whether to print progress messages.
#' @param seed Random seed.
#' @param ... Additional arguments passed to forest fitting helpers.
#'
#' @details `mrf3_init()` now performs initialization only (filtering, forest fitting,
#' and connection selection). IMD weights are computed downstream via
#' `get_multi_weights()` (for example inside `mrf3_fit(..., run_imd = TRUE)`).
#'
#' When `sub_mrf = TRUE`, the final model fit uses a sub-sampling ensemble
#' strategy: for each connection, features are randomly sub-sampled from both
#' response and predictor blocks, smaller MRFs are fitted, and the n x n
#' forest-weight and proximity matrices are averaged.  This trades a small
#' amount of matrix approximation accuracy for substantial speed gains when
#' p is large (e.g., 6-8x faster at p = 2000).  Connection selection (when
#' `connect_list = NULL`) still uses full MRF with OOB to ensure proper
#' direction determination.
#' @return mrf3 object
#' @export
#'
#'
# ---------------------------------------------------------------------------------
# Wrapper function for MRF initialization
# ---------------------------------------------------------------------------------
mrf3_init <- function(dat.list,
                 # Number of trees
                 ntree = 300,
                 scale = TRUE,
                 ytry = NULL,
                 samptype = c("swor", "swr"),
                 # Find connections
                 connect_list = NULL,
                 filter_mode = c("auto", "none", "manual"),
                 filter_method = c("mad", "variance"),
                 top_n_by_type = NULL,
                 top_n_manual = NULL,
                 filter_verbose = TRUE,
                 return_data = FALSE,
                 # Sub-MRF ensemble
                 sub_mrf = FALSE,
                 sub_mrf_args = list(),
                 verbose = TRUE,
                 seed = 529,
                 ...){

  dots <- list(...)
  removed_imd_args <- c(
    "calc_weights", "robust", "normalized", "use_depth", "calc",
    "parallel", "cores", "permute", "nperm", "alpha"
  )
  removed_hits <- intersect(names(dots), removed_imd_args)
  if (length(removed_hits) > 0) {
    stop(
      "`mrf3_init()` no longer accepts IMD arguments: ",
      paste(removed_hits, collapse = ", "),
      ". Run IMD in downstream stage via `mrf3_fit(..., run_imd = TRUE, imd_args = list(...))` ",
      "or call `get_multi_weights()` explicitly."
    )
  }

  if (length(dat.list) == 1) type = "unsupervised"
  if (length(dat.list) > 1) type = "regression"

  if (is.null(names(dat.list)) || any(names(dat.list) == "")) {
    names(dat.list) <- paste0("omics", seq_along(dat.list))
  }

  samptype <- match.arg(samptype)
  filter_mode <- match.arg(filter_mode)
  filter_method <- match.arg(filter_method)

  if (verbose) message("Filtering data..")
  prep_info <- filter_omics(
    dat.list = dat.list,
    filter_mode = filter_mode,
    filter_method = filter_method,
    top_n_by_type = top_n_by_type,
    top_n_manual = top_n_manual,
    return_summary = TRUE,
    verbose = filter_verbose
  )

  dat_fit <- prep_info$dat_filtered

  if (length(dat_fit) == 1) {
    connect_list <- list(c(names(dat_fit)))
  } else if (!is.null(connect_list)) {
    connect_list <- normalize_connect_list(
      connect_list = connect_list,
      n_blocks = length(dat.list),
      valid_names = names(dat_fit)
    )
  }

  if (scale) {
    new_dat <- purrr::map(dat_fit, ~as.data.frame(base::scale(.), check.names = FALSE))
  } else {
    new_dat <- dat_fit
  }

  if (verbose) message("Fitting models..")
  connection_score <- NULL
  connection_top_v_used <- NULL

  if (sub_mrf && length(new_dat) > 1) {
    ## ==== Sub-MRF path ====
    ## Fit all connections ONCE; reuse selected ones directly (no refit).
    sub_defaults <- list(
      n_sub = 15L,
      frac_response = 0.2,
      frac_predictor = 0.2,
      ntree_per_sub = 25L,
      min_response_for_sub = 500L,
      min_predictor_for_sub = 500L,
      ntree_full = ntree,
      enhanced = FALSE,
      compute_imd = FALSE,
      imd_args = list(),
      ytry = ytry,
      seed = seed,
      verbose = TRUE
    )
    sub_params <- utils::modifyList(sub_defaults, sub_mrf_args)

    if (is.null(connect_list)) {
      ## Generate all directed pairwise connections
      block_names <- names(new_dat)
      all_connections <- list()
      for (r in block_names) {
        for (p in block_names) {
          if (r != p) all_connections <- c(all_connections, list(c(r, p)))
        }
      }

      ## Fit ALL pairwise connections with sub-MRF
      if (verbose) message("  Fitting all pairwise (sub-MRF)..")
      sub_args <- sub_params
      sub_args$dat.list <- new_dat
      sub_args$connect_list <- all_connections
      mod_all <- do.call(fit_sub_multi_rfsrc, sub_args)

      ## Select best directions using OOB forest.wt for quality scores
      ## find_connection reads $forest.wt, so temporarily swap in the OOB version
      if (verbose) message("  Finding connections (OOB)..")
      mod_all_oob <- lapply(mod_all, function(m) {
        m$forest.wt <- m$forest.wt.oob
        m
      })
      conn_out <- find_connection(mod_all_oob, return_score = TRUE, drop_bottom_q = 0.2)
      connect_list <- conn_out$connect_list
      connection_score <- conn_out$score
      connection_top_v_used <- conn_out$top_v_used

      ## Reuse already-fitted models for the selected connections
      selected_names <- vapply(connect_list, paste0, collapse = "_",
                               FUN.VALUE = character(1))
      mod_list <- mod_all[selected_names]

    } else {
      ## connect_list provided — fit selected connections only
      if (verbose) message("  Fitting selected connections (sub-MRF)..")
      sub_args <- sub_params
      sub_args$dat.list <- new_dat
      sub_args$connect_list <- connect_list
      mod_list <- do.call(fit_sub_multi_rfsrc, sub_args)
    }

    oob_err <- NULL

  } else {
    ## ==== Standard full-MRF path ====
    if (length(new_dat) > 1 && is.null(connect_list)) {
      ## Fit ALL pairwise connections once with forest.wt = "all"
      mod_all <- fit_multi_forest(
        new_dat,
        connect_list = NULL,
        ntree = ntree,
        type = type,
        ytry = ytry,
        samptype = samptype,
        seed = seed,
        forest.wt = "all",
        ...
      )

      if (verbose) message("  Finding connections..")
      conn_out <- find_connection(mod_all, return_score = TRUE, drop_bottom_q = 0.2)
      connect_list <- conn_out$connect_list
      connection_score <- conn_out$score
      connection_top_v_used <- conn_out$top_v_used

      ## Reuse already-fitted models — no refit
      selected_names <- vapply(connect_list, paste0, collapse = "_",
                               FUN.VALUE = character(1))
      mod_list <- mod_all[selected_names]
    } else {
      mod_list <- fit_multi_forest(
        new_dat,
        connect_list = connect_list,
        ntree = ntree,
        type = type,
        ytry = ytry,
        samptype = samptype,
        seed = seed,
        forest.wt = "all",
        ...
      )
    }

    oob_err <- purrr::map(mod_list, ~get_r_sq(.))
    oob_err <- Reduce("+", oob_err)
  }

  if (!return_data) dat_fit <- NULL
  out <- list(
    mod = mod_list,
    oob_err = oob_err,
    type = type,
    connection = connect_list,
    connection_score = connection_score,
    connection_top_v_used = connection_top_v_used,
    ntree = ntree,
    ytry = ytry,
    sub_mrf = sub_mrf,
    dat.list = dat_fit,
    filter_summary = prep_info$filter_summary
  )

  attr(out, "class") <- "mrf3"

  if (verbose) message("Done!")

  return(out)
}
