# ============================================================================
# Accessor functions for multiRF pipeline objects
#
# Provides a clean, flat interface to commonly needed results without
# requiring users to know internal nested list structures.
# ============================================================================


# ---------- Cluster labels -------------------------------------------------

#' Get cluster labels
#'
#' Extract cluster assignments from any multiRF pipeline object.
#'
#' @param x An `mrf3_fit`, `mrf3`, `reconstr`, or `prox` object.
#' @param which For `mrf3_fit`: `"main"` (default) or `"robust"`.
#' @param ... Additional arguments (unused).
#'
#' @return A named vector of cluster assignments, or `NULL` if not available.
#' @export
get_clusters <- function(x, which = "main", ...) {
  if (inherits(x, "mrf3_fit")) {
    which <- match.arg(which, c("main", "robust"))
    if (identical(which, "robust")) {
      return(x$robust_clusters)
    }
    return(x$clusters)
  }
  if (inherits(x, "reconstr") || inherits(x, "prox")) {
    return(x$cl)
  }
  if (inherits(x, "mrf3") && !is.null(x$cl)) {
    return(x$cl)
  }
  NULL
}


# ---------- IMD weights ----------------------------------------------------

#' Get IMD weights
#'
#' Extract per-block feature importance weights.
#'
#' @param x An `mrf3_fit`, `mrf3`, or `cluster_imd` object.
#' @param cluster For `cluster_imd` objects, which cluster to extract.
#'   If `NULL`, returns all clusters as a named list of weight lists.
#' @param ... Additional arguments (unused).
#'
#' @return A named list of numeric vectors (one per block), or `NULL`.
#' @export
get_weights <- function(x, cluster = NULL, ...) {
  if (inherits(x, "cluster_imd")) {
    if (is.null(cluster)) {
      out <- lapply(x$by_cluster, function(cl) cl$imd$weight_list)
      return(out)
    }
    cl <- as.character(cluster)
    if (!cl %in% names(x$by_cluster)) {
      stop("Cluster `", cl, "` not found. Available: ",
           paste(names(x$by_cluster), collapse = ", "))
    }
    return(x$by_cluster[[cl]]$imd$weight_list)
  }
  if (inherits(x, "mrf3_fit")) {
    return(x$weights)
  }
  if (inherits(x, "mrf3")) {
    return(x$weights)
  }
  NULL
}


# ---------- Selected variables ---------------------------------------------

#' Get selected variable names
#'
#' Extract variable names selected by `mrf3_vs()`.
#'
#' @param x An `mrf3_fit`, `mrf3` (with class `vs`), or `cluster_imd` object.
#' @param cluster For `cluster_imd` objects, which cluster.  If `NULL`,
#'   returns all clusters as a named list.
#' @param ... Additional arguments (unused).
#'
#' @return A named list of character vectors (variable names per block),
#'   or `NULL` if variable selection has not been run.
#' @export
get_selected_vars <- function(x, cluster = NULL, ...) {
  if (inherits(x, "cluster_imd")) {
    extract_vs_vars <- function(cl) {
      if (is.null(cl$vs) || inherits(cl$vs, "error")) return(NULL)
      if (!is.null(cl$vs$dat.list)) {
        return(lapply(cl$vs$dat.list, colnames))
      }
      NULL
    }
    if (is.null(cluster)) {
      out <- lapply(x$by_cluster, extract_vs_vars)
      return(out)
    }
    cl <- as.character(cluster)
    if (!cl %in% names(x$by_cluster)) {
      stop("Cluster `", cl, "` not found. Available: ",
           paste(names(x$by_cluster), collapse = ", "))
    }
    return(extract_vs_vars(x$by_cluster[[cl]]))
  }
  if (inherits(x, "mrf3_fit")) {
    return(x$selected_vars)
  }
  if (inherits(x, "vs")) {
    if (!is.null(x$dat.list)) {
      return(lapply(x$dat.list, colnames))
    }
  }
  NULL
}


# ---------- Data -----------------------------------------------------------

#' Get data from a multiRF object
#'
#' Extract the data list attached to a pipeline object.
#'
#' @param x An `mrf3_fit`, `mrf3`, or `cluster_imd` object.
#' @param which For `mrf3_fit`: `"filtered"` (default) or `"selected"`
#'   (post-variable-selection).
#' @param cluster For `cluster_imd`, which cluster.
#' @param ... Additional arguments (unused).
#'
#' @return A named list of data frames, or `NULL`.
#' @export
get_data <- function(x, which = "filtered", cluster = NULL, ...) {
  if (inherits(x, "cluster_imd")) {
    if (is.null(cluster)) {
      return(lapply(x$by_cluster, function(cl) cl$data))
    }
    cl <- as.character(cluster)
    if (!cl %in% names(x$by_cluster)) {
      stop("Cluster `", cl, "` not found. Available: ",
           paste(names(x$by_cluster), collapse = ", "))
    }
    return(x$by_cluster[[cl]]$data)
  }
  if (inherits(x, "mrf3_fit")) {
    which <- match.arg(which, c("filtered", "selected"))
    if (identical(which, "selected")) {
      return(x$selected_data)
    }
    return(x$data)
  }
  if (inherits(x, "mrf3")) {
    return(x$dat.list)
  }
  NULL
}


# ---------- Connection -----------------------------------------------------

#' Get connection list
#'
#' Extract the connection (block pairing) structure.
#'
#' @param x An `mrf3_fit` or `mrf3` object.
#' @param ... Additional arguments (unused).
#'
#' @return A list of character vectors (each a block pair), or `NULL`.
#' @export
get_connection <- function(x, ...) {
  if (inherits(x, "mrf3_fit")) {
    return(x$connection)
  }
  if (inherits(x, "mrf3")) {
    return(x$connection)
  }
  NULL
}


# ---------- Models ---------------------------------------------------------

#' Get fitted RF models
#'
#' Extract the list of fitted random forest models.
#'
#' @param x An `mrf3_fit` or `mrf3` object.
#' @param which For `mrf3_fit`: `"init"` (default) or `"refit"`
#'   (post-variable-selection refit, if available).
#' @param ... Additional arguments (unused).
#'
#' @return A named list of rfsrc model objects, or `NULL`.
#' @export
get_models <- function(x, which = "init", ...) {
  if (inherits(x, "mrf3_fit")) {
    which <- match.arg(which, c("init", "refit"))
    if (identical(which, "refit")) {
      return(x$vs_detail$mod_list)
    }
    return(x$models)
  }
  if (inherits(x, "mrf3")) {
    return(x$mod)
  }
  NULL
}


# ---------- Shared / Specific ----------------------------------------------

#' Extract shared fraction table from a fitted pipeline
#'
#' Convenience accessor for the per-omics shared fraction stored
#' in an `mrf3_fit` object.  For computing shared fractions from raw
#' data, use `get_shared_frac()` instead.
#'
#' @param x An `mrf3_fit` object.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `data`, `shared_frac`, etc., or `NULL`.
extract_shared_frac <- function(x, ...) {
  if (inherits(x, "mrf3_fit")) {
    return(x$shared$frac)
  }
  NULL
}


# ---------- Top variables --------------------------------------------------

#' Get top weighted variables per block
#'
#' A convenience function that returns the top-n variables by IMD weight
#' for each omics block.
#'
#' @param x Any object accepted by `get_weights()`.
#' @param n Maximum number of top variables per block.
#' @param cluster Passed to `get_weights()`.
#' @param ... Additional arguments passed to `get_weights()`.
#'
#' @return A named list of data frames (one per block), each with columns
#'   `variable` and `weight`, sorted by weight descending.
#' @export
get_top_vars <- function(x, n = 20L, cluster = NULL, ...) {
  wl <- get_weights(x, cluster = cluster, ...)
  if (is.null(wl)) return(NULL)
  lapply(wl, function(w) {
    w <- sort(w[w != 0], decreasing = TRUE)
    w <- utils::head(w, n)
    data.frame(
      variable = names(w),
      weight = unname(w),
      stringsAsFactors = FALSE
    )
  })
}


# ---------- Variable selection summary -------------------------------------

#' Get variable selection summary
#'
#' Returns a tidy data frame summarizing how many features were selected
#' in each block.
#'
#' @param x An `mrf3_fit` or `cluster_imd` object.
#' @param ... Additional arguments (unused).
#'
#' @return A data frame with columns `block`, `p_total`, `p_selected`,
#'   `selected_ratio`.  For `cluster_imd`, an extra `cluster` column.
#' @export
get_vs_summary <- function(x, ...) {
  if (inherits(x, "mrf3_fit")) {
    return(x$vs_summary)
  }
  if (inherits(x, "cluster_imd")) {
    rows <- list()
    for (nm in names(x$by_cluster)) {
      cl <- x$by_cluster[[nm]]
      if (is.null(cl$vs) || inherits(cl$vs, "error")) next
      wl_orig <- cl$imd$weight_list
      vs_dat <- cl$vs$dat.list
      if (is.null(wl_orig) || is.null(vs_dat)) next
      for (blk in names(wl_orig)) {
        p_total <- length(wl_orig[[blk]])
        p_sel <- if (blk %in% names(vs_dat)) ncol(vs_dat[[blk]]) else 0L
        rows[[length(rows) + 1L]] <- data.frame(
          cluster = nm,
          block = blk,
          p_total = p_total,
          p_selected = p_sel,
          selected_ratio = if (p_total > 0) p_sel / p_total else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(rows) == 0L) return(NULL)
    return(do.call(rbind, rows))
  }
  NULL
}
