#' Compute Cluster-Specific IMD Without Pairwise IMD
#'
#' @param x A `mrf3` or `mrf3_fit` object.
#' @param cluster Optional cluster labels for samples.
#' @param dat.list Optional named omics list (samples in rows).
#' @param connect_list Optional global connection list used by all clusters.
#' @param min_cluster_size Minimum sample size required to run a cluster-specific IMD.
#' @param ntree Deprecated compatibility argument. Ignored in no-refit mode.
#' @param ytry `ytry` used for per-cluster IMD. If `NULL`, reuses object setting.
#' @param parallel Logical; whether IMD computation inside each cluster is parallelized.
#' @param imd_normalized_weights Logical; passed to `get_multi_weights()` as `normalized`.
#' @param fit_args Deprecated compatibility argument. Ignored in no-refit mode.
#' @param imd_args Optional named list merged into each per-cluster `get_multi_weights()` call.
#' @param run_vs Logical; whether to run `mrf3_vs()` variable selection on each
#'   cluster after IMD.
#' @param vs_args A named list of additional arguments passed to `mrf3_vs()`
#'   for each cluster.
#' @param keep_model Logical; whether to keep per-cluster mrf3 objects in output.
#' @param keep_data Logical; whether to keep per-cluster subset data in output.
#' @param seed Base seed; cluster `i` uses `seed + i - 1`.
#'
#' @return A list with cluster-level summary, per-cluster IMD outputs, and params.
#' @export
cluster_imd <- function(x,
                        cluster = NULL,
                        dat.list = NULL,
                        connect_list = NULL,
                        min_cluster_size = 20L,
                        ntree = NULL,
                        ytry = NULL,
                        parallel = TRUE,
                        imd_normalized_weights = TRUE,
                        run_vs = FALSE,
                        vs_args = list(),
                        fit_args = list(),
                        imd_args = list(),
                        keep_model = FALSE,
                        keep_data = FALSE,
                        seed = 529) {

  if (!is.list(fit_args)) stop("`fit_args` must be a named list.")
  if (!is.null(ntree)) {
    warning("`ntree` is ignored in no-refit `cluster_imd()`.", call. = FALSE)
  }
  if (length(fit_args) > 0L) {
    warning("`fit_args` is ignored in no-refit `cluster_imd()`.", call. = FALSE)
  }
  if (!is.list(imd_args)) stop("`imd_args` must be a named list.")
  if (!is.list(vs_args)) stop("`vs_args` must be a named list.")

  base <- list(mod_list = NULL, connect_list = NULL, ntree = NULL, ytry = NULL, dat = dat.list)
  if (inherits(x, "mrf3_fit")) {
    base$mod_list <- x$models
    base$connect_list <- x$connection
    base$ntree <- x$config$ntree
    base$ytry <- x$config$ytry
    if (is.null(base$dat)) base$dat <- x$data
    if (is.null(cluster)) cluster <- x$clusters
  } else if (inherits(x, "mrf3")) {
    base$mod_list <- x$mod
    base$connect_list <- x$connection
    base$ntree <- x$ntree
    base$ytry <- x$ytry
    if (is.null(base$dat)) base$dat <- x$dat.list
    if (is.null(cluster)) cluster <- x$cl
  } else {
    stop("`x` must be a `mrf3` or `mrf3_fit` object.")
  }

  if (!is.list(base$mod_list) || length(base$mod_list) == 0L) {
    stop("Cannot find fitted model list in `x`.")
  }
  if (!is.list(base$dat) || length(base$dat) == 0L) {
    stop("Cannot find `dat.list`. Provide `dat.list` explicitly or run `mrf3_fit(..., return_data = TRUE)`.")
  }
  if (is.null(cluster)) {
    stop("Cannot infer cluster labels from `x`. Please provide `cluster` explicitly.")
  }

  if (is.null(names(base$dat)) || any(names(base$dat) == "")) {
    names(base$dat) <- paste0("omics", seq_along(base$dat))
  }
  dat_names <- names(base$dat)
  n <- nrow(base$dat[[1]])
  if (any(vapply(base$dat, nrow, integer(1)) != n)) {
    stop("All matrices in `dat.list` must share the same number of rows.")
  }
  if (!is.null(names(cluster)) && !is.null(rownames(base$dat[[1]]))) {
    if (all(rownames(base$dat[[1]]) %in% names(cluster))) {
      cluster <- cluster[rownames(base$dat[[1]])]
    }
  }
  if (length(cluster) != n) {
    stop("`cluster` length does not match sample size.")
  }
  cluster <- as.character(cluster)
  labs <- unique(cluster[!is.na(cluster)])
  if (length(labs) == 0L) {
    stop("`cluster` has no non-NA labels.")
  }

  make_conn_name <- function(conn) paste(as.character(conn), collapse = "_")

  if (is.null(names(base$mod_list)) || any(names(base$mod_list) == "")) {
    if (!is.null(base$connect_list) && length(base$connect_list) == length(base$mod_list)) {
      names(base$mod_list) <- vapply(base$connect_list, make_conn_name, character(1))
    } else {
      stop("`base$mod_list` has no names; cannot map models to connections.")
    }
  }

  conn_global <- normalize_connect_list(connect_list, n_blocks = length(dat_names), valid_names = dat_names)
  if (is.null(conn_global)) {
    conn_global <- normalize_connect_list(base$connect_list, n_blocks = length(dat_names), valid_names = dat_names)
  }
  if (is.null(conn_global) && length(base$dat) == 1L) {
    conn_global <- list(c(dat_names))
  }
  if (is.null(conn_global)) {
    conn_out <- tryCatch(find_connection(base$mod_list, return_score = TRUE), error = function(e) e)
    if (inherits(conn_out, "error")) {
      stop("Failed to infer global connect_list from existing models via find_connection: ", conditionMessage(conn_out))
    }
    conn_global <- conn_out$connect_list
  }

  conn_names <- vapply(conn_global, make_conn_name, character(1))
  mod_template <- base$mod_list[conn_names]
  use_refit <- any(vapply(mod_template, function(m) !is.null(m$sub_mrf_info), logical(1)))
  ytry_use <- if (is.null(ytry)) base$ytry else ytry

  subset_model_for_cluster <- function(mod, conn, dat_sub, idx) {
    ms <- mod
    if (!is.null(mod$membership)) {
      ms$membership <- mod$membership[idx, , drop = FALSE]
    }
    if (!is.null(mod$forest.wt)) {
      ms$forest.wt <- mod$forest.wt[idx, idx, drop = FALSE]
    }
    if (!is.null(mod$proximity)) {
      ms$proximity <- mod$proximity[idx, idx, drop = FALSE]
    }
    x_name <- if (length(conn) == 1L) conn[[1]] else conn[[2]]
    x_cols <- intersect(colnames(mod$xvar), colnames(dat_sub[[x_name]]))
    ms$xvar <- as.data.frame(dat_sub[[x_name]][, x_cols, drop = FALSE], check.names = FALSE)
    if (length(conn) == 2L) {
      y_name <- conn[[1]]
      y_cols <- intersect(colnames(mod$yvar), colnames(dat_sub[[y_name]]))
      ms$yvar <- as.data.frame(dat_sub[[y_name]][, y_cols, drop = FALSE], check.names = FALSE)
    }
    ms
  }

  summary_tb <- data.frame(
    cluster = labs,
    n = integer(length(labs)),
    status = rep("ok", length(labs)),
    message = rep(NA_character_, length(labs)),
    stringsAsFactors = FALSE
  )
  by_cluster <- vector("list", length(labs))
  names(by_cluster) <- labs

  for (i in seq_along(labs)) {
    lb <- labs[[i]]
    idx <- which(cluster == lb)
    summary_tb$n[i] <- length(idx)
    if (length(idx) < min_cluster_size) {
      summary_tb$status[i] <- "skipped"
      summary_tb$message[i] <- "cluster too small"
      next
    }

    dat_sub <- lapply(base$dat, function(d) d[idx, , drop = FALSE])

    # ── Fast path: cluster-weighted IMD from pre-computed per-node scores ──
    # Uses tree_info$imd_x_score + membership to weight each node's
    # contribution by the fraction of this cluster's samples passing through it.
    has_node_imd <- all(vapply(mod_template, function(m) {
      !is.null(m$tree_info) &&
        !is.null(m$tree_info[[1]]$imd_x_score)
    }, logical(1)))

    if (has_node_imd && !use_refit) {
      # Compute cluster-weighted IMD for this cluster across all connections
      cluster_labels <- setNames(rep("other", length(cluster)), names(cluster))
      cluster_labels[idx] <- lb

      per_conn_imd <- tryCatch({
        lapply(seq_along(conn_global), function(j) {
          cw <- cluster_weighted_imd(
            mod = mod_template[[j]],
            cluster = cluster_labels,
            normalized = FALSE  # normalize after averaging across connections
          )
          cw[[lb]]  # list(X = vec, Y = vec)
        })
      }, error = function(e) e)

      if (inherits(per_conn_imd, "error")) {
        summary_tb$status[i] <- "error"
        summary_tb$message[i] <- paste0("Cluster IMD failed: ", conditionMessage(per_conn_imd))
        next
      }

      # Aggregate across connections: average by block
      block_names <- names(base$dat)
      weight_list <- lapply(block_names, function(bn) {
        ww <- list()
        for (j in seq_along(conn_global)) {
          conn <- conn_global[[j]]
          m_name_sep <- rev(as.character(conn))
          if (m_name_sep[1] == bn) ww <- c(ww, list(per_conn_imd[[j]]$X))
          if (length(m_name_sep) > 1 && m_name_sep[2] == bn) ww <- c(ww, list(per_conn_imd[[j]]$Y))
        }
        if (length(ww) == 0L) return(setNames(numeric(ncol(base$dat[[bn]])), colnames(base$dat[[bn]])))
        w_out <- Reduce("+", ww) / length(ww)
        if (imd_normalized_weights) {
          denom <- sqrt(sum(w_out^2))
          if (is.finite(denom) && denom > 0) w_out <- w_out / denom
        }
        w_out
      })
      names(weight_list) <- block_names
      imd_out <- list(weight_list = weight_list, weight_list_init = NULL, net = NULL)

    } else {
      # ── Slow path: subset model + recompute via get_multi_weights ──
      mod_sub <- tryCatch({
        if (use_refit) {
          fit_multi_forest(
            dat.list = dat_sub,
            connect_list = conn_global,
            ntree = if (is.null(base$ntree)) 200L else base$ntree,
            ytry = ytry_use,
            seed = seed + i - 1L
          )
        } else {
          out <- lapply(seq_along(conn_global), function(j) {
            ms <- subset_model_for_cluster(mod_template[[j]], conn_global[[j]], dat_sub, idx)
            ms$imd_weights <- NULL           # don't use global IMD
            ms$imd_weights_per_tree <- NULL
            ms
          })
          names(out) <- conn_names
          out
        }
      }, error = function(e) e)

      if (inherits(mod_sub, "error")) {
        summary_tb$status[i] <- "error"
        summary_tb$message[i] <- conditionMessage(mod_sub)
        next
      }

      imd_call <- c(
        list(
          mod_list = mod_sub,
          dat.list = dat_sub,
          ytry = ytry_use,
          parallel = parallel,
          normalized = imd_normalized_weights,
          seed = seed + i - 1L
        ),
        imd_args
      )
      imd_out <- tryCatch(do.call(get_multi_weights, imd_call), error = function(e) e)
      if (inherits(imd_out, "error")) {
        summary_tb$status[i] <- "error"
        summary_tb$message[i] <- paste0("IMD failed: ", conditionMessage(imd_out))
        next
      }
    }

    mrf_sub <- structure(
      list(
        mod = mod_sub,
        connection = conn_global,
        weights = imd_out$weight_list,
        net = imd_out$net,
        dat.list = dat_sub,
        ntree = if (is.null(base$ntree)) NA_integer_ else base$ntree,
        ytry = ytry_use
      ),
      class = "mrf3"
    )

    vs_out <- NULL
    if (isTRUE(run_vs)) {
      final_vs <- c(list(mod = mrf_sub, dat.list = dat_sub), vs_args)
      vs_out <- tryCatch(do.call(mrf3_vs, final_vs), error = function(e) e)
    }

    by_cluster[[lb]] <- list(
      imd = imd_out,
      vs = vs_out,
      model = if (isTRUE(keep_model)) mrf_sub else NULL,
      data = if (isTRUE(keep_data)) dat_sub else NULL
    )
  }

  out <- list(
    summary = summary_tb,
    by_cluster = by_cluster,
    params = list(
      connect_list = conn_global,
      min_cluster_size = min_cluster_size,
      ytry = ytry_use,
      parallel = parallel,
      run_vs = run_vs
    )
  )
  class(out) <- c("cluster_imd", "list")
out
}

#' Print `cluster_imd`
#'
#' @method print cluster_imd
#' @param x Output from `cluster_imd()`.
#' @param ... Unused.
#'
#' @return Input object invisibly.
#' @export
print.cluster_imd <- function(x, ...) {
  cat("cluster_imd summary\n")
  if (is.data.frame(x$summary)) {
    print(x$summary, row.names = FALSE)
  }
  cat("\nAccess: x$by_cluster[[cluster_name]]$imd / $vs / $model / $data\n")
  invisible(x)
}
