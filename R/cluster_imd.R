#' Compute Cluster-Specific IMD Without Pairwise IMD
#'
#' @param x A `mrf3` or `mrf3_fit` object.
#' @param cluster Optional cluster labels for samples.
#' @param dat.list Optional named omics list (samples in rows).
#' @param connect_list Optional global connection list used by all clusters.
#' @param min_cluster_size Minimum sample size required to run a cluster-specific IMD.
#' @param ntree Deprecated compatibility argument. Ignored in no-refit mode.
#' @param ytry `ytry` used for per-cluster IMD. If `NULL`, reuses object setting.
#'   Only used by the depth-based path.
#' @param parallel Logical; whether IMD computation inside each cluster is
#'   parallelized. Only used by the depth-based path (the node-score fast path
#'   is single-pass and needs no parallelism).
#' @param imd_normalized_weights Logical; passed to `get_multi_weights()` as
#'   `normalized`. The default is `FALSE` so cluster-level variable selection
#'   receives raw forest IMD on its 0-to-1 scale.
#' @param fit_args Deprecated compatibility argument. Ignored in no-refit mode.
#' @param imd_args Optional named list merged into each per-cluster `get_multi_weights()` call.
#'   Supplying any `imd_args` forces the depth-based path, because the
#'   node-score fast path cannot honor `get_multi_weights()` options.
#' @param run_vs Logical; whether to run `mrf3_vs()` variable selection on each
#'   cluster after IMD. `run_vs = TRUE` forces the depth-based path, which
#'   builds the per-cluster subset models that `mrf3_vs()` requires.
#' @param vs_args A named list of additional arguments passed to `mrf3_vs()`
#'   for each cluster.
#' @param use_node_imd Controls the estimator used for the per-cluster weights.
#'   `FALSE` (the default) computes Eq. 6-8 inverse minimal depth via the
#'   per-tree network traversal of `get_multi_weights()` on each cluster's
#'   subset model. `TRUE` uses the node-score fast path
#'   (`cluster_weighted_imd()`): the per-node split statistics stored by the
#'   native engine are averaged per variable, weighted by the fraction of the
#'   cluster's samples whose root-to-leaf paths visit each node. The fast
#'   path is deterministic and typically 100-900x faster, but it is a
#'   *different, split-score-weighted* estimator whose weights stay close to
#'   the global importance ranking: in benchmarks its per-cluster weights
#'   correlate at Spearman 0.9-0.99 *across* clusters (depth path: 0.6-0.8),
#'   and agree with the depth-based weights at Spearman ~0.4-0.6 (top-20
#'   overlap ~50-75%). Use it as a quick descriptive screen, not as a
#'   substitute for the depth-based cluster IMD (see
#'   `cluster_weighted_imd()` for the definition). `NULL` selects the fast
#'   path automatically whenever the models carry the required `tree_info`
#'   statistics and no blocking option (`imd_args`, `run_vs`) is requested,
#'   falling back to the depth-based path otherwise.
#' @param keep_model Logical; whether to keep per-cluster mrf3 objects in output.
#' @param keep_data Logical; whether to keep per-cluster subset data in output.
#' @param seed Base seed; cluster `i` uses `seed + i - 1`.
#'
#' @return A list with cluster-level summary, per-cluster IMD outputs, and
#'   params. `params$imd_path` records which estimator produced the weights:
#'   `"node_score"` (fast path) or `"depth_traversal"` (Eq. 6-8 IMD).
#' @export
cluster_imd <- function(x,
                        cluster = NULL,
                        dat.list = NULL,
                        connect_list = NULL,
                        min_cluster_size = 20L,
                        ntree = NULL,
                        ytry = NULL,
                        parallel = TRUE,
                        imd_normalized_weights = FALSE,
                        run_vs = FALSE,
                        vs_args = list(),
                        fit_args = list(),
                        imd_args = list(),
                        use_node_imd = FALSE,
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
    if (is.null(cluster)) {
      cl <- x$clusters
      cluster <- if (is.list(cl) && !is.null(names(cl))) cl$shared else cl
    }
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
  ref_samples <- rownames(base$dat[[1L]])
  if (is.null(ref_samples) || anyNA(ref_samples) || any(!nzchar(ref_samples)) || anyDuplicated(ref_samples)) {
    stop("Every block must have unique, non-missing sample names.")
  }
  for (j in seq_along(base$dat)[-1L]) {
    block_samples <- rownames(base$dat[[j]])
    if (!identical(block_samples, ref_samples)) {
      stop(
        "Sample identities and row order differ between blocks `", dat_names[[1L]],
        "` and `", dat_names[[j]], "`."
      )
    }
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

  if (is.null(names(base$mod_list)) || any(names(base$mod_list) == "")) {
    if (!is.null(base$connect_list) && length(base$connect_list) == length(base$mod_list)) {
      names(base$mod_list) <- connection_display_names(base$connect_list)
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

  conn_idx <- model_indices_for_connections(base$mod_list, conn_global)
  mod_template <- base$mod_list[conn_idx]
  use_refit <- any(vapply(mod_template, function(m) !is.null(m$sub_mrf_info), logical(1)))
  ytry_use <- if (is.null(ytry)) base$ytry else ytry

  # ── Node-score fast path availability ──────────────────────────────────
  # The native engine stores per-node split statistics (`imd_x_score`,
  # `imd_y_stats`) in `tree_info`; `cluster_weighted_imd()` turns them into
  # cluster-specific weights in a single deterministic pass, avoiding the
  # per-cluster per-tree network traversal entirely. Note this is a
  # split-score-weighted estimator, not Eq. 6-8 inverse minimal depth; see
  # `use_node_imd` in the roxygen above.
  has_node_imd <- !use_refit && all(vapply(mod_template, function(m) {
    !is.null(m$membership) && is.list(m$tree_info) &&
      length(m$tree_info) > 0L && !is.null(m$tree_info[[1L]]$imd_x_score)
  }, logical(1)))

  fast_blockers <- character(0)
  if (length(imd_args) > 0L) fast_blockers <- c(fast_blockers, "`imd_args`")
  if (isTRUE(run_vs)) fast_blockers <- c(fast_blockers, "`run_vs = TRUE`")

  if (is.null(use_node_imd)) {
    use_fast <- has_node_imd && length(fast_blockers) == 0L
    if (has_node_imd && length(fast_blockers) > 0L) {
      message(
        "cluster_imd: node-score fast path cannot honor ",
        paste(fast_blockers, collapse = " and "),
        "; using the depth-based traversal instead."
      )
    }
  } else if (isTRUE(use_node_imd)) {
    if (!has_node_imd) {
      stop(
        "`use_node_imd = TRUE` requires native-engine models carrying ",
        "`tree_info` with `imd_x_score` (and `membership`). Refit with the ",
        "native engine or use `use_node_imd = FALSE`."
      )
    }
    if (length(fast_blockers) > 0L) {
      stop(
        "`use_node_imd = TRUE` cannot honor ",
        paste(fast_blockers, collapse = " and "),
        "; drop those options or use `use_node_imd = FALSE`."
      )
    }
    use_fast <- TRUE
  } else {
    use_fast <- FALSE
  }

  # The per-label accumulators inside cluster_weighted_imd() are independent,
  # so one pass over the full label vector serves every cluster at once;
  # hoist it out of the per-cluster loop.
  cw_all <- NULL
  if (use_fast) {
    cluster_named <- stats::setNames(cluster, ref_samples)
    cw_all <- tryCatch(
      lapply(mod_template, function(m) {
        cluster_weighted_imd(m, cluster_named, normalized = FALSE)
      }),
      error = function(e) e
    )
    if (inherits(cw_all, "error")) {
      if (isTRUE(use_node_imd)) {
        stop("Node-score fast path failed: ", conditionMessage(cw_all))
      }
      message(
        "cluster_imd: node-score fast path failed (",
        conditionMessage(cw_all),
        "); falling back to the depth-based traversal."
      )
      use_fast <- FALSE
      cw_all <- NULL
    }
  }

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
    mod_sub <- NULL

    if (use_fast) {
      # ── Fast path: cluster-weighted IMD from pre-computed node scores ──
      # `cw_all[[j]][[lb]]` holds this cluster's X/Y weights for connection j
      # (computed once for all clusters before the loop). Aggregate across
      # connections per block with the same strict name alignment the
      # depth-based path uses.
      imd_out <- tryCatch({
        per_conn_imd <- lapply(cw_all, function(cw) cw[[lb]])
        if (any(vapply(per_conn_imd, is.null, logical(1)))) {
          stop("cluster label `", lb, "` missing from node-score weights.")
        }
        block_names <- names(base$dat)
        conn_labels <- names(mod_template)
        weight_list <- lapply(block_names, function(bn) {
          ww <- list()
          w_models <- character(0)
          for (j in seq_along(conn_global)) {
            conn <- as.character(conn_global[[j]])
            x_block <- if (length(conn) == 1L) conn[[1L]] else conn[[2L]]
            if (x_block == bn) {
              ww <- c(ww, list(per_conn_imd[[j]]$X))
              w_models <- c(w_models, conn_labels[[j]])
            }
            if (length(conn) == 2L && conn[[1L]] == bn) {
              ww <- c(ww, list(per_conn_imd[[j]]$Y))
              w_models <- c(w_models, conn_labels[[j]])
            }
          }
          .aggregate_imd_block(
            ww,
            feature_names = colnames(base$dat[[bn]]),
            block = bn,
            model_names = w_models,
            normalized = imd_normalized_weights
          )
        })
        names(weight_list) <- block_names
        list(weight_list = weight_list, weight_list_init = NULL, net = NULL)
      }, error = function(e) e)

      if (inherits(imd_out, "error")) {
        summary_tb$status[i] <- "error"
        summary_tb$message[i] <- paste0("Cluster IMD failed: ", conditionMessage(imd_out))
        next
      }

      if (isTRUE(keep_model)) {
        mod_sub <- tryCatch({
          out <- lapply(seq_along(conn_global), function(j) {
            ms <- subset_model_for_cluster(mod_template[[j]], conn_global[[j]], dat_sub, idx)
            ms$imd_weights <- NULL           # global IMD does not describe the cluster
            ms$imd_weights_per_tree <- NULL
            ms
          })
          names(out) <- names(mod_template)
          out
        }, error = function(e) NULL)
      }

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
          names(out) <- names(mod_template)
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
        imd = imd_out$weight_list,
        imd_ls = imd_out$weight_list_init,
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
      if (inherits(vs_out, "error")) {
        summary_tb$status[i] <- "vs_error"
        summary_tb$message[i] <- paste0("VS failed: ", conditionMessage(vs_out))
      }
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
      run_vs = run_vs,
      use_node_imd = use_node_imd,
      imd_path = if (use_fast) "node_score" else "depth_traversal"
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
