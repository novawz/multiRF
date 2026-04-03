#' Pairwise IMD Analysis After Model Fitting
#'
#' Run pairwise IMD as a standalone post-fit analysis. By default, this uses
#' the selected features from the variable-selection / robust-clustering stage
#' when available, rather than all variables.
#'
#' @param x An `mrf3_fit` or `mrf3` object containing `net`, `connection`, and
#'   IMD weights.
#' @param feature_source Which feature set to use:
#'   `"selected"` (default), `"weights_nonzero"`, or `"all"`.
#' @param normalized Logical; whether to return degree-normalized adjacency.
#'
#' @return A list with `adj_var_mat`, `adj_dat_mat`, `var_use`, and
#'   `feature_source`.
#' @export
pairwise_imd <- function(x,
                         feature_source = c("selected", "weights_nonzero", "all"),
                         normalized = FALSE) {
  feature_source <- match.arg(feature_source)

  mod <- pairwise_imd_extract_object(x)

  ## ---- Fast path: native engine with pre-computed pairwise_xy ----
  has_precomputed <- inherits(x, "mrf3_fit") && !is.null(x$models) &&
    any(vapply(x$models, function(m) !is.null(m$pairwise_xy), logical(1)))

  if (has_precomputed) {
    message("  Using pre-computed pairwise_xy from native engine.")
    imd_wts <- x$imd
    dat_names <- names(imd_wts)
    var_names_all <- lapply(imd_wts, names)
    var_use <- pairwise_imd_feature_set(x, mod, feature_source = feature_source)

    # Build block-level adjacency
    conn_mat <- do.call(rbind, x$connection)
    adj_dat_mat <- matrix(0, nrow = length(dat_names), ncol = length(dat_names),
                          dimnames = list(dat_names, dat_names))
    for (i in seq_len(nrow(conn_mat))) {
      adj_dat_mat[conn_mat[i, 1], conn_mat[i, 2]] <- 1
      adj_dat_mat[conn_mat[i, 2], conn_mat[i, 1]] <- 1
    }
    length_count <- table(conn_mat)[dat_names]
    length_count[is.na(length_count)] <- 0
    diag(adj_dat_mat) <- as.numeric(length_count)

    # Build variable-level adjacency from pairwise_xy matrices
    large_names <- unlist(var_names_all, use.names = FALSE)
    adj_var_mat <- matrix(0, nrow = length(large_names), ncol = length(large_names),
                          dimnames = list(large_names, large_names))

    for (m_name in names(x$models)) {
      pw <- x$models[[m_name]]$pairwise_xy
      if (is.null(pw)) next
      x_names <- rownames(pw)
      y_names <- colnames(pw)
      common_x <- intersect(x_names, large_names)
      common_y <- intersect(y_names, large_names)
      if (length(common_x) > 0 && length(common_y) > 0) {
        adj_var_mat[common_x, common_y] <- adj_var_mat[common_x, common_y] +
          pw[common_x, common_y, drop = FALSE]
        adj_var_mat[common_y, common_x] <- adj_var_mat[common_y, common_x] +
          t(pw[common_x, common_y, drop = FALSE])
      }
    }

    # Normalize by connection count
    for (i in seq_len(ncol(adj_dat_mat))) {
      for (j in i:ncol(adj_dat_mat)) {
        m1 <- colnames(adj_dat_mat)[i]
        m2 <- colnames(adj_dat_mat)[j]
        if (adj_dat_mat[i, j] != 0) {
          adj_var_mat[var_names_all[[m1]], var_names_all[[m2]]] <-
            adj_var_mat[var_names_all[[m1]], var_names_all[[m2]], drop = FALSE] / adj_dat_mat[i, j]
          if (i != j) {
            adj_var_mat[var_names_all[[m2]], var_names_all[[m1]]] <-
              adj_var_mat[var_names_all[[m2]], var_names_all[[m1]], drop = FALSE] / adj_dat_mat[i, j]
          }
        }
      }
    }

    keep_vars <- unlist(var_use, use.names = FALSE)
    keep_vars <- unique(keep_vars[keep_vars %in% rownames(adj_var_mat)])
    adj_var_mat <- adj_var_mat[keep_vars, keep_vars, drop = FALSE]

  } else {

    ## ---- Slow path: post-hoc tree traversal ----
    if (is.null(mod$net) || is.null(mod$connection) || is.null(mod$imd)) {
      stop("`pairwise_imd()` requires `net`, `connection`, and `imd`, ",
           "or models with pre-computed `pairwise_xy` (native engine).",
           call. = FALSE)
    }

    dat_names <- names(mod$imd)
    if (is.null(dat_names) || any(dat_names == "")) {
      stop("`imd` must be a named list.")
    }

    var_use <- pairwise_imd_feature_set(x, mod, feature_source = feature_source)
    var_names_all <- lapply(mod$imd, names)

    nt <- mod$ntree
    if (is.null(nt) || !is.finite(nt) || nt <= 0) nt <- 1L

    adj_dat_mat <- matrix(0, nrow = length(dat_names), ncol = length(dat_names),
                          dimnames = list(dat_names, dat_names))
    conn_mat <- do.call(rbind, mod$connection)
    length_count <- table(conn_mat)[dat_names]
    length_count[is.na(length_count)] <- 0
    for (i in seq_len(nrow(conn_mat))) {
      adj_dat_mat[conn_mat[i, 1], conn_mat[i, 2]] <- 1
      adj_dat_mat[conn_mat[i, 2], conn_mat[i, 1]] <- 1
    }
    diag(adj_dat_mat) <- as.numeric(length_count)

    large_names <- unlist(var_names_all, use.names = FALSE)
    adj_var_mat <- matrix(0, nrow = length(large_names), ncol = length(large_names),
                          dimnames = list(large_names, large_names))

    plyr::l_ply(
      seq_along(mod$net),
      .fun = function(i) {
        ne <- mod$net[[i]]
        vn <- var_names_all[mod$connection[[i]]]
        names(vn) <- c("yvar", "xvar")
        mat <- pairwise_imd_model(ne, vn)
        mat <- mat[rownames(mat) %in% rownames(adj_var_mat), colnames(mat) %in% colnames(adj_var_mat), drop = FALSE]
        adj_var_mat[rownames(mat), colnames(mat)] <<- adj_var_mat[rownames(mat), colnames(mat), drop = FALSE] + mat
      }
    )

    for (i in seq_len(ncol(adj_dat_mat))) {
      for (j in i:ncol(adj_dat_mat)) {
        m1 <- colnames(adj_dat_mat)[i]
        m2 <- colnames(adj_dat_mat)[j]
        if (adj_dat_mat[i, j] != 0) {
          adj_var_mat[var_names_all[[m1]], var_names_all[[m2]]] <-
            adj_var_mat[var_names_all[[m1]], var_names_all[[m2]], drop = FALSE] / adj_dat_mat[i, j]
          if (i != j) {
            adj_var_mat[var_names_all[[m2]], var_names_all[[m1]]] <-
              adj_var_mat[var_names_all[[m2]], var_names_all[[m1]], drop = FALSE] / adj_dat_mat[i, j]
          }
        }
      }
    }

    keep_vars <- unlist(var_use, use.names = FALSE)
    keep_vars <- unique(keep_vars[keep_vars %in% rownames(adj_var_mat)])
    adj_var_mat <- adj_var_mat[keep_vars, keep_vars, drop = FALSE]
    adj_var_mat <- adj_var_mat / nt
  }

  if (isTRUE(normalized) && nrow(adj_var_mat) > 0L) {
    d <- rowSums(adj_var_mat)
    ok <- d > 0
    inv_sqrt <- rep(0, length(d))
    inv_sqrt[ok] <- d[ok]^(-1 / 2)
    l <- diag(inv_sqrt) %*% adj_var_mat %*% diag(inv_sqrt)
    dimnames(l) <- dimnames(adj_var_mat)
    adj_var_mat <- l
  }

  out <- list(
    adj_var_mat = adj_var_mat,
    adj_dat_mat = adj_dat_mat,
    var_use = var_use,
    feature_source = feature_source
  )
  class(out) <- c("pairwise_imd_analysis", "list")
  out
}


#' Print `pairwise_imd_analysis`
#'
#' @method print pairwise_imd_analysis
#' @param x Output from `pairwise_imd()`.
#' @param ... Unused.
#'
#' @return Input object invisibly.
#' @export
print.pairwise_imd_analysis <- function(x, ...) {
  cat("pairwise_imd_analysis\n")
  cat("  feature_source:", x$feature_source, "\n")
  if (is.list(x$var_use) && length(x$var_use) > 0L) {
    cat("  selected features by block:\n")
    for (nm in names(x$var_use)) {
      cat("   - ", nm, ": ", length(x$var_use[[nm]]), "\n", sep = "")
    }
  }
  if (is.matrix(x$adj_var_mat)) {
    cat("  adj_var_mat:", nrow(x$adj_var_mat), "x", ncol(x$adj_var_mat), "\n")
  }
  if (is.matrix(x$adj_dat_mat)) {
    cat("  adj_dat_mat:", nrow(x$adj_dat_mat), "x", ncol(x$adj_dat_mat), "\n")
  }
  invisible(x)
}


pairwise_imd_extract_object <- function(x) {
  if (inherits(x, "mrf3_fit")) {
    return(list(
      net = x$imd_net,
      connection = x$connection,
      weights = x$imd,
      ntree = x$config$ntree
    ))
  }
  if (inherits(x, "mrf3")) {
    return(list(
      net = x$net,
      connection = x$connection,
      weights = x$imd,
      ntree = x$ntree
    ))
  }
  stop("`x` must be an `mrf3_fit` or `mrf3` object.")
}


pairwise_imd_feature_set <- function(x, mod, feature_source = "selected") {
  if (identical(feature_source, "selected")) {
    if (inherits(x, "mrf3_fit") && is.list(x$selected_data) && length(x$selected_data) > 0L) {
      return(lapply(x$selected_data, colnames))
    }
    if (inherits(x, "mrf3_fit") && is.list(x$selected_vars) && length(x$selected_vars) > 0L) {
      return(x$selected_vars)
    }
    warning("Selected features unavailable; falling back to non-zero IMD weights.", call. = FALSE)
    feature_source <- "weights_nonzero"
  }

  if (identical(feature_source, "weights_nonzero")) {
    return(lapply(mod$imd, function(w) names(w)[is.finite(w) & (w != 0)]))
  }

  lapply(mod$imd, names)
}


pairwise_imd_model <- function(net, var_name) {
  var_mat <- matrix(0,
                    nrow = length(var_name$xvar) + length(var_name$yvar),
                    ncol = length(var_name$xvar) + length(var_name$yvar))
  dimnames(var_mat) <- list(c(var_name$xvar, var_name$yvar), c(var_name$xvar, var_name$yvar))

  plyr::l_ply(
    net,
    .fun = function(i) {
      netdf <- data.frame(
        y = i$Y_id,
        x = gsub("^(.*)_.*", "\\1", i$from),
        imd = i$inv_d
      )
      netdf <- stats::na.omit(netdf)
      netdf <- dplyr::filter(netdf, y %in% var_name$yvar & x %in% var_name$xvar)
      if (nrow(netdf) == 0L) {
        return(NULL)
      }
      plyr::l_ply(seq_len(nrow(netdf)), .fun = function(j) {
        var_mat[netdf[j, "x"], netdf[j, "y"]] <<- var_mat[netdf[j, "x"], netdf[j, "y"]] + netdf[j, "imd"]
        var_mat[netdf[j, "y"], netdf[j, "x"]] <<- var_mat[netdf[j, "y"], netdf[j, "x"]] + netdf[j, "imd"]
      })
      NULL
    }
  )

  var_mat
}
