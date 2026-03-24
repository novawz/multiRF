#' Find optimal directional connections from fitted RF models
#'
#' Scores each directional RF model using modularity of the raw forest weight
#' matrix and OOB normalized MSE.  A two-step selection filters connections:
#' first an absolute quality gate removes clearly uninformative connections,
#' then rank-sum selection picks the best among survivors.
#'
#' @param mod.list A named list of fitted RF models. Model names must follow
#' `response_predictor` convention.
#' @param return_score Logical; whether to return the directional score matrix.
#' @param drop_bottom_q Proportion of models to drop based on quality rank-sum.
#' Must be in `[0, 1)`. Default is `0.2`.
#' @param select_one_per_pair Logical; whether to keep at most one direction per
#' omics pair among quality-filtered models.
#' @param quality_gate Character or NULL. Method for absolute quality filtering
#' before rank-sum selection. \code{"adaptive"} uses k-means (k=2) on
#' quality_score to separate good and bad connections when a clear gap exists.
#' \code{NULL} disables absolute filtering (backward compatible).
#' Default is \code{"adaptive"}.
#' @param min_modularity Numeric or NULL. If set, connections with modularity
#' below this value are excluded. Applied in addition to \code{quality_gate}.
#' @param max_oob_nmse Numeric or NULL. If set, connections with OOB nMSE
#' above this value are excluded. Applied in addition to \code{quality_gate}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A character vector of selected model names (`response_predictor`).
#' If `return_score = TRUE`, returns a list with `model_connection`,
#' `connect_list`, `score`, and model-level quality metrics.
# ---------------------------------------------------------------------------------
# Find directional connections from fitted models
# ---------------------------------------------------------------------------------
find_connection <- function(mod.list,
                            return_score = FALSE,
                            drop_bottom_q = 0.2,
                            select_one_per_pair = TRUE,
                            quality_gate = "adaptive",
                            min_modularity = NULL,
                            max_oob_nmse = NULL,
                            ...){

  if (inherits(mod.list, "mrf3")) {
    mod.list <- mod.list$mod
  }

  if (!is.list(mod.list) || length(mod.list) == 0) {
    stop("`mod.list` must be a non-empty list of fitted RF models.")
  }
  if (is.null(names(mod.list)) || any(names(mod.list) == "")) {
    stop("`mod.list` must be named using `response_predictor` format.")
  }
  if (!is.numeric(drop_bottom_q) || length(drop_bottom_q) != 1L ||
      !is.finite(drop_bottom_q) || drop_bottom_q < 0 || drop_bottom_q >= 1) {
    stop("`drop_bottom_q` must be a single numeric in [0, 1).")
  }

  model_names <- names(mod.list)

  quality_tbl <- data.frame(
    model = model_names,
    modularity = NA_real_,
    oob_nmse = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(model_names)) {
    mod_i <- mod.list[[i]]
    fw <- mod_i$forest.wt
    if (is.null(fw)) {
      stop("Model `", model_names[i], "` does not contain `forest.wt`.")
    }
    quality_tbl$modularity[i] <- calc_modularity(fw)
    quality_tbl$oob_nmse[i]   <- get_oob_nmse(mod_i)
  }

  # modularity: higher = better → rank descending
  # oob_nmse:   lower = better  → rank ascending
  quality_tbl$rank_modularity <- rank(-quality_tbl$modularity, ties.method = "average")
  quality_tbl$rank_oob        <- rank( quality_tbl$oob_nmse,   ties.method = "average")
  quality_tbl$quality_score   <- quality_tbl$modularity - quality_tbl$oob_nmse
  quality_tbl$rank_sum        <- quality_tbl$rank_modularity + quality_tbl$rank_oob

  # --- Step 1: Absolute quality gate ---
  # Filter out connections that are clearly uninformative before rank-based selection.
  quality_tbl$pass_gate <- TRUE
  n_models <- nrow(quality_tbl)

  # Hard thresholds (if specified)
  if (!is.null(min_modularity)) {
    quality_tbl$pass_gate <- quality_tbl$pass_gate & (quality_tbl$modularity >= min_modularity)
  }
  if (!is.null(max_oob_nmse)) {
    quality_tbl$pass_gate <- quality_tbl$pass_gate & (quality_tbl$oob_nmse <= max_oob_nmse)
  }

  # Adaptive gate: k-means(k=2) on quality_score
  if (!is.null(quality_gate) && identical(quality_gate, "adaptive") && n_models >= 4L) {
    qs <- quality_tbl$quality_score
    if (length(unique(qs)) >= 2L) {
      km <- stats::kmeans(qs, centers = 2L, nstart = 10L)
      # Identify the "good" cluster (higher quality_score = higher mean)
      good_cluster <- which.max(km$centers[, 1])
      bad_cluster  <- which.min(km$centers[, 1])
      gap <- km$centers[good_cluster, 1] - km$centers[bad_cluster, 1]
      # Only apply gate if there's a meaningful gap between clusters.
      # Use pooled within-cluster SD as reference; gate if gap > 1 SD.
      wss_good <- sum((qs[km$cluster == good_cluster] - km$centers[good_cluster, 1])^2)
      wss_bad  <- sum((qs[km$cluster == bad_cluster] - km$centers[bad_cluster, 1])^2)
      n_good <- sum(km$cluster == good_cluster)
      n_bad  <- sum(km$cluster == bad_cluster)
      pooled_sd <- sqrt((wss_good + wss_bad) / max(1, n_good + n_bad - 2L))
      # If pooled_sd == 0, separation is perfect → always apply gate.
      # Otherwise, require gap > 1 pooled SD.
      if (pooled_sd == 0 && gap > 0) {
        quality_tbl$pass_gate <- quality_tbl$pass_gate & (km$cluster == good_cluster)
      } else if (pooled_sd > 0 && gap / pooled_sd > 1.0) {
        quality_tbl$pass_gate <- quality_tbl$pass_gate & (km$cluster == good_cluster)
      }
    }
  }

  # --- Step 2: Rank-based filtering among gate-passed connections ---
  idx_passed <- which(quality_tbl$pass_gate)
  if (length(idx_passed) == 0L) {
    # All connections failed the gate → fall back to best single connection
    warning("All connections failed quality gate. Selecting the single best connection.")
    idx_passed <- which.min(quality_tbl$rank_sum)
  }

  # Within passed connections, apply drop_bottom_q
  n_passed <- length(idx_passed)
  n_drop <- floor(n_passed * drop_bottom_q)
  n_keep <- max(1L, n_passed - n_drop)
  keep_order <- order(
    quality_tbl$rank_sum[idx_passed],
    -quality_tbl$quality_score[idx_passed],
    quality_tbl$model[idx_passed],
    na.last = TRUE
  )
  idx_after_quality <- idx_passed[keep_order[seq_len(n_keep)]]

  keep_flag <- rep(FALSE, n_models)
  keep_flag[idx_after_quality] <- TRUE
  quality_tbl$keep_quality <- keep_flag

  split_names <- stringr::str_split(model_names, "_")
  lens <- lengths(split_names)
  if (any(lens < 2)) {
    stop("All model names must contain at least one underscore in `response_predictor` format.")
  }

  response <- vapply(split_names, `[`, FUN.VALUE = character(1), 1)
  predictor <- vapply(split_names, `[`, FUN.VALUE = character(1), 2)
  pair_id <- vapply(
    seq_along(model_names),
    function(i) paste(sort(c(response[i], predictor[i])), collapse = "__"),
    FUN.VALUE = character(1)
  )

  if (select_one_per_pair) {
    idx_selected <- unlist(
      lapply(split(idx_after_quality, pair_id[idx_after_quality]), function(idx) {
        o <- order(
          quality_tbl$rank_sum[idx],
          -quality_tbl$quality_score[idx],
          quality_tbl$model[idx]
        )
        idx[o[1]]
      }),
      use.names = FALSE
    )
  } else {
    idx_selected <- idx_after_quality
  }
  idx_selected <- sort(idx_selected)

  model_connection <- model_names[idx_selected]
  connect_list <- lapply(
    idx_selected,
    function(i) c(response[i], predictor[i])
  )

  if (!return_score) {
    return(model_connection)
  }

  dat_names <- unique(c(response, predictor))
  score <- matrix(NA_real_, nrow = length(dat_names), ncol = length(dat_names))
  dimnames(score) <- list(dat_names, dat_names)
  score[cbind(match(response, dat_names), match(predictor, dat_names))] <- as.numeric(quality_tbl$quality_score)
  diag(score) <- NA_real_

  list(
    model_connection = model_connection,
    connect_list = connect_list,
    score = score,
    top_v_used = NULL,
    quality = quality_tbl
  )
}

# Backward-compatible aliases for older scripts that still use camelCase names.
findConnection <- function(mod.list, ...) {
  find_connection(mod.list, ...)
}

#' @rdname find_connection
#' @param dat.list A named list of omics data blocks with samples in rows and
#' features in columns.

full_connect <- function(dat.list, ...){

  dat_names <- names(dat.list)
  mod_l <- plyr::llply(
    dat_names,
    .fun = function(d){
      response_d <- dat.list[[d]]
      predict_d <- dat.list[!names(dat.list) %in% d]

      mod <- plyr::llply(predict_d, .fun = function(pred){
        fit_forest(X = pred, Y = response_d, ...)
      })

      mod_names <- paste0(d, "_", names(predict_d))
      names(mod) <- mod_names

      return(mod)
    }
  )

  mod_l <- Reduce(c, mod_l)

  return(mod_l)
}

fullConnect <- function(dat.list, ...) {
  full_connect(dat.list, ...)
}
get_r_sq <- function(mod){
  
  fw <- mod$forest.wt
  
  ex <- mean(colMeans((fw %*% as.matrix(mod$xvar) - as.matrix(mod$xvar))^2)/matrixStats::colVars(as.matrix(mod$xvar)))
  if(is.null(mod$yvar)){
    return(ex)
  } else {
    if(!is.null(class(mod)[3]) && identical(class(mod)[3], "class+")){
      ey <- na.omit(mod$err.rate)
    } else {
      ey <- mean(colMeans((fw %*% as.matrix(mod$yvar) - as.matrix(mod$yvar))^2)/matrixStats::colVars(as.matrix(mod$yvar)))
    }
  }
  
  ey+ex
  
}

#' Compute OOB forest weight matrix
#'
#' For each sample, only trees where it was out-of-bag contribute to the
#' weight row.  Within each OOB tree, co-leaf samples contribute weight
#' proportional to their bootstrap frequency, normalized by the total
#' bootstrap mass in the leaf.
#'
#' @param mod A fitted model object from \code{fit_forest} (must contain
#'   \code{$membership} and \code{$inbag}).
#' @return An n x n numeric matrix of OOB forest weights.
#' @export
compute_oob_forest_wt <- function(mod) {
  if (is.null(mod$membership) || is.null(mod$inbag)) {
    stop("Model must contain $membership and $inbag matrices.")
  }
  W <- compute_oob_forest_wt_cpp(mod$membership, mod$inbag)
  diag(W) <- 0
  W
}

#' OOB normalized MSE for a fitted forest model
#'
#' Uses OOB forest weights to predict both xvar and yvar, then computes
#' the mean column-wise normalized MSE.  Lower is better.
#'
#' @param mod A fitted model from \code{fit_forest}.
#' @return A single numeric value (normalized MSE).
#' @export
get_oob_nmse <- function(mod) {
  W <- compute_oob_forest_wt(mod)
  X <- as.matrix(mod$xvar)
  pred_x <- W %*% X
  col_var_x <- apply(X, 2, var)
  col_var_x[col_var_x < 1e-12] <- 1e-12
  ex <- mean(colMeans((pred_x - X)^2) / col_var_x)

  if (is.null(mod$yvar)) return(ex)
  Y <- as.matrix(mod$yvar)
  pred_y <- W %*% Y
  col_var_y <- apply(Y, 2, var)
  col_var_y[col_var_y < 1e-12] <- 1e-12
  ey <- mean(colMeans((pred_y - Y)^2) / col_var_y)
  ex + ey
}

calc_weight_concentration <- function(W, eps = 1e-12) {
  W <- as.matrix(W)
  W[!is.finite(W)] <- 0
  W <- pmax(W, 0)

  rs <- rowSums(W)
  p <- W
  ok <- rs > eps
  if (any(ok)) {
    p[ok, ] <- W[ok, , drop = FALSE] / rs[ok]
  }
  if (any(!ok)) {
    p[!ok, ] <- 0
  }

  row_entropy <- apply(
    p,
    1,
    function(v) {
      v <- v[v > eps]
      k <- length(v)
      if (k <= 1L) {
        return(0)
      }
      -sum(v * log(v)) / log(k)
    }
  )
  c0 <- 1 - mean(row_entropy)
  max(0, min(1, c0))
}


calc_gcc_ratio <- function(W, edge_threshold = 0, symm = TRUE) {
  W <- as.matrix(W)
  W[!is.finite(W)] <- 0
  W <- pmax(W, 0)
  if (symm) {
    W <- pmax(W, t(W))
  }
  diag(W) <- 0
  A <- ifelse(W > edge_threshold, 1, 0)

  n <- nrow(A)
  if (n <= 1L) {
    return(1)
  }

  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  comp <- igraph::components(g)
  max(comp$csize) / n
}


#' Modularity of forest weight matrix
#'
#' Symmetrizes the raw forest weight matrix and computes weighted modularity
#' via Louvain community detection.  Higher modularity indicates clearer
#' block / cluster structure in the proximity graph.
#'
#' @param fw A raw forest weight matrix (n x n).
#' @return A single numeric modularity value (typically 0 to ~0.8).
#' @keywords internal
calc_modularity <- function(fw) {
  W <- as.matrix(fw)
  W[!is.finite(W)] <- 0
  W <- pmax(W, 0)
  W <- pmax(W, t(W))
  diag(W) <- 0
  n <- nrow(W)
  if (n <= 1L) return(0)
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected",
                                            weighted = TRUE, diag = FALSE)
  cl <- igraph::cluster_louvain(g)
  igraph::modularity(cl)
}
