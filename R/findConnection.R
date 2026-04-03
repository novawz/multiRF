#' Find optimal directional connections from fitted RF models
#'
#' Scores each directional RF model using modularity of the raw forest weight
#' matrix (and optionally OOB normalized MSE).  Two-step selection:
#' (1) drop the bottom \code{drop_bottom_q} fraction of connections by
#' rank-sum of quality metrics;
#' (2) \code{select_one_per_pair} keeps at most one direction per omics pair
#' (lowest rank-sum).
#'
#' @param mod.list A named list of fitted RF models. Model names must follow
#' `response_predictor` convention.
#' @param return_score Logical; whether to return the directional score matrix.
#' @param drop_bottom_q Proportion of models to drop based on quality rank-sum.
#' Must be in `[0, 1)`. Default is `0.2`.
#' @param select_one_per_pair Logical; whether to keep at most one direction per
#' omics pair among quality-filtered models.
#' @param compute_oob Logical; whether to compute OOB normalized MSE in
#' addition to modularity.  When \code{TRUE},
#' \code{quality_score = modularity - oob_nmse} and rank-sum uses both ranks;
#' when \code{FALSE}, \code{quality_score = modularity} and
#' \code{rank_sum = rank(-modularity)}.
#' Default is \code{FALSE}.
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
                            compute_oob = FALSE,
                            seed = 529L,
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
    quality_tbl$modularity[i] <- calc_modularity(fw, seed = seed)
    if (compute_oob) {
      quality_tbl$oob_nmse[i] <- get_oob_nmse(mod_i)
    }
  }

  # Scoring
  if (compute_oob) {
    quality_tbl$rank_modularity <- rank(-quality_tbl$modularity, ties.method = "average")
    quality_tbl$rank_oob        <- rank( quality_tbl$oob_nmse,   ties.method = "average")
    quality_tbl$quality_score   <- quality_tbl$modularity - quality_tbl$oob_nmse
    quality_tbl$rank_sum        <- quality_tbl$rank_modularity + quality_tbl$rank_oob
  } else {
    quality_tbl$quality_score <- quality_tbl$modularity
    quality_tbl$rank_sum      <- rank(-quality_tbl$modularity, ties.method = "average")
  }

  n_models <- nrow(quality_tbl)

  # --- Step 1: drop bottom fraction by rank_sum ---
  n_drop <- floor(n_models * drop_bottom_q)
  n_keep <- max(1L, n_models - n_drop)
  keep_order <- order(
    quality_tbl$rank_sum,
    -quality_tbl$quality_score,
    quality_tbl$model,
    na.last = TRUE
  )[seq_len(n_keep)]

  keep_flag <- rep(FALSE, n_models)
  keep_flag[keep_order] <- TRUE
  quality_tbl$keep_quality <- keep_flag

  # Parse model names into response / predictor
  split_names <- stringr::str_split(model_names, "_")
  lens <- lengths(split_names)
  if (any(lens < 2)) {
    stop("All model names must contain at least one underscore in `response_predictor` format.")
  }

  response  <- vapply(split_names, `[`, FUN.VALUE = character(1), 1)
  predictor <- vapply(split_names, `[`, FUN.VALUE = character(1), 2)
  pair_id   <- vapply(
    seq_along(model_names),
    function(i) paste(sort(c(response[i], predictor[i])), collapse = "__"),
    FUN.VALUE = character(1)
  )

  # --- Step 2: select one per pair from quality-filtered connections ---
  idx_after_quality <- which(quality_tbl$keep_quality)
  if (length(idx_after_quality) == 0L) {
    idx_after_quality <- which.min(quality_tbl$rank_sum)
  }

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

  quality_tbl$selected <- seq_len(n_models) %in% idx_selected

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
calc_modularity <- function(fw, seed = 529L) {
  W <- as.matrix(fw)
  W[!is.finite(W)] <- 0
  W <- pmax(W, 0)
  W <- pmax(W, t(W))
  diag(W) <- 0
  n <- nrow(W)
  if (n <= 1L) return(0)
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected",
                                            weighted = TRUE, diag = FALSE)
  # Louvain uses randomised node ordering internally; fix RNG state so the
  # modularity score is reproducible across runs.
  old_seed <- .Random.seed
  on.exit({ .Random.seed <<- old_seed }, add = TRUE)
  set.seed(seed)
  cl <- igraph::cluster_louvain(g)
  igraph::modularity(cl)
}
