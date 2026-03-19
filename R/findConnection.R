#' Find optimal directional connections from fitted RF models
#'
#' @param mod.list A named list of fitted RF models. Model names must follow
#' `response_predictor` convention.
#' @param return_score Logical; whether to return the directional score matrix.
#' @param drop_bottom_q Proportion of models to drop based on quality rank-sum.
#' Must be in `[0, 1)`. Default is `0.2`.
#' @param top_v Optional integer; if set, each row of forest weights keeps top-v
#' entries before scoring. If `NULL`, uses fixed defaults by sample size:
#' `n < 200 -> 15`, `200 <= n <= 800 -> 30`, `n > 800 -> 50`.
#' @param adjust_weights Logical; whether to apply adjusted-weight scaling.
#' @param row_normalize Logical; whether to row-normalize preprocessed weights.
#' @param keep_ties Logical; whether top-v truncation keeps ties at cutoff.
#' @param select_one_per_pair Logical; whether to keep at most one direction per
#' omics pair among quality-filtered models.
#' @param edge_threshold Non-negative threshold for binarizing weights when
#' computing GCC.
#' @param gcc_symm Logical; whether to symmetrize adjacency by `pmax(W, t(W))`
#' before GCC.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A character vector of selected model names (`response_predictor`).
#' If `return_score = TRUE`, returns a list with `model_connection`,
#' `connect_list`, `score`, `top_v_used`, and model-level quality metrics.
# ---------------------------------------------------------------------------------
# Find directional connections from fitted models
# ---------------------------------------------------------------------------------
find_connection <- function(mod.list,
                            return_score = FALSE,
                            drop_bottom_q = 0.2,
                            top_v = NULL,
                            adjust_weights = TRUE,
                            row_normalize = TRUE,
                            keep_ties = TRUE,
                            select_one_per_pair = TRUE,
                            edge_threshold = 0,
                            gcc_symm = TRUE,
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
  if (!is.numeric(edge_threshold) || length(edge_threshold) != 1L ||
      !is.finite(edge_threshold) || edge_threshold < 0) {
    stop("`edge_threshold` must be a single non-negative number.")
  }

  top_v_used <- top_v
  if (is.null(top_v_used)) {
    fw0 <- mod.list[[1]]$forest.wt
    if (is.null(fw0)) {
      stop("Model `", names(mod.list)[1], "` does not contain `forest.wt`.")
    }
    n <- nrow(fw0)
    if (n < 200) {
      top_v_used <- 15L
    } else if (n <= 800) {
      top_v_used <- 30L
    } else {
      top_v_used <- 50L
    }
  } else {
    if (!is.numeric(top_v_used) || length(top_v_used) != 1L || !is.finite(top_v_used)) {
      stop("`top_v` must be NULL or a single finite numeric value.")
    }
    top_v_used <- as.integer(top_v_used)
  }

  model_names <- names(mod.list)

  quality_tbl <- data.frame(
    model = model_names,
    entropy_concentration = NA_real_,
    gcc_ratio = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(model_names)) {
    fw <- mod.list[[i]]$forest.wt
    if (is.null(fw)) {
      stop("Model `", model_names[i], "` does not contain `forest.wt`.")
    }
    W <- prepare_weight_matrix(
      fw,
      adjust = adjust_weights,
      top_v = top_v_used,
      row_normalize = row_normalize,
      zero_diag = TRUE,
      keep_ties = keep_ties
    )
    quality_tbl$entropy_concentration[i] <- calc_weight_concentration(W)
    quality_tbl$gcc_ratio[i] <- calc_gcc_ratio(
      W,
      edge_threshold = edge_threshold,
      symm = gcc_symm
    )
  }

  quality_tbl$rank_entropy <- rank(-quality_tbl$entropy_concentration, ties.method = "average")
  quality_tbl$rank_gcc <- rank(-quality_tbl$gcc_ratio, ties.method = "average")
  quality_tbl$quality_score <- quality_tbl$entropy_concentration + quality_tbl$gcc_ratio
  quality_tbl$rank_sum <- quality_tbl$rank_entropy + quality_tbl$rank_gcc

  n_models <- nrow(quality_tbl)
  n_drop <- floor(n_models * drop_bottom_q)
  n_keep <- max(1L, n_models - n_drop)
  keep_idx <- order(
    quality_tbl$rank_sum,
    -quality_tbl$quality_score,
    quality_tbl$model,
    na.last = TRUE
  )[seq_len(n_keep)]
  keep_flag <- rep(FALSE, n_models)
  keep_flag[keep_idx] <- TRUE
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
    top_v_used = top_v_used,
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
        fit_rfsrc(X = pred, Y = response_d, ...)
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
