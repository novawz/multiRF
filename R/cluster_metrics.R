#' Clustering Metrics
#'
#' Utility functions for comparing two partitions (predicted vs reference),
#' including ARI, partition-level Jaccard, NMI, and purity.
#'
#' @param pred Predicted cluster labels.
#' @param ref Reference cluster labels.
#' @param na.rm Logical; whether to remove pairs with missing labels before computing metrics.
#'
#' @return For scalar metric functions, a numeric scalar (`NA` when undefined).
#' @export
cluster_ari <- function(pred, ref, na.rm = TRUE) {
  labs <- prepare_cluster_pair(pred, ref, na.rm = na.rm)
  if (is.null(labs)) {
    return(NA_real_)
  }
  adjusted_rand_index_internal(labs$pred, labs$ref)
}

#' @rdname cluster_ari
#' @export
cluster_jaccard <- function(pred, ref, na.rm = TRUE) {
  labs <- prepare_cluster_pair(pred, ref, na.rm = na.rm)
  if (is.null(labs)) {
    return(NA_real_)
  }
  partition_jaccard_index_internal(labs$pred, labs$ref)
}

#' @rdname cluster_ari
#' @export
cluster_nmi <- function(pred, ref, na.rm = TRUE) {
  labs <- prepare_cluster_pair(pred, ref, na.rm = na.rm)
  if (is.null(labs)) {
    return(NA_real_)
  }
  normalized_mutual_information_internal(labs$pred, labs$ref)
}

#' @rdname cluster_ari
#' @export
cluster_purity <- function(pred, ref, na.rm = TRUE) {
  labs <- prepare_cluster_pair(pred, ref, na.rm = na.rm)
  if (is.null(labs)) {
    return(NA_real_)
  }
  tb <- table(labs$pred, labs$ref)
  if (sum(tb) == 0L) {
    return(NA_real_)
  }
  sum(apply(tb, 1, max)) / sum(tb)
}


#' Compute multiple clustering metrics at once
#'
#' @param pred Predicted cluster labels.
#' @param ref Reference cluster labels.
#' @param na.rm Logical; whether to remove pairs with missing labels.
#' @param as_tibble Logical; whether to return a tibble.
#'
#' @return A one-row data frame/tibble with columns:
#' `n`, `ari`, `jaccard`, `nmi`, `purity`.
#' @export
cluster_metrics <- function(pred, ref, na.rm = TRUE, as_tibble = FALSE) {
  labs <- prepare_cluster_pair(pred, ref, na.rm = na.rm)
  n_use <- if (is.null(labs)) 0L else length(labs$pred)

  out <- data.frame(
    n = n_use,
    ari = if (is.null(labs)) NA_real_ else adjusted_rand_index_internal(labs$pred, labs$ref),
    jaccard = if (is.null(labs)) NA_real_ else partition_jaccard_index_internal(labs$pred, labs$ref),
    nmi = if (is.null(labs)) NA_real_ else normalized_mutual_information_internal(labs$pred, labs$ref),
    purity = if (is.null(labs)) NA_real_ else {
      tb <- table(labs$pred, labs$ref)
      if (sum(tb) == 0L) NA_real_ else sum(apply(tb, 1, max)) / sum(tb)
    },
    stringsAsFactors = FALSE
  )

  if (isTRUE(as_tibble)) {
    return(tibble::as_tibble(out))
  }
  out
}


#' Pairwise metric matrix across multiple clusterings
#'
#' @param cluster_list Named list of cluster-label vectors.
#' All vectors must be alignable to the same samples.
#' @param metric Metric to compute: `"ari"`, `"jaccard"`, `"nmi"`, or `"purity"`.
#' @param na.rm Logical; whether to remove missing labels pairwise.
#'
#' @return A symmetric matrix of pairwise metric values.
#' @export
cluster_metric_matrix <- function(cluster_list,
                                  metric = c("ari", "jaccard", "nmi", "purity"),
                                  na.rm = TRUE) {
  metric <- match.arg(metric)
  if (!is.list(cluster_list) || length(cluster_list) < 2L) {
    stop("`cluster_list` must be a list with at least 2 label vectors.")
  }
  if (is.null(names(cluster_list)) || any(names(cluster_list) == "")) {
    names(cluster_list) <- paste0("cluster_", seq_along(cluster_list))
  }

  aligned <- align_cluster_list(cluster_list)
  k <- length(aligned)
  out <- matrix(
    NA_real_,
    nrow = k,
    ncol = k,
    dimnames = list(names(aligned), names(aligned))
  )

  metric_fun <- switch(
    metric,
    ari = cluster_ari,
    jaccard = cluster_jaccard,
    nmi = cluster_nmi,
    purity = cluster_purity
  )

  for (i in seq_len(k)) {
    out[i, i] <- 1
    if (i == k) {
      next
    }
    for (j in seq.int(i + 1L, k)) {
      val <- metric_fun(aligned[[i]], aligned[[j]], na.rm = na.rm)
      out[i, j] <- val
      out[j, i] <- val
    }
  }
  out
}


# ---- helpers ----
prepare_cluster_pair <- function(pred, ref, na.rm = TRUE) {
  if (length(pred) != length(ref)) {
    stop("`pred` and `ref` must have the same length.")
  }
  pred <- as.character(pred)
  ref <- as.character(ref)
  if (isTRUE(na.rm)) {
    ok <- !is.na(pred) & !is.na(ref)
    pred <- pred[ok]
    ref <- ref[ok]
  }
  if (length(pred) < 2L) {
    return(NULL)
  }
  list(pred = pred, ref = ref)
}

align_cluster_list <- function(cluster_list) {
  all_named <- all(vapply(cluster_list, function(x) !is.null(names(x)), logical(1)))

  if (all_named) {
    common_ids <- Reduce(intersect, lapply(cluster_list, names))
    if (length(common_ids) < 2L) {
      stop("Named cluster vectors do not share enough common sample ids.")
    }
    out <- lapply(cluster_list, function(x) as.character(x[common_ids]))
    names(out) <- names(cluster_list)
    return(out)
  }

  n <- unique(vapply(cluster_list, length, integer(1)))
  if (length(n) != 1L) {
    stop("Unnamed cluster vectors must have the same length.")
  }
  lapply(cluster_list, as.character)
}

adjusted_rand_index_internal <- function(x, y) {
  tb <- table(x, y)
  n <- sum(tb)
  if (n <= 1L) {
    return(NA_real_)
  }
  comb2 <- function(v) sum(v * (v - 1) / 2)
  nij <- comb2(tb)
  ai <- comb2(rowSums(tb))
  bj <- comb2(colSums(tb))
  expected <- ai * bj / (n * (n - 1) / 2)
  max_index <- 0.5 * (ai + bj)
  den <- max_index - expected
  if (!is.finite(den) || den == 0) {
    return(NA_real_)
  }
  (nij - expected) / den
}

partition_jaccard_index_internal <- function(x, y) {
  n <- length(x)
  if (length(y) != n || n < 2L) {
    return(NA_real_)
  }
  sx <- outer(x, x, "==")
  sy <- outer(y, y, "==")
  idx <- upper.tri(sx, diag = FALSE)
  a <- sx[idx]
  b <- sy[idx]
  union <- sum(a | b)
  if (union == 0L) {
    return(NA_real_)
  }
  sum(a & b) / union
}

normalized_mutual_information_internal <- function(x, y) {
  tb <- table(x, y)
  n <- sum(tb)
  if (n <= 0L) {
    return(NA_real_)
  }
  pxy <- tb / n
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- pxy > 0
  mi <- sum(pxy[nz] * log(pxy[nz] / (px[row(pxy)[nz]] * py[col(pxy)[nz]])))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  den <- sqrt(hx * hy)
  if (!is.finite(den) || den == 0) {
    return(NA_real_)
  }
  mi / den
}
