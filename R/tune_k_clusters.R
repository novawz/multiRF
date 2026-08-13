# -------------------------------------------------------------------------------------------------------------
#' k clusters tuning function
#' @param x A proximity/similarity matrix (or embedding input) for clustering.
#' @param return_cluster A logical parameter that determine whether the clustering results to be returned
#' @param plot_k A logical parameter that determine whether the plot of tuning should be showed
#' @param method Clustering backend: `"PAM"` or `"Spectral"`.
#' @param tune_method Tuning criterion for PAM (`"silhouette"` or `"ratio"`).
#' @param gap_w Weighting scheme for spectral eigengap (`"uniform"` or `"log"`).
#' @param prox Logical; whether `x` is a proximity matrix (converted to dissimilarity for PAM).
#' @param k_tune Candidate cluster counts to evaluate.
#' @param d Optional degree vector used by spectral clustering.
#' @param ... Additional arguments passed to lower-level clustering functions.
#'
#' @return optimal number of k
#' @inheritParams cluster::pam

tune_k_clusters <- function(x,...){
  UseMethod("tune_k_clusters")
}

#' @rdname tune_k_clusters
tune_k_clusters.default <- function(x, return_cluster = FALSE, plot_k = FALSE,
                                    method = "Spectral", tune_method = "silhouette", gap_w = "uniform", prox = FALSE,...){

  n_obs <- if (is.matrix(x) || is.data.frame(x)) nrow(x) else length(x)
  k_max <- min(12L, as.integer(n_obs) - 1L)
  if (!is.finite(k_max) || k_max < 2L) {
    stop("At least three samples are required to tune the number of clusters.")
  }
  k_tune <- seq.int(2L, k_max)

  if(method == "Spectral"){
    cl <- spectral_cl(x, k_tune = k_tune, gap_w = gap_w,  ...)
  }

  if(method == "PAM"){
    if(prox){
      x <- 1-x
      diss <- TRUE
    } else {diss <- FALSE}
    cl <- pam_cl(x, k_tune = k_tune, diss = diss, tune_method = tune_method, ...)
  }

  if(plot_k){
    if(method == "PAM"){
      if(tune_method == "silhouette"){
        m <- 2:9
        ylab <- "Silhouette"
      } 
      if(tune_method == "ratio"){
        ylab <- "Ratio"
        m <- 2:10
        par(mfrow = c(1,2))
        col <- rep(1,8)
        col[cl$best_k - 1] <- 2
        plot(2:9, cl$diff, pch = 19, ylab = "Ratio Diff", xlab = "k", col = col)
      } 
      col <- rep(1,length(m))
      col[cl$best_k - 1] <- 2
      plot(m, cl$sil, pch = 19, ylab = ylab, xlab = "k", col = col)
    } else {
      par(mfrow = c(1,2))
      col <- rep(1,8)
      col[cl$best_k - 1] <- 2
      plot(2:9, cl$diff_e, pch = 19, ylab = "Diff eigenvalue", xlab = "k", col = col)
      col <- rep(1,9)
      col[cl$best_k - 1] <- 2
      plot(2:10, cl$ev, pch = 19, ylab = "Eigenvalues", xlab = "k", col = col)
    } 
      
  
  }


  if(return_cluster){
    # cl <- cluster::pam(p, k = k, diss = T, cluster.only = T, ...)
    k <- cl
  } else {k <- cl$best_k}

  return(k)
}

#' @rdname tune_k_clusters
tune_k_clusters.mrf3 <- function(x, return_cluster = FALSE, plot_k = FALSE, method = "Spectral", tune_method = "silhouette", gap_w = "uniform", prox = FALSE, ...){

  if(class(x)[2] %in% "cl"){
    x <- x$dat
  } else {
    stop("Must be a mrf cl object")
  }

  n_obs <- nrow(x)
  k_max <- min(12L, as.integer(n_obs) - 1L)
  if (!is.finite(k_max) || k_max < 2L) {
    stop("At least three samples are required to tune the number of clusters.")
  }
  k_tune <- seq.int(2L, k_max)

  if(method == "Spectral"){
    cl <- spectral_cl(x, k_tune = k_tune, gap_w = gap_w, ...)
  }

  if(method == "PAM"){
    if(prox){
      x <- 1-x
      diss <- TRUE
    } else {diss <- FALSE}
    cl <- pam_cl(x, k_tune = k_tune, diss = diss, tune_method = tune_method, ...)
  }
  
  if(plot_k){
    if(method == "PAM"){
      if(tune_method == "silhouette"){
        m <- 2:9
        ylab <- "Silhouette"
      } 
      if(tune_method == "ratio"){
        ylab <- "Ratio"
        m <- 2:10
        par(mfrow = c(1,2))
        col <- rep(1,8)
        col[cl$best_k - 1] <- 2
        plot(2:9, cl$diff, pch = 19, ylab = "Ratio Diff", xlab = "k", col = col)
      } 
      col <- rep(1,length(m))
      col[cl$best_k - 1] <- 2
      plot(m, cl$sil, pch = 19, ylab = ylab, xlab = "k", col = col)
    } else {
      par(mfrow = c(1,2))
      col <- rep(1,8)
      col[cl$best_k - 1] <- 2
      plot(2:9, cl$diff_e, pch = 19, ylab = "Diff eigenvalue", xlab = "k", col = col)
      col <- rep(1,9)
      col[cl$best_k - 1] <- 2
      plot(2:10, cl$ev, pch = 19, ylab = "Eigenvalues", xlab = "k", col = col)
    } 
    
    
  }
  
  
  if(return_cluster){
    # cl <- cluster::pam(p, k = k, diss = T, cluster.only = T, ...)
    k <- cl
  } else {k <- cl$best_k}
  return(k)
}

#' @rdname tune_k_clusters
#' @export
spectral_cl <- function(x, k_tune = seq(2,12,by = 1), gap_w = "uniform", d = NULL, ...){
  x <- as.matrix(x)
  if (!is.numeric(x) || nrow(x) != ncol(x)) {
    stop("`x` must be a square numeric matrix.")
  }
  x[!is.finite(x)] <- 0
  x <- (x + t(x)) / 2
  x[x < 0] <- 0
  n <- nrow(x)
  k_grid <- unique(as.integer(k_tune[is.finite(k_tune)]))
  k_grid <- k_grid[k_grid >= 2L & k_grid < n]
  if (length(k_grid) == 0L) {
    stop("`k_tune` must contain at least one integer in [2, n - 1].")
  }

  if (is.null(d)) {
    d <- rowSums(x)
  }
  d <- as.numeric(d)
  if (length(d) != nrow(x)) {
    d <- rowSums(x)
  }
  eps <- 1e-8
  d[!is.finite(d)] <- eps
  d[d <= eps] <- eps
  inv_sqrt_d <- 1 / sqrt(d)
  l <- diag(1, nrow(x)) - diag(inv_sqrt_d) %*% x %*% diag(inv_sqrt_d)
  l[!is.finite(l)] <- 0
  l <- (l + t(l)) / 2

  e <- tryCatch(eigen(l, symmetric = TRUE), error = function(e) NULL)
  if (is.null(e)) {
    pam_fit <- pam_cl(1 - x, k_tune = k_grid, diss = TRUE, tune_method = "silhouette", ...)
    return(list(
      best_k = pam_fit$best_k,
      cl = pam_fit$cl,
      diff_e = NA_real_,
      ev = NA_real_,
      embed = NULL,
      obj = pam_fit$obj,
      sil = NA_real_
    ))
  }
  eigenvectors <- e$vectors
  lambda <- rev(e$values)
  gap_w <- match.arg(gap_w, c("uniform", "log"))
  w <- if (identical(gap_w, "log")) log(k_grid) else rep(1, length(k_grid))
  raw_gap <- lambda[k_grid + 1L] - lambda[k_grid]
  diff_e <- raw_gap * w
  k <- k_grid[which.max(diff_e)]
  eigenvalues <- lambda[seq_len(min(n, max(k_grid) + 1L))]
  
  mat <- as.matrix(eigenvectors[,(n-k+1):(n)])
  mat_norm <- mat / pmax(sqrt(rowSums(mat^2)), .Machine$double.eps)
 
  cl <- cluster::pam(mat, k = k, ...)

  embed <- mat_norm
  
  return(list(
    best_k = k, cl = cl$cluster, diff_e = diff_e, ev = eigenvalues, embed = embed, obj = cl$objective[2], sil = cl$silinfo$avg.width
  ))
}

#' @rdname tune_k_clusters
#' @export
pam_cl <- function(x, k_tune = seq(2,9,by = 1), diss = TRUE, tune_method = "silhouette", ...){
  tune_method <- match.arg(tune_method, c("silhouette", "ratio"))
  n <- if (is.matrix(x) || is.data.frame(x)) nrow(x) else attr(x, "Size")
  if (is.null(n) || !is.finite(n)) n <- length(x)
  k_tune <- unique(as.integer(k_tune[is.finite(k_tune)]))
  k_tune <- k_tune[k_tune >= 2L & k_tune < n]
  if (length(k_tune) == 0L) {
    stop("`k_tune` must contain at least one integer in [2, n - 1].")
  }

  sil <- numeric(0)
  diss_use <- isTRUE(diss)
  if(length(k_tune) > 1){
    if(!diss_use){
      x <- as.matrix(dist(x))
      diss_use <- TRUE
    }
    
    if(tune_method == "ratio"){
      extra_k <- max(k_tune) + 1L
      k_eval <- if (extra_k < n) c(k_tune, extra_k) else k_tune
      m <- 1-x
      dist_all <- mean(m)
    } else {
      k_eval <- k_tune
    }
    sil <- numeric(length(k_eval))
    for (ii in seq_along(k_eval)){
      k <- k_eval[ii]
      cl <- cluster::pam(x, k = k, diss = diss_use, ...)
      if(tune_method == "silhouette"){
        sil[ii] <- cl$objective[1] - cl$objective[2]
      }
      if(tune_method == "ratio"){
        cl_unique <- unique(cl$cluster)
        
        class_dist <- sapply(cl_unique, function(i){
          mean(m[cl$cluster == i,cl$cluster != i])
        })
        
        sil[ii] <- (mean(class_dist))/dist_all
      }
     
    }
    if(tune_method == "silhouette"){
      k <- k_eval[which.max(sil)]
      diff_S <- NULL
    }
    if(tune_method == "ratio"){
      diff_S <- abs(diff(sil)) * log(k_eval[-length(k_eval)])
      if (length(diff_S) == 0L || !any(is.finite(diff_S))) {
        k <- k_tune[1L]
      } else {
        k <- k_eval[which.max(diff_S)]
      }
    }
  } else {
    k <- k_tune[1L]
    diff_S <- NULL
  }
  
  cl <- cluster::pam(x, k = k, diss = diss_use, ...)
  
  return(
    list(best_k = k, cl = cl$cluster, sil = sil, diff = diff_S, obj = cl$objective[2])
  )
  
}
