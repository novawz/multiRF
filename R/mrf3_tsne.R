#' compute tSNE embedding for MRF similarity matrix
#' @param mod mrf3 model
#' @param dims Number of dims
#' @param max_iter integer; Number of iterations (default: 1000)
#' @param learning_rate numeric; Learning rate (default: 200.0)
#' @param verbose Logical; whether to print optimization progress.
#' @param tol Convergence tolerance.
#' @param patience Early-stopping patience measured in iterations.
#' @return Two dimensional tSNE embedding
#'
#' @export mrf3_tsne
#' 

mrf3_tsne <- function(mod,  dims = 2,  max_iter = 1000, learning_rate = 200, verbose = FALSE, tol = 1e-05, patience = 10) {
  
  X <- t(apply(mod$dat, 1, function(row) row / sum(row)))
  X <- X/sum(X)
  
  tsne(X = X, 
       dims = dims, 
       max_iter = max_iter, 
       learning_rate = learning_rate, 
       verbose = verbose, 
       patience = patience,
       tol = tol)
}

tsne <- function(X, dims = 2,  max_iter = 1000, learning_rate = 200, verbose = FALSE, patience = 10, tol = 1e-05) {
  n <- nrow(X)
  Y <- matrix(rnorm(n * dims), n, dims)
  cost_history <- rep(NA_real_, max_iter)
  patience_count <- 0
  # Faster squared Euclidean distance matrix via matrix algebra.
  sq_dist_mat <- function(A) {
    rs <- rowSums(A * A)
    D <- outer(rs, rs, "+") - 2 * tcrossprod(A)
    pmax(D, 0)
  }
  
  for (iter in 1:max_iter) {
    grad_Y <- tsne_cost_gradient_cpp(Y, X)
    Y <- Y - learning_rate * grad_Y
    Q <- 1 / (1 + sq_dist_mat(Y))
    diag(Q) <- 0
    Q <- Q / sum(Q)
    
    cost <- sum(X * log((X + .Machine$double.eps) / (Q + .Machine$double.eps)))
    cost_history[iter] <- cost

    
    # Check convergence every patience iterations
    if (iter > patience) {
      recent_costs <- cost_history[(iter - patience + 1):iter]

      max_diff <- abs(max(recent_costs[-1] - recent_costs[-patience]))

      if (max_diff < tol) {
        patience_count <- patience_count + 1
      } else {
        patience_count <- 0
      }
      
      if (patience_count >= patience) {
        converged <- TRUE
        break
      }
    }
    if(verbose) {
      if (iter %% 10 == 0) {
        cat("Iteration", iter, "Cost", cost, "\n")
      }
    }
   
  }
  
  colnames(Y) <- paste0("tSNE", 1:ncol(Y))
  return(Y)
}
