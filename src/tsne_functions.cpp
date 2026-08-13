
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix tsne_cost_gradient_cpp(NumericMatrix Y, NumericMatrix P) {
  int n = Y.nrow();
  int d = Y.ncol();

  if (n < 2) {
    stop("`Y` must contain at least two observations.");
  }
  if (d < 1) {
    stop("`Y` must have at least one embedding dimension.");
  }
  if (P.nrow() != n || P.ncol() != n) {
    stop("`P` must be a square matrix with one row and column per row of `Y`.");
  }

  double off_diagonal_mass = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < d; ++k) {
      if (!std::isfinite(Y(i, k))) {
        stop("`Y` must contain only finite values.");
      }
    }
    for (int j = 0; j < n; ++j) {
      const double p = P(i, j);
      if (!std::isfinite(p)) {
        stop("`P` must contain only finite values.");
      }
      if (p < 0.0) {
        stop("`P` must contain only non-negative probabilities.");
      }
      if (i != j) {
        off_diagonal_mass += p;
      }
    }
  }
  if (!std::isfinite(off_diagonal_mass) || off_diagonal_mass <= 0.0) {
    stop("`P` must have positive off-diagonal probability mass.");
  }
  
  NumericMatrix dist_matrix_Y(n, n);
  NumericMatrix Q(n, n);
  NumericMatrix grad_Y(n, d);
  
  // Compute pairwise Euclidean distances in the low-dimensional space
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dist = 0;
      for (int k = 0; k < d; ++k) {
        const double delta = Y(i, k) - Y(j, k);
        dist += delta * delta;
      }
      if (!std::isfinite(dist)) {
        stop("Pairwise squared distances in `Y` must be finite.");
      }
      dist_matrix_Y(i, j) = dist;
      dist_matrix_Y(j, i) = dist;
      Q(i, j) = 1 / (1 + dist);
      Q(j, i) = Q(i, j);
    }
  }
  
  // Normalize Q
  double sum_Q = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i != j) sum_Q += Q(i, j);
    }
  }
  if (!std::isfinite(sum_Q) || sum_Q <= 0.0) {
    stop("The low-dimensional probability normalizer must be finite and positive.");
  }
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      Q(i, j) /= sum_Q;
    }
  }
  
  // Compute the gradient
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i != j) {
        double PQ_diff = P(i, j) - Q(i, j);
        for (int k = 0; k < d; ++k) {
          grad_Y(i, k) += 4 * PQ_diff * (Y(i, k) - Y(j, k)) / (1 + dist_matrix_Y(i, j));
        }
      }
    }
  }
  
  return grad_Y;
}
