
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix tsne_cost_gradient_cpp(NumericMatrix Y, NumericMatrix P) {
  int n = Y.nrow();
  int d = Y.ncol();
  
  NumericMatrix dist_matrix_Y(n, n);
  NumericMatrix Q(n, n);
  NumericMatrix grad_Y(n, d);
  
  // Compute pairwise Euclidean distances in the low-dimensional space
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double dist = 0;
      for (int k = 0; k < d; ++k) {
        dist += pow(Y(i, k) - Y(j, k), 2);
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
