// oob_forest_wt.cpp — Compute OOB forest weight matrix from membership + inbag
//
// For sample i, only trees where inbag(i,t) == 0 contribute.
// In those trees, all samples j sharing the same leaf as i get weight
// proportional to their bootstrap frequency (inbag(j,t)), normalized
// by the total bootstrap mass in the leaf.
//
// This mirrors the inbag forest weight logic in mv_forest.cpp but restricts
// accumulation to OOB trees for each row sample.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_oob_forest_wt_cpp(IntegerMatrix membership,
                                         IntegerMatrix inbag) {
  int n = membership.nrow();
  int ntree = membership.ncol();

  // Output: n x n weight matrix + per-row OOB tree count
  std::vector<double> W(n * n, 0.0);
  std::vector<double> denom(n, 0.0);

  for (int t = 0; t < ntree; t++) {
    // Build leaf -> sample list for this tree
    // Also compute bootstrap leaf sizes
    std::unordered_map<int, std::vector<int>> leaf_samples;
    for (int i = 0; i < n; i++) {
      leaf_samples[membership(i, t)].push_back(i);
    }

    // For each leaf, accumulate OOB forest weights
    for (auto& kv : leaf_samples) {
      const std::vector<int>& group = kv.second;
      int g = (int)group.size();

      // Compute bootstrap leaf size (sum of inbag frequencies)
      double boot_leaf_size = 0.0;
      for (int a = 0; a < g; a++) {
        int freq = inbag(group[a], t);
        if (freq > 0) boot_leaf_size += (double)freq;
      }
      if (boot_leaf_size < 1.0) continue;
      double inv_bls = 1.0 / boot_leaf_size;

      // For OOB samples in this leaf: accumulate weights
      for (int a = 0; a < g; a++) {
        int ia = group[a];
        if (inbag(ia, t) != 0) continue;  // skip inbag samples

        denom[ia] += 1.0;
        double* row = &W[ia * n];
        for (int b = 0; b < g; b++) {
          int jb = group[b];
          int freq = inbag(jb, t);
          if (freq > 0) {
            row[jb] += (double)freq * inv_bls;
          }
        }
      }
    }
  }

  // Normalize and copy to R matrix
  NumericMatrix out(n, n);
  for (int i = 0; i < n; i++) {
    if (denom[i] > 0.0) {
      double inv_d = 1.0 / denom[i];
      for (int j = 0; j < n; j++) {
        out(i, j) = W[i * n + j] * inv_d;
      }
    }
  }

  // Copy dimnames from membership rownames
  if (membership.hasAttribute("dimnames")) {
    List dn = membership.attr("dimnames");
    if (!Rf_isNull(dn[0])) {
      CharacterVector rn = dn[0];
      out.attr("dimnames") = List::create(rn, clone(rn));
    }
  }

  return out;
}
