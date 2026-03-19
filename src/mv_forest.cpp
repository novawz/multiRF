// mv_forest.cpp — Multivariate regression random forest engine for multiRF
//
// Implements a lightweight multivariate forest that outputs:
//   - forest_wt   : n x n  (forest weight matrix)
//   - proximity   : n x n  (terminal-node co-occurrence)
//
// OpenMP is used to parallelize tree building across threads.
//   - membership  : n x ntree (terminal node IDs)
//   - tree_info   : per-tree structure for downstream IMD
//
// Split criterion: multivariate between-group sum of squares on Y
//   score = sum_j [ (sum_L_j)^2/nL + (sum_R_j)^2/nR ]
//   where Y is column-centered within the node.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ──────────────── Node structure ────────────────
struct Node {
  int id;
  int left;       // child index (-1 = leaf)
  int right;
  int split_var;  // column index in X (-1 = leaf)
  double split_val;
  int depth;
  std::vector<int> samples;  // row indices in this node
  // IMD: per-Y-column split contribution at this node
  // imd_y_stats[j] = (sum_L_j^2/nL + sum_R_j^2/nR) for full Y
  // Only populated for internal nodes; empty for leaves.
  std::vector<double> imd_y_stats;
  // IMD: which X variable was used to split (same as split_var but
  // stored here for clarity; the importance goes to this X var)
  double imd_x_score;  // total split score at this node
};

// Lightweight column-major matrix view backed by raw memory.
// This avoids constructing any R/Rcpp objects inside OpenMP regions.
struct MatrixView {
  const double* data;
  int nrow_;
  int ncol_;

  MatrixView(const double* data_, int nrow, int ncol)
    : data(data_), nrow_(nrow), ncol_(ncol) {}

  inline int nrow() const { return nrow_; }
  inline int ncol() const { return ncol_; }
  inline double operator()(int i, int j) const { return data[i + j * nrow_]; }
};

// ──────────────── Single tree ────────────────

// Find the best split for a node
// X: n x px, Y: n x qy (full matrices, use sample indices)
// mtry: number of candidate X vars, ytry: number of candidate Y vars
template <typename XMat, typename YMat>
static bool find_best_split(
    const XMat& X, const YMat& Y,
    const std::vector<int>& samples,
    int mtry, int ytry, int nodesize_min,
    std::mt19937& rng,
    int& best_var, double& best_val, double& best_score,
    std::vector<int>& left_samples, std::vector<int>& right_samples,
    std::vector<double>& best_y_stats)  // IMD: per-Y split stats for best split
{
  int n_node = (int)samples.size();
  int px = X.ncol();
  int qy = Y.ncol();

  if (n_node < 2 * nodesize_min) return false;

  // Random subset of X columns (mtry)
  std::vector<int> x_candidates(px);
  std::iota(x_candidates.begin(), x_candidates.end(), 0);
  std::shuffle(x_candidates.begin(), x_candidates.end(), rng);
  int n_x_try = std::min(mtry, px);

  // Random subset of Y columns (ytry)
  std::vector<int> y_candidates(qy);
  std::iota(y_candidates.begin(), y_candidates.end(), 0);
  std::shuffle(y_candidates.begin(), y_candidates.end(), rng);
  int n_y_try = std::min(ytry, qy);

  // Pre-compute column means for selected Y within this node
  std::vector<double> y_means(n_y_try, 0.0);
  for (int jj = 0; jj < n_y_try; jj++) {
    int j = y_candidates[jj];
    for (int k = 0; k < n_node; k++) {
      y_means[jj] += Y(samples[k], j);
    }
    y_means[jj] /= n_node;
  }

  best_score = -1.0;
  best_var = -1;
  best_val = 0.0;

  for (int xi = 0; xi < n_x_try; xi++) {
    int xvar = x_candidates[xi];

    // Get X values and sort
    std::vector<std::pair<double, int>> x_sorted(n_node);
    for (int k = 0; k < n_node; k++) {
      x_sorted[k] = {X(samples[k], xvar), k};
    }
    std::sort(x_sorted.begin(), x_sorted.end());

    // Running sums for left child (per Y column)
    std::vector<double> sum_L(n_y_try, 0.0);
    int nL = 0;

    // Try splits between consecutive sorted X values
    for (int s = 0; s < n_node - 1; s++) {
      int sample_idx = x_sorted[s].second;
      nL++;
      int nR = n_node - nL;

      // Update left sums
      for (int jj = 0; jj < n_y_try; jj++) {
        int j = y_candidates[jj];
        double val = Y(samples[sample_idx], j) - y_means[jj];
        sum_L[jj] += val;
      }

      // Skip if same X value as next
      if (x_sorted[s].first == x_sorted[s + 1].first) continue;

      // Skip if children too small
      if (nL < nodesize_min || nR < nodesize_min) continue;

      // Compute score: sum_j [ sum_L_j^2/nL + sum_R_j^2/nR ]
      // sum_R_j = total_j - sum_L_j, but total_j = 0 (centered)
      // so sum_R_j = -sum_L_j
      double score = 0.0;
      for (int jj = 0; jj < n_y_try; jj++) {
        double sL = sum_L[jj];
        score += (sL * sL) / nL + (sL * sL) / nR;
      }

      if (score > best_score) {
        best_score = score;
        // Record per-Y split stats for IMD (only for the ytry columns tried)
        best_y_stats.assign(qy, 0.0);
        for (int jj = 0; jj < n_y_try; jj++) {
          double sL = sum_L[jj];
          best_y_stats[y_candidates[jj]] = (sL * sL) / nL + (sL * sL) / nR;
        }
        best_var = xvar;
        best_val = (x_sorted[s].first + x_sorted[s + 1].first) / 2.0;
      }
    }
  }

  if (best_var < 0) return false;

  // Partition samples
  left_samples.clear();
  right_samples.clear();
  left_samples.reserve(n_node);
  right_samples.reserve(n_node);
  for (int k = 0; k < n_node; k++) {
    if (X(samples[k], best_var) <= best_val) {
      left_samples.push_back(samples[k]);
    } else {
      right_samples.push_back(samples[k]);
    }
  }

  return !left_samples.empty() && !right_samples.empty();
}

// Build a single tree
// Returns vector of Nodes
template <typename XMat, typename YMat>
static std::vector<Node> build_tree(
    const XMat& X, const YMat& Y,
    const std::vector<int>& bag_samples,
    int mtry, int ytry, int nodesize_min, int max_depth,
    std::mt19937& rng)
{
  std::vector<Node> nodes;
  nodes.reserve(256);

  // Root node
  Node root;
  root.id = 0;
  root.left = -1;
  root.right = -1;
  root.split_var = -1;
  root.split_val = 0.0;
  root.depth = 0;
  root.samples = bag_samples;
  nodes.push_back(root);

  // BFS-style tree building
  std::vector<int> to_split = {0};

  while (!to_split.empty()) {
    std::vector<int> next_split;

    for (int ni : to_split) {
      Node& node = nodes[ni];

      if ((int)node.samples.size() < 2 * nodesize_min) continue;
      if (max_depth > 0 && node.depth >= max_depth) continue;

      int bv;
      double bval, bscore;
      std::vector<int> lsamp, rsamp;
      std::vector<double> by_stats;

      if (!find_best_split(X, Y, node.samples, mtry, ytry, nodesize_min,
                           rng, bv, bval, bscore, lsamp, rsamp, by_stats)) {
        continue;
      }

      node.split_var = bv;
      node.split_val = bval;
      node.imd_y_stats = std::move(by_stats);
      node.imd_x_score = bscore;

      int left_id = (int)nodes.size();
      Node left_node;
      left_node.id = left_id;
      left_node.left = -1;
      left_node.right = -1;
      left_node.split_var = -1;
      left_node.split_val = 0.0;
      left_node.depth = node.depth + 1;
      left_node.samples = std::move(lsamp);
      nodes.push_back(left_node);

      int right_id = (int)nodes.size();
      Node right_node;
      right_node.id = right_id;
      right_node.left = -1;
      right_node.right = -1;
      right_node.split_var = -1;
      right_node.split_val = 0.0;
      right_node.depth = node.depth + 1;
      right_node.samples = std::move(rsamp);
      nodes.push_back(right_node);

      // Update parent pointers (re-reference since vector may reallocate)
      nodes[ni].left = left_id;
      nodes[ni].right = right_id;

      next_split.push_back(left_id);
      next_split.push_back(right_id);
    }

    to_split = std::move(next_split);
  }

  return nodes;
}

// Predict terminal node for a sample
template <typename XMat>
static int predict_leaf(const std::vector<Node>& tree,
                        const XMat& X, int sample_idx) {
  int node_idx = 0;
  while (tree[node_idx].split_var >= 0) {
    if (X(sample_idx, tree[node_idx].split_var) <= tree[node_idx].split_val) {
      node_idx = tree[node_idx].left;
    } else {
      node_idx = tree[node_idx].right;
    }
  }
  return node_idx;
}

// ──────────────── Forest + matrix accumulation ────────────────

// [[Rcpp::export]]
List fit_mv_forest_cpp(NumericMatrix X, NumericMatrix Y,
                       int ntree = 500,
                       int mtry = 0,
                       int ytry = 0,
                       int nodesize_min = 5,
                       int max_depth = 0,
                       int seed = -1,
                       int nthread = 0) {

  int n = X.nrow();
  int px = X.ncol();
  int qy = Y.ncol();

  // Defaults
  if (mtry <= 0) mtry = std::max(1, (int)std::sqrt((double)px));
  if (ytry <= 0) ytry = std::max(1, (int)std::ceil(std::sqrt((double)qy)));
  if (max_depth <= 0) max_depth = 150;

  #ifdef _OPENMP
  if (nthread > 0) omp_set_num_threads(nthread);
  #endif

  // Seed
  unsigned int actual_seed = (seed < 0) ? std::random_device{}() : (unsigned int)seed;

  // Output: plain C++ buffers for thread-safe accumulation
  std::vector<double> fw_buf(n * n, 0.0);
  std::vector<double> prox_buf(n * n, 0.0);

  // Per-tree results: tree structure + leaf membership + inbag mask
  struct TreeResult {
    std::vector<Node> nodes;
    std::vector<int> leaf_ids;
    std::vector<int> inbag;   // 0/1 per sample
  };
  std::vector<TreeResult> tree_results(ntree);

  // Parallel tree building with thread-local accumulation buffers.
  #ifdef _OPENMP
  #pragma omp parallel
  {
    std::vector<double> fw_local(n * n, 0.0);
    std::vector<double> prox_local(n * n, 0.0);
  #pragma omp for schedule(dynamic)
  #endif
  for (int t = 0; t < ntree; t++) {
    // Per-thread RNG (deterministic: seed + tree index)
    std::mt19937 rng_t(actual_seed + (unsigned int)t);
    std::uniform_int_distribution<int> boot_dist(0, n - 1);

    // Bootstrap sample WITH replacement (keep duplicates, like rfsrc)
    std::vector<int> bag(n);
    std::vector<int> inbag_vec(n, 0);
    for (int i = 0; i < n; i++) {
      int idx = boot_dist(rng_t);
      bag[i] = idx;
      inbag_vec[idx] = 1;
    }
    tree_results[t].inbag = std::move(inbag_vec);
    // Keep duplicates — a sample drawn twice has more influence on the split
    std::sort(bag.begin(), bag.end());

    // Build tree
    tree_results[t].nodes = build_tree(X, Y, bag, mtry, ytry,
                                       nodesize_min, max_depth, rng_t);

    // Predict leaf for ALL samples
    tree_results[t].leaf_ids.resize(n);
    for (int i = 0; i < n; i++) {
      tree_results[t].leaf_ids[i] = predict_leaf(tree_results[t].nodes, X, i);
    }

    // Accumulate proximity and forest weights into thread-local buffers.
    std::unordered_map<int, std::vector<int>> leaf_groups;
    for (int i = 0; i < n; i++) {
      leaf_groups[tree_results[t].leaf_ids[i]].push_back(i);
    }

    for (auto& kv : leaf_groups) {
      const std::vector<int>& group = kv.second;
      int g = (int)group.size();
      if (g == 0) continue;
      double wt = 1.0 / g;
      for (int a = 0; a < g; a++) {
        int ia = group[a];
        // Diagonal
        fw_local[ia * n + ia] += wt;
        prox_local[ia * n + ia] += 1.0;
        // Off-diagonal
        for (int b = a + 1; b < g; b++) {
          int ib = group[b];
          prox_local[ia * n + ib] += 1.0;
          prox_local[ib * n + ia] += 1.0;
          fw_local[ia * n + ib] += wt;
          fw_local[ib * n + ia] += wt;
        }
      }
    }
  }
  #ifdef _OPENMP
  #pragma omp critical
    {
      for (std::size_t idx = 0; idx < fw_buf.size(); ++idx) {
        fw_buf[idx] += fw_local[idx];
        prox_buf[idx] += prox_local[idx];
      }
    }
  }
  #endif

  // Copy to R matrices + normalize
  NumericMatrix forest_wt(n, n);
  NumericMatrix prox(n, n);
  IntegerMatrix membership(n, ntree);
  IntegerMatrix inbag_mat(n, ntree);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      forest_wt(i, j) = fw_buf[i * n + j] / ntree;
      prox(i, j) = prox_buf[i * n + j] / ntree;
    }
    prox(i, i) = 1.0;
  }

  // ──── IMD aggregation ────
  // Global accumulators
  std::vector<double> imd_x_sum(px, 0.0);
  std::vector<int>    imd_x_cnt(px, 0);
  std::vector<double> imd_y_sum(qy, 0.0);
  std::vector<int>    imd_y_cnt(qy, 0);

  // Per-tree IMD: ntree x (px + qy) stored as R matrix (for method="test")
  NumericMatrix imd_x_per_tree(px, ntree);
  NumericMatrix imd_y_per_tree(qy, ntree);

  // Pairwise X-Y co-occurrence matrix (px x qy) for pairwise_imd
  // pairwise_xy[x_var, y_var] = sum of y_stats at splits using x_var
  std::vector<double> pairwise_buf(px * qy, 0.0);

  for (int t = 0; t < ntree; t++) {
    // Per-tree accumulators
    std::vector<double> tx_sum(px, 0.0);
    std::vector<int>    tx_cnt(px, 0);
    std::vector<double> ty_sum(qy, 0.0);
    std::vector<int>    ty_cnt(qy, 0);

    for (const auto& node : tree_results[t].nodes) {
      if (node.split_var < 0) continue;
      int xv = node.split_var;
      // X importance
      tx_sum[xv] += node.imd_x_score;
      tx_cnt[xv]++;
      imd_x_sum[xv] += node.imd_x_score;
      imd_x_cnt[xv]++;
      // Y importance + pairwise
      for (int j = 0; j < (int)node.imd_y_stats.size() && j < qy; j++) {
        if (node.imd_y_stats[j] > 0.0) {
          ty_sum[j] += node.imd_y_stats[j];
          ty_cnt[j]++;
          imd_y_sum[j] += node.imd_y_stats[j];
          imd_y_cnt[j]++;
          // Pairwise: this split used X var xv and Y var j contributed
          pairwise_buf[xv * qy + j] += node.imd_y_stats[j];
        }
      }
    }

    // Store per-tree averages
    for (int j = 0; j < px; j++) {
      imd_x_per_tree(j, t) = (tx_cnt[j] > 0) ? tx_sum[j] / tx_cnt[j] : 0.0;
    }
    for (int j = 0; j < qy; j++) {
      imd_y_per_tree(j, t) = (ty_cnt[j] > 0) ? ty_sum[j] / ty_cnt[j] : 0.0;
    }
  }

  // Normalize pairwise matrix by ntree
  NumericMatrix pairwise_xy(px, qy);
  for (int i = 0; i < px; i++) {
    for (int j = 0; j < qy; j++) {
      pairwise_xy(i, j) = pairwise_buf[i * qy + j] / ntree;
    }
  }

  // Global averages + L2 normalize
  NumericVector imd_x(px);
  NumericVector imd_y(qy);
  for (int j = 0; j < px; j++) {
    imd_x[j] = (imd_x_cnt[j] > 0) ? imd_x_sum[j] / imd_x_cnt[j] : 0.0;
  }
  for (int j = 0; j < qy; j++) {
    imd_y[j] = (imd_y_cnt[j] > 0) ? imd_y_sum[j] / imd_y_cnt[j] : 0.0;
  }
  double norm_x = 0.0, norm_y = 0.0;
  for (int j = 0; j < px; j++) norm_x += imd_x[j] * imd_x[j];
  for (int j = 0; j < qy; j++) norm_y += imd_y[j] * imd_y[j];
  norm_x = std::sqrt(norm_x);
  norm_y = std::sqrt(norm_y);
  if (norm_x > 0) for (int j = 0; j < px; j++) imd_x[j] /= norm_x;
  if (norm_y > 0) for (int j = 0; j < qy; j++) imd_y[j] /= norm_y;

  // Free heavy per-node data no longer needed (samples + imd_y_stats)
  // This can reclaim ~80 MB for 500 trees.
  for (int t = 0; t < ntree; t++) {
    for (auto& node : tree_results[t].nodes) {
      std::vector<int>().swap(node.samples);
      std::vector<double>().swap(node.imd_y_stats);
    }
  }

  // Fill membership + inbag matrices, build tree_info R list
  List tree_info_list(ntree);
  for (int t = 0; t < ntree; t++) {
    for (int i = 0; i < n; i++) {
      membership(i, t) = tree_results[t].leaf_ids[i];
      inbag_mat(i, t) = tree_results[t].inbag[i];
    }

    int n_nodes = (int)tree_results[t].nodes.size();
    IntegerVector t_split_var(n_nodes);
    NumericVector t_split_val(n_nodes);
    IntegerVector t_left(n_nodes);
    IntegerVector t_right(n_nodes);
    IntegerVector t_depth(n_nodes);
    IntegerVector t_nodesize(n_nodes);
    LogicalVector t_is_leaf(n_nodes);

    for (int ni = 0; ni < n_nodes; ni++) {
      t_split_var[ni] = tree_results[t].nodes[ni].split_var;
      t_split_val[ni] = tree_results[t].nodes[ni].split_val;
      t_left[ni] = tree_results[t].nodes[ni].left;
      t_right[ni] = tree_results[t].nodes[ni].right;
      t_depth[ni] = tree_results[t].nodes[ni].depth;
      t_nodesize[ni] = (int)tree_results[t].nodes[ni].samples.size();
      t_is_leaf[ni] = (tree_results[t].nodes[ni].split_var < 0);
    }

    tree_info_list[t] = List::create(
      Named("split_var") = t_split_var,
      Named("split_val") = t_split_val,
      Named("left") = t_left,
      Named("right") = t_right,
      Named("depth") = t_depth,
      Named("nodesize") = t_nodesize,
      Named("is_leaf") = t_is_leaf
    );
  }

  return List::create(
    Named("forest.wt") = forest_wt,
    Named("proximity") = prox,
    Named("membership") = membership,
    Named("inbag") = inbag_mat,
    Named("imd_x") = imd_x,
    Named("imd_y") = imd_y,
    Named("imd_x_per_tree") = imd_x_per_tree,
    Named("imd_y_per_tree") = imd_y_per_tree,
    Named("pairwise_xy") = pairwise_xy,
    Named("ntree") = ntree,
    Named("n") = n,
    Named("px") = px,
    Named("qy") = qy,
    Named("tree_info") = tree_info_list
  );
}


// ──────────────── Unsupervised forest ────────────────
// Each tree randomly splits columns into pseudo-X and pseudo-Y,
// then fits a multivariate regression tree.  This emulates
// randomForestSRC's unsupervised mode entirely in C++.

// [[Rcpp::export]]
List fit_mv_forest_unsup_cpp(NumericMatrix data,
                              int ntree = 500,
                              int ytry = 10,
                              int nodesize_min = 5,
                              int max_depth = 0,
                              int seed = -1,
                              int nthread = 0) {

  int n = data.nrow();
  int p = data.ncol();

  if (max_depth <= 0) max_depth = 150;

  #ifdef _OPENMP
  if (nthread > 0) omp_set_num_threads(nthread);
  #endif

  unsigned int actual_seed = (seed < 0) ? std::random_device{}() : (unsigned int)seed;

  std::vector<double> fw_buf(n * n, 0.0);
  std::vector<double> prox_buf(n * n, 0.0);

  struct UnsupTreeResult {
    std::vector<Node> nodes;
    std::vector<int> leaf_ids;
    std::vector<int> col_perm;  // column permutation for this tree
    int n_y;
    int n_x;
  };
  std::vector<UnsupTreeResult> tree_results(ntree);

  #ifdef _OPENMP
  #pragma omp parallel
  {
    std::vector<double> fw_local(n * n, 0.0);
    std::vector<double> prox_local(n * n, 0.0);
  #pragma omp for schedule(dynamic)
  #endif
  for (int t = 0; t < ntree; t++) {
    std::mt19937 rng_t(actual_seed + (unsigned int)t);
    std::uniform_int_distribution<int> boot_dist(0, n - 1);

    // Random column split
    std::vector<int> col_perm(p);
    std::iota(col_perm.begin(), col_perm.end(), 0);
    std::shuffle(col_perm.begin(), col_perm.end(), rng_t);
    int n_y = std::max(1, p / 2);
    int n_x = p - n_y;
    if (n_x == 0) { n_x = 1; n_y = p - 1; }

    // Build sub-matrices (thread-local, no R allocations)
    // Use raw vectors to avoid Rcpp thread-safety issues
    std::vector<double> X_buf(n * n_x), Y_buf(n * n_y);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n_x; j++) X_buf[i + j * n] = data(i, col_perm[n_y + j]);
      for (int j = 0; j < n_y; j++) Y_buf[i + j * n] = data(i, col_perm[j]);
    }

    // Pure C++ matrix views; avoid creating Rcpp objects inside OpenMP threads.
    MatrixView X_sub(X_buf.data(), n, n_x);
    MatrixView Y_sub(Y_buf.data(), n, n_y);

    int mtry_t = std::max(1, (int)std::sqrt((double)n_x));
    int ytry_t = std::min(ytry, n_y);

    // Bootstrap WITH replacement (keep duplicates)
    std::vector<int> bag(n);
    for (int i = 0; i < n; i++) bag[i] = boot_dist(rng_t);
    std::sort(bag.begin(), bag.end());

    std::vector<Node> tree = build_tree(X_sub, Y_sub, bag,
                                        mtry_t, ytry_t, nodesize_min,
                                        max_depth, rng_t);

    std::vector<int> leaf_ids(n);
    for (int i = 0; i < n; i++) {
      leaf_ids[i] = predict_leaf(tree, X_sub, i);
    }

    // Accumulate into thread-local buffers.
    std::unordered_map<int, std::vector<int>> leaf_groups;
    for (int i = 0; i < n; i++) leaf_groups[leaf_ids[i]].push_back(i);

    for (auto& kv : leaf_groups) {
      const std::vector<int>& group = kv.second;
      int g = (int)group.size();
      if (g == 0) continue;
      double wt = 1.0 / g;
      for (int a = 0; a < g; a++) {
        int ia = group[a];
        // Diagonal
        fw_local[ia * n + ia] += wt;
        prox_local[ia * n + ia] += 1.0;
        // Off-diagonal
        for (int b = a + 1; b < g; b++) {
          int ib = group[b];
          prox_local[ia * n + ib] += 1.0;
          prox_local[ib * n + ia] += 1.0;
          fw_local[ia * n + ib] += wt;
          fw_local[ib * n + ia] += wt;
        }
      }
    }

    tree_results[t].nodes = std::move(tree);
    tree_results[t].leaf_ids = std::move(leaf_ids);
    tree_results[t].col_perm = std::move(col_perm);
    tree_results[t].n_y = n_y;
    tree_results[t].n_x = n_x;
  }
  #ifdef _OPENMP
  #pragma omp critical
    {
      for (std::size_t idx = 0; idx < fw_buf.size(); ++idx) {
        fw_buf[idx] += fw_local[idx];
        prox_buf[idx] += prox_local[idx];
      }
    }
  }
  #endif

  // Copy to R matrices + normalize
  NumericMatrix forest_wt(n, n);
  NumericMatrix prox(n, n);
  IntegerMatrix membership(n, ntree);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      forest_wt(i, j) = fw_buf[i * n + j] / ntree;
      prox(i, j) = prox_buf[i * n + j] / ntree;
    }
    prox(i, i) = 1.0;
  }

  List tree_info_list(ntree);
  for (int t = 0; t < ntree; t++) {
    for (int i = 0; i < n; i++) {
      membership(i, t) = tree_results[t].leaf_ids[i];
    }

    int n_nodes = (int)tree_results[t].nodes.size();
    int n_y = tree_results[t].n_y;
    int n_x = tree_results[t].n_x;
    const auto& col_perm = tree_results[t].col_perm;

    IntegerVector t_split_var(n_nodes);
    NumericVector t_split_val(n_nodes);
    IntegerVector t_left(n_nodes);
    IntegerVector t_right(n_nodes);
    IntegerVector t_depth(n_nodes);
    IntegerVector t_nodesize(n_nodes);
    LogicalVector t_is_leaf(n_nodes);

    // Remap split_var from local X_sub index to global column index
    for (int ni = 0; ni < n_nodes; ni++) {
      int sv = tree_results[t].nodes[ni].split_var;
      t_split_var[ni] = (sv >= 0) ? col_perm[n_y + sv] : -1;
      t_split_val[ni] = tree_results[t].nodes[ni].split_val;
      t_left[ni] = tree_results[t].nodes[ni].left;
      t_right[ni] = tree_results[t].nodes[ni].right;
      t_depth[ni] = tree_results[t].nodes[ni].depth;
      t_nodesize[ni] = (int)tree_results[t].nodes[ni].samples.size();
      t_is_leaf[ni] = (sv < 0);
    }

    IntegerVector y_cols(n_y);
    IntegerVector x_cols(n_x);
    for (int j = 0; j < n_y; j++) y_cols[j] = col_perm[j] + 1;
    for (int j = 0; j < n_x; j++) x_cols[j] = col_perm[n_y + j] + 1;

    tree_info_list[t] = List::create(
      Named("split_var") = t_split_var,
      Named("split_val") = t_split_val,
      Named("left") = t_left,
      Named("right") = t_right,
      Named("depth") = t_depth,
      Named("nodesize") = t_nodesize,
      Named("is_leaf") = t_is_leaf,
      Named("y_cols") = y_cols,
      Named("x_cols") = x_cols
    );
  }

  return List::create(
    Named("forest.wt") = forest_wt,
    Named("proximity") = prox,
    Named("membership") = membership,
    Named("ntree") = ntree,
    Named("n") = n,
    Named("p") = p,
    Named("tree_info") = tree_info_list
  );
}
