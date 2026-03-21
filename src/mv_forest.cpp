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
// Split criterion: normalized multivariate between-group sum of squares on Y
//   score = sum_j [ ((sum_L_j)^2/nL + (sum_R_j)^2/nR) / var_j ]
//   where Y is column-centered within the node and var_j is the
//   within-node variance of Y column j (matching randomForestSRC's
//   standardized composite splitting rule).

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <array>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

static constexpr double RF_EPSILON = 1.0e-9;

// ──────────────── Node structure ────────────────
struct Node {
  int id;
  int left;       // child index (-1 = leaf)
  int right;
  int split_var;  // column index in X (-1 = leaf)
  double split_val;
  int depth;
  int nodesize;
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

// rfsrc-style SWOR sampling:
// repeatedly draw one remaining position and remove it by swap-with-last.
// This matches the algorithmic structure in randomForestSRC/bootstrap.c
// more closely than std::shuffle + prefix slice.
struct RfsrcRan1;

template <typename RNG>
static double random_unit(RNG& rng) {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(rng);
}

template <typename RNG>
static std::vector<int> sample_swor_rfsrc_style(int n, int sample_size, RNG& rng) {
  sample_size = std::max(1, std::min(sample_size, n));
  std::vector<int> pool(n);
  std::iota(pool.begin(), pool.end(), 0);
  std::vector<int> out(sample_size);
  int remaining = n;
  for (int i = 0; i < sample_size; i++) {
    int pick = std::max(0, std::min(remaining - 1,
      (int)std::ceil(random_unit(rng) * remaining) - 1));
    out[i] = pool[pick];
    pool[pick] = pool[remaining - 1];
    remaining--;
  }
  return out;
}

struct RfsrcRan1 {
  static constexpr int IA = 16807;
  static constexpr int IM = 2147483647;
  static constexpr int IQ = 127773;
  static constexpr int IR = 2836;
  static constexpr int NTAB = 32;
  static constexpr int NDIV = 1 + (IM - 1) / NTAB;
  static constexpr double AM = 1.0 / IM;
  static constexpr double EPS = 1.2e-7;
  static constexpr double RNMX = 1.0 - EPS;

  int iy = 0;
  std::array<int, NTAB> iv{};
  int seed = -1;

  explicit RfsrcRan1(int init_seed = -1) : seed(init_seed) {}

  double next() {
    int j, k;
    double temp;
    if (seed <= 0 || iy == 0) {
      seed = (-seed < 1) ? 1 : -seed;
      for (j = NTAB + 7; j >= 0; j--) {
        k = seed / IQ;
        seed = IA * (seed - k * IQ) - IR * k;
        if (seed < 0) seed += IM;
        if (j < NTAB) iv[j] = seed;
      }
      iy = iv[0];
    }
    k = seed / IQ;
    seed = IA * (seed - k * IQ) - IR * k;
    if (seed < 0) seed += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = seed;
    temp = AM * iy;
    return (temp > RNMX) ? RNMX : temp;
  }
};

static double random_unit(RfsrcRan1& rng) {
  return rng.next();
}

static unsigned int lcg_next(unsigned int seed, bool reset) {
  constexpr unsigned int LCG_IM = 714025;
  constexpr unsigned int LCG_IA = 1366;
  constexpr unsigned int LCG_IC = 150889;
  if (reset) {
    return (seed >= LCG_IM) ? (seed % LCG_IM) : seed;
  }
  return (LCG_IA * seed + LCG_IC) % LCG_IM;
}

template <typename RNG>
static std::vector<int> sample_from_pool_rfsrc_style(const std::vector<int>& pool_in,
                                                     int sample_size,
                                                     RNG& rng) {
  if (pool_in.empty()) return {};
  sample_size = std::max(1, std::min(sample_size, (int)pool_in.size()));
  std::vector<int> pool = pool_in;
  std::vector<int> out;
  out.reserve(sample_size);
  int remaining = (int)pool.size();
  for (int i = 0; i < sample_size; i++) {
    int pick = std::max(0, std::min(remaining - 1,
      (int)std::ceil(random_unit(rng) * remaining) - 1));
    out.push_back(pool[pick]);
    pool[pick] = pool[remaining - 1];
    remaining--;
  }
  return out;
}

template <typename RNG>
static std::vector<int> sample_cutpoint_positions(const std::vector<int>& candidates,
                                                  int nsplit,
                                                  RNG& rng) {
  if (nsplit <= 0 || (int)candidates.size() <= nsplit) {
    return candidates;
  }
  std::vector<int> out = sample_from_pool_rfsrc_style(candidates, nsplit, rng);
  std::sort(out.begin(), out.end());
  return out;
}

// ──────────────── Pre-sorted split (fast path) ────────────────

// Pre-sort all X columns. Returns sort_order[j] = sample indices sorted by X[,j].
template <typename XMat>
static std::vector<std::vector<int>> presort_columns(
    const XMat& X, int n, int px)
{
  std::vector<std::vector<int>> sort_order(px);
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int j = 0; j < px; j++) {
    sort_order[j].resize(n);
    std::iota(sort_order[j].begin(), sort_order[j].end(), 0);
    std::sort(sort_order[j].begin(), sort_order[j].end(),
              [&](int a, int b){ return X(a, j) < X(b, j); });
  }
  return sort_order;
}

// ──────────────── Partition-based split (ranger-style) ────────────────
// Each node holds its own sorted index subsets (node_sorted[j] has only
// n_node entries, already in sorted order for variable j).  No in_node
// scanning needed — every entry belongs to this node.

template <typename XMat, typename YMat, typename RNG>
static bool find_best_split_part(
    const XMat& X, const YMat& Y,
    const std::vector<int>& samples,                    // samples in this node
    const std::vector<std::vector<int>>& node_sorted,   // [px][n_node] pre-sorted per var
    std::vector<int>& sample_pos,                       // reusable scratch [n_total]
    int mtry, int ytry, int nodesize_min, int nsplit,
    RNG& rng,
    int& best_var, double& best_val, double& best_score,
    std::vector<int>& left_samples, std::vector<int>& right_samples,
    std::vector<double>& best_y_stats)
{
  int n_node = (int)samples.size();
  int px = X.ncol();
  int qy = Y.ncol();

  if (n_node < 2 * nodesize_min) return false;

  // Random subset of X columns (mtry)
  std::vector<int> x_pool(px);
  std::iota(x_pool.begin(), x_pool.end(), 0);
  int n_x_try = std::min(mtry, px);
  std::vector<int> x_candidates = sample_from_pool_rfsrc_style(x_pool, n_x_try, rng);

  // Random subset of Y columns (ytry)
  std::vector<int> y_pool(qy);
  std::iota(y_pool.begin(), y_pool.end(), 0);
  int n_y_try = std::min(ytry, qy);
  std::vector<int> y_candidates = sample_from_pool_rfsrc_style(y_pool, n_y_try, rng);

  // Pre-standardize selected Y columns within this node
  std::vector<double> y_means(n_y_try, 0.0);
  std::vector<double> y_sds(n_y_try, 1.0);
  for (int jj = 0; jj < n_y_try; jj++) {
    int j = y_candidates[jj];
    double sum = 0.0;
    for (int k = 0; k < n_node; k++) sum += Y(samples[k], j);
    y_means[jj] = sum / n_node;
  }
  for (int jj = 0; jj < n_y_try; jj++) {
    int j = y_candidates[jj];
    double ss = 0.0;
    for (int k = 0; k < n_node; k++) {
      double d = Y(samples[k], j) - y_means[jj];
      ss += d * d;
    }
    double var = (n_node > 1) ? ss / n_node : 1.0;
    y_sds[jj] = (var > 1e-6) ? std::sqrt(var) : 0.0;
  }

  best_score = -1.0;
  best_var = -1;
  best_val = 0.0;

  // Pre-compute standardized Y; build sample_id -> local index map
  std::vector<double> y_std_flat(n_y_try * n_node);
  for (int jj = 0; jj < n_y_try; jj++) {
    int j = y_candidates[jj];
    double inv_sd = (y_sds[jj] > 0.0) ? 1.0 / y_sds[jj] : 0.0;
    for (int k = 0; k < n_node; k++) {
      y_std_flat[jj * n_node + k] =
        (Y(samples[k], j) - y_means[jj]) * inv_sd;
    }
  }
  for (int k = 0; k < n_node; k++) sample_pos[samples[k]] = k;

  for (int xi = 0; xi < n_x_try; xi++) {
    int xvar = x_candidates[xi];
    const std::vector<int>& order = node_sorted[xvar];
    // order has exactly n_node entries, all in this node — no skipping
    // Note: nodesize is enforced only at the parent level (n_node >= 2*nodesize_min)
    // to match rfsrc, which checks leftSize>0 && rghtSize>0 at the cutpoint level.
    // A child smaller than nodesize simply becomes a terminal leaf.
    std::vector<int> split_positions;
    split_positions.reserve(n_node);
    for (int s = 1; s < n_node; s++) {
      if (X(order[s - 1], xvar) == X(order[s], xvar)) continue;
      split_positions.push_back(s);
    }
    if (split_positions.empty()) continue;
    std::vector<int> eval_positions = sample_cutpoint_positions(split_positions, nsplit, rng);
    int eval_idx = 0;
    int next_eval = eval_positions[eval_idx];

    std::vector<double> sum_L(n_y_try, 0.0);
    int nL = 0;

    for (int s = 0; s < n_node; s++) {
      int si = order[s];
      double x_val = X(si, xvar);
      int local_k = sample_pos[si];

      // Check split before adding this sample to left
      if (eval_idx < (int)eval_positions.size() && s == next_eval) {
        double prev_x = X(order[s - 1], xvar);
        double score = 0.0;
        int deltaNorm = 0;
        for (int jj = 0; jj < n_y_try; jj++) {
          if (y_sds[jj] > 1e-6) {
            double sL = sum_L[jj];
            score += (sL * sL) / nL + (sL * sL) / (n_node - nL);
            deltaNorm++;
          }
        }
        if (deltaNorm > 0) {
          score /= deltaNorm;

          if ((score - best_score) > RF_EPSILON) {
            best_score = score;
            best_y_stats.assign(qy, 0.0);
            for (int jj = 0; jj < n_y_try; jj++) {
              double sL = sum_L[jj];
              best_y_stats[y_candidates[jj]] = (sL * sL) / nL + (sL * sL) / (n_node - nL);
            }
            best_var = xvar;
            best_val = prev_x;
          }
        }
        // Always advance eval_idx even when deltaNorm == 0;
        // the old `else continue` would skip nL++ and sum_L update,
        // corrupting all subsequent cutpoint scores for this variable.
        eval_idx++;
        if (eval_idx < (int)eval_positions.size()) {
          next_eval = eval_positions[eval_idx];
        } else {
          next_eval = n_node;
        }
      }

      nL++;
      for (int jj = 0; jj < n_y_try; jj++) {
        sum_L[jj] += y_std_flat[jj * n_node + local_k];
      }
    }
  }

  // Clean up sample_pos
  for (int k = 0; k < n_node; k++) sample_pos[samples[k]] = -1;

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

// Build tree with partition-based sorted indices (ranger-style).
// Each node owns its sorted index subsets; children get partitioned copies.
template <typename XMat, typename YMat, typename RNG>
static std::vector<Node> build_tree_part(
    const XMat& X, const YMat& Y,
    const std::vector<int>& bag_samples,
    const std::vector<std::vector<int>>& sort_order,  // global pre-sorted [px][n]
    int n_total, int px,
    int mtry, int ytry, int nodesize_min, int max_depth, int nsplit,
    RNG& rng)
{
  std::vector<Node> nodes;
  nodes.reserve(256);

  std::vector<int> sample_pos(n_total, -1);  // reusable scratch
  std::vector<char> left_flag(n_total, 0);   // reusable scratch for partitioning

  // BFS task: node id + its per-variable sorted indices
  struct SplitTask {
    int node_id;
    std::vector<std::vector<int>> sorted;  // [px][n_node]
  };

  // Create root node
  Node root;
  root.id = 0;
  root.left = -1;
  root.right = -1;
  root.split_var = -1;
  root.split_val = 0.0;
  root.depth = 0;
  root.samples = bag_samples;
  root.nodesize = (int)bag_samples.size();
  nodes.push_back(root);

  // Build root's sorted indices: filter global sort_order to bag samples
  // Mark bag membership
  std::vector<char> in_bag(n_total, 0);
  for (int si : bag_samples) in_bag[si] = 1;

  SplitTask root_task;
  root_task.node_id = 0;
  root_task.sorted.resize(px);
  for (int j = 0; j < px; j++) {
    root_task.sorted[j].reserve(bag_samples.size());
    for (int si : sort_order[j]) {
      if (in_bag[si]) root_task.sorted[j].push_back(si);
    }
  }

  std::vector<SplitTask> to_split;
  to_split.push_back(std::move(root_task));

  while (!to_split.empty()) {
    std::vector<SplitTask> next_split;

    for (auto& task : to_split) {
      Node& node = nodes[task.node_id];

      if ((int)node.samples.size() < 2 * nodesize_min) continue;
      if (max_depth > 0 && node.depth >= max_depth) continue;

      int bv;
      double bval, bscore;
      std::vector<int> lsamp, rsamp;
      std::vector<double> by_stats;

      bool found = find_best_split_part(
        X, Y, node.samples, task.sorted, sample_pos,
        mtry, ytry, nodesize_min, nsplit, rng,
        bv, bval, bscore, lsamp, rsamp, by_stats);

      if (!found) continue;

      node.split_var = bv;
      node.split_val = bval;
      node.imd_y_stats = std::move(by_stats);
      node.imd_x_score = bscore;

      // Mark left samples for partitioning
      for (int si : lsamp) left_flag[si] = 1;

      // Partition each variable's sorted indices into left/right
      std::vector<std::vector<int>> left_sorted(px), right_sorted(px);
      for (int j = 0; j < px; j++) {
        left_sorted[j].reserve(lsamp.size());
        right_sorted[j].reserve(rsamp.size());
        for (int si : task.sorted[j]) {
          if (left_flag[si]) left_sorted[j].push_back(si);
          else right_sorted[j].push_back(si);
        }
      }

      // Clear left_flag
      for (int si : lsamp) left_flag[si] = 0;

      // Free parent's sorted indices (no longer needed)
      task.sorted.clear();
      task.sorted.shrink_to_fit();

      int left_id = (int)nodes.size();
      Node left_node;
      left_node.id = left_id;
      left_node.left = -1;
      left_node.right = -1;
      left_node.split_var = -1;
      left_node.split_val = 0.0;
      left_node.depth = node.depth + 1;
      left_node.samples = std::move(lsamp);
      left_node.nodesize = (int)left_node.samples.size();
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
      right_node.nodesize = (int)right_node.samples.size();
      nodes.push_back(right_node);

      nodes[task.node_id].left = left_id;
      nodes[task.node_id].right = right_id;

      // Queue children with their partitioned sorted indices
      SplitTask left_task;
      left_task.node_id = left_id;
      left_task.sorted = std::move(left_sorted);
      next_split.push_back(std::move(left_task));

      SplitTask right_task;
      right_task.node_id = right_id;
      right_task.sorted = std::move(right_sorted);
      next_split.push_back(std::move(right_task));
    }

    to_split = std::move(next_split);
  }

  return nodes;
}

// ──────────────── Global-scan split (kept for fallback) ────────────────

// Fast split using pre-sorted indices.
// in_node[i] = true if sample i belongs to current node.
// sort_order[j] = global sorted indices for X column j.
template <typename XMat, typename YMat, typename RNG>
static bool find_best_split_fast(
    const XMat& X, const YMat& Y,
    const std::vector<int>& samples,          // samples in this node
    const std::vector<char>& in_node,         // n-length flag
    const std::vector<std::vector<int>>& sort_order,
    std::vector<int>& sample_pos,             // reusable n-length scratch buffer
    int mtry, int ytry, int nodesize_min, int nsplit,
    RNG& rng,
    int& best_var, double& best_val, double& best_score,
    std::vector<int>& left_samples, std::vector<int>& right_samples,
    std::vector<double>& best_y_stats)
{
  int n_node = (int)samples.size();
  int n_total = (int)sample_pos.size();
  int px = X.ncol();
  int qy = Y.ncol();

  if (n_node < 2 * nodesize_min) return false;

  // Random subset of X columns (mtry)
  std::vector<int> x_pool(px);
  std::iota(x_pool.begin(), x_pool.end(), 0);
  int n_x_try = std::min(mtry, px);
  std::vector<int> x_candidates = sample_from_pool_rfsrc_style(x_pool, n_x_try, rng);

  // Random subset of Y columns (ytry)
  std::vector<int> y_pool(qy);
  std::iota(y_pool.begin(), y_pool.end(), 0);
  int n_y_try = std::min(ytry, qy);
  std::vector<int> y_candidates = sample_from_pool_rfsrc_style(y_pool, n_y_try, rng);

  // Pre-standardize selected Y columns within this node
  std::vector<double> y_means(n_y_try, 0.0);
  std::vector<double> y_sds(n_y_try, 1.0);
  for (int jj = 0; jj < n_y_try; jj++) {
    int j = y_candidates[jj];
    double sum = 0.0;
    for (int k = 0; k < n_node; k++) sum += Y(samples[k], j);
    y_means[jj] = sum / n_node;
  }
  for (int jj = 0; jj < n_y_try; jj++) {
    int j = y_candidates[jj];
    double ss = 0.0;
    for (int k = 0; k < n_node; k++) {
      double d = Y(samples[k], j) - y_means[jj];
      ss += d * d;
    }
    double var = (n_node > 1) ? ss / n_node : 1.0;
    y_sds[jj] = (var > 1e-6) ? std::sqrt(var) : 0.0;
  }

  best_score = -1.0;
  best_var = -1;
  best_val = 0.0;

  // Pre-compute standardized Y for node samples (avoid repeated computation)
  // y_std[jj][k] = standardized Y for sample samples[k], Y-candidate jj
  // We need to map sample_id -> position for fast lookup during scanning
  std::vector<double> y_std_flat(n_y_try * n_node);
  for (int jj = 0; jj < n_y_try; jj++) {
    int j = y_candidates[jj];
    double inv_sd = (y_sds[jj] > 0.0) ? 1.0 / y_sds[jj] : 0.0;
    for (int k = 0; k < n_node; k++) {
      y_std_flat[jj * n_node + k] =
        (Y(samples[k], j) - y_means[jj]) * inv_sd;
    }
  }
  // Set sample_id -> local index map (buffer passed in, cleared on exit)
  for (int k = 0; k < n_node; k++) sample_pos[samples[k]] = k;

  for (int xi = 0; xi < n_x_try; xi++) {
    int xvar = x_candidates[xi];
    const std::vector<int>& order = sort_order[xvar];

    std::vector<int> split_positions;
    split_positions.reserve(n_node);
    int nL_scan = 0;
    double prev_scan = 0.0;
    bool first_scan = true;
    for (int s = 0; s < n_total; s++) {
      int si = order[s];
      if (!in_node[si]) continue;
      double x_val = X(si, xvar);
      // nodesize enforced at parent level only; match rfsrc cutpoint check (>0)
      if (!first_scan && x_val != prev_scan) {
        split_positions.push_back(s);
      }
      nL_scan++;
      prev_scan = x_val;
      first_scan = false;
    }
    if (split_positions.empty()) continue;
    std::vector<int> eval_positions = sample_cutpoint_positions(split_positions, nsplit, rng);
    int eval_idx = 0;
    int next_eval = eval_positions[eval_idx];

    // Scan pre-sorted indices; skip samples not in this node
    std::vector<double> sum_L(n_y_try, 0.0);
    int nL = 0;
    double prev_x = 0.0;
    bool first = true;

    for (int s = 0; s < n_total; s++) {
      int si = order[s];
      if (!in_node[si]) continue;

      double x_val = X(si, xvar);
      int local_k = sample_pos[si];

      if (eval_idx < (int)eval_positions.size() && s == next_eval) {
        // Evaluate split between prev_x and x_val
        int nR = n_node - nL;
        double score = 0.0;
        int deltaNorm = 0;
        for (int jj = 0; jj < n_y_try; jj++) {
          if (y_sds[jj] > 1e-6) {
            double sL = sum_L[jj];
            score += (sL * sL) / nL + (sL * sL) / nR;
            deltaNorm++;
          }
        }
        if (deltaNorm > 0) score /= deltaNorm;

        if (deltaNorm > 0 && (score - best_score) > RF_EPSILON) {
          best_score = score;
          best_y_stats.assign(qy, 0.0);
          for (int jj = 0; jj < n_y_try; jj++) {
            double sL = sum_L[jj];
            best_y_stats[y_candidates[jj]] = (sL * sL) / nL + (sL * sL) / nR;
          }
          best_var = xvar;
          best_val = prev_x;
        }
        eval_idx++;
        if (eval_idx < (int)eval_positions.size()) {
          next_eval = eval_positions[eval_idx];
        } else {
          next_eval = n_total;
        }
      }

      // Add this sample to left child
      nL++;
      for (int jj = 0; jj < n_y_try; jj++) {
        sum_L[jj] += y_std_flat[jj * n_node + local_k];
      }
      prev_x = x_val;
      first = false;
    }
  }

  // Clean up sample_pos
  for (int k = 0; k < n_node; k++) sample_pos[samples[k]] = -1;

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

// Build tree using pre-sorted indices (fast path)
template <typename XMat, typename YMat, typename RNG>
static std::vector<Node> build_tree_fast(
    const XMat& X, const YMat& Y,
    const std::vector<int>& bag_samples,
    const std::vector<std::vector<int>>& sort_order,
    int n_total,
    int mtry, int ytry, int nodesize_min, int max_depth, int nsplit,
    RNG& rng)
{
  std::vector<Node> nodes;
  nodes.reserve(256);

  // Reusable scratch buffers (allocated once, cleared per-node)
  std::vector<char> in_node(n_total, 0);
  std::vector<int> sample_pos(n_total, -1);

  Node root;
  root.id = 0;
  root.left = -1;
  root.right = -1;
  root.split_var = -1;
  root.split_val = 0.0;
  root.depth = 0;
  root.samples = bag_samples;
  root.nodesize = (int)bag_samples.size();
  nodes.push_back(root);

  std::vector<int> to_split = {0};

  while (!to_split.empty()) {
    std::vector<int> next_split;

    for (int ni : to_split) {
      Node& node = nodes[ni];

      if ((int)node.samples.size() < 2 * nodesize_min) continue;
      if (max_depth > 0 && node.depth >= max_depth) continue;

      // Set in_node flags for current node
      for (int si : node.samples) in_node[si] = true;

      int bv;
      double bval, bscore;
      std::vector<int> lsamp, rsamp;
      std::vector<double> by_stats;

      bool found = find_best_split_fast(
        X, Y, node.samples, in_node, sort_order, sample_pos,
        mtry, ytry, nodesize_min, nsplit, rng,
        bv, bval, bscore, lsamp, rsamp, by_stats);

      // Clear in_node flags
      for (int si : node.samples) in_node[si] = false;

      if (!found) continue;

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
      left_node.nodesize = (int)left_node.samples.size();
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
      right_node.nodesize = (int)right_node.samples.size();
      nodes.push_back(right_node);

      nodes[ni].left = left_id;
      nodes[ni].right = right_id;

      next_split.push_back(left_id);
      next_split.push_back(right_id);
    }

    to_split = std::move(next_split);
  }

  return nodes;
}

// ──────────────── Original split (kept for reference) ────────────────

// Find the best split for a node
// X: n x px, Y: n x qy (full matrices, use sample indices)
// mtry: number of candidate X vars, ytry: number of candidate Y vars
template <typename XMat, typename YMat, typename RNG>
static bool find_best_split(
    const XMat& X, const YMat& Y,
    const std::vector<int>& samples,
    int mtry, int ytry, int nodesize_min, int nsplit,
    RNG& rng,
    int& best_var, double& best_val, double& best_score,
    std::vector<int>& left_samples, std::vector<int>& right_samples,
    std::vector<double>& best_y_stats)  // IMD: per-Y split stats for best split
{
  int n_node = (int)samples.size();
  int px = X.ncol();
  int qy = Y.ncol();

  if (n_node < 2 * nodesize_min) return false;

  // Random subset of X columns (mtry)
  std::vector<int> x_pool(px);
  std::iota(x_pool.begin(), x_pool.end(), 0);
  int n_x_try = std::min(mtry, px);
  std::vector<int> x_candidates = sample_from_pool_rfsrc_style(x_pool, n_x_try, rng);

  // Random subset of Y columns (ytry)
  std::vector<int> y_pool(qy);
  std::iota(y_pool.begin(), y_pool.end(), 0);
  int n_y_try = std::min(ytry, qy);
  std::vector<int> y_candidates = sample_from_pool_rfsrc_style(y_pool, n_y_try, rng);

  // Pre-standardize selected Y columns within this node: Y* = (Y - mean) / sd
  // This matches rfsrc's normalized composite splitting rule.
  // n_node here is the bootstrap bag size (including duplicates).
  std::vector<double> y_means(n_y_try, 0.0);
  std::vector<double> y_sds(n_y_try, 1.0);
  for (int jj = 0; jj < n_y_try; jj++) {
    int j = y_candidates[jj];
    double sum = 0.0;
    for (int k = 0; k < n_node; k++) {
      sum += Y(samples[k], j);
    }
    y_means[jj] = sum / n_node;
  }
  for (int jj = 0; jj < n_y_try; jj++) {
    int j = y_candidates[jj];
    double ss = 0.0;
    for (int k = 0; k < n_node; k++) {
      double d = Y(samples[k], j) - y_means[jj];
      ss += d * d;
    }
    // Standard deviation with n denominator (population sd within node)
    double var = (n_node > 1) ? ss / n_node : 1.0;
    y_sds[jj] = (var > 1e-6) ? std::sqrt(var) : 0.0;
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

    // nodesize enforced at parent level only; match rfsrc cutpoint check (>0)
    std::vector<int> split_positions;
    split_positions.reserve(n_node);
    for (int s = 0; s < n_node - 1; s++) {
      if (x_sorted[s].first == x_sorted[s + 1].first) continue;
      split_positions.push_back(s);
    }
    if (split_positions.empty()) continue;
    std::vector<int> eval_positions = sample_cutpoint_positions(split_positions, nsplit, rng);
    int eval_idx = 0;
    int next_eval = eval_positions[eval_idx];

    // Running sums for left child (per Y column)
    std::vector<double> sum_L(n_y_try, 0.0);
    int nL = 0;

    // Try splits between consecutive sorted X values
    for (int s = 0; s < n_node - 1; s++) {
      int sample_idx = x_sorted[s].second;
      nL++;
      int nR = n_node - nL;

      // Update left sums of standardized Y*
      for (int jj = 0; jj < n_y_try; jj++) {
        int j = y_candidates[jj];
        double y_star = (y_sds[jj] > 0.0) ? (Y(samples[sample_idx], j) - y_means[jj]) / y_sds[jj] : 0.0;
        sum_L[jj] += y_star;
      }

      if (s != next_eval) continue;

      // Score = average over informative Y columns (matching rfsrc deltaNorm)
      double score = 0.0;
      int deltaNorm = 0;
      for (int jj = 0; jj < n_y_try; jj++) {
        if (y_sds[jj] > 1e-6) {
          double sL = sum_L[jj];
          score += (sL * sL) / nL + (sL * sL) / nR;
          deltaNorm++;
        }
      }
      if (deltaNorm > 0) {
        score /= deltaNorm;

        if ((score - best_score) > RF_EPSILON) {
          best_score = score;
          // Record per-Y split stats for IMD (only for the ytry columns tried)
          // Y* is already standardized, so no separate variance division needed
          best_y_stats.assign(qy, 0.0);
          for (int jj = 0; jj < n_y_try; jj++) {
            double sL = sum_L[jj];
            best_y_stats[y_candidates[jj]] = (sL * sL) / nL + (sL * sL) / nR;
          }
          best_var = xvar;
          best_val = x_sorted[s].first;
        }
      }
      // Always advance eval_idx even when deltaNorm == 0;
      // the old `else continue` would freeze next_eval, causing all
      // remaining cutpoints for this X variable to be skipped.
      eval_idx++;
      if (eval_idx < (int)eval_positions.size()) {
        next_eval = eval_positions[eval_idx];
      } else {
        next_eval = n_node;
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

// ──────────────── Unsupervised split ────────────────
// Matches rfsrc's unsupervised split: X = Y = same data matrix.
// For each candidate split variable, exclude it from pseudo-Y,
// then randomly sample ytry pseudo-responses from the remaining columns.
template <typename Mat>
static bool find_best_split_unsup(
    const Mat& D,                     // n x p data (used as both X and Y)
    const std::vector<int>& samples,
    const std::vector<char>& in_node,
    const std::vector<std::vector<int>>& sort_order,
    std::vector<int>& sample_pos,     // reusable n-length scratch buffer
    int n_total,
    int mtry, int ytry, int nodesize_min,
    std::mt19937& rng,
    int& best_var, double& best_val, double& best_score,
    std::vector<int>& left_samples, std::vector<int>& right_samples)
{
  int n_node = (int)samples.size();
  int p = D.ncol();

  if (n_node < 2 * nodesize_min) return false;

  // Random subset of candidate split variables (mtry)
  std::vector<int> x_candidates(p);
  std::iota(x_candidates.begin(), x_candidates.end(), 0);
  std::shuffle(x_candidates.begin(), x_candidates.end(), rng);
  int n_x_try = std::min(mtry, p);

  best_score = -1.0;
  best_var = -1;
  best_val = 0.0;

  // Set sample_id -> local index map (buffer passed in, cleared on exit)
  for (int k = 0; k < n_node; k++) sample_pos[samples[k]] = k;

  for (int xi = 0; xi < n_x_try; xi++) {
    int xvar = x_candidates[xi];

    // Select ytry pseudo-responses from columns OTHER than xvar
    std::vector<int> y_pool;
    y_pool.reserve(p - 1);
    for (int j = 0; j < p; j++) {
      if (j != xvar) y_pool.push_back(j);
    }
    std::shuffle(y_pool.begin(), y_pool.end(), rng);
    int n_y_try = std::min(ytry, (int)y_pool.size());

    // Pre-standardize selected pseudo-Y columns within this node
    std::vector<double> y_means(n_y_try, 0.0);
    std::vector<double> y_sds(n_y_try, 1.0);
    for (int jj = 0; jj < n_y_try; jj++) {
      int j = y_pool[jj];
      double sum = 0.0;
      for (int k = 0; k < n_node; k++) sum += D(samples[k], j);
      y_means[jj] = sum / n_node;
    }
    for (int jj = 0; jj < n_y_try; jj++) {
      int j = y_pool[jj];
      double ss = 0.0;
      for (int k = 0; k < n_node; k++) {
        double d = D(samples[k], j) - y_means[jj];
        ss += d * d;
      }
      double var = (n_node > 1) ? ss / n_node : 1.0;
      y_sds[jj] = (var > 1e-6) ? std::sqrt(var) : 0.0;
    }

    bool any_informative = false;
    for (int jj = 0; jj < n_y_try; jj++) {
      if (y_sds[jj] > 1e-6) { any_informative = true; break; }
    }
    if (!any_informative) continue;

    // Pre-compute standardized pseudo-Y for node samples
    std::vector<double> y_std_flat(n_y_try * n_node);
    for (int jj = 0; jj < n_y_try; jj++) {
      int j = y_pool[jj];
      double inv_sd = (y_sds[jj] > 0.0) ? 1.0 / y_sds[jj] : 0.0;
      for (int k = 0; k < n_node; k++) {
        y_std_flat[jj * n_node + k] = (D(samples[k], j) - y_means[jj]) * inv_sd;
      }
    }

    // Scan pre-sorted indices for this split variable
    const std::vector<int>& order = sort_order[xvar];
    std::vector<double> sum_L(n_y_try, 0.0);
    int nL = 0;
    double prev_x = 0.0;
    bool first = true;

    for (int s = 0; s < n_total; s++) {
      int si = order[s];
      if (!in_node[si]) continue;

      double x_val = D(si, xvar);
      int local_k = sample_pos[si];

      // nodesize enforced at parent level only; match rfsrc cutpoint check (>0)
      if (!first && x_val != prev_x) {
        int nR = n_node - nL;
        double score = 0.0;
        int deltaNorm = 0;
        for (int jj = 0; jj < n_y_try; jj++) {
          if (y_sds[jj] > 1e-6) {
            double sL = sum_L[jj];
            score += (sL * sL) / nL + (sL * sL) / nR;
            deltaNorm++;
          }
        }
        if (deltaNorm > 0) score /= deltaNorm;

        if (deltaNorm > 0 && (score - best_score) > RF_EPSILON) {
          best_score = score;
          best_var = xvar;
          best_val = (prev_x + x_val) / 2.0;
        }
      }

      nL++;
      for (int jj = 0; jj < n_y_try; jj++) {
        sum_L[jj] += y_std_flat[jj * n_node + local_k];
      }
      prev_x = x_val;
      first = false;
    }
  }

  // Clean up sample_pos
  for (int k = 0; k < n_node; k++) sample_pos[samples[k]] = -1;

  if (best_var < 0) return false;

  left_samples.clear();
  right_samples.clear();
  left_samples.reserve(n_node);
  right_samples.reserve(n_node);
  for (int k = 0; k < n_node; k++) {
    if (D(samples[k], best_var) <= best_val) {
      left_samples.push_back(samples[k]);
    } else {
      right_samples.push_back(samples[k]);
    }
  }

  return !left_samples.empty() && !right_samples.empty();
}

// Build an unsupervised tree using pre-sorted indices
template <typename Mat>
static std::vector<Node> build_tree_unsup(
    const Mat& D,
    const std::vector<int>& bag_samples,
    const std::vector<std::vector<int>>& sort_order,
    int n_total,
    int mtry, int ytry, int nodesize_min, int max_depth,
    std::mt19937& rng)
{
  std::vector<Node> nodes;
  nodes.reserve(256);

  std::vector<char> in_node(n_total, 0);
  std::vector<int> sample_pos(n_total, -1);

  Node root;
  root.id = 0;
  root.left = -1;
  root.right = -1;
  root.split_var = -1;
  root.split_val = 0.0;
  root.depth = 0;
  root.samples = bag_samples;
  root.nodesize = (int)bag_samples.size();
  nodes.push_back(root);

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

      // Set in_node flags for current node
      for (int si : node.samples) in_node[si] = true;

      if (!find_best_split_unsup(D, node.samples, in_node, sort_order, sample_pos, n_total,
                                  mtry, ytry, nodesize_min,
                                  rng, bv, bval, bscore, lsamp, rsamp)) {
        for (int si : node.samples) in_node[si] = false;
        continue;
      }

      // Clear in_node flags after split
      for (int si : node.samples) in_node[si] = false;

      node.split_var = bv;
      node.split_val = bval;

      int left_id = (int)nodes.size();
      Node left_node;
      left_node.id = left_id;
      left_node.left = -1;
      left_node.right = -1;
      left_node.split_var = -1;
      left_node.split_val = 0.0;
      left_node.depth = node.depth + 1;
      left_node.samples = std::move(lsamp);
      left_node.nodesize = (int)left_node.samples.size();
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
      right_node.nodesize = (int)right_node.samples.size();
      nodes.push_back(right_node);

      nodes[ni].left = left_id;
      nodes[ni].right = right_id;

      next_split.push_back(left_id);
      next_split.push_back(right_id);
    }

    to_split = std::move(next_split);
  }

  return nodes;
}

// Build a single tree (supervised)
// Returns vector of Nodes
template <typename XMat, typename YMat, typename RNG>
static std::vector<Node> build_tree(
    const XMat& X, const YMat& Y,
    const std::vector<int>& bag_samples,
    int mtry, int ytry, int nodesize_min, int max_depth, int nsplit,
    RNG& rng)
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
  root.nodesize = (int)bag_samples.size();
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
                           nsplit, rng, bv, bval, bscore, lsamp, rsamp, by_stats)) {
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
      left_node.nodesize = (int)left_node.samples.size();
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
      right_node.nodesize = (int)right_node.samples.size();
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

// ──────────────── Spearman correlation helper ────────────────

// Spearman correlation between two vectors of length d.
// For small d (≈10) this is trivially fast.
static double spearman_corr(const double* a, const double* b, int d) {
  if (d < 2) return 0.0;

  // Rank vectors using simple insertion-sort index (d is tiny)
  std::vector<int> ord_a(d), ord_b(d);
  std::iota(ord_a.begin(), ord_a.end(), 0);
  std::iota(ord_b.begin(), ord_b.end(), 0);
  std::sort(ord_a.begin(), ord_a.end(), [&](int i, int j){ return a[i] < a[j]; });
  std::sort(ord_b.begin(), ord_b.end(), [&](int i, int j){ return b[i] < b[j]; });

  // Assign ranks (average rank for ties)
  auto assign_ranks = [&](const double* v, const std::vector<int>& ord, std::vector<double>& ranks) {
    ranks.resize(d);
    int i = 0;
    while (i < d) {
      int j = i + 1;
      while (j < d && v[ord[j]] == v[ord[i]]) j++;
      double avg_rank = 0.5 * (i + j - 1);  // 0-based average
      for (int k = i; k < j; k++) ranks[ord[k]] = avg_rank;
      i = j;
    }
  };

  std::vector<double> ra, rb;
  assign_ranks(a, ord_a, ra);
  assign_ranks(b, ord_b, rb);

  // Pearson on ranks
  double mean_a = 0, mean_b = 0;
  for (int i = 0; i < d; i++) { mean_a += ra[i]; mean_b += rb[i]; }
  mean_a /= d; mean_b /= d;

  double cov = 0, va = 0, vb = 0;
  for (int i = 0; i < d; i++) {
    double da = ra[i] - mean_a, db = rb[i] - mean_b;
    cov += da * db;
    va += da * da;
    vb += db * db;
  }
  double denom = std::sqrt(va * vb);
  return (denom > 1e-12) ? cov / denom : 0.0;
}

// ──────────────── Forest + matrix accumulation ────────────────

// [[Rcpp::export]]
List fit_mv_forest_cpp(NumericMatrix X, NumericMatrix Y,
                       int ntree = 500,
                       int mtry = 0,
                       int ytry = 0,
                       int nsplit = 10,
                       int nodesize_min = 5,
                       int max_depth = 0,
                       int seed = -1,
                       int nthread = 0,
                       int samptype = 0,
                       int prox_mode = 0,
                       Nullable<NumericMatrix> embed = R_NilValue,
                       double sibling_gamma = 0.5,
                       int enhanced_prox_mode = 0) {
  // enhanced_prox_mode: 0 = off, 1 = compute enhanced proximity

  int n = X.nrow();
  int px = X.ncol();
  int qy = Y.ncol();

  // Defaults: ceiling(p/3) for mtry, min(qy, ceiling(p/3)) for ytry,
  // and nsplit = 10 to match randomForestSRC's default randomized cut search.
  if (mtry <= 0) mtry = std::max(1, (int)std::ceil((double)px / 3.0));
  if (ytry <= 0) ytry = std::min(qy, std::max(1, (int)std::ceil((double)px / 3.0)));
  if (nsplit < 0) nsplit = 10;
  // max_depth <= 0 means unlimited (grow until nodesize constraint only)

  #ifdef _OPENMP
  if (nthread > 0) omp_set_num_threads(nthread);
  #endif

  // Copy to thread-safe MatrixView (avoid Rcpp operator() inside OpenMP)
  std::vector<double> X_buf(n * px), Y_buf(n * qy);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < px; j++) X_buf[i + j * n] = X(i, j);
    for (int j = 0; j < qy; j++) Y_buf[i + j * n] = Y(i, j);
  }
  MatrixView Xv(X_buf.data(), n, px);
  MatrixView Yv(Y_buf.data(), n, qy);

  // Seed
  unsigned int actual_seed = (seed < 0) ? std::random_device{}() : (unsigned int)seed;
  unsigned int seed_lc = lcg_next(actual_seed, true);
  std::vector<int> chain_seed_a(ntree), chain_seed_b(ntree);
  for (int t = 0; t < ntree; t++) {
    do {
      seed_lc = lcg_next(seed_lc, false);
      seed_lc = lcg_next(seed_lc, false);
    } while (seed_lc == 0);
    chain_seed_a[t] = -(int)seed_lc;
  }
  for (int t = 0; t < ntree; t++) {
    do {
      seed_lc = lcg_next(seed_lc, false);
      seed_lc = lcg_next(seed_lc, false);
    } while (seed_lc == 0);
    chain_seed_b[t] = -(int)seed_lc;
  }

  // ──── Enhanced proximity setup ────
  bool compute_enhanced = (enhanced_prox_mode > 0) && embed.isNotNull();
  int embed_dim = 0;
  std::vector<double> embed_buf;
  if (compute_enhanced) {
    NumericMatrix embed_mat(embed.get());
    embed_dim = embed_mat.ncol();
    embed_buf.resize(n * embed_dim);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < embed_dim; j++)
        embed_buf[i + j * n] = embed_mat(i, j);
  }

  // Output: plain C++ buffers for thread-safe accumulation
  // prox_mode: -1 = skip proximity entirely, 0 = all, 1 = inbag, 2 = oob
  bool compute_prox = (prox_mode >= 0);
  std::vector<double> fw_buf(n * n, 0.0);
  std::vector<double> fw_denom_buf(n, 0.0);
  std::vector<double> prox_buf(compute_prox ? n * n : 0, 0.0);
  std::vector<double> prox_denom_buf(prox_mode > 0 ? n * n : 0, 0.0);
  std::vector<double> eprox_buf(compute_enhanced ? n * n : 0, 0.0);

  // Per-tree results: tree structure + leaf membership + bootstrap frequency
  struct TreeResult {
    std::vector<Node> nodes;
    std::vector<int> leaf_ids;
    std::vector<int> inbag;   // bootstrap frequency n_{i,b} per sample
  };
  std::vector<TreeResult> tree_results(ntree);

  // Pre-sort all X columns once (shared across trees, read-only)
  auto sort_order = presort_columns(Xv, n, px);

  // Parallel tree building with thread-local accumulation buffers.
  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<double> fw_local(n * n, 0.0);
    std::vector<double> fw_denom_local(n, 0.0);
    std::vector<double> prox_local(compute_prox ? n * n : 0, 0.0);
    std::vector<double> prox_denom_local(prox_mode > 0 ? n * n : 0, 0.0);
    std::vector<double> eprox_local(compute_enhanced ? n * n : 0, 0.0);
    // Scratch buffers for enhanced proximity centroid computation
    std::vector<double> centroid_scratch(compute_enhanced ? 2 * embed_dim : 0, 0.0);
    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
  for (int t = 0; t < ntree; t++) {
    RfsrcRan1 rng_boot(chain_seed_a[t]);
    RfsrcRan1 rng_split(chain_seed_b[t]);

    // Bootstrap sampling
    std::vector<int> bag;
    std::vector<int> inbag_freq(n, 0);

    if (samptype == 1) {
      // SWR: sample WITH replacement (n draws from {0, ..., n-1})
      // Matches rfsrc samptype="swr": same index can appear multiple times
      bag.resize(n);
      for (int i = 0; i < n; i++) {
        int idx = std::max(0, std::min(n - 1,
          (int)std::ceil(rng_boot.next() * n) - 1));
        bag[i] = idx;
        inbag_freq[idx]++;
      }
    } else {
      // SWOR: sample WITHOUT replacement (default, matches rfsrc samptype="swor")
      // Draw ~63.2% of samples (matching expected unique count of standard bootstrap)
      int samp_size = std::max(1, (int)std::round(n * (1.0 - std::exp(-1.0))));
      bag = sample_swor_rfsrc_style(n, samp_size, rng_boot);
      for (int i = 0; i < samp_size; i++) {
        inbag_freq[bag[i]] = 1;
      }
    }
    tree_results[t].inbag = inbag_freq;
    std::sort(bag.begin(), bag.end());

    // Build tree with partition-based sorted indices (ranger-style)
    tree_results[t].nodes = build_tree_part(Xv, Yv, bag, sort_order, n, px,
                                             mtry, ytry, nodesize_min,
                                             max_depth, nsplit, rng_split);

    // Release sample vectors from all nodes — they are only needed during
    // tree construction; leaf membership is obtained via predict_leaf.
    for (auto& node : tree_results[t].nodes) {
      node.samples.clear();
      node.samples.shrink_to_fit();
    }

    // Predict leaf for ALL samples (both IB and OOB)
    tree_results[t].leaf_ids.resize(n);
    for (int i = 0; i < n; i++) {
      tree_results[t].leaf_ids[i] = predict_leaf(tree_results[t].nodes, Xv, i);
    }

    // Accumulate forest weights and proximity.
    // Forest-weight semantics are currently inbag-style:
    //   - rows are all samples predicted through the tree
    //   - columns are inbag donors only
    //   - each row is normalized by the number of trees contributing to it
    std::unordered_map<int, std::vector<int>> leaf_groups;
    for (int i = 0; i < n; i++) {
      leaf_groups[tree_results[t].leaf_ids[i]].push_back(i);
    }

    for (auto& kv : leaf_groups) {
      const std::vector<int>& group = kv.second;
      int g = (int)group.size();
      if (g == 0) continue;

      // Pre-extract inbag donors and their weights in this leaf.
      // This avoids the branch `if (inbag_freq[ib] > 0)` in the inner loop.
      std::vector<int> ib_idx;
      std::vector<double> ib_wt;
      double boot_leaf_size = 0.0;
      ib_idx.reserve(g);
      ib_wt.reserve(g);
      for (int a = 0; a < g; a++) {
        int freq = inbag_freq[group[a]];
        if (freq > 0) {
          ib_idx.push_back(group[a]);
          ib_wt.push_back((double)freq);
          boot_leaf_size += freq;
        }
      }
      if (boot_leaf_size < 1.0) continue;

      // Normalize weights once
      double inv_bls = 1.0 / boot_leaf_size;
      int n_ib = (int)ib_idx.size();

      for (int a = 0; a < g; a++) {
        int ia = group[a];
        fw_denom_local[ia] += 1.0;
        double* row = &fw_local[ia * n];
        for (int b = 0; b < n_ib; b++) {
          row[ib_idx[b]] += ib_wt[b] * inv_bls;
        }
      }
    }

    // Proximity accumulation (skipped when prox_mode == -1).
    if (compute_prox) {
      if (prox_mode == 0) {
        // All samples — reuse leaf_groups; denom = ntree (handled at finalize)
        for (auto& kv : leaf_groups) {
          const std::vector<int>& group = kv.second;
          int g = (int)group.size();
          for (int a = 0; a < g; a++) {
            int ia = group[a];
            prox_local[ia * n + ia] += 1.0;
            for (int b = a + 1; b < g; b++) {
              int ib = group[b];
              prox_local[ia * n + ib] += 1.0;
              prox_local[ib * n + ia] += 1.0;
            }
          }
        }
      } else {
        // inbag or oob mode: subset of samples, need per-pair denom
        std::vector<int> prox_members;
        prox_members.reserve(n);
        if (prox_mode == 1) {
          for (int i = 0; i < n; i++) {
            if (inbag_freq[i] > 0) prox_members.push_back(i);
          }
        } else {
          for (int i = 0; i < n; i++) {
            if (inbag_freq[i] == 0) prox_members.push_back(i);
          }
        }

        std::unordered_map<int, std::vector<int>> prox_leaf_groups;
        for (int idx : prox_members) {
          prox_leaf_groups[tree_results[t].leaf_ids[idx]].push_back(idx);
        }
        for (int a = 0; a < (int)prox_members.size(); a++) {
          int ia = prox_members[a];
          prox_denom_local[ia * n + ia] += 1.0;
          for (int b = a + 1; b < (int)prox_members.size(); b++) {
            int ib = prox_members[b];
            prox_denom_local[ia * n + ib] += 1.0;
            prox_denom_local[ib * n + ia] += 1.0;
          }
        }
        for (auto& kv : prox_leaf_groups) {
          const std::vector<int>& group = kv.second;
          int g = (int)group.size();
          for (int a = 0; a < g; a++) {
            int ia = group[a];
            prox_local[ia * n + ia] += 1.0;
            for (int b = a + 1; b < g; b++) {
              int ib = group[b];
              prox_local[ia * n + ib] += 1.0;
              prox_local[ib * n + ia] += 1.0;
            }
          }
        }
      }
    }

    // ──── Enhanced proximity accumulation ────
    // For each sibling-leaf pair (two leaves sharing the same parent),
    // compute Spearman correlation of leaf centroids in embedding space.
    // Same-leaf pairs get weight 1; sibling-leaf pairs get gamma * max(corr, 0).
    if (compute_enhanced) {
      const auto& nodes = tree_results[t].nodes;
      const auto& lids = tree_results[t].leaf_ids;

      // 1. Same-leaf contribution (weight = 1)
      for (auto& kv : leaf_groups) {
        const std::vector<int>& group = kv.second;
        int g = (int)group.size();
        for (int a = 0; a < g; a++) {
          int ia = group[a];
          eprox_local[ia * n + ia] += 1.0;
          for (int b = a + 1; b < g; b++) {
            int ib = group[b];
            eprox_local[ia * n + ib] += 1.0;
            eprox_local[ib * n + ia] += 1.0;
          }
        }
      }

      // 2. Sibling-leaf contribution
      //    Find internal nodes whose left AND right children are both leaves.
      for (const auto& nd : nodes) {
        if (nd.split_var < 0) continue;  // leaf node, skip
        int li = nd.left, ri = nd.right;
        if (li < 0 || ri < 0) continue;
        if (nodes[li].split_var >= 0 || nodes[ri].split_var >= 0) continue;
        // Both children are leaves — this is a sibling pair.

        // Find samples in each sibling leaf via leaf_groups
        auto it_l = leaf_groups.find(li);
        auto it_r = leaf_groups.find(ri);
        if (it_l == leaf_groups.end() || it_r == leaf_groups.end()) continue;
        const std::vector<int>& grp_l = it_l->second;
        const std::vector<int>& grp_r = it_r->second;
        if (grp_l.empty() || grp_r.empty()) continue;

        // Compute centroids in embedding space
        double* cent_l = centroid_scratch.data();
        double* cent_r = centroid_scratch.data() + embed_dim;
        std::fill(cent_l, cent_l + embed_dim, 0.0);
        std::fill(cent_r, cent_r + embed_dim, 0.0);
        for (int si : grp_l)
          for (int d = 0; d < embed_dim; d++)
            cent_l[d] += embed_buf[si + d * n];
        for (int si : grp_r)
          for (int d = 0; d < embed_dim; d++)
            cent_r[d] += embed_buf[si + d * n];
        double inv_l = 1.0 / grp_l.size();
        double inv_r = 1.0 / grp_r.size();
        for (int d = 0; d < embed_dim; d++) { cent_l[d] *= inv_l; cent_r[d] *= inv_r; }

        // Spearman correlation of centroids
        double corr = spearman_corr(cent_l, cent_r, embed_dim);
        double w = sibling_gamma * std::max(corr, 0.0);
        if (w <= 0.0) continue;
        if (w > 1.0) w = 1.0;

        // Add sibling proximity
        for (int a : grp_l) {
          for (int b : grp_r) {
            eprox_local[a * n + b] += w;
            eprox_local[b * n + a] += w;
          }
        }
      }
    }
  }
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
      for (std::size_t i = 0; i < fw_denom_buf.size(); ++i) {
        fw_denom_buf[i] += fw_denom_local[i];
      }
      if (prox_mode > 0) {
        for (std::size_t idx = 0; idx < prox_denom_buf.size(); ++idx) {
          prox_denom_buf[idx] += prox_denom_local[idx];
        }
      }
      for (std::size_t idx = 0; idx < fw_buf.size(); ++idx) {
        fw_buf[idx] += fw_local[idx];
      }
      if (compute_prox) {
        for (std::size_t idx = 0; idx < prox_buf.size(); ++idx) {
          prox_buf[idx] += prox_local[idx];
        }
      }
      if (compute_enhanced) {
        for (std::size_t idx = 0; idx < eprox_buf.size(); ++idx) {
          eprox_buf[idx] += eprox_local[idx];
        }
      }
    }
  }

  // Copy to R matrices + normalize
  NumericMatrix forest_wt(n, n);
  NumericMatrix prox(n, n);
  IntegerMatrix membership(n, ntree);
  IntegerMatrix inbag_mat(n, ntree);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      forest_wt(i, j) = (fw_denom_buf[i] > 0.0) ? fw_buf[i * n + j] / fw_denom_buf[i] : NA_REAL;
    }
  }
  if (compute_prox) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (prox_mode == 0) {
          prox(i, j) = prox_buf[i * n + j] / ntree;
        } else {
          prox(i, j) = (prox_denom_buf[i * n + j] > 0.0) ? prox_buf[i * n + j] / prox_denom_buf[i * n + j] : NA_REAL;
        }
      }
    }
  }

  // ──── Enhanced proximity normalization ────
  NumericMatrix enhanced_prox(compute_enhanced ? n : 0, compute_enhanced ? n : 0);
  if (compute_enhanced) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        enhanced_prox(i, j) = eprox_buf[i * n + j] / ntree;
      }
    }
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
      t_nodesize[ni] = tree_results[t].nodes[ni].nodesize;
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
    Named("enhanced_prox") = enhanced_prox,
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
// Matches rfsrc unsupervised mode: at each split, the current covariate
// is excluded and ytry pseudo-responses are drawn from the remaining
// columns.  All p columns serve as both candidate split variables and
// candidate pseudo-responses.

// [[Rcpp::export]]
List fit_mv_forest_unsup_cpp(NumericMatrix data,
                              int ntree = 500,
                              int ytry = 10,
                              int nodesize_min = 5,
                              int max_depth = 0,
                              int seed = -1,
                              int nthread = 0,
                              int samptype = 0,
                              int prox_mode = 0) {  // 0 = all, 1 = inbag, 2 = oob

  int n = data.nrow();
  int p = data.ncol();

  // max_depth <= 0 means unlimited (grow until nodesize constraint only)

  #ifdef _OPENMP
  if (nthread > 0) omp_set_num_threads(nthread);
  #endif

  unsigned int actual_seed = (seed < 0) ? std::random_device{}() : (unsigned int)seed;

  bool compute_prox = (prox_mode >= 0);
  std::vector<double> fw_buf(n * n, 0.0);
  std::vector<double> fw_denom_buf(n, 0.0);
  std::vector<double> prox_buf(compute_prox ? n * n : 0, 0.0);
  std::vector<double> prox_denom_buf(prox_mode > 0 ? n * n : 0, 0.0);

  struct UnsupTreeResult {
    std::vector<Node> nodes;
    std::vector<int> leaf_ids;
  };
  std::vector<UnsupTreeResult> tree_results(ntree);

  // Copy data to a thread-safe MatrixView (avoid Rcpp inside OpenMP)
  std::vector<double> data_buf(n * p);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < p; j++)
      data_buf[i + j * n] = data(i, j);
  MatrixView D(data_buf.data(), n, p);

  // mtry default: ceiling(sqrt(p)), matching rfsrc unsupervised
  int mtry_default = std::max(1, (int)std::ceil(std::sqrt((double)p)));
  // ytry default: min(ytry_param, p-1); ytry <= 0 means use p-1
  int ytry_use = (ytry <= 0) ? std::max(1, p - 1) : std::min(ytry, p - 1);

  // Pre-sort all columns once (shared across trees, read-only)
  auto sort_order_unsup = presort_columns(D, n, p);

  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    std::vector<double> fw_local(n * n, 0.0);
    std::vector<double> fw_denom_local(n, 0.0);
    std::vector<double> prox_local(compute_prox ? n * n : 0, 0.0);
    std::vector<double> prox_denom_local(prox_mode > 0 ? n * n : 0, 0.0);
    #ifdef _OPENMP
    #pragma omp for schedule(dynamic)
    #endif
  for (int t = 0; t < ntree; t++) {
    std::mt19937 rng_t(actual_seed + (unsigned int)t);
    std::uniform_int_distribution<int> boot_dist(0, n - 1);

    // Bootstrap sampling (same logic as supervised engine)
    std::vector<int> bag;
    std::vector<int> inbag_freq(n, 0);

    if (samptype == 1) {
      // SWR: sample WITH replacement
      bag.resize(n);
      for (int i = 0; i < n; i++) {
        int idx = boot_dist(rng_t);
        bag[i] = idx;
        inbag_freq[idx]++;
      }
    } else {
      // SWOR: sample WITHOUT replacement (default)
      int samp_size_u = std::max(1, (int)std::round(n * (1.0 - std::exp(-1.0))));
      bag = sample_swor_rfsrc_style(n, samp_size_u, rng_t);
      for (int i = 0; i < samp_size_u; i++) {
        inbag_freq[bag[i]] = 1;
      }
    }
    std::sort(bag.begin(), bag.end());

    // Build unsupervised tree: all p columns as both X and Y,
    // pseudo-Y selected dynamically per split (excluding split var)
    std::vector<Node> tree = build_tree_unsup(D, bag, sort_order_unsup, n,
                                               mtry_default, ytry_use,
                                               nodesize_min, max_depth, rng_t);

    // Release sample vectors (only needed during tree construction)
    for (auto& node : tree) {
      node.samples.clear();
      node.samples.shrink_to_fit();
    }

    std::vector<int> leaf_ids(n);
    for (int i = 0; i < n; i++) {
      leaf_ids[i] = predict_leaf(tree, D, i);
    }

    // Accumulate forest weights and proximity with the same inbag-style
    // semantics used by the supervised engine: rows are all samples,
    // columns are inbag donors, and rows are normalized by their
    // effective tree count.
    std::unordered_map<int, std::vector<int>> leaf_groups;
    for (int i = 0; i < n; i++) leaf_groups[leaf_ids[i]].push_back(i);

    for (auto& kv : leaf_groups) {
      const std::vector<int>& group = kv.second;
      int g = (int)group.size();
      if (g == 0) continue;

      // Pre-extract inbag donors (same optimization as supervised engine)
      std::vector<int> ib_idx;
      std::vector<double> ib_wt;
      double boot_leaf_size = 0.0;
      ib_idx.reserve(g);
      ib_wt.reserve(g);
      for (int a = 0; a < g; a++) {
        int freq = inbag_freq[group[a]];
        if (freq > 0) {
          ib_idx.push_back(group[a]);
          ib_wt.push_back((double)freq);
          boot_leaf_size += freq;
        }
      }
      if (boot_leaf_size < 1.0) continue;

      double inv_bls = 1.0 / boot_leaf_size;
      int n_ib = (int)ib_idx.size();

      for (int a = 0; a < g; a++) {
        int ia = group[a];
        fw_denom_local[ia] += 1.0;
        double* row = &fw_local[ia * n];
        for (int b = 0; b < n_ib; b++) {
          row[ib_idx[b]] += ib_wt[b] * inv_bls;
        }
      }
    }

    if (compute_prox) {
      if (prox_mode == 0) {
        for (auto& kv : leaf_groups) {
          const std::vector<int>& group = kv.second;
          int g = (int)group.size();
          for (int a = 0; a < g; a++) {
            int ia = group[a];
            prox_local[ia * n + ia] += 1.0;
            for (int b = a + 1; b < g; b++) {
              int ib = group[b];
              prox_local[ia * n + ib] += 1.0;
              prox_local[ib * n + ia] += 1.0;
            }
          }
        }
      } else {
        std::vector<int> prox_members;
        prox_members.reserve(n);
        if (prox_mode == 1) {
          for (int i = 0; i < n; i++) {
            if (inbag_freq[i] > 0) prox_members.push_back(i);
          }
        } else {
          for (int i = 0; i < n; i++) {
            if (inbag_freq[i] == 0) prox_members.push_back(i);
          }
        }

        std::unordered_map<int, std::vector<int>> prox_leaf_groups;
        for (int idx : prox_members) {
          prox_leaf_groups[leaf_ids[idx]].push_back(idx);
        }
        for (int a = 0; a < (int)prox_members.size(); a++) {
          int ia = prox_members[a];
          prox_denom_local[ia * n + ia] += 1.0;
          for (int b = a + 1; b < (int)prox_members.size(); b++) {
            int ib = prox_members[b];
            prox_denom_local[ia * n + ib] += 1.0;
            prox_denom_local[ib * n + ia] += 1.0;
          }
        }
        for (auto& kv : prox_leaf_groups) {
          const std::vector<int>& group = kv.second;
          int g = (int)group.size();
          for (int a = 0; a < g; a++) {
            int ia = group[a];
            prox_local[ia * n + ia] += 1.0;
            for (int b = a + 1; b < g; b++) {
              int ib = group[b];
              prox_local[ia * n + ib] += 1.0;
              prox_local[ib * n + ia] += 1.0;
            }
          }
        }
      }
    }

    tree_results[t].nodes = std::move(tree);
    tree_results[t].leaf_ids = std::move(leaf_ids);
  }
    #ifdef _OPENMP
    #pragma omp critical
    #endif
    {
      for (std::size_t i = 0; i < fw_denom_buf.size(); ++i) {
        fw_denom_buf[i] += fw_denom_local[i];
      }
      if (prox_mode > 0) {
        for (std::size_t idx = 0; idx < prox_denom_buf.size(); ++idx) {
          prox_denom_buf[idx] += prox_denom_local[idx];
        }
      }
      for (std::size_t idx = 0; idx < fw_buf.size(); ++idx) {
        fw_buf[idx] += fw_local[idx];
      }
      if (compute_prox) {
        for (std::size_t idx = 0; idx < prox_buf.size(); ++idx) {
          prox_buf[idx] += prox_local[idx];
        }
      }
    }
  }

  // Copy to R matrices + normalize
  NumericMatrix forest_wt(n, n);
  NumericMatrix prox(n, n);
  IntegerMatrix membership(n, ntree);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      forest_wt(i, j) = (fw_denom_buf[i] > 0.0) ? fw_buf[i * n + j] / fw_denom_buf[i] : NA_REAL;
      if (compute_prox) {
        if (prox_mode == 0) {
          prox(i, j) = prox_buf[i * n + j] / ntree;
        } else {
          prox(i, j) = (prox_denom_buf[i * n + j] > 0.0) ? prox_buf[i * n + j] / prox_denom_buf[i * n + j] : NA_REAL;
        }
      }
    }
  }

  List tree_info_list(ntree);
  for (int t = 0; t < ntree; t++) {
    for (int i = 0; i < n; i++) {
      membership(i, t) = tree_results[t].leaf_ids[i];
    }

    int n_nodes = (int)tree_results[t].nodes.size();

    IntegerVector t_split_var(n_nodes);
    NumericVector t_split_val(n_nodes);
    IntegerVector t_left(n_nodes);
    IntegerVector t_right(n_nodes);
    IntegerVector t_depth(n_nodes);
    IntegerVector t_nodesize(n_nodes);
    LogicalVector t_is_leaf(n_nodes);

    // split_var is already a global column index (no remap needed)
    for (int ni = 0; ni < n_nodes; ni++) {
      int sv = tree_results[t].nodes[ni].split_var;
      t_split_var[ni] = sv;
      t_split_val[ni] = tree_results[t].nodes[ni].split_val;
      t_left[ni] = tree_results[t].nodes[ni].left;
      t_right[ni] = tree_results[t].nodes[ni].right;
      t_depth[ni] = tree_results[t].nodes[ni].depth;
      t_nodesize[ni] = tree_results[t].nodes[ni].nodesize;
      t_is_leaf[ni] = (sv < 0);
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
    Named("ntree") = ntree,
    Named("n") = n,
    Named("p") = p,
    Named("tree_info") = tree_info_list
  );
}
