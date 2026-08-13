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
#include <map>

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
  // Multivariate splitting response variable (MSRV): the response with the
  // largest component split statistic at this node.  A leaf (or an
  // unsupervised node) has no MSRV and keeps the sentinel -1.
  int msrv = -1;
  // IMD: which X variable was used to split (same as split_var but
  // stored here for clarity; the importance goes to this X var)
  double imd_x_score = 0.0;  // total split score at this node
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
  inline double operator()(int i, int j) const {
    return data[(std::size_t)i + (std::size_t)j * nrow_];
  }
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

// Weighted sampling without replacement for candidate variables.  Variables
// with zero weight are not eligible while positive mass exists.  The uniform
// helper above remains the exact path when no weights are supplied, preserving
// the established positive-seed stream for ordinary fits.
template <typename RNG>
static std::vector<int> sample_from_pool_weighted(
    const std::vector<int>& pool_in,
    const std::vector<double>& weights,
    int sample_size,
    RNG& rng)
{
  std::vector<int> pool;
  std::vector<double> wt;
  pool.reserve(pool_in.size());
  wt.reserve(pool_in.size());
  for (int idx : pool_in) {
    double w = (idx >= 0 && idx < (int)weights.size()) ? weights[idx] : 0.0;
    if (std::isfinite(w) && w > 0.0) {
      pool.push_back(idx);
      wt.push_back(w);
    }
  }
  if (pool.empty()) {
    return sample_from_pool_rfsrc_style(pool_in, sample_size, rng);
  }

  sample_size = std::max(1, std::min(sample_size, (int)pool.size()));
  std::vector<int> out;
  out.reserve(sample_size);
  for (int draw = 0; draw < sample_size; draw++) {
    double total = std::accumulate(wt.begin(), wt.end(), 0.0);
    double target = random_unit(rng) * total;
    double cumulative = 0.0;
    int pick = (int)pool.size() - 1;
    for (int k = 0; k < (int)pool.size(); k++) {
      cumulative += wt[k];
      if (target < cumulative) { pick = k; break; }
    }
    out.push_back(pool[pick]);
    pool.erase(pool.begin() + pick);
    wt.erase(wt.begin() + pick);
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
    const std::vector<double>* xvar_wt,
    const std::vector<double>* yvar_wt,
    RNG& rng,
    int& best_var, double& best_val, double& best_score,
    std::vector<int>& left_samples, std::vector<int>& right_samples,
    std::vector<double>& best_y_stats, int& best_msrv)
{
  int n_node = (int)samples.size();
  int px = X.ncol();
  int qy = Y.ncol();

  if (n_node < 2 * nodesize_min) return false;

  // Random subset of X columns (mtry)
  std::vector<int> x_pool(px);
  std::iota(x_pool.begin(), x_pool.end(), 0);
  std::vector<int> x_candidates = xvar_wt
    ? sample_from_pool_weighted(x_pool, *xvar_wt, std::min(mtry, px), rng)
    : sample_from_pool_rfsrc_style(x_pool, std::min(mtry, px), rng);
  int n_x_try = (int)x_candidates.size();

  // Random subset of Y columns (ytry)
  std::vector<int> y_pool(qy);
  std::iota(y_pool.begin(), y_pool.end(), 0);
  std::vector<int> y_candidates = yvar_wt
    ? sample_from_pool_weighted(y_pool, *yvar_wt, std::min(ytry, qy), rng)
    : sample_from_pool_rfsrc_style(y_pool, std::min(ytry, qy), rng);
  int n_y_try = (int)y_candidates.size();

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
    double var = (n_node > 1) ? ss / n_node : 0.0;
    y_sds[jj] = (var > 0.0) ? std::sqrt(var) : 0.0;
  }

  best_score = -1.0;
  best_var = -1;
  best_val = 0.0;
  best_msrv = -1;

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
      int local_k = sample_pos[si];

      // Check split before adding this sample to left
      if (eval_idx < (int)eval_positions.size() && s == next_eval) {
        double prev_x = X(order[s - 1], xvar);
        double score = 0.0;
        int deltaNorm = 0;
        for (int jj = 0; jj < n_y_try; jj++) {
          if (y_sds[jj] > 0.0) {
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
            double max_component = -1.0;
            for (int jj = 0; jj < n_y_try; jj++) {
              double sL = sum_L[jj];
              double component = (sL * sL) / nL +
                (sL * sL) / (n_node - nL);
              int yv = y_candidates[jj];
              best_y_stats[yv] = component;
              if (y_sds[jj] > 0.0 && component > max_component) {
                max_component = component;
                best_msrv = yv;
              }
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
    const std::vector<double>* xvar_wt,
    const std::vector<double>* yvar_wt,
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

  // Build root's sorted indices with one entry per bootstrap draw.  Keeping
  // multiplicity here is essential for SWR: node.samples contains duplicate
  // row IDs and every split statistic must count those draws independently.
  std::vector<int> bag_frequency(n_total, 0);
  for (int si : bag_samples) bag_frequency[si]++;

  SplitTask root_task;
  root_task.node_id = 0;
  root_task.sorted.resize(px);
  for (int j = 0; j < px; j++) {
    root_task.sorted[j].reserve(bag_samples.size());
    for (int si : sort_order[j]) {
      for (int k = 0; k < bag_frequency[si]; k++) {
        root_task.sorted[j].push_back(si);
      }
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
      int bmsrv = -1;

      bool found = find_best_split_part(
        X, Y, node.samples, task.sorted, sample_pos,
        mtry, ytry, nodesize_min, nsplit, xvar_wt, yvar_wt, rng,
        bv, bval, bscore, lsamp, rsamp, by_stats, bmsrv);

      if (!found) continue;

      node.split_var = bv;
      node.split_val = bval;
      node.imd_y_stats = std::move(by_stats);
      node.imd_x_score = bscore;
      node.msrv = bmsrv;

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
    double var = (n_node > 1) ? ss / n_node : 0.0;
    y_sds[jj] = (var > 0.0) ? std::sqrt(var) : 0.0;
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
          if (y_sds[jj] > 0.0) {
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
    double var = (n_node > 1) ? ss / n_node : 0.0;
    y_sds[jj] = (var > 0.0) ? std::sqrt(var) : 0.0;
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
        if (y_sds[jj] > 0.0) {
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
    const std::vector<int>& node_frequency,
    const std::vector<std::vector<int>>& sort_order,
    std::vector<int>& sample_pos,     // reusable n-length scratch buffer
    int n_total,
    int mtry, int ytry, int nodesize_min, int nsplit,
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
      double var = (n_node > 1) ? ss / n_node : 0.0;
      y_sds[jj] = (var > 0.0) ? std::sqrt(var) : 0.0;
    }

    bool any_informative = false;
    for (int jj = 0; jj < n_y_try; jj++) {
      if (y_sds[jj] > 0.0) { any_informative = true; break; }
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

    // Identify all distinct-value boundaries in the node.  Positions are
    // measured in bootstrap draws (not unique row IDs), so SWR multiplicity
    // contributes to both child sizes and the pseudo-response sums.
    const std::vector<int>& order = sort_order[xvar];
    std::vector<int> boundary_rows;
    boundary_rows.reserve(n_node);
    double boundary_prev_x = 0.0;
    bool boundary_first = true;
    for (int s = 0; s < n_total; s++) {
      int si = order[s];
      if (node_frequency[si] <= 0) continue;
      double x_val = D(si, xvar);
      if (!boundary_first && x_val != boundary_prev_x) boundary_rows.push_back(s);
      boundary_prev_x = x_val;
      boundary_first = false;
    }
    if (boundary_rows.empty()) continue;
    std::vector<int> eval_rows = sample_cutpoint_positions(boundary_rows, nsplit, rng);
    int eval_idx = 0;

    // Scan pre-sorted unique row IDs, adding each row's bootstrap frequency.
    std::vector<double> sum_L(n_y_try, 0.0);
    int nL = 0;
    double prev_x = 0.0;
    bool first = true;

    for (int s = 0; s < n_total; s++) {
      int si = order[s];
      int multiplicity = node_frequency[si];
      if (multiplicity <= 0) continue;

      double x_val = D(si, xvar);
      int local_k = sample_pos[si];

      // nodesize enforced at parent level only; match rfsrc cutpoint check (>0)
      if (!first && x_val != prev_x &&
          eval_idx < (int)eval_rows.size() && s == eval_rows[eval_idx]) {
        int nR = n_node - nL;
        double score = 0.0;
        int deltaNorm = 0;
        for (int jj = 0; jj < n_y_try; jj++) {
          if (y_sds[jj] > 0.0) {
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
        eval_idx++;
      }

      nL += multiplicity;
      for (int jj = 0; jj < n_y_try; jj++) {
        sum_L[jj] += multiplicity * y_std_flat[jj * n_node + local_k];
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
    int mtry, int ytry, int nodesize_min, int max_depth, int nsplit,
    std::mt19937& rng)
{
  std::vector<Node> nodes;
  nodes.reserve(256);

  std::vector<int> node_frequency(n_total, 0);
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

      // Count bootstrap multiplicity for the current node.
      for (int si : node.samples) node_frequency[si]++;

      if (!find_best_split_unsup(D, node.samples, node_frequency, sort_order,
                                  sample_pos, n_total,
                                  mtry, ytry, nodesize_min, nsplit,
                                  rng, bv, bval, bscore, lsamp, rsamp)) {
        for (int si : node.samples) node_frequency[si] = 0;
        continue;
      }

      // Clear multiplicities after split.
      for (int si : node.samples) node_frequency[si] = 0;

      node.split_var = bv;
      node.split_val = bval;
      node.imd_x_score = bscore;  // IMD: store split score for X importance

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

// For supervised enhanced proximity, the legacy R implementation builds
// predictor and response PCA embeddings separately, computes one Spearman
// correlation in each coordinate system, then averages the two correlations.
// `split` is the number of leading predictor-embedding columns. A non-positive
// split keeps the one-embedding behavior used by unsupervised forests.
static double grouped_spearman_corr(const double* a, const double* b,
                                    int d, int split) {
  if (split <= 0 || split >= d) return spearman_corr(a, b, d);
  double corr_x = spearman_corr(a, b, split);
  double corr_y = spearman_corr(a + split, b + split, d - split);
  return 0.5 * (corr_x + corr_y);
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
                       int embed_split = 0,
                       double sibling_gamma = 0.5,
                       int enhanced_prox_mode = 0,
                       int forest_wt_mode = 0,
                       Nullable<NumericVector> xvar_wt = R_NilValue,
                       Nullable<NumericVector> yvar_wt = R_NilValue) {
  // enhanced_prox_mode: 0 = off, 1 = compute enhanced proximity

  int n = X.nrow();
  int px = X.ncol();
  int qy = Y.ncol();

  if (n < 2 || px < 1) {
    stop("X must contain at least two rows and one column");
  }
  if (Y.nrow() != n || qy < 1) {
    stop("Y must have nrow(X) rows and at least one column");
  }
  if (ntree < 1) stop("ntree must be a positive integer");
  if (seed == NA_INTEGER) stop("seed must be a finite integer");
  if (nodesize_min < 1) stop("nodesize_min must be a positive integer");
  if (max_depth < 0) stop("max_depth must be a non-negative integer");
  if (nthread < 0) stop("nthread must be a non-negative integer");
  if (samptype < 0 || samptype > 1) {
    stop("samptype must be 0 ('swor') or 1 ('swr')");
  }
  if (prox_mode < -1 || prox_mode > 2) {
    stop("prox_mode must be -1 ('none'), 0 ('all'), 1 ('inbag'), or 2 ('oob')");
  }
  if (enhanced_prox_mode < 0 || enhanced_prox_mode > 1) {
    stop("enhanced_prox_mode must be 0 or 1");
  }
  if (!std::isfinite(sibling_gamma) || sibling_gamma < 0.0) {
    stop("sibling_gamma must be finite and non-negative");
  }
  if (embed_split < 0) stop("embed_split must be a non-negative integer");

  // Defaults: ceiling(px/3) for mtry and ceiling(qy/3) for ytry,
  // and nsplit = 10 to match randomForestSRC's default randomized cut search.
  if (mtry <= 0) mtry = std::max(1, (int)std::ceil((double)px / 3.0));
  if (ytry <= 0) ytry = std::max(1, (int)std::ceil((double)qy / 3.0));
  mtry = std::min(mtry, px);
  ytry = std::min(ytry, qy);
  if (nsplit < 0) stop("nsplit must be a non-negative integer");
  if (forest_wt_mode < 0 || forest_wt_mode > 2) {
    stop("forest_wt_mode must be 0 ('all'), 1 ('inbag'), or 2 ('oob')");
  }
  // max_depth <= 0 means unlimited (grow until nodesize constraint only)

  #ifdef _OPENMP
  if (nthread > 0) omp_set_num_threads(nthread);
  #endif

  // Copy to thread-safe MatrixView (avoid Rcpp operator() inside OpenMP)
  std::vector<double> X_buf((std::size_t)n * px);
  std::vector<double> Y_buf((std::size_t)n * qy);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < px; j++) {
      double value = X(i, j);
      if (!std::isfinite(value)) stop("X must contain only finite values");
      X_buf[(std::size_t)i + (std::size_t)j * n] = value;
    }
    for (int j = 0; j < qy; j++) {
      double value = Y(i, j);
      if (!std::isfinite(value)) stop("Y must contain only finite values");
      Y_buf[(std::size_t)i + (std::size_t)j * n] = value;
    }
  }
  MatrixView Xv(X_buf.data(), n, px);
  MatrixView Yv(Y_buf.data(), n, qy);

  // Optional variable-selection weights.  Copy before entering OpenMP so no
  // R API is accessed by worker threads.
  std::vector<double> xvar_wt_buf, yvar_wt_buf;
  const std::vector<double>* xvar_wt_ptr = nullptr;
  const std::vector<double>* yvar_wt_ptr = nullptr;
  if (xvar_wt.isNotNull()) {
    NumericVector w(xvar_wt.get());
    if (w.size() != px) stop("xvar_wt must have length ncol(X)");
    xvar_wt_buf.assign(w.begin(), w.end());
    double total = 0.0;
    for (double value : xvar_wt_buf) {
      if (!std::isfinite(value) || value < 0.0) {
        stop("xvar_wt must be finite and non-negative");
      }
      total += value;
    }
    if (!(total > 0.0) || !std::isfinite(total)) {
      stop("xvar_wt must contain positive finite mass");
    }
    xvar_wt_ptr = &xvar_wt_buf;
  }
  if (yvar_wt.isNotNull()) {
    NumericVector w(yvar_wt.get());
    if (w.size() != qy) stop("yvar_wt must have length ncol(Y)");
    yvar_wt_buf.assign(w.begin(), w.end());
    double total = 0.0;
    for (double value : yvar_wt_buf) {
      if (!std::isfinite(value) || value < 0.0) {
        stop("yvar_wt must be finite and non-negative");
      }
      total += value;
    }
    if (!(total > 0.0) || !std::isfinite(total)) {
      stop("yvar_wt must contain positive finite mass");
    }
    yvar_wt_ptr = &yvar_wt_buf;
  }

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
    if (embed_mat.nrow() != n || embed_mat.ncol() < 1) {
      stop("embed must have nrow(X) rows and at least one column");
    }
    embed_dim = embed_mat.ncol();
    if (embed_split >= embed_dim) {
      stop("embed_split must be smaller than ncol(embed)");
    }
    embed_buf.resize((std::size_t)n * embed_dim);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < embed_dim; j++) {
        double value = embed_mat(i, j);
        if (!std::isfinite(value)) stop("embed must contain only finite values");
        embed_buf[(std::size_t)i + (std::size_t)j * n] = value;
      }
  } else if (enhanced_prox_mode > 0) {
    stop("embed must be supplied when enhanced_prox_mode is enabled");
  }

  // Output: plain C++ buffers for thread-safe accumulation
  // prox_mode: -1 = skip proximity entirely, 0 = all, 1 = inbag, 2 = oob
  bool compute_prox = (prox_mode >= 0);
  const std::size_t nn = (std::size_t)n * (std::size_t)n;
  std::vector<double> fw_buf(nn, 0.0);
  std::vector<double> fw_denom_buf(n, 0.0);
  std::vector<double> prox_buf(compute_prox ? nn : 0, 0.0);
  std::vector<double> prox_denom_buf(prox_mode > 0 ? nn : 0, 0.0);
  std::vector<double> eprox_buf(compute_enhanced ? nn : 0, 0.0);

  // Per-tree results: tree structure + leaf membership + bootstrap frequency
  struct TreeResult {
    std::vector<Node> nodes;
    std::vector<int> leaf_ids;
    std::vector<int> inbag;   // bootstrap frequency n_{i,b} per sample
  };
  std::vector<TreeResult> tree_results(ntree);

  // Pre-sort all X columns once (shared across trees, read-only)
  auto sort_order = presort_columns(Xv, n, px);

  // Phase 1: build trees in parallel.  Each tree owns its RNG streams and its
  // result slot, so scheduling cannot change the fitted forest.
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
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
      int samp_size = std::max(1, (int)std::round(0.632 * n));
      bag = sample_swor_rfsrc_style(n, samp_size, rng_boot);
      for (int i = 0; i < samp_size; i++) {
        inbag_freq[bag[i]] = 1;
      }
    }
    tree_results[t].inbag = std::move(inbag_freq);
    std::sort(bag.begin(), bag.end());

    // Build tree with partition-based sorted indices (ranger-style)
    tree_results[t].nodes = build_tree_part(Xv, Yv, bag, sort_order, n, px,
                                             mtry, ytry, nodesize_min,
                                             max_depth, nsplit, xvar_wt_ptr,
                                             yvar_wt_ptr, rng_split);

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
  }

  // Phase 2: accumulate strictly in tree-index order.  Besides making floating
  // point results bit-identical across thread counts, this avoids allocating
  // multiple n-by-n matrices per OpenMP worker.
  std::vector<double> centroid_scratch(compute_enhanced ? 2 * embed_dim : 0, 0.0);

  for (int t = 0; t < ntree; t++) {
    const auto& inbag_freq = tree_results[t].inbag;
    const auto& leaf_ids = tree_results[t].leaf_ids;

    // Accumulate forest weights.
    //   all   (0): all query rows, all observed rows as uniform donors (Eq. 3)
    //   inbag (1): all query rows, bootstrap-frequency-weighted inbag donors
    //   oob   (2): OOB query rows only, bootstrap-frequency-weighted donors
    std::map<int, std::vector<int>> leaf_groups;
    for (int i = 0; i < n; i++) {
      leaf_groups[leaf_ids[i]].push_back(i);
    }

    for (auto& kv : leaf_groups) {
      const std::vector<int>& group = kv.second;
      int g = (int)group.size();
      if (g == 0) continue;

      std::vector<int> donor_idx;
      std::vector<double> donor_wt;
      double donor_mass = 0.0;
      donor_idx.reserve(g);
      donor_wt.reserve(g);
      for (int idx : group) {
        double wt = (forest_wt_mode == 0) ? 1.0 : (double)inbag_freq[idx];
        if (wt > 0.0) {
          donor_idx.push_back(idx);
          donor_wt.push_back(wt);
          donor_mass += wt;
        }
      }
      if (donor_mass < 1.0) continue;
      double inv_mass = 1.0 / donor_mass;

      for (int a = 0; a < g; a++) {
        int ia = group[a];
        if (forest_wt_mode == 2 && inbag_freq[ia] > 0) continue;
        fw_denom_buf[ia] += 1.0;
        double* row = &fw_buf[(std::size_t)ia * n];
        for (int b = 0; b < (int)donor_idx.size(); b++) {
          row[donor_idx[b]] += donor_wt[b] * inv_mass;
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
            prox_buf[(std::size_t)ia * n + ia] += 1.0;
            for (int b = a + 1; b < g; b++) {
              int ib = group[b];
              prox_buf[(std::size_t)ia * n + ib] += 1.0;
              prox_buf[(std::size_t)ib * n + ia] += 1.0;
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

        std::map<int, std::vector<int>> prox_leaf_groups;
        for (int idx : prox_members) {
          prox_leaf_groups[leaf_ids[idx]].push_back(idx);
        }
        for (int a = 0; a < (int)prox_members.size(); a++) {
          int ia = prox_members[a];
          prox_denom_buf[(std::size_t)ia * n + ia] += 1.0;
          for (int b = a + 1; b < (int)prox_members.size(); b++) {
            int ib = prox_members[b];
            prox_denom_buf[(std::size_t)ia * n + ib] += 1.0;
            prox_denom_buf[(std::size_t)ib * n + ia] += 1.0;
          }
        }
        for (auto& kv : prox_leaf_groups) {
          const std::vector<int>& group = kv.second;
          int g = (int)group.size();
          for (int a = 0; a < g; a++) {
            int ia = group[a];
            prox_buf[(std::size_t)ia * n + ia] += 1.0;
            for (int b = a + 1; b < g; b++) {
              int ib = group[b];
              prox_buf[(std::size_t)ia * n + ib] += 1.0;
              prox_buf[(std::size_t)ib * n + ia] += 1.0;
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
      // 1. Same-leaf contribution (weight = 1)
      for (auto& kv : leaf_groups) {
        const std::vector<int>& group = kv.second;
        int g = (int)group.size();
        for (int a = 0; a < g; a++) {
          int ia = group[a];
          eprox_buf[(std::size_t)ia * n + ia] += 1.0;
          for (int b = a + 1; b < g; b++) {
            int ib = group[b];
            eprox_buf[(std::size_t)ia * n + ib] += 1.0;
            eprox_buf[(std::size_t)ib * n + ia] += 1.0;
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
            cent_l[d] += embed_buf[(std::size_t)si + (std::size_t)d * n];
        for (int si : grp_r)
          for (int d = 0; d < embed_dim; d++)
            cent_r[d] += embed_buf[(std::size_t)si + (std::size_t)d * n];
        double inv_l = 1.0 / grp_l.size();
        double inv_r = 1.0 / grp_r.size();
        for (int d = 0; d < embed_dim; d++) { cent_l[d] *= inv_l; cent_r[d] *= inv_r; }

        // Spearman correlation of centroids
        double corr = grouped_spearman_corr(
          cent_l, cent_r, embed_dim, embed_split
        );
        double w = sibling_gamma * std::max(corr, 0.0);
        if (w <= 0.0) continue;
        if (w > 1.0) w = 1.0;

        // Add sibling proximity
        for (int a : grp_l) {
          for (int b : grp_r) {
            eprox_buf[(std::size_t)a * n + b] += w;
            eprox_buf[(std::size_t)b * n + a] += w;
          }
        }
      }
    }
  }

  // Copy one dense buffer at a time and release it before allocating the next
  // output matrix, keeping peak memory close to the size of the returned data.
  NumericMatrix forest_wt(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      forest_wt(i, j) = (fw_denom_buf[i] > 0.0) ? fw_buf[(std::size_t)i * n + j] / fw_denom_buf[i] : NA_REAL;
    }
  }
  std::vector<double>().swap(fw_buf);
  std::vector<double>().swap(fw_denom_buf);

  NumericMatrix prox(n, n);
  if (compute_prox) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (prox_mode == 0) {
          prox(i, j) = prox_buf[(std::size_t)i * n + j] / ntree;
        } else {
          const std::size_t idx = (std::size_t)i * n + j;
          prox(i, j) = (prox_denom_buf[idx] > 0.0) ? prox_buf[idx] / prox_denom_buf[idx] : NA_REAL;
        }
      }
    }
  }
  std::vector<double>().swap(prox_buf);
  std::vector<double>().swap(prox_denom_buf);

  // ──── Enhanced proximity normalization ────
  NumericMatrix enhanced_prox(compute_enhanced ? n : 0, compute_enhanced ? n : 0);
  if (compute_enhanced) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        enhanced_prox(i, j) = eprox_buf[(std::size_t)i * n + j] / ntree;
      }
    }
  }
  std::vector<double>().swap(eprox_buf);

  IntegerMatrix membership(n, ntree);
  IntegerMatrix inbag_mat(n, ntree);

  // ──── Inverse minimal depth (IMD) aggregation ────
  // For each tree and variable, IMD = 1 / (minimum split depth + 1), or zero
  // when absent.  Y variables are considered only when selected as the node's
  // MSRV (largest component G_j among the ytry responses at that split).
  NumericMatrix imd_x_per_tree(px, ntree);
  NumericMatrix imd_y_per_tree(qy, ntree);

  // Pairwise semantics: every split contributes its depth weight to the one
  // X/MSRV pair assigned to that node, then the matrix is averaged over trees.
  std::vector<double> pairwise_buf((std::size_t)px * qy, 0.0);

  for (int t = 0; t < ntree; t++) {
    std::vector<int> min_x_depth(px, -1);
    std::vector<int> min_y_depth(qy, -1);

    for (const auto& node : tree_results[t].nodes) {
      if (node.split_var < 0) continue;
      int xv = node.split_var;
      if (min_x_depth[xv] < 0 || node.depth < min_x_depth[xv]) {
        min_x_depth[xv] = node.depth;
      }
      int yv = node.msrv;
      if (yv >= 0 && yv < qy) {
        if (min_y_depth[yv] < 0 || node.depth < min_y_depth[yv]) {
          min_y_depth[yv] = node.depth;
        }
        pairwise_buf[(std::size_t)xv * qy + yv] += 1.0 / (node.depth + 1.0);
      }
    }

    for (int j = 0; j < px; j++) {
      imd_x_per_tree(j, t) = (min_x_depth[j] >= 0)
        ? 1.0 / (min_x_depth[j] + 1.0) : 0.0;
    }
    for (int j = 0; j < qy; j++) {
      imd_y_per_tree(j, t) = (min_y_depth[j] >= 0)
        ? 1.0 / (min_y_depth[j] + 1.0) : 0.0;
    }
  }

  // Normalize pairwise matrix by ntree
  NumericMatrix pairwise_xy(px, qy);
  for (int i = 0; i < px; i++) {
    for (int j = 0; j < qy; j++) {
      pairwise_xy(i, j) = pairwise_buf[(std::size_t)i * qy + j] / ntree;
    }
  }

  // Forest IMD is the simple mean of per-tree inverse minimal depths.
  NumericVector imd_x(px);
  NumericVector imd_y(qy);
  for (int j = 0; j < px; j++) {
    double sum = 0.0;
    for (int t = 0; t < ntree; t++) sum += imd_x_per_tree(j, t);
    imd_x[j] = sum / ntree;
  }
  for (int j = 0; j < qy; j++) {
    double sum = 0.0;
    for (int t = 0; t < ntree; t++) sum += imd_y_per_tree(j, t);
    imd_y[j] = sum / ntree;
  }

  // Free heavy per-node data no longer needed (samples only).
  // Keep imd_x_score and imd_y_stats for cluster-weighted IMD.
  for (int t = 0; t < ntree; t++) {
    for (auto& node : tree_results[t].nodes) {
      std::vector<int>().swap(node.samples);
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
    NumericVector t_imd_x_score(n_nodes);
    IntegerVector t_msrv(n_nodes, NA_INTEGER);

    for (int ni = 0; ni < n_nodes; ni++) {
      t_split_var[ni] = tree_results[t].nodes[ni].split_var;
      t_split_val[ni] = tree_results[t].nodes[ni].split_val;
      t_left[ni] = tree_results[t].nodes[ni].left;
      t_right[ni] = tree_results[t].nodes[ni].right;
      t_depth[ni] = tree_results[t].nodes[ni].depth;
      t_nodesize[ni] = tree_results[t].nodes[ni].nodesize;
      t_is_leaf[ni] = (tree_results[t].nodes[ni].split_var < 0);
      t_imd_x_score[ni] = tree_results[t].nodes[ni].imd_x_score;
      if (tree_results[t].nodes[ni].msrv >= 0) {
        t_msrv[ni] = tree_results[t].nodes[ni].msrv;
      }
    }

    // Per-node Y IMD stats: n_internal_nodes x qy matrix (only for internal nodes)
    int n_internal = 0;
    for (int ni = 0; ni < n_nodes; ni++) {
      if (tree_results[t].nodes[ni].split_var >= 0) n_internal++;
    }
    NumericMatrix t_imd_y_stats(n_internal > 0 ? qy : 0, n_internal);
    {
      int col = 0;
      for (int ni = 0; ni < n_nodes; ni++) {
        if (tree_results[t].nodes[ni].split_var < 0) continue;
        const auto& ys = tree_results[t].nodes[ni].imd_y_stats;
        for (int j = 0; j < qy && j < (int)ys.size(); j++) {
          t_imd_y_stats(j, col) = ys[j];
        }
        col++;
      }
    }

    tree_info_list[t] = List::create(
      Named("split_var") = t_split_var,
      Named("split_val") = t_split_val,
      Named("left") = t_left,
      Named("right") = t_right,
      Named("depth") = t_depth,
      Named("nodesize") = t_nodesize,
      Named("is_leaf") = t_is_leaf,
      Named("imd_x_score") = t_imd_x_score,
      Named("imd_y_stats") = t_imd_y_stats,
      Named("msrv") = t_msrv
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
    Named("forest_wt_mode") = forest_wt_mode,
    Named("ntree") = ntree,
    Named("mtry") = mtry,
    Named("ytry") = ytry,
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
                              int prox_mode = 0,
                              Nullable<NumericMatrix> embed = R_NilValue,
                              double sibling_gamma = 0.5,
                              int enhanced_prox_mode = 0,
                              int forest_wt_mode = 0,
                              int nsplit = 10) {
  // enhanced_prox_mode: 0 = off, 1 = compute enhanced proximity

  int n = data.nrow();
  int p = data.ncol();
  if (n < 2 || p < 1) {
    stop("data must contain at least two rows and one column");
  }
  if (ntree < 1) stop("ntree must be a positive integer");
  if (seed == NA_INTEGER) stop("seed must be a finite integer");
  if (nodesize_min < 1) stop("nodesize_min must be a positive integer");
  if (max_depth < 0) stop("max_depth must be a non-negative integer");
  if (nthread < 0) stop("nthread must be a non-negative integer");
  if (samptype < 0 || samptype > 1) {
    stop("samptype must be 0 ('swor') or 1 ('swr')");
  }
  if (prox_mode < -1 || prox_mode > 2) {
    stop("prox_mode must be -1 ('none'), 0 ('all'), 1 ('inbag'), or 2 ('oob')");
  }
  if (enhanced_prox_mode < 0 || enhanced_prox_mode > 1) {
    stop("enhanced_prox_mode must be 0 or 1");
  }
  if (!std::isfinite(sibling_gamma) || sibling_gamma < 0.0) {
    stop("sibling_gamma must be finite and non-negative");
  }
  if (forest_wt_mode < 0 || forest_wt_mode > 2) {
    stop("forest_wt_mode must be 0 ('all'), 1 ('inbag'), or 2 ('oob')");
  }
  if (nsplit < 0) stop("nsplit must be a non-negative integer");

  // max_depth <= 0 means unlimited (grow until nodesize constraint only)

  #ifdef _OPENMP
  if (nthread > 0) omp_set_num_threads(nthread);
  #endif

  unsigned int actual_seed = (seed < 0) ? std::random_device{}() : (unsigned int)seed;

  // ──── Enhanced proximity setup ────
  bool compute_enhanced = (enhanced_prox_mode > 0) && embed.isNotNull();
  int embed_dim = 0;
  std::vector<double> embed_buf;
  if (compute_enhanced) {
    NumericMatrix embed_mat(embed.get());
    if (embed_mat.nrow() != n || embed_mat.ncol() < 1) {
      stop("embed must have nrow(data) rows and at least one column");
    }
    embed_dim = embed_mat.ncol();
    embed_buf.resize((std::size_t)n * embed_dim);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < embed_dim; j++) {
        double value = embed_mat(i, j);
        if (!std::isfinite(value)) stop("embed must contain only finite values");
        embed_buf[(std::size_t)i + (std::size_t)j * n] = value;
      }
  } else if (enhanced_prox_mode > 0) {
    stop("embed must be supplied when enhanced_prox_mode is enabled");
  }

  bool compute_prox = (prox_mode >= 0);
  const std::size_t nn = (std::size_t)n * (std::size_t)n;
  std::vector<double> fw_buf(nn, 0.0);
  std::vector<double> fw_denom_buf(n, 0.0);
  std::vector<double> prox_buf(compute_prox ? nn : 0, 0.0);
  std::vector<double> prox_denom_buf(prox_mode > 0 ? nn : 0, 0.0);
  std::vector<double> eprox_buf(compute_enhanced ? nn : 0, 0.0);

  struct UnsupTreeResult {
    std::vector<Node> nodes;
    std::vector<int> leaf_ids;
    std::vector<int> inbag;
  };
  std::vector<UnsupTreeResult> tree_results(ntree);

  // Copy data to a thread-safe MatrixView (avoid Rcpp inside OpenMP)
  std::vector<double> data_buf((std::size_t)n * p);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < p; j++) {
      double value = data(i, j);
      if (!std::isfinite(value)) stop("data must contain only finite values");
      data_buf[(std::size_t)i + (std::size_t)j * n] = value;
    }
  MatrixView D(data_buf.data(), n, p);

  // mtry default: ceiling(sqrt(p)), matching rfsrc unsupervised
  int mtry_default = std::max(1, (int)std::ceil(std::sqrt((double)p)));
  // ytry default: min(ytry_param, p-1); ytry <= 0 means use p-1
  int ytry_use = (ytry <= 0) ? std::max(1, p - 1) : std::min(ytry, p - 1);

  // Pre-sort all columns once (shared across trees, read-only)
  auto sort_order_unsup = presort_columns(D, n, p);

  // Phase 1: Parallel tree building (unsupervised).
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (int t = 0; t < ntree; t++) {
    std::mt19937 rng_t(actual_seed + (unsigned int)t);
    std::uniform_int_distribution<int> boot_dist(0, n - 1);

    std::vector<int> bag;
    std::vector<int> inbag_freq(n, 0);

    if (samptype == 1) {
      bag.resize(n);
      for (int i = 0; i < n; i++) {
        int idx = boot_dist(rng_t);
        bag[i] = idx;
        inbag_freq[idx]++;
      }
    } else {
      int samp_size_u = std::max(1, (int)std::round(0.632 * n));
      bag = sample_swor_rfsrc_style(n, samp_size_u, rng_t);
      for (int i = 0; i < samp_size_u; i++) {
        inbag_freq[bag[i]] = 1;
      }
    }
    std::sort(bag.begin(), bag.end());

    std::vector<Node> tree = build_tree_unsup(D, bag, sort_order_unsup, n,
                                               mtry_default, ytry_use,
                                               nodesize_min, max_depth,
                                               nsplit, rng_t);

    for (auto& node : tree) {
      node.samples.clear();
      node.samples.shrink_to_fit();
    }

    std::vector<int> leaf_ids(n);
    for (int i = 0; i < n; i++) {
      leaf_ids[i] = predict_leaf(tree, D, i);
    }

    tree_results[t].nodes = std::move(tree);
    tree_results[t].leaf_ids = std::move(leaf_ids);
    tree_results[t].inbag = std::move(inbag_freq);
  }

  // Phase 2: Deterministic serial accumulation (unsupervised).
  std::vector<double> centroid_scratch(compute_enhanced ? 2 * embed_dim : 0, 0.0);

  for (int t = 0; t < ntree; t++) {
    const auto& inbag_freq = tree_results[t].inbag;
    const auto& leaf_ids = tree_results[t].leaf_ids;

    std::map<int, std::vector<int>> leaf_groups;
    for (int i = 0; i < n; i++) leaf_groups[leaf_ids[i]].push_back(i);

    for (auto& kv : leaf_groups) {
      const std::vector<int>& group = kv.second;
      int g = (int)group.size();
      if (g == 0) continue;

      std::vector<int> donor_idx;
      std::vector<double> donor_wt;
      double donor_mass = 0.0;
      donor_idx.reserve(g);
      donor_wt.reserve(g);
      for (int idx : group) {
        double wt = (forest_wt_mode == 0) ? 1.0 : (double)inbag_freq[idx];
        if (wt > 0.0) {
          donor_idx.push_back(idx);
          donor_wt.push_back(wt);
          donor_mass += wt;
        }
      }
      if (donor_mass < 1.0) continue;
      double inv_mass = 1.0 / donor_mass;

      for (int a = 0; a < g; a++) {
        int ia = group[a];
        if (forest_wt_mode == 2 && inbag_freq[ia] > 0) continue;
        fw_denom_buf[ia] += 1.0;
        double* row = &fw_buf[ia * n];
        for (int b = 0; b < (int)donor_idx.size(); b++) {
          row[donor_idx[b]] += donor_wt[b] * inv_mass;
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
            prox_buf[ia * n + ia] += 1.0;
            for (int b = a + 1; b < g; b++) {
              int ib = group[b];
              prox_buf[ia * n + ib] += 1.0;
              prox_buf[ib * n + ia] += 1.0;
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

        std::map<int, std::vector<int>> prox_leaf_groups;
        for (int idx : prox_members) {
          prox_leaf_groups[leaf_ids[idx]].push_back(idx);
        }
        for (int a = 0; a < (int)prox_members.size(); a++) {
          int ia = prox_members[a];
          prox_denom_buf[ia * n + ia] += 1.0;
          for (int b = a + 1; b < (int)prox_members.size(); b++) {
            int ib = prox_members[b];
            prox_denom_buf[ia * n + ib] += 1.0;
            prox_denom_buf[ib * n + ia] += 1.0;
          }
        }
        for (auto& kv : prox_leaf_groups) {
          const std::vector<int>& group = kv.second;
          int g = (int)group.size();
          for (int a = 0; a < g; a++) {
            int ia = group[a];
            prox_buf[ia * n + ia] += 1.0;
            for (int b = a + 1; b < g; b++) {
              int ib = group[b];
              prox_buf[ia * n + ib] += 1.0;
              prox_buf[ib * n + ia] += 1.0;
            }
          }
        }
      }
    }

    // ──── Enhanced proximity accumulation (unsupervised) ────
    if (compute_enhanced) {
      // 1. Same-leaf contribution (weight = 1)
      for (auto& kv : leaf_groups) {
        const std::vector<int>& group = kv.second;
        int g = (int)group.size();
        for (int a = 0; a < g; a++) {
          int ia = group[a];
          eprox_buf[ia * n + ia] += 1.0;
          for (int b = a + 1; b < g; b++) {
            int ib = group[b];
            eprox_buf[ia * n + ib] += 1.0;
            eprox_buf[ib * n + ia] += 1.0;
          }
        }
      }
      // 2. Sibling-leaf contribution
      const auto& tree_nodes = tree_results[t].nodes;
      for (const auto& nd : tree_nodes) {
        if (nd.split_var < 0) continue;
        int li = nd.left, ri = nd.right;
        if (li < 0 || ri < 0) continue;
        if (tree_nodes[li].split_var >= 0 || tree_nodes[ri].split_var >= 0) continue;
        auto it_l = leaf_groups.find(li);
        auto it_r = leaf_groups.find(ri);
        if (it_l == leaf_groups.end() || it_r == leaf_groups.end()) continue;
        const std::vector<int>& grp_l = it_l->second;
        const std::vector<int>& grp_r = it_r->second;
        if (grp_l.empty() || grp_r.empty()) continue;
        double* cent_l = centroid_scratch.data();
        double* cent_r = centroid_scratch.data() + embed_dim;
        std::fill(cent_l, cent_l + embed_dim, 0.0);
        std::fill(cent_r, cent_r + embed_dim, 0.0);
        for (int si : grp_l)
          for (int d = 0; d < embed_dim; d++)
            cent_l[d] += embed_buf[(std::size_t)si + (std::size_t)d * n];
        for (int si : grp_r)
          for (int d = 0; d < embed_dim; d++)
            cent_r[d] += embed_buf[(std::size_t)si + (std::size_t)d * n];
        double inv_l = 1.0 / grp_l.size();
        double inv_r = 1.0 / grp_r.size();
        for (int d = 0; d < embed_dim; d++) { cent_l[d] *= inv_l; cent_r[d] *= inv_r; }
        double corr = spearman_corr(cent_l, cent_r, embed_dim);
        double w = sibling_gamma * std::max(corr, 0.0);
        if (w <= 0.0) continue;
        if (w > 1.0) w = 1.0;
        for (int a : grp_l) {
          for (int b : grp_r) {
            eprox_buf[a * n + b] += w;
            eprox_buf[b * n + a] += w;
          }
        }
      }
    }
  }

  // ──── Enhanced proximity normalization (unsupervised) ────
  NumericMatrix enhanced_prox(compute_enhanced ? n : 0, compute_enhanced ? n : 0);
  if (compute_enhanced) {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        enhanced_prox(i, j) = eprox_buf[i * n + j] / ntree;
  }

  // Copy to R matrices + normalize
  NumericMatrix forest_wt(n, n);
  NumericMatrix prox(n, n);
  IntegerMatrix membership(n, ntree);
  IntegerMatrix inbag_mat(n, ntree);

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

  // ──── Inverse minimal depth (unsupervised: X-side only) ────
  NumericMatrix imd_x_per_tree(p, ntree);

  for (int t = 0; t < ntree; t++) {
    std::vector<int> min_depth(p, -1);

    for (const auto& node : tree_results[t].nodes) {
      if (node.split_var < 0) continue;
      int xv = node.split_var;
      if (min_depth[xv] < 0 || node.depth < min_depth[xv]) {
        min_depth[xv] = node.depth;
      }
    }

    for (int j = 0; j < p; j++) {
      imd_x_per_tree(j, t) = (min_depth[j] >= 0)
        ? 1.0 / (min_depth[j] + 1.0) : 0.0;
    }
  }

  // Simple forest mean; no across-variable normalization.
  NumericVector imd_x(p);
  for (int j = 0; j < p; j++) {
    double sum = 0.0;
    for (int t = 0; t < ntree; t++) sum += imd_x_per_tree(j, t);
    imd_x[j] = sum / ntree;
  }

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
    NumericVector t_imd_x_score(n_nodes);

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
      t_imd_x_score[ni] = tree_results[t].nodes[ni].imd_x_score;
    }

    tree_info_list[t] = List::create(
      Named("split_var") = t_split_var,
      Named("split_val") = t_split_val,
      Named("left") = t_left,
      Named("right") = t_right,
      Named("depth") = t_depth,
      Named("nodesize") = t_nodesize,
      Named("is_leaf") = t_is_leaf,
      Named("imd_x_score") = t_imd_x_score
    );
  }

  return List::create(
    Named("forest.wt") = forest_wt,
    Named("proximity") = prox,
    Named("enhanced_prox") = enhanced_prox,
    Named("membership") = membership,
    Named("inbag") = inbag_mat,
    Named("forest_wt_mode") = forest_wt_mode,
    Named("ntree") = ntree,
    Named("n") = n,
    Named("p") = p,
    Named("tree_info") = tree_info_list,
    Named("imd_x") = imd_x,
    Named("imd_x_per_tree") = imd_x_per_tree
  );
}
