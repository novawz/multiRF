#include <Rcpp.h>
using namespace Rcpp;

// ---------- build_tree_network_cpp ----------
// Replaces the R for-loop in get_tree_net that builds the edge network.
//
// Inputs:
//   var_conc   : character vector of concatenated var names (length = n_nodes)
//   var_tip_id : integer vector of tip IDs (length = n_nodes)
//   nodeSZ     : integer vector of node sizes (length = n_nodes)
//   dpthST     : numeric vector of depth stats (length = n_nodes)
//   is_leaf    : logical vector, TRUE if node is a leaf (length = n_nodes)
//
// Returns a List with from, to, from_id, to_id, inv_d, edge, nodesize vectors.
// [[Rcpp::export]]
List build_tree_network_cpp(CharacterVector var_conc,
                            IntegerVector var_tip_id,
                            IntegerVector nodeSZ,
                            NumericVector dpthST,
                            LogicalVector is_leaf) {

  int n = var_conc.size();
  if (n < 2) {
    return List::create(
      Named("from") = CharacterVector(0),
      Named("to") = CharacterVector(0),
      Named("from_id") = IntegerVector(0),
      Named("to_id") = IntegerVector(0),
      Named("inv_d") = NumericVector(0),
      Named("edge") = NumericVector(0),
      Named("nodesize") = IntegerVector(0)
    );
  }

  int max_edges = n - 1;
  CharacterVector from_vec(max_edges);
  CharacterVector to_vec(max_edges);
  IntegerVector from_id_vec(max_edges);
  IntegerVector to_id_vec(max_edges);
  NumericVector inv_d_vec(max_edges);
  NumericVector edge_vec(max_edges);
  IntegerVector nodesize_vec(max_edges);

  // Track children count per internal node
  // Use parallel vectors for name -> count mapping
  std::vector<std::string> nc_names;
  std::vector<int> nc_counts;
  nc_names.reserve(n);
  nc_counts.reserve(n);

  // Pre-populate with internal nodes
  for (int i = 0; i < n; i++) {
    if (!is_leaf[i]) {
      std::string nm = as<std::string>(var_conc[i]);
      // Check for duplicate
      bool found = false;
      for (size_t k = 0; k < nc_names.size(); k++) {
        if (nc_names[k] == nm) { found = true; break; }
      }
      if (!found) {
        nc_names.push_back(nm);
        nc_counts.push_back(0);
      }
    }
  }

  // Helper lambda: find count index by name
  auto find_nc = [&](const std::string& nm) -> int {
    for (size_t k = 0; k < nc_names.size(); k++) {
      if (nc_names[k] == nm) return (int)k;
    }
    return -1;
  };

  // Helper lambda: find parent in edges built so far
  // network$from[network$to == from_node]
  auto find_parent = [&](const std::string& child, int edge_idx) -> std::string {
    for (int k = edge_idx - 1; k >= 0; k--) {
      if (as<std::string>(to_vec[k]) == child) {
        return as<std::string>(from_vec[k]);
      }
    }
    return child; // shouldn't happen
  };

  auto find_parent_id = [&](int child_id, int edge_idx) -> int {
    for (int k = edge_idx - 1; k >= 0; k--) {
      if (to_id_vec[k] == child_id) {
        return from_id_vec[k];
      }
    }
    return child_id;
  };

  std::string from_node = as<std::string>(var_conc[0]);
  int from_id = var_tip_id[0];
  double ns_all = (double)nodeSZ[0];
  int edge_idx = 0;

  for (int i = 1; i < n; i++) {
    std::string to_node = as<std::string>(var_conc[i]);
    int to_id = var_tip_id[i];
    double ns_part = (double)nodeSZ[i];
    double depth_val = dpthST[i];
    double inv_d = (depth_val != 0.0) ? 1.0 / depth_val : 1.0;

    from_vec[edge_idx] = from_node;
    to_vec[edge_idx] = to_node;
    from_id_vec[edge_idx] = from_id;
    to_id_vec[edge_idx] = to_id;
    inv_d_vec[edge_idx] = inv_d;
    edge_vec[edge_idx] = (ns_part != 0.0) ? ns_all / ns_part : 1.0;
    nodesize_vec[edge_idx] = (int)ns_part;

    int nc_idx = find_nc(from_node);
    if (nc_idx >= 0) nc_counts[nc_idx]++;

    edge_idx++;

    if (!is_leaf[i]) {
      from_node = to_node;
      from_id = to_id;
      ns_all = ns_part;
    } else if (i != n - 1) {
      // Walk back up until we find a node with < 2 children
      int cur_nc = find_nc(from_node);
      while (cur_nc >= 0 && nc_counts[cur_nc] == 2) {
        from_node = find_parent(from_node, edge_idx);
        from_id = find_parent_id(from_id, edge_idx);
        ns_all = ns_part;
        cur_nc = find_nc(from_node);
      }
    }
  }

  // Trim to actual number of edges
  if (edge_idx < max_edges) {
    from_vec = from_vec[Range(0, edge_idx - 1)];
    to_vec = to_vec[Range(0, edge_idx - 1)];
    from_id_vec = from_id_vec[Range(0, edge_idx - 1)];
    to_id_vec = to_id_vec[Range(0, edge_idx - 1)];
    inv_d_vec = inv_d_vec[Range(0, edge_idx - 1)];
    edge_vec = edge_vec[Range(0, edge_idx - 1)];
    nodesize_vec = nodesize_vec[Range(0, edge_idx - 1)];
  }

  return List::create(
    Named("from") = from_vec,
    Named("to") = to_vec,
    Named("from_id") = from_id_vec,
    Named("to_id") = to_id_vec,
    Named("inv_d") = inv_d_vec,
    Named("edge") = edge_vec,
    Named("nodesize") = nodesize_vec
  );
}


// ---------- compute_split_stats_cpp ----------
// Replaces the inner loop of get_Y_imp: for each node, compute
// the between-split statistic for every response column.
//
// For a given node with left/right membership:
//   stat_j = (sum(Yt[L, j]))^2 / nL + (sum(Yt[R, j]))^2 / nR
//
// where Yt is the column-centered and scaled Y matrix for the node samples.
//
// Inputs:
//   Y         : full response matrix (n_samples x q)
//   membership: tree membership vector (length n_samples), integer node IDs
//   mem_left  : integer, membership ID for left child
//   mem_right : integer, membership ID for right child
//
// Returns a NumericVector of length q with the split statistics.
// [[Rcpp::export]]
NumericVector compute_split_stats_cpp(NumericMatrix Y,
                                      IntegerVector membership,
                                      int mem_left,
                                      int mem_right) {

  int n = Y.nrow();
  int q = Y.ncol();

  // Find indices belonging to this node (left or right child)
  std::vector<int> idx_all;
  std::vector<bool> is_left_flag;
  idx_all.reserve(n);
  is_left_flag.reserve(n);

  for (int i = 0; i < n; i++) {
    if (membership[i] == mem_left || membership[i] == mem_right) {
      idx_all.push_back(i);
      is_left_flag.push_back(membership[i] == mem_left);
    }
  }

  int n_node = (int)idx_all.size();
  if (n_node < 2) {
    return NumericVector(q, 0.0);
  }

  int nL = 0, nR = 0;
  for (size_t i = 0; i < is_left_flag.size(); i++) {
    if (is_left_flag[i]) nL++; else nR++;
  }
  if (nL == 0 || nR == 0) {
    return NumericVector(q, 0.0);
  }

  NumericVector stats(q);

  for (int j = 0; j < q; j++) {
    // Column-center and scale (like R's scale())
    double col_mean = 0.0;
    for (int k = 0; k < n_node; k++) {
      col_mean += Y(idx_all[k], j);
    }
    col_mean /= n_node;

    double col_ss = 0.0;
    for (int k = 0; k < n_node; k++) {
      double v = Y(idx_all[k], j) - col_mean;
      col_ss += v * v;
    }
    double col_sd = (n_node > 1 && col_ss > 0.0) ? std::sqrt(col_ss / (n_node - 1)) : 1.0;

    // Compute left/right sums on scaled data
    double sum_L = 0.0, sum_R = 0.0;
    for (int k = 0; k < n_node; k++) {
      double scaled = (Y(idx_all[k], j) - col_mean) / col_sd;
      if (is_left_flag[k]) {
        sum_L += scaled;
      } else {
        sum_R += scaled;
      }
    }

    stats[j] = (sum_L * sum_L) / nL + (sum_R * sum_R) / nR;
  }

  return stats;
}
