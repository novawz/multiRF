# Compute Cluster-Specific IMD Without Pairwise IMD

Compute Cluster-Specific IMD Without Pairwise IMD

## Usage

``` r
cluster_imd(
  x,
  cluster = NULL,
  dat.list = NULL,
  connect_list = NULL,
  min_cluster_size = 20L,
  ntree = NULL,
  ytry = NULL,
  parallel = TRUE,
  imd_normalized_weights = FALSE,
  run_vs = FALSE,
  vs_args = list(),
  fit_args = list(),
  imd_args = list(),
  use_node_imd = FALSE,
  keep_model = FALSE,
  keep_data = FALSE,
  seed = 529
)
```

## Arguments

- x:

  A `mrf3` or `mrf3_fit` object.

- cluster:

  Optional cluster labels for samples.

- dat.list:

  Optional named omics list (samples in rows).

- connect_list:

  Optional global connection list used by all clusters.

- min_cluster_size:

  Minimum sample size required to run a cluster-specific IMD.

- ntree:

  Deprecated compatibility argument. Ignored in no-refit mode.

- ytry:

  `ytry` used for per-cluster IMD. If `NULL`, reuses object setting.
  Only used by the depth-based path.

- parallel:

  Logical; whether IMD computation inside each cluster is parallelized.
  Only used by the depth-based path (the node-score fast path is
  single-pass and needs no parallelism).

- imd_normalized_weights:

  Logical; passed to [`get_multi_weights()`](get_multi_weights.md) as
  `normalized`. The default is `FALSE` so cluster-level variable
  selection receives raw forest IMD on its 0-to-1 scale.

- run_vs:

  Logical; whether to run [`mrf3_vs()`](mrf3_vs.md) variable selection
  on each cluster after IMD. `run_vs = TRUE` forces the depth-based
  path, which builds the per-cluster subset models that
  [`mrf3_vs()`](mrf3_vs.md) requires.

- vs_args:

  A named list of additional arguments passed to
  [`mrf3_vs()`](mrf3_vs.md) for each cluster.

- fit_args:

  Deprecated compatibility argument. Ignored in no-refit mode.

- imd_args:

  Optional named list merged into each per-cluster
  [`get_multi_weights()`](get_multi_weights.md) call. Supplying any
  `imd_args` forces the depth-based path, because the node-score fast
  path cannot honor [`get_multi_weights()`](get_multi_weights.md)
  options.

- use_node_imd:

  Controls the estimator used for the per-cluster weights. `FALSE` (the
  default) computes Eq. 6-8 inverse minimal depth via the per-tree
  network traversal of [`get_multi_weights()`](get_multi_weights.md) on
  each cluster's subset model. `TRUE` uses the node-score fast path
  ([`cluster_weighted_imd()`](cluster_weighted_imd.md)): the per-node
  split statistics stored by the native engine are averaged per
  variable, weighted by the fraction of the cluster's samples whose
  root-to-leaf paths visit each node. The fast path is deterministic and
  typically 100-900x faster, but it is a *different,
  split-score-weighted* estimator whose weights stay close to the global
  importance ranking: in benchmarks its per-cluster weights correlate at
  Spearman 0.9-0.99 *across* clusters (depth path: 0.6-0.8), and agree
  with the depth-based weights at Spearman ~0.4-0.6 (top-20 overlap
  ~50-75%). Use it as a quick descriptive screen, not as a substitute
  for the depth-based cluster IMD (see
  [`cluster_weighted_imd()`](cluster_weighted_imd.md) for the
  definition). `NULL` selects the fast path automatically whenever the
  models carry the required `tree_info` statistics and no blocking
  option (`imd_args`, `run_vs`) is requested, falling back to the
  depth-based path otherwise.

- keep_model:

  Logical; whether to keep per-cluster mrf3 objects in output.

- keep_data:

  Logical; whether to keep per-cluster subset data in output.

- seed:

  Base seed; cluster `i` uses `seed + i - 1`.

## Value

A list with cluster-level summary, per-cluster IMD outputs, and params.
`params$imd_path` records which estimator produced the weights:
`"node_score"` (fast path) or `"depth_traversal"` (Eq. 6-8 IMD).
