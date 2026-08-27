# Compute cluster-weighted IMD from pre-computed per-node split scores

Aggregates the per-node split statistics that the native C++ engine
stores on every fitted forest (`tree_info[[t]]$imd_x_score` and
`tree_info[[t]]$imd_y_stats`) into cluster-specific variable weights,
without refitting the forest and without re-traversing trees per sample
in R.

## Usage

``` r
cluster_weighted_imd(mod, cluster, normalized = TRUE)
```

## Arguments

- mod:

  A single fitted forest model from the native engine (with
  `$tree_info`, `$membership`, `$xvar`, `$yvar`). `$membership` is
  expected in the remapped form stored by
  [`fit_mv_forest()`](fit_mv_forest.md): sequential 1-based DFS leaf IDs
  per tree (this function inverts that remap internally to recover raw
  node indices).

- cluster:

  Named character/factor vector of cluster labels for ALL samples. When
  named, it is aligned to `rownames(mod$xvar)`.

- normalized:

  Logical; L2-normalize the output weights per side.

## Value

A named list: one element per cluster label, each containing
`list(X = named_vector, Y = named_vector)` of cluster-weighted scores.
`Y` has length zero for unsupervised (single-block) models.

## Estimator

For cluster \\c\\ and predictor \\j\\, the weight is the
cluster-occupancy-weighted mean of node split scores: \$\$w_j(c) =
\frac{\sum\_{v : \mathrm{split}(v) = j} s_v \\ n_c(v)} {\sum\_{v :
\mathrm{split}(v) = j} n_c(v)}\$\$ where \\s_v\\ is the (standardized
multivariate) split score recorded at internal node \\v\\ during fitting
and \\n_c(v)\\ is the number of cluster-\\c\\ samples whose root-to-leaf
path passes through \\v\\. The response-side weight replaces \\s_v\\
with the per-response component `imd_y_stats[j, v]` (nodes where that
component is zero are excluded from both numerator and denominator).

This is a *split-score-weighted* descriptive estimator. It is **not**
the Eq. 6-8 inverse-minimal-depth (IMD) statistic computed by
[`get_multi_weights()`](get_multi_weights.md): depth-based IMD scores a
variable by how close to the root it splits, while this estimator scores
it by the size of the variance reduction at its split nodes, localized
to the cluster through path occupancy. Because the node scores \\s_v\\
are global statistics (computed from all inbag samples during fitting)
and occupancy of the high-score shallow nodes is similar for every
cluster, the weights stay close to the global importance ranking: in
benchmarks (synthetic with planted cluster-specific drivers, and
TCGA-BRCA subtypes) per-cluster weights correlate at Spearman 0.9-0.99
*across* clusters, versus 0.6-0.8 for the depth-based path, and agree
with the depth-based per-cluster weights at Spearman ~0.4-0.6 (top-20
overlap ~50-75%). Treat the output as a fast, deterministic descriptive
screen with a mild cluster tilt – not as a substitute for the
depth-based cluster IMD.
