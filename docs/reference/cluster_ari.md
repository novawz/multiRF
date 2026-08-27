# Clustering Metrics

Utility functions for comparing two partitions (predicted vs reference),
including ARI, partition-level Jaccard, NMI, and purity.

## Usage

``` r
cluster_ari(pred, ref, na.rm = TRUE)

cluster_jaccard(pred, ref, na.rm = TRUE)

cluster_nmi(pred, ref, na.rm = TRUE)

cluster_purity(pred, ref, na.rm = TRUE)
```

## Arguments

- pred:

  Predicted cluster labels.

- ref:

  Reference cluster labels.

- na.rm:

  Logical; whether to remove pairs with missing labels before computing
  metrics. When `FALSE`, the table-based metrics (ARI, NMI, purity)
  still silently drop pairs with missing labels, whereas the partition
  Jaccard index propagates them and returns `NA`.

## Value

For scalar metric functions, a numeric scalar (`NA` when undefined).
