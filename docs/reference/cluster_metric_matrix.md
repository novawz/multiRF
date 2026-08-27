# Pairwise metric matrix across multiple clusterings

Pairwise metric matrix across multiple clusterings

## Usage

``` r
cluster_metric_matrix(
  cluster_list,
  metric = c("ari", "jaccard", "nmi", "purity"),
  na.rm = TRUE
)
```

## Arguments

- cluster_list:

  Named list of cluster-label vectors. All vectors must be alignable to
  the same samples.

- metric:

  Metric to compute: `"ari"`, `"jaccard"`, `"nmi"`, or `"purity"`.

- na.rm:

  Logical; whether to remove missing labels pairwise.

## Value

A matrix of pairwise metric values. The matrix is symmetric for `"ari"`,
`"jaccard"`, and `"nmi"`; for `"purity"`, entry `[i, j]` is the purity
of clustering `i` evaluated against clustering `j` as reference, so the
matrix may be asymmetric.
