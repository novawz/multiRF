# Tune the number of clusters

Tune the number of clusters

## Usage

``` r
tune_k_clusters(x, ...)

# Default S3 method
tune_k_clusters(
  x,
  return_cluster = FALSE,
  plot_k = FALSE,
  method = "Spectral",
  tune_method = "silhouette",
  gap_w = "uniform",
  prox = FALSE,
  ...
)

# S3 method for class 'mrf3'
tune_k_clusters(
  x,
  return_cluster = FALSE,
  plot_k = FALSE,
  method = "Spectral",
  tune_method = "silhouette",
  gap_w = "uniform",
  prox = FALSE,
  ...
)

spectral_cl(x, k_tune = seq(2, 12, by = 1), gap_w = "uniform", d = NULL, ...)

pam_cl(
  x,
  k_tune = seq(2, 9, by = 1),
  diss = TRUE,
  tune_method = "silhouette",
  ...
)
```

## Arguments

- x:

  A proximity/similarity matrix (or embedding input) for clustering.

- ...:

  Additional arguments passed to lower-level clustering functions.

- return_cluster:

  Logical; return the complete tuning result instead of only the
  selected number of clusters.

- plot_k:

  Logical; draw the tuning criterion over the evaluated values of `k`.

- method:

  Clustering backend: `"PAM"` or `"Spectral"`.

- tune_method:

  Tuning criterion for PAM (`"silhouette"` or `"ratio"`). The `"ratio"`
  criterion requires dissimilarities scaled to `[0, 1]`.

- gap_w:

  Weighting scheme for spectral eigengap (`"uniform"` or `"log"`).

- prox:

  Logical; whether `x` is a proximity matrix (converted to a
  dissimilarity matrix for PAM).

- k_tune:

  Candidate cluster counts to evaluate.

- d:

  Optional degree vector used by spectral clustering.

- diss:

  logical flag: if TRUE (default for `dist` or `dissimilarity` objects),
  then `x` will be considered as a dissimilarity matrix. If FALSE, then
  `x` will be considered as a matrix of observations by variables.

## Value

The selected number of clusters, or the complete clustering result when
`return_cluster = TRUE`.
