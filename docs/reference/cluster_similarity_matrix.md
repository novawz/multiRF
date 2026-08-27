# Cluster A Similarity Matrix

Cluster A Similarity Matrix

## Usage

``` r
cluster_similarity_matrix(
  S,
  k = NULL,
  method = c("Spectral", "PAM"),
  tune_method = "silhouette",
  gap_w = "uniform",
  ...
)
```

## Arguments

- S:

  A square similarity matrix.

- k:

  Optional integer. If `NULL`, the function tunes `k`.

- method:

  Clustering backend: `"Spectral"` or `"PAM"`.

- tune_method:

  Tuning criterion for PAM when `k` is `NULL`.

- gap_w:

  Weighting scheme for spectral eigengap when `k` is `NULL`.

- ...:

  Additional arguments passed to clustering backends.

## Value

A list with `cl`, `cl_mod`, `k`, `method`, and `embed`.
