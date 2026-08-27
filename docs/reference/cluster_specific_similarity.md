# Cluster Specific Similarities Per Omics Block

Cluster Specific Similarities Per Omics Block

## Usage

``` r
cluster_specific_similarity(
  shared_specific,
  k = NULL,
  method = c("Spectral", "PAM", "Proximity", "Enhanced_Proximity"),
  prox_method_cl = c("PAM", "Spectral"),
  tune_method = "silhouette",
  gap_w = "uniform",
  similarity_type = c("second", "first"),
  ...
)
```

## Arguments

- shared_specific:

  Output from
  [`get_shared_specific_weights()`](get_shared_specific_weights.md).

- k:

  Optional specific-clustering `k` configuration. `NULL` tunes each
  block; single integer uses same `k` for all blocks; named numeric/list
  allows block-specific `k`.

- method:

  Clustering backend for each specific block. Choose from `"Spectral"`,
  `"PAM"`, `"Proximity"`, and `"Enhanced_Proximity"`.

- prox_method_cl:

  Proximity clustering backend used when `method` is `"Proximity"` or
  `"Enhanced_Proximity"`.

- tune_method:

  Tuning criterion for PAM when block-level `k` is `NULL`.

- gap_w:

  Weighting scheme for spectral eigengap when block-level `k` is `NULL`.

- similarity_type:

  Similarity structure for each specific block.

- ...:

  Additional arguments passed to clustering backends.

## Value

A list with per-omics specific similarities and clustering outputs.
