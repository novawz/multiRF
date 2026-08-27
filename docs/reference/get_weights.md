# Get IMD weights

Extract per-block feature importance weights.

## Usage

``` r
get_weights(x, cluster = NULL, ...)
```

## Arguments

- x:

  An `mrf3_fit`, `mrf3`, or `cluster_imd` object.

- cluster:

  For `cluster_imd` objects, which cluster to extract. If `NULL`,
  returns all clusters as a named list of weight lists.

- ...:

  Additional arguments (unused).

## Value

A named list of numeric vectors (one per block), or `NULL`.
