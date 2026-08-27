# Get data from a multiRF object

Extract the data list attached to a pipeline object.

## Usage

``` r
get_data(x, which = "filtered", cluster = NULL, ...)
```

## Arguments

- x:

  An `mrf3_fit`, `mrf3`, or `cluster_imd` object.

- which:

  For `mrf3_fit`: `"filtered"` (default) or `"selected"`
  (post-variable-selection).

- cluster:

  For `cluster_imd`, which cluster.

- ...:

  Additional arguments (unused).

## Value

A named list of data frames, or `NULL`.
