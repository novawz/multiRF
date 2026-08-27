# Get cluster labels

Extract cluster assignments from any multiRF pipeline object.

## Usage

``` r
get_clusters(x, which = "main", block = NULL, ...)
```

## Arguments

- x:

  An `mrf3_fit`, `mrf3`, `reconstr`, or `prox` object.

- which:

  For `mrf3_fit`: `"main"` (default) or `"robust"`.

- block:

  Optional block name. For an `mrf3_fit` main or robust result, returns
  the corresponding block-specific labels instead of the shared labels.

- ...:

  Additional arguments (unused).

## Value

A named vector of cluster assignments, or `NULL` if not available.
