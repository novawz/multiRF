# Get selected variable names

Extract variable names selected by
[`mrf3_vs()`](https://novawz.github.io/multiRF/reference/mrf3_vs.md).

## Usage

``` r
get_selected_vars(x, cluster = NULL, ...)
```

## Arguments

- x:

  An `mrf3_fit`, `mrf3` (with class `vs`), or `cluster_imd` object.

- cluster:

  For `cluster_imd` objects, which cluster. If `NULL`, returns all
  clusters as a named list.

- ...:

  Additional arguments (unused).

## Value

A named list of character vectors (variable names per block), or `NULL`
if variable selection has not been run.
