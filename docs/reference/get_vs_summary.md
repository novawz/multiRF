# Get variable selection summary

Returns a tidy data frame summarizing how many features were selected in
each block.

## Usage

``` r
get_vs_summary(x, ...)
```

## Arguments

- x:

  An `mrf3_fit` or `cluster_imd` object.

- ...:

  Additional arguments (unused).

## Value

A data frame with columns `block`, `p_total`, `p_selected`,
`selected_ratio`. For `cluster_imd`, an extra `cluster` column.
