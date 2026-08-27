# Get top weighted variables per block

A convenience function that returns the top-n variables by IMD weight
for each omics block.

## Usage

``` r
get_top_vars(x, n = 20L, cluster = NULL, ...)
```

## Arguments

- x:

  Any object accepted by [`get_weights()`](get_weights.md).

- n:

  Maximum number of top variables per block.

- cluster:

  Passed to [`get_weights()`](get_weights.md).

- ...:

  Additional arguments passed to [`get_weights()`](get_weights.md).

## Value

A named list of data frames (one per block), each with columns
`variable` and `weight`, sorted by weight descending.
