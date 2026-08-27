# Compute multiple clustering metrics at once

Compute multiple clustering metrics at once

## Usage

``` r
cluster_metrics(pred, ref, na.rm = TRUE, as_tibble = FALSE)
```

## Arguments

- pred:

  Predicted cluster labels.

- ref:

  Reference cluster labels.

- na.rm:

  Logical; whether to remove pairs with missing labels. When `FALSE`,
  the table-based metrics (ARI, NMI, purity) still silently drop pairs
  with missing labels, whereas the partition Jaccard index propagates
  them and returns `NA`.

- as_tibble:

  Logical; whether to return a tibble.

## Value

A one-row data frame/tibble with columns: `n`, `ari`, `jaccard`, `nmi`,
`purity`.
