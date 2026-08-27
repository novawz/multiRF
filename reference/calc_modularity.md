# Modularity of forest weight matrix

Removes self influence, row-normalizes, symmetrizes, and computes
weighted modularity. Top-v truncation is intentionally not applied here
because the modularity scores are needed before the shared top-v value
is tuned. via Louvain community detection. Higher modularity indicates
clearer block / cluster structure in the proximity graph.

## Usage

``` r
calc_modularity(fw, seed = 529L)
```

## Arguments

- fw:

  A raw forest weight matrix (n x n).

## Value

A single numeric modularity value (typically 0 to ~0.8).
