# Compute effective neighbourhood size for each row of a weight matrix

For each row \\w_i\\, the effective neighbourhood size is defined as
\\n\_{\mathrm{eff},i} = \exp(H(w_i))\\, where \\H(w_i) = -\sum_j w\_{ij}
\log w\_{ij}\\ is the Shannon entropy. This equals the number of
equally-weighted neighbours that would produce the same entropy.

## Usage

``` r
effective_neighbourhood_size(W, eps = 1e-12)
```

## Arguments

- W:

  A square row-stochastic weight matrix (or will be row-normalised).

- eps:

  Small positive value to avoid log(0).

## Value

A numeric vector of length `nrow(W)` containing \\n\_{\mathrm{eff},i}\\.
