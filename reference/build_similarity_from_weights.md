# Build Similarity Matrix From Forest Weights

Build Similarity Matrix From Forest Weights

## Usage

``` r
build_similarity_from_weights(
  W,
  similarity_type = c("second", "first"),
  zero_diag = TRUE,
  symm = TRUE
)
```

## Arguments

- W:

  A square forest-weight matrix.

- similarity_type:

  Similarity structure: `"second"` uses `W %*% t(W)`; `"first"` uses
  `(W + t(W))/2`.

- zero_diag:

  Logical; whether to zero out diagonal entries in similarity.

- symm:

  Logical; whether to force symmetry by averaging with transpose.

## Value

A sample-by-sample similarity matrix.
