# Preprocess a weight matrix for downstream steps

Preprocess a weight matrix for downstream steps

## Usage

``` r
prepare_weight_matrix(
  W,
  adjust = TRUE,
  top_v = NULL,
  row_normalize = TRUE,
  zero_diag = TRUE,
  eps = 1e-08,
  keep_ties = TRUE
)

adjust_weight_matrix(W, zero_diag = TRUE, eps = 1e-08)

truncate_top_v_rows(W, top_v = NULL, keep_ties = TRUE, sparse = FALSE)

row_normalize_weights(W, eps = 1e-12)
```

## Arguments

- W:

  A square numeric weight matrix.

- adjust:

  Logical; whether to apply adjusted-weight scaling.

- top_v:

  Optional integer. If set, each row keeps only the top-v entries.

- row_normalize:

  Logical; whether to row-normalize after adjustment/truncation.

- zero_diag:

  Logical; whether to set diagonal to zero.

- eps:

  Small positive value used to avoid division by zero.

- keep_ties:

  Logical; whether top-v truncation keeps ties at the cutoff.

- sparse:

  Logical; whether `truncate_top_v_rows()` should return a sparse
  matrix.

## Value

A processed square numeric matrix.

Adjusted weight matrix using row-wise scaling by `1 - diag(W)`.

Top-v truncated matrix (row-wise).

Row-normalized weight matrix.
