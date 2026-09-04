# Compute Shared Fraction Per Omics Block

Shared fraction is computed by the signal-based definition:
`1 - (||R||_F^2 / ||X||_F^2)`, where `R = X - X_pred`.

## Usage

``` r
get_shared_frac(
  dat.list = NULL,
  residual = NULL,
  shared_specific = NULL,
  eps = 1e-12
)
```

## Arguments

- dat.list:

  Optional named list of original omics matrices (`X`).

- residual:

  Optional named list of residual matrices (`R`).

- shared_specific:

  Optional output from
  [`get_shared_specific_weights()`](https://novawz.github.io/multiRF/reference/get_shared_specific_weights.md).
  If provided, `residual` and (when needed) reconstructed `dat.list` are
  extracted from it.

- eps:

  Small positive threshold to avoid division by zero.

## Value

A data.frame with one row per omics block: `data`, `frob_shared`,
`frob_specific`, `specific_ratio`, `shared_frac`, and `method`.
`frob_shared` is the Frobenius norm of the reconstructed (shared) matrix
`X_pred = X - R`. For backward compatibility, it also includes
`frob_data`, `frob_residual`, and `residual_ratio`.
