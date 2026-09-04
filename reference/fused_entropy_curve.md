# Entropy of row-wise top-v truncated weights for every truncation level

Computes, in closed form, the mean normalized row entropy that
`calc_fused_weight_entropy(postprocess_fused_weight(W, top_v = v))`
would return for every `v` from 1 to `ncol(W)`, plus the no-truncation
value. Each row is sorted once; the entropy of its top-`v` renormalized
weights is `log(S_v) - T_v / S_v` with `S_v` the cumulative weight and
`T_v` the cumulative `w * log(w)`, so the whole curve costs one sort per
row.

## Usage

``` r
fused_entropy_curve(W, keep_ties = TRUE, eps = 1e-12)
```

## Arguments

- W:

  Square numeric weight matrix (rows are samples).

- keep_ties:

  Logical; whether truncation keeps ties at the cutoff, as
  [`truncate_top_v_rows()`](prepare_weight_matrix.md) does by default.

- eps:

  Weights at or below this value are treated as zero.

## Value

A list with `v` (integer vector `1:ncol(W)`), `entropy` (mean normalized
row entropy at each `v`), and `entropy_inf` (no truncation).
