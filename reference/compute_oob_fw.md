# Compute OOB forest weights from membership and inbag info

Given an rfsrc model fitted with `membership = TRUE`, reconstructs the
out-of-bag forest weight matrix. For each sample i, only trees where i
is OOB contribute to its weight row. The target sample is excluded from
its own donor row (zero diagonal) and rows are renormalized, matching
`compute_oob_forest_wt_cpp()` semantics.

## Usage

``` r
compute_oob_fw(mod)
```

## Arguments

- mod:

  A fitted rfsrc model with `membership` and `inbag` slots.

## Value

An n x n OOB forest-weight matrix.
