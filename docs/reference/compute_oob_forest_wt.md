# Compute OOB forest weight matrix

For each sample, only trees where it was out-of-bag contribute to the
weight row. Within each OOB tree, co-leaf samples contribute weight
proportional to their bootstrap frequency, normalized by the total
bootstrap mass in the leaf.

## Usage

``` r
compute_oob_forest_wt(mod)
```

## Arguments

- mod:

  A fitted model object from `fit_forest` (must contain `$membership`
  and `$inbag`), or any model carrying a pre-computed `$forest.wt.oob`
  matrix (e.g. from `fit_sub_mrf`), which is returned directly.

## Value

An n x n numeric matrix of OOB forest weights. Rows are `NA` for samples
that were not out-of-bag in any tree.
