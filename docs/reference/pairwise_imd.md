# Pairwise IMD Analysis After Model Fitting

Run pairwise IMD as a standalone post-fit analysis. By default, this
uses the selected features from the variable-selection /
robust-clustering stage when available, rather than all variables.

## Usage

``` r
pairwise_imd(
  x,
  feature_source = c("selected", "weights_nonzero", "all"),
  normalized = FALSE
)
```

## Arguments

- x:

  An `mrf3_fit` or `mrf3` object containing `net`, `connection`, and IMD
  weights.

- feature_source:

  Which feature set to use: `"selected"` (default), `"weights_nonzero"`,
  or `"all"`.

- normalized:

  Logical; whether to return degree-normalized adjacency.

## Value

A list with `adj_var_mat`, `adj_dat_mat`, `var_use`, and
`feature_source`.

## Details

Pairwise IMD is only defined for two-block connections; fits with
single-block (self) connections are rejected with an error. Feature
names must also be unique across omics blocks.
