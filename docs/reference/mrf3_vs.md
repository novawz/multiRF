# MRF variable selection

MRF variable selection

## Usage

``` r
mrf3_vs(
  mod,
  dat.list = NULL,
  method = "filter",
  signal = c("shared", "specific", "all"),
  se = 1,
  c1 = "normal",
  c2 = "normal",
  level = 0.05,
  tscore = FALSE,
  use_distribution = TRUE,
  re_weights = FALSE,
  re_fit = TRUE,
  ntree = NULL,
  scale = FALSE,
  k = 3,
  tol = 0.01,
  iter = 1000,
  eps = 1e-05,
  normalized = FALSE,
  select = "ALL",
  ...
)
```

## Arguments

- mod:

  A fitted `mrf3_fit` object, or an `mrf3`-like object that already
  contains raw forest IMD in `mod$imd`.

- dat.list:

  A named list of omics matrices used for feature selection and optional
  refitting.

- method:

  Feature-selection rule. `"filter"` adaptively tunes the cutoff
  `tau * sd(IMD)` using OOB error; `"mixture"` fits a
  point-mass/two-component model; and `"test"` (alias
  `"transformation"`) selects by the Eq. 16 t-score IMD: each feature's
  forest IMD is standardized against the mean forest IMD of all features
  in the block, using the feature's across-tree standard error, and
  features in the upper tail of a Student-t reference with `ntree - 1`
  degrees of freedom (`p < level`) are kept; see Details. `"thres"`
  applies a fixed `se * sd(IMD)` cutoff.

- signal:

  Which signal component to select on: `"shared"` (cross-modal),
  `"specific"` (per-block residual), or `"all"` (the union of both).

- se:

  The fixed value of tau used by `method = "thres"`.

- c1:

  Distribution family for the lower-IMD (noise) component in mixture
  mode: `"normal"` or `"truncn"`.

- c2:

  Distribution family for the higher-IMD component in mixture mode:
  `"normal"` or `"gamma"`.

- level:

  One-sided significance level for the transformation, or the maximum
  posterior noise probability for mixture selection. A single numeric
  value is applied to every block; a named vector/list supplies
  block-specific values. `"auto"` retains the optional Bayesian-FDR rule
  in mixture mode.

- tscore:

  Deprecated compatibility switch. Use `method = "transformation"`
  instead.

- use_distribution:

  Deprecated compatibility argument.

- re_weights:

  Logical; whether selected IMD values are supplied as variable-sampling
  weights during the final refit.

- re_fit:

  Logical; whether to refit forests after feature selection.

- ntree:

  Number of trees used in the optional final refit. `NULL` reuses the
  fitted object's tree count, falling back to 300.

- scale:

  Logical; whether to standardize selected data before refitting.

- k:

  Number of repeated forest fits at each candidate filtering cutoff.

- tol:

  Tolerable adjacent change in mean OOB normalized MSE used to choose
  the filtering cutoff.

- iter:

  Maximum EM iterations in mixture mode.

- eps:

  EM convergence tolerance in mixture mode.

- normalized:

  Logical; whether to L2-normalize the retained weights after selection.
  These methods are always evaluated on raw IMD; therefore the default
  is `FALSE`.

- select:

  Data block(s) to which selection is applied; default `"ALL"`.

- ...:

  Additional arguments passed to forest fitting.

## Value

An object inheriting from `mrf3` and `vs`. For `signal = "all"`,
`dat.list` and `imd` contain the union selection and the two component
results are retained in `signal_results`.

## Details

For `method = "transformation"`, every connected forest supplies a
per-tree IMD matrix (features x trees) per block. Feature `v` is
standardized as the t-score IMD of Eq. 16, `t_v = (M_v - mu) / SE(M_v)`,
where `M_v` is the feature's forest IMD (its per-tree IMD averaged over
the `B` trees), `mu` is the mean forest IMD over all features of the
block, and `SE(M_v)` is the feature's across-tree standard error
(`sd(per-tree IMD) / sqrt(B)`). Features with an upper-tail
`p = P(T_{B-1} > t_v) < level` are selected, i.e. features whose forest
IMD lies significantly above the block-wide mean. When several connected
forests cover the same block, each forest is evaluated separately and a
feature is retained by strict majority vote across forests; the reported
t-scores and p-values are averaged. Because per-tree IMDs within a
forest are correlated (trees share the training data) and `mu` is
treated as a fixed reference, the p-values are screening scores rather
than exact test levels.
