# Filter multi-omics feature blocks before RF fitting

Filter multi-omics feature blocks before RF fitting

## Usage

``` r
filter_omics(
  dat.list,
  filter_mode = c("auto", "none", "manual"),
  filter_method = c("mad", "variance"),
  top_n_by_type = NULL,
  top_n_manual = NULL,
  return_summary = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- dat.list:

  A named list containing multi-omics datasets with samples in rows and
  features in columns.

- filter_mode:

  Feature filtering mode. `"auto"` applies omics-aware defaults,
  `"none"` keeps all features, and `"manual"` uses `top_n_manual`.

- filter_method:

  Feature dispersion metric used for ranking. Choose from `"mad"` or
  `"variance"`.

- top_n_by_type:

  Optional named numeric/list object used in `"auto"` mode. Supported
  names are `rna`, `mirna`, `methylation`, `protein`, and `unknown`.
  Defaults are RNA = 5000, miRNA = adaptive 500-1000, methylation =
  adaptive 5000-20000, protein = all, unknown = all.

- top_n_manual:

  Optional named or positional numeric/list object used in `"manual"`
  mode. Names should match `dat.list` block names.

- return_summary:

  Logical; whether to return filtering summary metadata.

- verbose:

  Logical; whether to print per-block filtering messages.

- ...:

  Additional arguments (currently ignored).

## Value

If `return_summary = FALSE`, returns filtered `dat.list`. Otherwise
returns a list with `dat_filtered` and `filter_summary`.

## Details

Before feature ranking, all blocks are checked to have the same samples
in the same order, unique sample names, numeric finite values, and at
least two samples. Zero-variance features are removed before optional
scaling or RF fitting because they cannot define a split and would
become non-finite under z-standardization. Column names are preserved
verbatim.
