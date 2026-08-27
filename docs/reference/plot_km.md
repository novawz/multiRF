# Kaplan Meier plot

Kaplan Meier plot

## Usage

``` r
plot_km(
  test_var,
  time_var,
  event_var,
  pheno_mat,
  cut = c("median", "mean", "maxstat"),
  na.rm = TRUE,
  base_size = 8,
  base_family = "sans",
  ...
)
```

## Arguments

- test_var:

  Sample-level grouping variable or continuous score. Numeric vectors
  with more than five unique values are treated as continuous;
  low-cardinality numeric vectors are treated as group labels.

- time_var:

  Column name containing survival or follow-up time.

- event_var:

  Column name containing a binary 0/1 event indicator.

- pheno_mat:

  Data frame containing the time and event columns.

- cut:

  For a continuous `test_var`, split at the `"median"`, `"mean"`, or a
  maximally selected rank-statistic (`"maxstat"`) cutoff.

- na.rm:

  Logical; whether to remove rows with missing values.

- base_size:

  Base font size.

- base_family:

  Font family.

- ...:

  Additional arguments passed to
  [`survminer::ggsurvplot()`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html).

## Value

A `ggsurvplot` object.
