# Resolve a tuning parameter (mtry or ytry) from integer, formula string, or NULL

Accepts an integer (used as-is), a formula string like `"sqrt(p)"`,
`"p/3"`, `"p/2"`, or `NULL` (returns `default`). In the formula, `p`
refers to the number of columns (px for mtry, qy for ytry).

## Usage

``` r
resolve_param(value, p, default, name = "param")
```

## Arguments

- value:

  Integer, character formula, or `NULL`.

- p:

  Number of columns to substitute into the formula.

- default:

  Value to return when `value` is `NULL`.

- name:

  Parameter name for error messages (e.g. "mtry", "ytry").

## Value

An integer value (\>= 1).
