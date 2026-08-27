# Select top-v from effective neighbourhood size

Sets \\v = \lceil \mathrm{quantile}\_q(n\_{\mathrm{eff},i}) \rceil\\
where \\n\_{\mathrm{eff},i} = \exp(H(w_i))\\.

## Usage

``` r
select_top_v_neff(W, quantile_prob = 0.5, min_v = 10L, eps = 1e-12)
```

## Arguments

- W:

  A square weight matrix.

- quantile_prob:

  Quantile probability. Default 0.5 (median).

- min_v:

  Floor on the returned value.

- eps:

  Small positive value for entropy computation.

## Value

A single integer: the selected top-v.
