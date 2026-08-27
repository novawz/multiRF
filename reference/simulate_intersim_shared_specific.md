# Simulate Shared + Specific Multi-Omics Data via InterSIM

Generate multi-omics data where shared signal is provided by
[`InterSIM::InterSIM()`](https://rdrr.io/pkg/InterSIM/man/InterSIM.html)
and additional omics-specific signal is injected by per-omics latent
labels `U`.

## Usage

``` r
simulate_intersim_shared_specific(
  n = 500,
  Z_prop = c(0.3, 0.3, 0.4),
  delta_methyl = 2,
  delta_expr = 2,
  delta_protein = 2,
  p_DMP = 0.2,
  p_DEG = NULL,
  p_DEP = NULL,
  p_spec_feat = 0.05,
  delta_spec = c(methyl = 1, gene = 1, protein = 1),
  U_nclass = 2,
  U_within_Z = TRUE,
  max_features_per_block = 400,
  sample_prefix = "S",
  return_workflow_data = TRUE,
  keep_sim_raw = TRUE,
  seed = NULL
)
```

## Arguments

- n:

  Number of samples.

- Z_prop:

  Cluster proportion for InterSIM shared clusters.

- delta_methyl:

  Shared-signal strength for methylation in InterSIM.

- delta_expr:

  Shared-signal strength for expression in InterSIM.

- delta_protein:

  Shared-signal strength for protein in InterSIM.

- p_DMP:

  Proportion of differential methylation features in InterSIM. Default
  matches
  [`InterSIM::InterSIM()`](https://rdrr.io/pkg/InterSIM/man/InterSIM.html)
  (`0.2`).

- p_DEG:

  Proportion of differential expression features in InterSIM. Default
  matches
  [`InterSIM::InterSIM()`](https://rdrr.io/pkg/InterSIM/man/InterSIM.html)
  (`NULL`).

- p_DEP:

  Proportion of differential protein features in InterSIM. Default
  matches
  [`InterSIM::InterSIM()`](https://rdrr.io/pkg/InterSIM/man/InterSIM.html)
  (`NULL`).

- p_spec_feat:

  Proportion of features used as specific set `T` per omics.

- delta_spec:

  Named numeric vector with specific shift size for `methyl`, `gene`,
  and `protein`.

- U_nclass:

  Number of categories for each specific latent label `U`.

- U_within_Z:

  Logical; whether each `U` is generated within each shared cluster `Z`
  (recommended).

- max_features_per_block:

  Maximum number of features kept per omics block in returned `dat.list`
  when `return_workflow_data = TRUE`.

- sample_prefix:

  Prefix for generated sample IDs in `dat.list`.

- return_workflow_data:

  Logical; whether to return workflow-ready `dat.list`, `truth`, and
  `k_ref`.

- keep_sim_raw:

  Logical; whether to keep original InterSIM output as `sim_raw`.

- seed:

  Optional random seed. If `NULL`, current RNG state is used.

## Value

A list containing:

- `X_list`: simulated matrices (`methyl`, `gene`, `protein`)

- `Z`: shared labels from InterSIM

- `U`: specific labels per omics

- `S`: shared feature index sets per omics

- `T`: specific feature index sets per omics

- `sim_raw`: raw InterSIM output (optional)

- `dat.list`, `truth`, `k_ref`: workflow-ready outputs (optional)

## Details

This function wraps the notebook simulation workflow into a reusable
package API and optionally returns workflow-ready `dat.list` and truth
labels.
