# Changelog

## multiRF 0.3.0

- Added permutation-null filtering for residual-specific IMD; shared
  filtering continues to tune its cutoff with OOB error.
- Made saturation-based entropy selection the default for model and
  fused top-v, with closed-form fused entropy and the earlier methods
  retained.
- Replaced fork-based parallel paths with PSOCK workers and safer serial
  defaults for tuning grids.
- Added reconstruction fallbacks for missing proximity and response
  directions.

## multiRF 0.2.3

- Added the unified `plot(fit, type = ...)` interface and pkgdown site.
- Fixed clustering, OOB weights, IMD, stability, and related edge cases.
- Improved forest-fitting speed, parallel safety, and cross-platform
  reproducibility; seeded unsupervised results change once from 0.2.2.
- Added an optional fast
  [`cluster_imd()`](https://novawz.github.io/multiRF/reference/cluster_imd.md)
  path and simplified internals, dependencies, tests, and documentation.
- Relicensed under GPL (\>= 3) with RF-SRC provenance and citations.

## multiRF 0.2.2

- Added
  [`plot_cluster_composition()`](https://novawz.github.io/multiRF/reference/plot_cluster_composition.md)
  for auditable cluster-by-annotation composition heatmaps with counts,
  normalized percentages, and consistent package styling.

## multiRF 0.2.1

- **Fix**: Corrected the default `ytry` for supervised multivariate
  regression.

## multiRF 0.2.0

- Added multivariate regression and unsupervised forests with parallel
  fitting.
- Partition-based ranger-style tree building for faster splits.
- IMD importance computed during tree building (zero-cost).

## multiRF 0.1.1

- Prepared package metadata for publishing (`DESCRIPTION`, license, and
  dependency declarations).
- Added a practical `README.md` with installation and usage
  instructions.
- Added `NEWS.md` to track user-facing changes.
- Hardened build settings to exclude development artifacts.
- Removed compiled artifacts from `src/` for clean source distribution.
