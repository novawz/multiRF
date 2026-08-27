# Package index

## Fitting and engines

Fit directed multivariate random forests across omics blocks, discover
optimal connections, and compute out-of-bag quantities.

- [`fit_forest()`](fit_forest.md) [`fit_multi_forest()`](fit_forest.md)
  [`fit_rfsrc()`](fit_forest.md) [`fit_multi_rfsrc()`](fit_forest.md) :
  Fit random forest model and create model lists
- [`fit_sub_mrf()`](fit_sub_mrf.md) : Fit an ensemble of sub-sampled MRF
  models for a single connection
- [`fit_sub_multi_rfsrc()`](fit_sub_multi_rfsrc.md) : Fit sub-MRF
  ensemble for a multi-omics connection list
- [`find_connection()`](find_connection.md) : Find optimal directional
  connections from fitted RF models
- [`filter_omics()`](filter_omics.md) : Filter multi-omics feature
  blocks before RF fitting
- [`compute_oob_forest_wt()`](compute_oob_forest_wt.md) : Compute OOB
  forest weight matrix
- [`get_oob_nmse()`](get_oob_nmse.md) : OOB normalized MSE for a fitted
  forest model

## The mrf3 workflow

High-level entry points for the staged mrf3 pipeline, from initial fit
to the full analysis.

- [`mrf3()`](mrf3.md) : Simplified mrf3 Entry Point

- [`mrf3_init()`](mrf3_init.md) : Fit an initial MRF model

- [`mrf3_fit()`](mrf3_fit.md) : Stage-based multiRF fitting pipeline

- [`print(`*`<mrf3>`*`)`](print.mrf3.md) :

  Print method for `mrf3` objects

- [`print(`*`<mrf3_fit>`*`)`](print.mrf3_fit.md) :

  Print `mrf3_fit` object in compact form

## Clustering

Cluster samples from learned similarity, tune the number of clusters,
and score agreement between clusterings.

- [`mrf3_cl_prox()`](mrf3_cl_prox.md) : MRF unsupervised clustering –
  Proximity method
- [`mrf3_specific_clustering()`](mrf3_specific_clustering.md) : Joint
  Shared/Specific Similarity Clustering Wrapper
- [`run_cluster_pipeline()`](run_cluster_pipeline.md) : Run
  Shared/Specific Clustering Pipeline
- [`tune_k_clusters()`](tune_k_clusters.md)
  [`spectral_cl()`](tune_k_clusters.md) [`pam_cl()`](tune_k_clusters.md)
  : Tune the number of clusters
- [`get_clusters()`](get_clusters.md) : Get cluster labels
- [`cluster_similarity_matrix()`](cluster_similarity_matrix.md) :
  Cluster A Similarity Matrix
- [`cluster_shared_similarity()`](cluster_shared_similarity.md) :
  Cluster Shared Similarity From Reconstruction Weights
- [`cluster_specific_similarity()`](cluster_specific_similarity.md) :
  Cluster Specific Similarities Per Omics Block
- [`cluster_metrics()`](cluster_metrics.md) : Compute multiple
  clustering metrics at once
- [`cluster_ari()`](cluster_ari.md)
  [`cluster_jaccard()`](cluster_ari.md)
  [`cluster_nmi()`](cluster_ari.md) [`cluster_purity()`](cluster_ari.md)
  : Clustering Metrics
- [`cluster_metric_matrix()`](cluster_metric_matrix.md) : Pairwise
  metric matrix across multiple clusterings

## Variable selection and tuning

Select informative features per omics block and tune top-v truncation
parameters.

- [`mrf3_vs()`](mrf3_vs.md) : MRF variable selection
- [`get_selected_vars()`](get_selected_vars.md) : Get selected variable
  names
- [`get_vs_summary()`](get_vs_summary.md) : Get variable selection
  summary
- [`get_top_vars()`](get_top_vars.md) : Get top weighted variables per
  block
- [`select_top_v_neff()`](select_top_v_neff.md) : Select top-v from
  effective neighbourhood size
- [`tune_model_top_v()`](tune_model_top_v.md) : Tune model-level top-v
- [`tune_fused_top_v()`](tune_fused_top_v.md) : Tune fused weight top-v
  truncation

## Weights, IMD and reconstruction

Forest weights, intersample measure of dissimilarity (IMD),
reconstruction of shared and omics-specific structure, and related print
methods.

- [`get_weights()`](get_weights.md) : Get IMD weights

- [`get_multi_weights()`](get_multi_weights.md) : Get multi-omics
  weights

- [`get_imp_forest()`](get_imp_forest.md) : Get forest importance

- [`prepare_weight_matrix()`](prepare_weight_matrix.md)
  [`adjust_weight_matrix()`](prepare_weight_matrix.md)
  [`truncate_top_v_rows()`](prepare_weight_matrix.md)
  [`row_normalize_weights()`](prepare_weight_matrix.md) : Preprocess a
  weight matrix for downstream steps

- [`build_similarity_from_weights()`](build_similarity_from_weights.md)
  : Build Similarity Matrix From Forest Weights

- [`mrf3_reconstr()`](mrf3_reconstr.md)
  [`get_reconstr_matrix()`](mrf3_reconstr.md) : MRF unsupervised
  clustering – Reconstruction method

- [`get_shared_specific_weights()`](get_shared_specific_weights.md) :
  Build Shared/Specific Weights From Reconstruction

- [`get_shared_frac()`](get_shared_frac.md) : Compute Shared Fraction
  Per Omics Block

- [`pairwise_imd()`](pairwise_imd.md) : Pairwise IMD Analysis After
  Model Fitting

- [`cluster_imd()`](cluster_imd.md) : Compute Cluster-Specific IMD
  Without Pairwise IMD

- [`print(`*`<pairwise_imd_analysis>`*`)`](print.pairwise_imd_analysis.md)
  :

  Print `pairwise_imd_analysis`

- [`print(`*`<cluster_imd>`*`)`](print.cluster_imd.md) :

  Print `cluster_imd`

## Stability and summaries

Stability evaluation of a fitted pipeline, summaries, and tidy accessors
for downstream analysis.

- [`mrf3_stability()`](mrf3_stability.md) : mrf3 Stability Evaluation
  (No Refit)

- [`print(`*`<mrf3_stability>`*`)`](print.mrf3_stability.md) :

  Print `mrf3_stability`

- [`summary(`*`<mrf3_fit>`*`)`](summary.mrf3_fit.md) :

  Summarize `mrf3_fit` output

- [`print(`*`<summary.mrf3_fit>`*`)`](print.summary.mrf3_fit.md) : Print
  mrf3 summary

## Visualization

Embeddings and plots for similarity matrices, clusters, weights, and
survival outcomes.

- [`plot(`*`<mrf3_fit>`*`)`](plot.mrf3_fit.md) :

  Plot an `mrf3_fit` object

- [`mrf3_tsne()`](mrf3_tsne.md) : compute tSNE embedding for MRF
  similarity matrix

- [`plot_tsne()`](plot_tsne.md) [`plot_network()`](plot_tsne.md)
  [`plot_weights()`](plot_tsne.md) [`plot_circos()`](plot_tsne.md)
  [`plot_embed()`](plot_tsne.md) [`plot_umap()`](plot_tsne.md) :
  Multiple plot functions

- [`plot_cluster_composition()`](plot_cluster_composition.md) : Plot
  cluster composition as a heatmap

- [`plot_km()`](plot_km.md) : Kaplan Meier plot

## Data

Bundled TCGA BRCA example data and multi-omics simulation.

- [`tcga_brca`](tcga_brca.md) : TCGA BRCA Expression Data
- [`tcga_brca_clinical`](tcga_brca_clinical.md) : TCGA BRCA Clinical
  Data
- [`simulate_intersim_shared_specific()`](simulate_intersim_shared_specific.md)
  : Simulate Shared + Specific Multi-Omics Data via InterSIM

## Utilities

Accessors for components of fitted multiRF objects.

- [`get_data()`](get_data.md) : Get data from a multiRF object
- [`get_models()`](get_models.md) : Get fitted RF models
- [`get_connection()`](get_connection.md) : Get connection list

## Deprecated

Backward-compatible aliases retained from earlier releases; use the
replacements noted on each page.

- [`simulate_shared_specific_intersim()`](simulate_shared_specific_intersim.md)
  :

  Backward-compatible alias of
  [`simulate_intersim_shared_specific()`](../reference/simulate_intersim_shared_specific.md)

- [`fit_forest()`](fit_forest.md) [`fit_multi_forest()`](fit_forest.md)
  [`fit_rfsrc()`](fit_forest.md) [`fit_multi_rfsrc()`](fit_forest.md) :
  Fit random forest model and create model lists
