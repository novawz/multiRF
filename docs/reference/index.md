# Package index

## Forest fitting

Fit directed multivariate random forests across omics blocks, discover
optimal connections, and compute out-of-bag quantities.

- [`fit_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
  [`fit_multi_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
  [`fit_rfsrc()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
  [`fit_multi_rfsrc()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
  : Fit random forest model and create model lists
- [`fit_sub_mrf()`](https://novawz.github.io/multiRF/reference/fit_sub_mrf.md)
  : Fit an ensemble of sub-sampled MRF models for a single connection
- [`fit_sub_multi_rfsrc()`](https://novawz.github.io/multiRF/reference/fit_sub_multi_rfsrc.md)
  : Fit sub-MRF ensemble for a multi-omics connection list
- [`find_connection()`](https://novawz.github.io/multiRF/reference/find_connection.md)
  : Find optimal directional connections from fitted RF models
- [`filter_omics()`](https://novawz.github.io/multiRF/reference/filter_omics.md)
  : Filter multi-omics feature blocks before RF fitting
- [`compute_oob_forest_wt()`](https://novawz.github.io/multiRF/reference/compute_oob_forest_wt.md)
  : Compute OOB forest weight matrix
- [`get_oob_nmse()`](https://novawz.github.io/multiRF/reference/get_oob_nmse.md)
  : OOB normalized MSE for a fitted forest model

## The mrf3 workflow

High-level entry points for the staged mrf3 pipeline, from initial fit
to the full analysis.

- [`mrf3()`](https://novawz.github.io/multiRF/reference/mrf3.md) :
  Simplified mrf3 Entry Point

- [`mrf3_init()`](https://novawz.github.io/multiRF/reference/mrf3_init.md)
  : Fit an initial MRF model

- [`mrf3_fit()`](https://novawz.github.io/multiRF/reference/mrf3_fit.md)
  : Stage-based multiRF fitting pipeline

- [`print(`*`<mrf3>`*`)`](https://novawz.github.io/multiRF/reference/print.mrf3.md)
  :

  Print method for `mrf3` objects

- [`print(`*`<mrf3_fit>`*`)`](https://novawz.github.io/multiRF/reference/print.mrf3_fit.md)
  :

  Print `mrf3_fit` object in compact form

## Clustering

Cluster samples from learned similarity, tune the number of clusters,
and score agreement between clusterings.

- [`mrf3_cl_prox()`](https://novawz.github.io/multiRF/reference/mrf3_cl_prox.md)
  : MRF unsupervised clustering – Proximity method
- [`mrf3_specific_clustering()`](https://novawz.github.io/multiRF/reference/mrf3_specific_clustering.md)
  : Joint Shared/Specific Similarity Clustering Wrapper
- [`run_cluster_pipeline()`](https://novawz.github.io/multiRF/reference/run_cluster_pipeline.md)
  : Run Shared/Specific Clustering Pipeline
- [`tune_k_clusters()`](https://novawz.github.io/multiRF/reference/tune_k_clusters.md)
  [`spectral_cl()`](https://novawz.github.io/multiRF/reference/tune_k_clusters.md)
  [`pam_cl()`](https://novawz.github.io/multiRF/reference/tune_k_clusters.md)
  : Tune the number of clusters
- [`get_clusters()`](https://novawz.github.io/multiRF/reference/get_clusters.md)
  : Get cluster labels
- [`cluster_similarity_matrix()`](https://novawz.github.io/multiRF/reference/cluster_similarity_matrix.md)
  : Cluster A Similarity Matrix
- [`cluster_shared_similarity()`](https://novawz.github.io/multiRF/reference/cluster_shared_similarity.md)
  : Cluster Shared Similarity From Reconstruction Weights
- [`cluster_specific_similarity()`](https://novawz.github.io/multiRF/reference/cluster_specific_similarity.md)
  : Cluster Specific Similarities Per Omics Block
- [`cluster_metrics()`](https://novawz.github.io/multiRF/reference/cluster_metrics.md)
  : Compute multiple clustering metrics at once
- [`cluster_ari()`](https://novawz.github.io/multiRF/reference/cluster_ari.md)
  [`cluster_jaccard()`](https://novawz.github.io/multiRF/reference/cluster_ari.md)
  [`cluster_nmi()`](https://novawz.github.io/multiRF/reference/cluster_ari.md)
  [`cluster_purity()`](https://novawz.github.io/multiRF/reference/cluster_ari.md)
  : Clustering Metrics
- [`cluster_metric_matrix()`](https://novawz.github.io/multiRF/reference/cluster_metric_matrix.md)
  : Pairwise metric matrix across multiple clusterings

## Variable selection and tuning

Select informative features per omics block and tune top-v truncation
parameters.

- [`mrf3_vs()`](https://novawz.github.io/multiRF/reference/mrf3_vs.md) :
  MRF variable selection
- [`get_selected_vars()`](https://novawz.github.io/multiRF/reference/get_selected_vars.md)
  : Get selected variable names
- [`get_vs_summary()`](https://novawz.github.io/multiRF/reference/get_vs_summary.md)
  : Get variable selection summary
- [`get_top_vars()`](https://novawz.github.io/multiRF/reference/get_top_vars.md)
  : Get top weighted variables per block
- [`select_top_v_neff()`](https://novawz.github.io/multiRF/reference/select_top_v_neff.md)
  : Select top-v from effective neighbourhood size
- [`tune_model_top_v()`](https://novawz.github.io/multiRF/reference/tune_model_top_v.md)
  : Tune model-level top-v
- [`tune_fused_top_v()`](https://novawz.github.io/multiRF/reference/tune_fused_top_v.md)
  : Tune fused weight top-v truncation

## Weights, IMD and reconstruction

Forest weights, intersample measure of dissimilarity (IMD),
reconstruction of shared and omics-specific structure, and related print
methods.

- [`get_weights()`](https://novawz.github.io/multiRF/reference/get_weights.md)
  : Get IMD weights

- [`get_multi_weights()`](https://novawz.github.io/multiRF/reference/get_multi_weights.md)
  : Get multi-omics weights

- [`get_imp_forest()`](https://novawz.github.io/multiRF/reference/get_imp_forest.md)
  : Get forest importance

- [`prepare_weight_matrix()`](https://novawz.github.io/multiRF/reference/prepare_weight_matrix.md)
  [`adjust_weight_matrix()`](https://novawz.github.io/multiRF/reference/prepare_weight_matrix.md)
  [`truncate_top_v_rows()`](https://novawz.github.io/multiRF/reference/prepare_weight_matrix.md)
  [`row_normalize_weights()`](https://novawz.github.io/multiRF/reference/prepare_weight_matrix.md)
  : Preprocess a weight matrix for downstream steps

- [`build_similarity_from_weights()`](https://novawz.github.io/multiRF/reference/build_similarity_from_weights.md)
  : Build Similarity Matrix From Forest Weights

- [`mrf3_reconstr()`](https://novawz.github.io/multiRF/reference/mrf3_reconstr.md)
  [`get_reconstr_matrix()`](https://novawz.github.io/multiRF/reference/mrf3_reconstr.md)
  : MRF unsupervised clustering – Reconstruction method

- [`get_shared_specific_weights()`](https://novawz.github.io/multiRF/reference/get_shared_specific_weights.md)
  : Build Shared/Specific Weights From Reconstruction

- [`get_shared_frac()`](https://novawz.github.io/multiRF/reference/get_shared_frac.md)
  : Compute Shared Fraction Per Omics Block

- [`pairwise_imd()`](https://novawz.github.io/multiRF/reference/pairwise_imd.md)
  : Pairwise IMD Analysis After Model Fitting

- [`cluster_imd()`](https://novawz.github.io/multiRF/reference/cluster_imd.md)
  : Compute Cluster-Specific IMD Without Pairwise IMD

- [`print(`*`<pairwise_imd_analysis>`*`)`](https://novawz.github.io/multiRF/reference/print.pairwise_imd_analysis.md)
  :

  Print `pairwise_imd_analysis`

- [`print(`*`<cluster_imd>`*`)`](https://novawz.github.io/multiRF/reference/print.cluster_imd.md)
  :

  Print `cluster_imd`

## Stability and summaries

Stability evaluation of a fitted pipeline, summaries, and tidy accessors
for downstream analysis.

- [`mrf3_stability()`](https://novawz.github.io/multiRF/reference/mrf3_stability.md)
  : mrf3 Stability Evaluation (No Refit)

- [`print(`*`<mrf3_stability>`*`)`](https://novawz.github.io/multiRF/reference/print.mrf3_stability.md)
  :

  Print `mrf3_stability`

- [`summary(`*`<mrf3_fit>`*`)`](https://novawz.github.io/multiRF/reference/summary.mrf3_fit.md)
  :

  Summarize `mrf3_fit` output

- [`print(`*`<summary.mrf3_fit>`*`)`](https://novawz.github.io/multiRF/reference/print.summary.mrf3_fit.md)
  : Print mrf3 summary

## Visualization

Embeddings and plots for similarity matrices, clusters, weights, and
survival outcomes.

- [`plot(`*`<mrf3_fit>`*`)`](https://novawz.github.io/multiRF/reference/plot.mrf3_fit.md)
  :

  Plot an `mrf3_fit` object

- [`mrf3_tsne()`](https://novawz.github.io/multiRF/reference/mrf3_tsne.md)
  : compute tSNE embedding for MRF similarity matrix

- [`plot_tsne()`](https://novawz.github.io/multiRF/reference/plot_tsne.md)
  [`plot_network()`](https://novawz.github.io/multiRF/reference/plot_tsne.md)
  [`plot_weights()`](https://novawz.github.io/multiRF/reference/plot_tsne.md)
  [`plot_circos()`](https://novawz.github.io/multiRF/reference/plot_tsne.md)
  [`plot_embed()`](https://novawz.github.io/multiRF/reference/plot_tsne.md)
  [`plot_umap()`](https://novawz.github.io/multiRF/reference/plot_tsne.md)
  : Multiple plot functions

- [`plot_cluster_composition()`](https://novawz.github.io/multiRF/reference/plot_cluster_composition.md)
  : Plot cluster composition as a heatmap

- [`plot_km()`](https://novawz.github.io/multiRF/reference/plot_km.md) :
  Kaplan Meier plot

## Data

Bundled TCGA BRCA example data and multi-omics simulation.

- [`tcga_brca`](https://novawz.github.io/multiRF/reference/tcga_brca.md)
  : TCGA BRCA Expression Data
- [`tcga_brca_clinical`](https://novawz.github.io/multiRF/reference/tcga_brca_clinical.md)
  : TCGA BRCA Clinical Data
- [`simulate_intersim_shared_specific()`](https://novawz.github.io/multiRF/reference/simulate_intersim_shared_specific.md)
  : Simulate Shared + Specific Multi-Omics Data via InterSIM

## Utilities

Accessors for components of fitted multiRF objects.

- [`get_data()`](https://novawz.github.io/multiRF/reference/get_data.md)
  : Get data from a multiRF object
- [`get_models()`](https://novawz.github.io/multiRF/reference/get_models.md)
  : Get fitted RF models
- [`get_connection()`](https://novawz.github.io/multiRF/reference/get_connection.md)
  : Get connection list

## Deprecated

Backward-compatible aliases retained from earlier releases; use the
replacements noted on each page.

- [`simulate_shared_specific_intersim()`](https://novawz.github.io/multiRF/reference/simulate_shared_specific_intersim.md)
  :

  Backward-compatible alias of
  [`simulate_intersim_shared_specific()`](https://novawz.github.io/multiRF/reference/simulate_intersim_shared_specific.md)

- [`fit_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
  [`fit_multi_forest()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
  [`fit_rfsrc()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
  [`fit_multi_rfsrc()`](https://novawz.github.io/multiRF/reference/fit_forest.md)
  : Fit random forest model and create model lists
