library(pkgload)

load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
load("data/tcga_brca_data.rda")

dir.create("docs/assets/img", recursive = TRUE, showWarnings = FALSE)

fit <- mrf3(
  tcga_brca,
  k = 4,
  ntree = 50,
  filter_mode = "none",
  run_imd = TRUE,
  seed = 529
)

fit_similarity <- mrf3(
  tcga_brca,
  k = 4,
  main_clustering = "similarity",
  ntree = 50,
  seed = 529
)

fit_proximity <- mrf3(
  tcga_brca,
  k = 4,
  main_clustering = "proximity",
  ntree = 50,
  seed = 529
)

fit_enhanced <- mrf3(
  tcga_brca,
  k = 4,
  main_clustering = "enhanced_proximity",
  ntree = 50,
  seed = 529
)

fit_robust <- mrf3_fit(
  dat.list = tcga_brca,
  k = 4,
  ntree = 50,
  run_imd = TRUE,
  run_variable_selection = TRUE,
  run_cluster_imd = FALSE,
  run_robust_clustering = TRUE,
  variable_selection_args = list(method = "mixture"),
  seed = 529
)

fit_sub <- mrf3_fit(
  dat.list = tcga_brca,
  k = 4,
  ntree = 200,
  sub_mrf = TRUE,
  sub_mrf_args = list(
    n_sub = 15,
    frac_response = 0.2,
    frac_predictor = 0.2,
    ntree_per_sub = 10
  ),
  run_imd = TRUE,
  seed = 529
)

common_ids <- intersect(names(get_clusters(fit)), tcga_brca_clinical$sampleID)
pam50 <- tcga_brca_clinical$BRCA_Subtype_PAM50[
  match(common_ids, tcga_brca_clinical$sampleID)
]
keep <- !is.na(pam50) & nzchar(pam50)
common_ids <- common_ids[keep]
pam50 <- pam50[keep]
robust_clusters <- get_clusters(fit_robust, which = "robust")

metrics_tbl <- rbind(
  similarity = cluster_metrics(get_clusters(fit_similarity)[common_ids], pam50),
  proximity = cluster_metrics(get_clusters(fit_proximity)[common_ids], pam50),
  enhanced_proximity = cluster_metrics(get_clusters(fit_enhanced)[common_ids], pam50),
  robust = cluster_metrics(robust_clusters[common_ids], pam50),
  sub_mrf = cluster_metrics(get_clusters(fit_sub)[common_ids], pam50)
)

metrics_out <- data.frame(method = rownames(metrics_tbl), metrics_tbl, row.names = NULL)
write.csv(metrics_out, "docs/assets/img/metrics_summary.csv", row.names = FALSE)
capture.output(metrics_tbl, file = "docs/assets/img/method_metrics.txt")
capture.output(summary(fit_robust), file = "docs/assets/img/robust_summary.txt")

cat("Generated docs/assets/img/metrics_summary.csv\n")
cat("Generated docs/assets/img/method_metrics.txt\n")
cat("Generated docs/assets/img/robust_summary.txt\n")
