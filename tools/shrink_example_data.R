#!/usr/bin/env Rscript
# Shrink tcga_brca to gene=200, methy=200, mirna=100 (top-variance features)
# and prune tcga_brca_clinical to the columns used by the package docs,
# vignette, and survival examples (keep list below).
# Run from package root: Rscript tools/shrink_example_data.R

rda_path <- "data/tcga_brca_data.rda"
load(rda_path)

cat("Before:\n")
cat("  gene: ", paste(dim(tcga_brca$gene), collapse = " x "), "\n")
cat("  methy:", paste(dim(tcga_brca$methy), collapse = " x "), "\n")
cat("  mirna:", paste(dim(tcga_brca$mirna), collapse = " x "), "\n")

top_var <- function(mat, n) {
  v <- apply(mat, 2, var, na.rm = TRUE)
  keep <- names(sort(v, decreasing = TRUE))[seq_len(min(n, ncol(mat)))]
  mat[, keep, drop = FALSE]
}

tcga_brca$gene  <- top_var(tcga_brca$gene,  200)
tcga_brca$methy <- top_var(tcga_brca$methy, 200)
tcga_brca$mirna <- top_var(tcga_brca$mirna, 100)

cat("\nAfter:\n")
cat("  gene: ", paste(dim(tcga_brca$gene), collapse = " x "), "\n")
cat("  methy:", paste(dim(tcga_brca$methy), collapse = " x "), "\n")
cat("  mirna:", paste(dim(tcga_brca$mirna), collapse = " x "), "\n")

# Prune the clinical table to the columns the package actually uses:
# - sampleID, sample_type, BRCA_Subtype_PAM50: vignette + docs examples
# - OS, OS.time, PFI, PFI.time: survival columns for plot_km() examples
# - age_at_initial_pathologic_diagnosis, pathologic_stage: common covariates
# Keep this list in sync with the @format block in R/data.R.
clinical_keep <- c(
  "sampleID", "sample_type", "BRCA_Subtype_PAM50",
  "OS", "OS.time", "PFI", "PFI.time",
  "age_at_initial_pathologic_diagnosis", "pathologic_stage"
)
cat("\nClinical before:", paste(dim(tcga_brca_clinical), collapse = " x "), "\n")
stopifnot(all(clinical_keep %in% colnames(tcga_brca_clinical)))
tcga_brca_clinical <- tcga_brca_clinical[, clinical_keep, drop = FALSE]
cat("Clinical after: ", paste(dim(tcga_brca_clinical), collapse = " x "), "\n")

save(tcga_brca, tcga_brca_clinical, file = rda_path, compress = "xz")

fsize <- file.info(rda_path)$size / 1024 / 1024
cat("\nSaved:", rda_path, sprintf("(%.1f MB)\n", fsize))
