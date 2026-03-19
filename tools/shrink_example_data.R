#!/usr/bin/env Rscript
# Shrink tcga_brca to gene=200, methy=200, mirna=100 (top-variance features)
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

save(tcga_brca, tcga_brca_clinical, file = rda_path, compress = "xz")

fsize <- file.info(rda_path)$size / 1024 / 1024
cat("\nSaved:", rda_path, sprintf("(%.1f MB)\n", fsize))
