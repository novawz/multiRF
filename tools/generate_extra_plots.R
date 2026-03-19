library(pkgload)

load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
load("data/tcga_brca_data.rda")
dir.create("docs/assets/img", recursive = TRUE, showWarnings = FALSE)

fit <- mrf3(
  tcga_brca,
  k = 4,
  ntree = 30,
  filter_mode = "none",
  run_imd = TRUE,
  seed = 529
)

png("docs/assets/img/circos_wall.png", width = 1400, height = 1400, res = 160)
plot_circos(mat = fit$reconstruction$W$W_all)
dev.off()

png("docs/assets/img/heatmap_fused_gene.png", width = 1200, height = 1000, res = 160)
plot_heatmap(dat = fit, source = "reconstruction_fused_block", omics = "gene")
dev.off()

cat("Generated circos_wall.png and heatmap_fused_gene.png\n")
