if (!"package:multiRF" %in% search() &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(export_all = TRUE, helpers = FALSE, quiet = TRUE)
}
