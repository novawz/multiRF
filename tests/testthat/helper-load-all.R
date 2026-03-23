if (!"package:multiRF" %in% search()) {
  loaded <- FALSE

  if (requireNamespace("pkgload", quietly = TRUE)) {
    loaded <- tryCatch({
      pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
  }

  if (!loaded) {
    pkg_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
    r_files <- list.files(file.path(pkg_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)
    for (path in r_files) {
      source(path, local = FALSE, chdir = FALSE)
    }
  }
}
