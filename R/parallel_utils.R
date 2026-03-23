sanitize_mc_cores <- function(cores = NULL, fallback = 1L) {
  detect <- suppressWarnings(parallel::detectCores())
  if (length(detect) != 1L || is.na(detect) || !is.finite(detect) || detect < 1L) {
    detect <- fallback
  }
  if (length(cores) != 1L || is.na(cores) || !is.finite(cores) || cores < 1L) {
    cores <- fallback
  }
  as.integer(max(1L, min(as.integer(cores), as.integer(detect))))
}
