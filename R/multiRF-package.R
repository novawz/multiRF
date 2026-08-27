#' multiRF: A package for multiRF methods
#'
#' @useDynLib multiRF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom dplyr filter group_by n slice_max slice_min
#' @importFrom stats ave density na.omit rnorm setNames var
#' @importFrom utils tail
#'
#' @keywords internal
"_PACKAGE"
