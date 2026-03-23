#' multiRF: A package for multiRF methods
#'
#' @useDynLib multiRF, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr %>%
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom plyr l_ply llply laply
#' @importFrom dplyr filter group_by mutate n select slice_max slice_min summarize_all summarise_at
#' @importFrom graphics legend pairs par plot
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.formula ave cor density dgamma dist dnorm median model.matrix na.omit pnorm prcomp quantile reorder rnorm sd setNames var
#' @importFrom utils combn tail
#' @import circlize
#' @import cluster
#' @import tibble
#' @import tidyr
#' @import truncnorm
#'
#' @keywords internal
"_PACKAGE"
