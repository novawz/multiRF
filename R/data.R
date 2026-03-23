#' TCGA BRCA Expression Data
#'
#' A dataset containing matched mRNA, miRNA, and DNA methylation features from
#' the TCGA BRCA cohort.
#'
#' @name tcga_brca
#' @format A named list with three data frames:
#' \describe{
#'   \item{gene}{mRNA expression features.}
#'   \item{miRNA}{miRNA expression features.}
#'   \item{methy}{DNA methylation features.}
#' }
#' @source \url{https://www.cancer.gov/tcga}
NULL

#' TCGA BRCA Clinical Data
#'
#' Clinical metadata matched to the bundled TCGA BRCA omics cohort.
#'
#' @name tcga_brca_clinical
#' @format A data frame with one row per sample.
#' @source \url{https://www.cancer.gov/tcga}
NULL
