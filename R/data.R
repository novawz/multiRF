#' TCGA BRCA Expression Data
#'
#' A dataset containing matched mRNA, miRNA, and DNA methylation features from
#' the TCGA BRCA cohort.
#'
#' @name tcga_brca
#' @format A named list with three data frames:
#' \describe{
#'   \item{gene}{mRNA expression features.}
#'   \item{methy}{DNA methylation features.}
#'   \item{mirna}{miRNA expression features.}
#' }
#' @source \url{https://www.cancer.gov/tcga}
NULL

#' TCGA BRCA Clinical Data
#'
#' Clinical metadata matched to the bundled TCGA BRCA omics cohort.
#'
#' @name tcga_brca_clinical
#' @format A data frame with 674 rows (one per sample) and 9 columns:
#' \describe{
#'   \item{sampleID}{TCGA sample barcode matching the row names of the
#'     \code{tcga_brca} omics blocks.}
#'   \item{sample_type}{Sample type: \code{"Primary Tumor"},
#'     \code{"Solid Tissue Normal"}, or \code{"Metastatic"}.}
#'   \item{BRCA_Subtype_PAM50}{PAM50 intrinsic subtype call (\code{"Basal"},
#'     \code{"Her2"}, \code{"LumA"}, \code{"LumB"}, \code{"Normal"}, or
#'     \code{NA}).}
#'   \item{OS}{Overall survival event indicator (1 = death, 0 = censored).}
#'   \item{OS.time}{Overall survival time in days.}
#'   \item{PFI}{Progression-free interval event indicator (1 = event,
#'     0 = censored).}
#'   \item{PFI.time}{Progression-free interval time in days.}
#'   \item{age_at_initial_pathologic_diagnosis}{Age in years at initial
#'     pathologic diagnosis.}
#'   \item{pathologic_stage}{AJCC pathologic tumor stage (e.g.
#'     \code{"Stage IIA"}).}
#' }
#' @source \url{https://www.cancer.gov/tcga}
NULL
