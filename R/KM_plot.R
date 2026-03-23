
#' Kaplan Meier plot
#'
#' @export
#' @param test_var Vector for test. Can be numeric or character vector. If the vector is numeric, than a cut will be used to divide the vector into two groups.
#' @param time_var Variable name that indicates the survival time
#' @param event_var Variable name that indicates the survival event
#' @param pheno_mat A data frame that contains the survival time and event variable
#' @param cut If test_var is numeric, the cut will be applied to divide the vector into two groups. The default is "median" cut. Can select from "median", "mean" and "maxstat".
#' If select "maxstat", the cut will determine by the Maximally Selected Rank statistics.
#'
#' @return Return a ggsurvplot.
#' @param ... Additional arguments passed to `survminer::ggsurvplot()`.
#'
#' @import survival
#' @import survminer

# plot KM
plot_km <- function(test_var, time_var, event_var, pheno_mat, cut = "median", ...){

  if(is.numeric(test_var)){

    if(cut == "median"){
      m <- median(test_var)
    }

    if(cut == "mean"){
      m <- mean(test_var)
    }

    if(cut == "maxstat"){

      df <- data.frame(cluster = test_var, time = pheno_mat[[time_var]], death = pheno_mat[[event_var]])
      m <- maxstat::maxstat.test(Surv(time, death) ~ cluster, data = df, smethod = "LogRank")$estimate

    }
    gene_cut <- ifelse(test_var < m, "low", "high")

  } else {

    gene_cut = test_var

  }

  df <- data.frame(cluster = gene_cut, time = pheno_mat[[time_var]], death = pheno_mat[[event_var]])

  fo <- as.formula(paste0("Surv(time, death) ~ cluster"))
  fit <- surv_fit(fo, data = df)
  p <- round(-log10(survdiff(fo, data = df)$p),4)

  survminer::ggsurvplot(
    fit,
    data = df,
    size = 1,
    palette = "jco",
    conf.int = FALSE,
    pval = TRUE,
    pval.coord = c(0, 0.15),
    legend = c(0.8, 0.85),
    risk.table = TRUE,
    ...
  )
}
