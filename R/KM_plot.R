
#' Kaplan Meier plot
#'
#' @export
#' @param test_var Sample-level grouping variable or continuous score. Numeric
#'   vectors with more than five unique values are treated as continuous;
#'   low-cardinality numeric vectors are treated as group labels.
#' @param time_var Column name containing survival or follow-up time.
#' @param event_var Column name containing a binary 0/1 event indicator.
#' @param pheno_mat Data frame containing the time and event columns.
#' @param cut For a continuous `test_var`, split at the `"median"`, `"mean"`,
#'   or a maximally selected rank-statistic (`"maxstat"`) cutoff.
#' @param base_size Base font size.
#' @param base_family Font family.
#'
#' @return A `ggsurvplot` object.
#' @param ... Additional arguments passed to `survminer::ggsurvplot()`.
#' @param na.rm Logical; whether to remove rows with missing values.

# plot KM
plot_km <- function(test_var, time_var, event_var, pheno_mat,
                    cut = c("median", "mean", "maxstat"), na.rm = TRUE,
                    base_size = 8, base_family = "sans", ...){
  if (!requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("survminer", quietly = TRUE)) {
    stop("Packages `survival` and `survminer` are required for `plot_km()`.")
  }
  if (!is.data.frame(pheno_mat)) pheno_mat <- as.data.frame(pheno_mat)
  if (!all(c(time_var, event_var) %in% names(pheno_mat))) {
    stop("`time_var` and `event_var` must name columns in `pheno_mat`.")
  }
  if (length(test_var) != nrow(pheno_mat)) {
    stop("`test_var` must have one value per row of `pheno_mat`.")
  }
  dots <- list(...)
  cut <- match.arg(cut)
  df <- data.frame(
    test = test_var,
    time = pheno_mat[[time_var]],
    event = pheno_mat[[event_var]],
    stringsAsFactors = FALSE
  )
  complete <- stats::complete.cases(df)
  if (any(!complete)) {
    if (!isTRUE(na.rm)) stop("Missing values are present; set `na.rm = TRUE` to omit them.")
    df <- df[complete, , drop = FALSE]
  }
  if (nrow(df) < 3L || any(!is.finite(df$time)) || any(df$time < 0) ||
      any(!df$event %in% c(0, 1))) {
    stop("Survival data must contain at least three complete rows, non-negative finite times, and 0/1 events.")
  }
  is_continuous <- is.numeric(df$test) && length(unique(df$test)) > 5L
  if (is_continuous) {
    if (cut == "maxstat") {
      if (!requireNamespace("maxstat", quietly = TRUE)) {
        stop("Package `maxstat` is required when `cut = 'maxstat'`.")
      }
      ms <- maxstat::maxstat.test(
        survival::Surv(time, event) ~ test,
        data = df,
        smethod = "LogRank"
      )
      threshold <- as.numeric(ms$estimate)
      if (is.null(dots$pval)) {
        warning(
          "The ordinary log-rank p-value is hidden after a data-selected maxstat cutoff; ",
          "report an adjusted maxstat result for inference.",
          call. = FALSE
        )
      }
    } else {
      threshold <- if (cut == "median") stats::median(df$test) else mean(df$test)
    }
    df$cluster <- factor(ifelse(df$test < threshold, "Low", "High"), levels = c("Low", "High"))
  } else {
    df$cluster <- factor(df$test)
  }
  if (nlevels(droplevels(df$cluster)) < 2L) stop("`test_var` must define at least two non-empty groups.")
  fit <- survival::survfit(survival::Surv(time, event) ~ cluster, data = df)
  defaults <- list(
    fit = fit,
    data = df,
    size = 0.65,
    palette = unname(mrf3_plot_palette(nlevels(df$cluster))),
    linetype = "strata",
    conf.int = FALSE,
    pval = !identical(cut, "maxstat"),
    risk.table = TRUE,
    legend.title = "",
    legend.labs = levels(df$cluster),
    ggtheme = theme_mrf3(base_size = base_size, base_family = base_family),
    tables.theme = theme_mrf3(base_size = base_size, base_family = base_family)
  )
  args <- utils::modifyList(defaults, dots)
  do.call(survminer::ggsurvplot, args)
}
