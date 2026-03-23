#' Print method for `mrf3` objects
#'
#' Displays a compact summary of fitted multiRF models including
#' block names, connection structure, weight statistics, and class info.
#'
#' @method print mrf3
#' @param x An `mrf3` object.
#' @param max_weights Maximum number of top variables to display per block.
#' @param ... Additional arguments (unused).
#'
#' @return The input object, invisibly.
#' @export
print.mrf3 <- function(x, max_weights = 5L, ...) {
  is_vs <- inherits(x, "vs")
  tag <- if (is_vs) "mrf3 (variable-selected)" else "mrf3"
  cat(tag, "\n", sep = "")
  cat(strrep("-", nchar(tag)), "\n", sep = "")

  ## Type & tree count
  type_str <- if (!is.null(x$type)) x$type else "unknown"
  ntree_str <- if (!is.null(x$ntree)) as.character(x$ntree) else "NA"
  cat("  Type  : ", type_str, "\n", sep = "")
  cat("  Trees : ", ntree_str, "\n", sep = "")

  ## Connection
  if (!is.null(x$connection) && is.list(x$connection)) {
    conn_strs <- vapply(x$connection, function(cc) paste(cc, collapse = " -> "), character(1))
    cat("  Connections (", length(conn_strs), "):\n", sep = "")
    for (cs in conn_strs) cat("    ", cs, "\n", sep = "")
  }

  ## OOB error
  if (!is.null(x$oob_err) && is.numeric(x$oob_err) && is.finite(x$oob_err)) {
    cat("  OOB error : ", round(x$oob_err, 4), "\n", sep = "")
  }

  ## Weights summary
  if (!is.null(x$weights) && is.list(x$weights)) {
    cat("\n  Weights per block:\n")
    for (nm in names(x$weights)) {
      w <- x$weights[[nm]]
      p_total <- length(w)
      p_nonzero <- sum(w != 0)
      cat("    ", nm, ": ", p_nonzero, "/", p_total, " non-zero", sep = "")
      if (p_nonzero > 0L) {
        top_idx <- utils::head(order(w, decreasing = TRUE), max_weights)
        top_nm <- names(w)[top_idx]
        top_val <- round(w[top_idx], 4)
        top_str <- paste0(top_nm, "(", top_val, ")")
        cat("  top: ", paste(top_str, collapse = ", "), sep = "")
      }
      cat("\n")
    }
  }

  ## Data attached?
  if (!is.null(x$dat.list) && is.list(x$dat.list)) {
    dims <- vapply(x$dat.list, function(d) paste0(nrow(d), "x", ncol(d)), character(1))
    cat("\n  Data: ", paste(paste0(names(x$dat.list), "[", dims, "]"), collapse = ", "), "\n", sep = "")
  }

  ## VS-specific info
  if (is_vs && !is.null(x$thres)) {
    cat("\n  Variable selection applied (class 'vs').\n")
  }

  invisible(x)
}


#' Print method for `mrf3_vs` output (vs class)
#'
#' Provides a concise overview of variable selection results.
#'
#' @method print vs
#' @param x An object with class `c("mrf3", "vs")` returned by `mrf3_vs()`.
#' @param ... Additional arguments passed to `print.mrf3()`.
#'
#' @return The input object, invisibly.
#' @keywords internal
print.vs <- function(x, ...) {
  ## Delegate to print.mrf3 which already handles the vs subclass
  print.mrf3(x, ...)
}
