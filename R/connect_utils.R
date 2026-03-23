normalize_connect_list <- function(connect_list,
                                   n_blocks = NULL,
                                   valid_names = NULL,
                                   arg_name = "connect_list") {
  if (is.null(connect_list)) {
    return(NULL)
  }

  if (is.character(connect_list) && !is.list(connect_list)) {
    connect_list <- stringr::str_split(connect_list, "_")
  }
  if (!is.list(connect_list)) {
    stop(
      "`", arg_name, "` must be a list of length-2 character vectors ",
      "or a character vector with `response_predictor` names."
    )
  }

  connect_list <- lapply(connect_list, as.character)

  if (!is.null(n_blocks) && is.finite(n_blocks) && as.integer(n_blocks) > 1L &&
      any(lengths(connect_list) != 2L)) {
    stop(
      "Each element of `", arg_name, "` must contain exactly two block names: ",
      "c(response, predictor)."
    )
  }

  if (!is.null(valid_names)) {
    valid_names <- as.character(valid_names)
    bad <- setdiff(unique(unlist(connect_list)), valid_names)
    if (length(bad) > 0L) {
      stop(
        "`", arg_name, "` contains unknown block name(s): ",
        paste(bad, collapse = ", ")
      )
    }
  }

  connect_list
}
