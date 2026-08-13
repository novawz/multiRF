normalize_connect_list <- function(connect_list,
                                   n_blocks = NULL,
                                   valid_names = NULL,
                                   arg_name = "connect_list") {
  if (is.null(connect_list)) {
    return(NULL)
  }

  if (is.character(connect_list) && !is.list(connect_list)) {
    labels <- as.character(connect_list)
    if (!is.null(valid_names)) {
      valid_names <- as.character(valid_names)
      candidates <- unlist(
        lapply(valid_names, function(response) {
          lapply(valid_names, function(predictor) c(response, predictor))
        }),
        recursive = FALSE
      )
      connect_list <- lapply(labels, function(label) {
        hit <- vapply(
          candidates,
          function(pair) identical(paste(pair, collapse = "_"), label),
          logical(1)
        )
        if (sum(hit) != 1L) {
          reason <- if (sum(hit) > 1L) "is ambiguous" else "does not match any block pair"
          stop(
            "Connection label `", label, "` ", reason,
            ". Use list(c(response, predictor)) when block names contain underscores."
          )
        }
        candidates[[which(hit)]]
      })
    } else {
      connect_list <- lapply(labels, function(label) {
        pair <- strsplit(label, "_", fixed = TRUE)[[1L]]
        if (length(pair) != 2L) {
          stop(
            "Connection label `", label, "` is not unambiguous. ",
            "Use list(c(response, predictor))."
          )
        }
        pair
      })
    }
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

# Return a stable, unique display name for each connection.  Downstream code
# must use the explicit response/predictor metadata for interpretation; the
# names are identifiers only and therefore may receive make.unique suffixes.
connection_display_names <- function(connect_list) {
  make.unique(vapply(
    connect_list,
    function(conn) paste(as.character(conn), collapse = "_"),
    character(1)
  ))
}

model_indices_for_connections <- function(mod_list, connect_list) {
  model_pairs <- lapply(seq_along(mod_list), function(i) {
    parse_model_pair(names(mod_list)[i], model = mod_list[[i]])
  })
  vapply(connect_list, function(conn) {
    conn <- as.character(conn)
    hit <- which(vapply(model_pairs, identical, logical(1), conn))
    if (length(hit) != 1L) {
      stop(
        "Expected exactly one fitted model for connection ",
        paste(conn, collapse = " <- "), "; found ", length(hit), "."
      )
    }
    hit
  }, integer(1))
}
