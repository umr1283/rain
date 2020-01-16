#' check_bin
#'
#' @inheritParams estimate_ethnicity
#'
#' @keywords internal
check_bin <- function(bin_path) {
  if (!all(sapply(bin_path, file.exists))) {
    message_prefix <- "[CARoT] "
    stop(
      message_prefix,
      paste0(
        "No binary found for the following tools: ",
        glue::glue_collapse(
          x = paste0(
            names(bin_path[!sapply(bin_path, file.exists)]),
            ' ("', bin_path[!sapply(bin_path, file.exists)], '")'
          ),
          sep = ", ",
          last = " and "
        ),
        "!"
      )
    )
  }
}
