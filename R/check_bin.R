#' check_bin
#'
#' @inheritParams estimate_ethnicity
#'
#' @keywords internal
check_bin <- function(bin_path) {
  if (!all(sapply(bin_path, file.exists))) {
    message_prefix <- "[rain] "
    text <- paste0(
      names(bin_path[!sapply(bin_path, file.exists)]),
      ' ("', bin_path[!sapply(bin_path, file.exists)], '")'
    )
    stop(
      message_prefix,
      paste0(
        "No binary found for the following tools: ",
        paste(text[-length(text)], collapse = ", "),
        " and ",
        text[length(text)],
        "!"
      )
    )
  }
}
