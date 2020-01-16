#' check_input
#'
#' @param input A `character`. A file path.
#'
#' @keywords internal
check_input <- function(input) {
  message_prefix <- "[CARoT] "

  name <- deparse(substitute(input))
  is_input_good <- (length(input) == 1 && (fs::is_dir(input) | fs::is_file(input))) |
    all(fs::is_file(input))

  if (!is_input_good) {
    stop(
      message_prefix, 'A valid "', name, '" must be provided, ',
      'either a directory (with VCF files) or vcf file(s)!'
    )
  }
  if (length(input) == 1 && fs::is_dir(input)) {
    list_input <- list.files(path = input, pattern = ".vcf.gz$", full.names = TRUE)
  } else {
    list_input <- input
  }
  if (!all(grepl(".vcf.gz$", list_input) & fs::is_file(list_input))) {
    stop(message_prefix, 'VCF files must be compressed using bgzip!')
  }
  if (length(list_input)==0) {
    stop(
      message_prefix, 'A valid "', name, '" must be provided, ',
      'either a directory (with VCF files) or a vcf file!'
    )
  }
  invisible(list_input)
}
