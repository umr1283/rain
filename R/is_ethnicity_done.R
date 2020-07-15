#' is_ethnicity_done
#'
#' @inheritParams estimate_ethnicity
#'
#' @keywords internal
is_ethnicity_done <- function(output_directory) {
  message_prefix <- "[rain] "

  plink_files_exists <- all(file.exists(file.path(output_directory, paste0("all", c(".bim", ".fam", ".bed")))))

  if (plink_files_exists & interactive()) {
    message(paste0(
      message_prefix,
      "bim/bed/fam files already exist!\n",
      "  (y) to format (again) the VCF files and to perform the PCA.\n",
      "  (n) to cancel and call `compute_pca`.\n",
      "  Please choose (y) or (n): "
    ))
    answer <- readline(prompt = "")

    out <- grepl("y", answer, ignore.case = TRUE)
  } else {
    out <- TRUE
  }
}
