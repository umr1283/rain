#' merge_vcf
#'
#' @inheritParams estimate_ethnicity
#'
#' @keywords internal
merge_vcf <- function(input_vcfs, bin_path) {
  if (length(input_vcfs) > 1) {
    output_temp <- paste0(tempdir(), "/samples_merged.vcf.gz")
    vcf_list <- paste0(tempdir(), "/samples_merged.txt")
    cat(input_vcfs, sep = "\n", file = vcf_list)
    system(
      intern = TRUE, wait = TRUE,
      command = paste(
        bin_path[["bcftools"]], "merge --merge none",
        " --file-list", vcf_list,
        "--output-type z",
        "--output", output_temp,
        "&&",
        bin_path[["tabix"]], "-p vcf", output_temp
      )
    )
    unlink(vcf_list)
  } else {
    output_temp <- input_vcfs
  }

  output_temp
}
