#' format_array_all
#'
#' @inheritParams estimate_ethnicity
#'
#' @keywords internal
format_array_all <- function(
  cohort_name,
  input_vcfs,
  output_directory,
  ref1kg_vcfs,
  ref1kg_maf,
  quality_tag,
  quality_threshold,
  recode,
  bin_path
) {
  format_vcf(
    input_vcfs = input_vcfs,
    ref1kg_vcfs = ref1kg_vcfs,
    ref1kg_maf = ref1kg_maf,
    ichr = "ALL",
    quality_tag = quality_tag,
    quality_threshold = quality_threshold,
    output_directory = output_directory,
    recode,
    bin_path = bin_path
  )

  system(
    intern = TRUE, wait = TRUE,
    command = paste(
      bin_path[["plink"]],
      "--vcf", list.files(path = output_directory, pattern = "_merged.vcf.gz$", full.names = TRUE),
      "--snps-only",
      "--maf", ref1kg_maf,
      "--hwe 0.0001",
      "--geno 0.1",
      "--make-bed",
      "--double-id",
      "--out", file.path(output_directory, "all")
    )
  )

  unlink(
    x = list.files(path = output_directory, pattern = "_merged.vcf.gz", full.names = TRUE),
    recursive = TRUE
  )

  invisible()
}
