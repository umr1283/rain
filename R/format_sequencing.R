#' format_sequencing
#'
#' @inheritParams estimate_ethnicity
#'
#' @keywords internal
format_sequencing <- function(
  cohort_name,
  input_vcfs,
  output_directory,
  ref1kg_vcfs,
  ref1kg_maf,
  recode,
  vcf_half_call,
  bin_path
) {
  merged_vcfs <- merge_vcf(
    input_vcfs = input_vcfs,
    bin_path = bin_path
  )

  format_vcf(
    input_vcfs = merged_vcfs,
    ref1kg_vcfs = ref1kg_vcfs,
    ref1kg_maf = ref1kg_maf,
    ichr = "ALL",
    quality_tag = NULL,
    quality_threshold = NULL,
    output_directory = output_directory,
    recode = recode,
    bin_path = bin_path
  )

  system(
    intern = TRUE, wait = TRUE,
    command = paste(
      bin_path[["plink"]],
      "--vcf", list.files(path = output_directory, pattern = "_merged.vcf.gz$", full.names = TRUE),
      "--snps-only",
      "--maf", ref1kg_maf,
      # "--hwe 0.0001",
      "--geno 0.1",
      "--make-bed",
      "--double-id",
      "--vcf-half-call", paste0("'", vcf_half_call, "'"),
      "--out", file.path(output_directory, "all")
    )
  )

  # unlink(
  #   x = list.files(path = output_directory, pattern = "_merged.vcf.gz", full.names = TRUE),
  #   recursive = TRUE
  # )

  invisible()
}
