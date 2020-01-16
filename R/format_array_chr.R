#' format_array_chr
#'
#' @inheritParams estimate_ethnicity
#'
#' @keywords internal
format_array_chr <- function(
  cohort_name,
  input_vcfs,
  output_directory,
  ref1kg_vcfs,
  ref1kg_maf,
  quality_tag,
  quality_threshold,
  recode,
  n_cores,
  bin_path
) {
  out <- parallel::mclapply(
    X = 1:22,
    mc.preschedule = FALSE,
    mc.cores = min(parallel::detectCores(), n_cores),
    mc_input_vcfs = input_vcfs,
    mc_ref1kg_vcfs = ref1kg_vcfs,
    mc_ref1kg_maf = ref1kg_maf,
    mc_quality_tag = quality_tag,
    mc_quality_threshold = quality_threshold,
    mc_output_directory = output_directory,
    FUN = function(
      ichr,
      mc_input_vcfs,
      mc_ref1kg_vcfs,
      mc_ref1kg_maf,
      mc_quality_tag,
      mc_quality_threshold,
      mc_output_directory
    ) {
      ipattern <- paste0("^[^0-9]*chr", ichr, "[^0-9]+.*vcf.gz$")
      iinput_vcfs <- mc_input_vcfs[grep(pattern = gsub("chr", "", ipattern), x = basename(mc_input_vcfs))]
      iref1kg_vcfs <- mc_ref1kg_vcfs[grep(pattern = ipattern, x = basename(mc_ref1kg_vcfs))]

      format_vcf(
        input_vcfs = iinput_vcfs,
        ref1kg_vcfs = iref1kg_vcfs,
        ref1kg_maf = mc_ref1kg_maf,
        ichr = ichr,
        quality_tag = mc_quality_tag,
        quality_threshold = mc_quality_threshold,
        output_directory = mc_output_directory,
        recode = recode,
        bin_path = bin_path
      )
  })


  temp_file <- tempfile(fileext = ".merge")
  cat(
    list.files(path = output_directory, pattern = "_merged.vcf.gz$", full.names = TRUE),
    sep = "\n",
    file = temp_file
  )
  system(
    intern = TRUE, wait = TRUE,
    command = paste(
      bin_path[["bcftools"]], "concat",
      "--file-list", temp_file,
      "--allow-overlaps",
      "--output-type z",
      "--output", paste0(output_directory, "/all.vcf.gz")
    )
  )
  unlink(x = temp_file, recursive = TRUE)


  system(
    intern = TRUE, wait = TRUE,
    command = paste(
      bin_path[["plink1.9"]],
      "--vcf", paste0(output_directory, "/all.vcf.gz"),
      "--snps-only",
      "--maf", ref1kg_maf,
      "--hwe 0.0001",
      "--geno 0.1",
      "--make-bed",
      "--double-id",
      "--out", paste0(output_directory, "/all")
    )
  )

  unlink(
    x = list.files(path = output_directory, pattern = "_merged.vcf.gz", full.names = TRUE),
    recursive = TRUE
  )
  unlink(
    x = list.files(path = output_directory, pattern = "all.vcf.gz", full.names = TRUE),
    recursive = TRUE
  )

  invisible()
}
