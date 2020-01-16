#' Format VCF files and compute the genomic components (and some figures) for ethnicity.
#'
#' @param cohort_name A `character`. A name to describe the studied population compared to 1,000 Genomes.
#' @param input_vcfs A `character`. A path to one or several VCFs file.
#' @param input_type A `character`. Either `"array"` or `"sequencing"`.
#' @param output_directory A `character`. The path where the data and figures is written.
#' @param ref1kg_vcfs A `character`. A path to the reference VCFs files (i.e., 1,000 Genomes sequencing data).
#' @param ref1kg_population A `character`. A file which describe samples and their ethnicity.
#' @param ref1kg_maf A `numeric`. MAF threshold for SNPs in 1,000 Genomes
#' @param splitted_by_chr A `logical`. Is the VCFs files splitted by chromosome?
#' @param quality_tag A `character`. Name of the imputation quality tag for `"array"`,
#'   *e.g.*, `"INFO"` or `"R2"`. Default is `NULL`.
#' @param quality_threshold A `numeric`. The threshold to keep/discard SNPs based on their imputation quality.
#' @param recode A `character`. Which VCF should be filtered and recode, either `"all"` or `"input"`.
#' @param vcf_half_call A `character`. The mode to handle half-call.
#'     + 'haploid'/'h': Treat half-calls as haploid/homozygous (the PLINK 1 file format does not distinguish between the two). This maximizes similarity between the VCF and BCF2 parsers.
#'     + 'missing'/'m': Treat half-calls as missing (default).
#'     + 'reference'/'r': Treat the missing part as reference.
#' @param n_cores An `integer`. The number of CPUs to use to estimate the ethnicity.
#' @param bin_path A `list(character)`. A list giving the binary path of
#'   `vcftools`, `bcftools`, `bgzip`, `tabix` and `plink1.9`.
#'
#' @return A `data.frame`.
#' @export
estimate_ethnicity <- function(
  cohort_name,
  input_vcfs,
  input_type,
  output_directory,
  ref1kg_vcfs,
  ref1kg_population,
  ref1kg_maf = 0.05,
  splitted_by_chr = TRUE,
  quality_tag = NULL,
  quality_threshold = 0.9,
  recode = "all",
  vcf_half_call = "missing",
  n_cores = 6,
  bin_path = list(
    vcftools = "/usr/bin/vcftools",
    bcftools = "/usr/bin/bcftools",
    bgzip = "/usr/bin/bgzip",
    tabix = "/usr/bin/tabix",
    plink1.9 = "/usr/bin/plink1.9"
  )
) {
  message_prefix <- "[CARoT] "

  check_bin(bin_path)

  if (!is_ethnicity_done(output_directory)) {
    stop(message_prefix, '"estimate_ethnicity" has been canceled by the user!')
  }

  list_input <- check_input(input = input_vcfs)
  list_ref <- check_input(input = ref1kg_vcfs)

  if (!input_type%in%c("array", "sequencing")) {
    stop(message_prefix, '"input_type" must be either "array" or "sequencing"!')
  }

  if (input_type=="sequencing" & length(list_ref)!=1) {
    stop(
      message_prefix, 'A unique vcf file ("ref1kg_vcfs") must be provided ',
      'with `input_type = "sequencing"`!'
    )
  }
  if (input_type=="array" & !splitted_by_chr & length(list_ref)!=1) {
    stop(
      message_prefix, 'A unique vcf file ("ref1kg_vcfs") must be provided ',
      'with `input_type = "array"` & `splitted_by_chr = FALSE`!'
    )
  }

  if (input_type=="sequencing") {
    quality_tag <- NULL
  }

  ######################
  ### Formating VCFs ###
  ######################
  message(message_prefix, "Formating VCFs ...")
  switch(
    EXPR = input_type,
    "array" = {
      if (splitted_by_chr) {
        format_array_chr(
          cohort_name = cohort_name,
          input_vcfs = list_input,
          output_directory = output_directory,
          ref1kg_vcfs = list_ref,
          ref1kg_maf = ref1kg_maf,
          quality_tag = quality_tag,
          quality_threshold = quality_threshold,
          recode = recode,
          n_cores = n_cores,
          bin_path = bin_path
        )
      } else {
        format_array_all(
          cohort_name = cohort_name,
          input_vcfs = list_input,
          output_directory = output_directory,
          ref1kg_vcfs = list_ref,
          ref1kg_maf = ref1kg_maf,
          quality_tag = quality_tag,
          quality_threshold = quality_threshold,
          bin_path = bin_path
        )
      }
    },
    "sequencing" = {
      format_sequencing(
        cohort_name = cohort_name,
        input_vcfs = list_input,
        output_directory = output_directory,
        ref1kg_vcfs = list_ref,
        ref1kg_maf = ref1kg_maf,
        recode = recode,
        vcf_half_call = vcf_half_call,
        bin_path = bin_path
      )
    }
  )


  ######################
  ### Performing PCA ###
  ######################
  compute_pca(
    cohort_name = cohort_name,
    input_plink = paste0(output_directory, "/all"),
    output_directory = output_directory,
    ref1kg_population = ref1kg_population
  )

}
