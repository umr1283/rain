#' format_vcf
#'
#' @inheritParams estimate_ethnicity
#' @param ichr A `character` or `numeric`. The chromosome identifier.
#'
#' @keywords internal
format_vcf <- function(
  input_vcfs,
  ref1kg_vcfs,
  ref1kg_maf,
  ichr,
  quality_tag,
  quality_threshold,
  output_directory,
  recode,
  bin_path
) {
  temp_directory <- paste0(file.path(tempdir(), "chr"), ichr)
  invisible(sapply(
    X = file.path(temp_directory, c("study", "ref", "isec")),
    FUN = dir.create,
    recursive = TRUE, showWarnings = FALSE, mode = '0777'
  ))
  output_study_temp <- file.path(temp_directory, "study", paste0("filtered_", basename(input_vcfs)))
  output_study_final <- file.path(temp_directory, "study", paste0("final_", basename(input_vcfs)))
  output_ref <- file.path(temp_directory, "ref", paste0("filtered_", basename(input_vcfs)))
  if (is.character(ichr)) {
    output_merge_temp <- file.path(temp_directory, paste0("chr", ichr, "_merged.vcf.gz"))
    output_merge <- file.path(output_directory, paste0("chr", ichr, "_merged.vcf.gz"))
  } else {
    output_merge_temp <- file.path(temp_directory, sprintf("chr%02d_merged.vcf.gz", ichr))
    output_merge <- file.path(output_directory, sprintf("chr%02d_merged.vcf.gz", ichr))
  }

  switch(
    EXPR = recode,
    "all" = {
      out_cmd <- system(
        intern = TRUE, wait = TRUE,
        command = paste(
          bin_path[["vcftools"]],
          "--gzvcf", ref1kg_vcfs,
          "--maf", ref1kg_maf,
          "--recode",
          "--stdout",
          "|", bin_path[["bgzip"]], "-c >", output_ref,
          "&&",
          bin_path[["tabix"]], "-p vcf", output_ref
        )
      )
    },
    "input" = { output_ref <- ref1kg_vcfs }
  )

  chr_conv_file <- tempfile(fileext = ".conv")
  data.table::fwrite(
    x = rbind(
      matrix(paste0(c("chr", ""), rep(c(1:22, "X", "Y"), each = 2)), byrow = TRUE, ncol = 2),
      matrix(rep(c(1:22, "X", "Y"), each = 2), byrow = TRUE, ncol = 2)
    ),
    file = chr_conv_file,
    quote = FALSE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
  )

  out_cmd <- system(
    intern = TRUE, wait = TRUE,
    command = paste(
      if (!is.null(quality_tag)) {
        paste(bin_path[["vcftools"]],
          "--gzvcf", input_vcfs,
          "--get-INFO", quality_tag,
          "--out", gsub("filtered_", "excluded_", output_study_temp),
          "&&",
          'awk \'{if($5<', quality_threshold, ') print $1"\t"$2}\'',
          paste0(gsub("filtered_", "excluded_", output_study_temp), ".INFO"),
          ">", paste0(gsub("filtered_", "excluded_", output_study_temp), ".exclude"),
          "&&"
        )
      },
      bin_path[["vcftools"]],
      "--gzvcf", input_vcfs,
      if (!is.null(quality_tag)) {
        paste0(
          "--exclude-positions ", gsub("filtered_", "excluded_", output_study_temp), ".exclude"
        )
      },
      "--remove-indels",
      "--remove-filtered-all",
      "--max-missing-count 1",
      "--recode",
      "--stdout",
      "|", bin_path[["bgzip"]], "-c >", output_study_temp,
      "&&",
      bin_path[["tabix"]], "-p vcf", output_study_temp,
      "&&",
      bin_path[["bcftools"]], "annotate",
      "--rename-chrs", chr_conv_file, output_study_temp,
      "|", bin_path[["bgzip"]], "-c >", output_study_final,
      "&&",
      bin_path[["tabix"]], "-p vcf", output_study_final,
      "&&",
      bin_path[["bcftools"]], "isec",
      "--collapse none",
      "--nfiles=2",
      output_study_final,
      output_ref,
      "--output-type z",
      "--prefix", paste0(temp_directory, "/isec"),
      "&&",
      bin_path[["bcftools"]], "merge --merge none",
      paste0(temp_directory, "/isec/0000.vcf.gz"),
      paste0(temp_directory, "/isec/0001.vcf.gz"),
      "--output-type z",
      "--output", output_merge_temp,
      "&&",
      bin_path[["tabix"]], "-p vcf", output_merge_temp,
      "&&",
      bin_path[["bcftools"]], "annotate", "-x INFO,^FORMAT/GT", output_merge_temp,
      "--output-type z",
      "--output", output_merge,
      "&&",
      bin_path[["tabix"]], "-p vcf", output_merge
    )
  )

  unlink(x = temp_directory, recursive = TRUE)

  invisible()
}
