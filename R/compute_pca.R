#' Compute the genomic components (and some figures) for ethnicity based on VCF files.
#'
#' @inheritParams estimate_ethnicity
#' @param input_plink A `character`. The path to plink format files (i.e., `.bed`, `.bim` and `.fam` files).
#'
#' @return A `data.frame`.
#' @export
compute_pca <- function(cohort_name, input_plink, output_directory, ref1kg_population) {
  message_prefix <- "[CARoT] "

  ######################
  ### Performing PCA ###
  ######################
  message(message_prefix, "Performing PCA ...")

  res_pca <- flashpcaR::flashpca(input_plink, ndim = 10)

  pca_gg <- as.data.frame(res_pca[["projection"]])
  colnames(pca_gg) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_gg))))

  pca_gg <- dplyr::as_tibble(pca_gg)

  fid_iid <- utils::read.table(paste0(input_plink, ".fam"), stringsAsFactors = FALSE)[, c(1, 2)]

  if (all.equal(fid_iid[[1]], fid_iid[[2]])) {
    pca_gg[["sample"]] <- fid_iid[[2]]
  } else {
    pca_gg[["sample"]] <- paste(fid_iid[[1]], fid_iid[[2]], sep = "/")
  }

  pca_gg <- dplyr::left_join(
    x = pca_gg,
    y = suppressWarnings(
      readr::read_tsv(
        file = ref1kg_population,
        col_types = readr::cols_only(
          sample = readr::col_character(),
          pop = readr::col_character(),
          super_pop = readr::col_character()
        )
      )
    ),
    by = "sample"
  )
  pca_gg[["cohort"]] <- factor(
    x = ifelse(is.na(pca_gg[["pop"]]), cohort_name, "1,000 Genomes"),
    levels = c(cohort_name, "1,000 Genomes")
  )
  pca_gg[["pop"]] <- factor(
    x = ifelse(is.na(pca_gg[["pop"]]), cohort_name, pca_gg[["pop"]]),
    levels = c(cohort_name, sort(unique(stats::na.omit(pca_gg[["pop"]]))))
  )
  pca_gg[["super_pop"]] <- factor(
    x = ifelse(is.na(pca_gg[["super_pop"]] ), cohort_name, pca_gg[["super_pop"]]),
    levels = c(cohort_name, "AFR", "AMR", "EAS", "SAS", "EUR")
  )
  pca_gg <- dplyr::select(.data = pca_gg, "sample", dplyr::everything())

  p_ethni <- ggplot2::ggplot(
    data = pca_gg,
    mapping = ggplot2::aes(x = .data[["PC01"]], y = .data[["PC02"]], colour = .data[["super_pop"]])
  ) +
    ggplot2::theme_light(base_size = 12) +
  	ggplot2::geom_hline(yintercept = 0, linetype = 1, size = 0.5, na.rm = TRUE) +
  	ggplot2::geom_vline(xintercept = 0, linetype = 1, size = 0.5, na.rm = TRUE) +
    ggforce::geom_mark_ellipse(mapping = ggplot2::aes_string(fill = "super_pop"), con.cap = 0) +
    ggplot2::geom_point(mapping = ggplot2::aes_string(shape = "super_pop"), na.rm = TRUE) +
  	ggplot2::scale_colour_viridis_d(na.translate = FALSE, drop = FALSE, end = 0.9) +
    ggplot2::scale_fill_viridis_d(na.translate = FALSE, drop = FALSE, end = 0.9) +
    ggplot2::scale_shape_manual(values = c(3, rep(1, 5))) +
    ggplot2::labs(
      shape = NULL,
      colour = NULL,
      fill = NULL,
      caption = paste(
        "SNPs:",
        scales::comma(R.utils::countLines(paste0(input_plink, ".bim")))
      )
    ) +
    ggforce::facet_zoom(
      xlim = range(dplyr::filter(pca_gg, .data[["cohort"]] == !!cohort_name)[["PC01"]]),
      ylim = range(dplyr::filter(pca_gg, .data[["cohort"]] == !!cohort_name)[["PC02"]]),
      zoom.size = 0.5,
      horizontal = FALSE
    )

  pop_centre <- purrr::map_df(c("super_pop", "pop") , function(ipop) {
    pca_gg %>%
      dplyr::filter(.data[["cohort"]] != !!cohort_name) %>%
      dplyr::select(-"cohort") %>%
      dplyr::group_by(.data[[ipop]]) %>%
      dplyr::summarise(
        PC01_centre = mean(.data[["PC01"]]),
        PC02_centre = mean(.data[["PC02"]])
      ) %>%
      dplyr::mutate(pop_type = ipop) %>%
      dplyr::rename("pop_closest" = !!ipop) %>%
      dplyr::mutate_if(.predicate = is.factor, .funs = as.character)
  })

  pca_gg_pred <- pca_gg %>%
    dplyr::select(-c(.data[["pop"]], .data[["super_pop"]])) %>%
    dplyr::filter(.data[["cohort"]] == !!cohort_name) %>%
    dplyr::mutate(pop_centre = list(.data[["pop_centre"]])) %>%
    tidyr::unnest("pop_centre") %>%
    dplyr::group_by(.data[["sample"]], .data[["pop_type"]]) %>%
    dplyr::mutate(
      pop_dist = purrr::pmap_dbl(
        .l = list(
          .x = .data[["PC01"]],
          .y = .data[["PC02"]],
          .x_centre = .data[["PC01_centre"]],
          .y_centre = .data[["PC02_centre"]],
          .pop = .data[["pop_closest"]]
        ),
        .f = function(.x, .y, .x_centre, .y_centre, .pop) {
          sqrt((.x - .x_centre)^2 + (.y - .y_centre)^2)
        }
      )
    ) %>%
    dplyr::ungroup()

  pca_gg_pred_best <- pca_gg_pred %>%
    dplyr::group_by(.data[["sample"]], .data[["pop_type"]]) %>%
    dplyr::slice(which.min(.data[["pop_dist"]])) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::ends_with("_centre"), -.data[["pop_dist"]]) %>%
    tidyr::pivot_wider(names_from = "pop_type", values_from = "pop_closest")

  #################
  ### Exporting ###
  #################
  message(message_prefix, "Exporting ...")
  ggplot2::ggsave(
    filename = paste0(output_directory, "/", cohort_name, "_ethnicity.png"),
    plot = p_ethni,
    width = 6.3,
    height = 4.7 * 1.5,
    units = "in",
    dpi = 300
  )

  invisible(
    readr::write_csv(
      x = pca_gg_pred_best,
      path = paste0(output_directory, "/", cohort_name, "_ethnicity.csv")
    )
  )
}
