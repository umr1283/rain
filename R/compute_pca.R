#' Compute the genomic components (and some figures) for ethnicity based on VCF files.
#'
#' @inheritParams estimate_ethnicity
#' @param input_plink A `character`. The path to plink format files (i.e., `.bed`, `.bim` and `.fam` files).
#'
#' @return A `data.frame`.
#'
#' @export
compute_pca <- function(cohort_name, input_plink, output_directory, ref1kg_population) {
  message_prefix <- "[rain] "

  ######################
  ### Performing PCA ###
  ######################
  message(message_prefix, "Performing PCA ...")

  fid_iid <- readr::read_delim(
    file = paste0(input_plink, ".fam"),
    delim = " ",
    col_types = readr::cols_only(
      X1 = readr::col_character(),
      X2 = readr::col_character()
    ),
    col_names = FALSE
  )

  pca_ndim <- 10

  pca_res <- flashpcaR::flashpca(input_plink, ndim = pca_ndim)
  pca_dfxy <- data.table::as.data.table(pca_res[["projection"]])
  data.table::setnames(
    x = pca_dfxy,
    old = setdiff(names(pca_dfxy), id_var),
    new = sprintf("PC%02d", as.numeric(gsub("V", "", setdiff(names(pca_dfxy), id_var))))
  )

  if (all.equal(fid_iid[[1]], fid_iid[[2]])) {
    pca_dfxy[["sample"]] <- fid_iid[[2]]
  } else {
    pca_dfxy[["sample"]] <- paste(fid_iid[[1]], fid_iid[[2]], sep = "/")
  }

  pca_dfxy <- pca_dfxy %>%
    dplyr::left_join(
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
    ) %>%
    dplyr::mutate(
      cohort = ifelse(is.na(.data[["pop"]]), cohort_name, "1,000 Genomes"),
      pop = ifelse(is.na(.data[["pop"]]), cohort_name, .data[["pop"]]),
      super_pop = ifelse(is.na(.data[["super_pop"]]), cohort_name, .data[["super_pop"]])
    ) %>%
    dplyr::mutate_at(
      .vars = dplyr::vars(.data[["cohort"]], .data[["pop"]], .data[["super_pop"]]),
      .funs = ~ stats::relevel(factor(.x), ref = cohort_name)
    ) %>%
    dplyr::select("sample", dplyr::everything())

  pop_centre <- purrr::map_df(c("super_pop", "pop") , function(ipop) {
    pca_dfxy %>%
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

  pca_gg_pred_best <- pca_dfxy %>%
    dplyr::select(-c(.data[["pop"]], .data[["super_pop"]])) %>%
    dplyr::filter(.data[["cohort"]] == !!cohort_name) %>%
    dplyr::mutate(pop_centre = list(pop_centre)) %>%
    tidyr::unnest("pop_centre") %>%
    dplyr::group_by(.data[["sample"]], .data[["pop_type"]]) %>%
    dplyr::mutate(
      pop_dist = purrr::pmap_dbl(
        .l = list(
          .x = .data[["PC01"]], .y = .data[["PC02"]],
          .x_centre = .data[["PC01_centre"]], .y_centre = .data[["PC02_centre"]]
        ),
        .f = function(.x, .y, .x_centre, .y_centre) {
          sqrt((.x - .x_centre)^2 + (.y - .y_centre)^2)
        }
      )
    ) %>%
    dplyr::slice(which.min(.data[["pop_dist"]])) %>%
    dplyr::ungroup() %>%
    dplyr::select(-dplyr::ends_with("_centre"), -.data[["pop_dist"]]) %>%
    tidyr::pivot_wider(names_from = "pop_type", values_from = "pop_closest")

  p <- purrr::map(.x = c("pop" = "pop", "super_pop" = "super_pop"), .f = function(ipop) {
    ggplot2::ggplot(
      data = pca_dfxy,
      mapping = ggplot2::aes(x = .data[["PC01"]], y = .data[["PC02"]], colour = .data[[ipop]])
    ) +
      ggplot2::theme_light(base_size = 12) +
    	ggplot2::geom_hline(yintercept = 0, linetype = 1, size = 0.5, na.rm = TRUE) +
    	ggplot2::geom_vline(xintercept = 0, linetype = 1, size = 0.5, na.rm = TRUE) +
      ggforce::geom_mark_ellipse(mapping = ggplot2::aes(fill = .data[[ipop]]), con.cap = 0, alpha = 0.1) +
      ggplot2::geom_point(mapping = ggplot2::aes(shape = .data[[ipop]]), na.rm = TRUE) +
    	ggplot2::scale_colour_viridis_d(na.translate = FALSE, drop = FALSE, end = 0.9) +
      ggplot2::scale_fill_viridis_d(na.translate = FALSE, drop = FALSE, end = 0.9) +
      ggplot2::scale_shape_manual(values = c(3, rep(1, length(unique(pca_dfxy[[ipop]]))))) +
      ggplot2::labs(
        shape = c("super_pop" = "Super Population", "pop" = "Population")[ipop],
        colour = c("super_pop" = "Super Population", "pop" = "Population")[ipop],
        fill = c("super_pop" = "Super Population", "pop" = "Population")[ipop]
      ) +
      ggplot2::facet_grid(cols = ggplot2::vars(.data[["cohort"]])) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          ncol = 2,
          label.theme = ggplot2::element_text(size = 8),
          keyheight = ggplot2::unit(9, "point")
        )
      )
  })

  p_zoom <- purrr::map(p, ~ .x + ggplot2::coord_cartesian(
    xlim = range(dplyr::filter(pca_dfxy, .data[["cohort"]] == !!cohort_name)[["PC01"]]),
    ylim = range(dplyr::filter(pca_dfxy, .data[["cohort"]] == !!cohort_name)[["PC02"]])
  ))

  p_ethni <- patchwork::wrap_plots(
    patchwork::wrap_plots(p[["super_pop"]], p_zoom[["super_pop"]]) +
      patchwork::plot_layout(tag_level = "new", guides = "collect"),
    patchwork::wrap_plots(p[["pop"]], p_zoom[["pop"]]) +
      patchwork::plot_layout(tag_level = "new", guides = "collect")
  ) +
    patchwork::plot_layout(nrow = 2, ncol = 1) +
    patchwork::plot_annotation(
      tag_levels = c("A", "1"),
      title = "Estimated Ethnicity",
      subtitle = paste(
        "Based on Principal Component Analysis using",
        scales::comma(nrow(data.table::fwrite(paste0(input_plink, ".bim")))),
        "SNPs."
      )
    )


  #################
  ### Exporting ###
  #################
  message(message_prefix, "Exporting ...")
  ggplot2::ggsave(
    filename = file.path(output_directory, paste0(cohort_name, "_ethnicity.pdf")),
    plot = p_ethni,
    width = 29.7 - 5,
    height = 21 - 5,
    units = "cm",
    dpi = 120
  )

  invisible(
    readr::write_csv(
      x = pca_gg_pred_best,
      path = file.path(output_directory, paste0(cohort_name, "_ethnicity.csv"))
    )
  )
}
