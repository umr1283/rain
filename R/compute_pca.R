#' Compute the genomic components (and some figures) for ethnicity based on VCF files.
#'
#' @inheritParams estimate_ethnicity
#' @param input_plink A `character`. The path to plink format files (i.e., `.bed`, `.bim` and `.fam` files).
#'
#' @return A `data.frame`.
#'
#' @import data.table
#' @import ggplot2
#' @import gt
#' @import patchwork
#'
#' @export
compute_pca <- function(cohort_name, input_plink, output_directory, ref1kg_population, n_comp = 10) {
  message_prefix <- "[rain] "

  message(message_prefix, "Performing PCA ...", appendLF = TRUE)

  fid_iid <- data.table::fread(
    file = paste0(input_plink, ".fam"),
    colClasses = c("V1" = "character", "V2" = "character"),
    header = FALSE
  )

  pca_res <- flashpcaR::flashpca(input_plink, ndim = n_comp)

  pca_contrib <- sprintf(
    fmt = "PC%02d (%s %%)",
    1:n_comp,
    format(pca_res$values / sum(pca_res$values) * 100, digits = 2, nsmall = 2, trim = TRUE)
  )
  names(pca_contrib) <- sprintf("PC%02d", 1:n_comp)

  pca_dfxy <- data.table::as.data.table(pca_res[["projection"]], keep.rownames = "iid")
  data.table::setnames(
    x = pca_dfxy,
    old = setdiff(names(pca_dfxy), "iid"),
    new = sprintf("PC%02d", as.numeric(gsub("V", "", setdiff(names(pca_dfxy), "iid"))))
  )

  if (all.equal(fid_iid[[1]], fid_iid[[2]])) {
    pca_dfxy[, "iid" := fid_iid[[2]]]
  } else {
    pca_dfxy[, "iid" := paste(fid_iid[[1]], fid_iid[[2]], sep = "/")]
  }

  ref_pop_table <- data.table::fread(
    file = ref1kg_population,
    colClasses = "character",
    header = FALSE,
    col.names = c("iid", "pop", "super_pop", "sex")
  )

  pca_dfxy <- merge(
    x = pca_dfxy,
    y = ref_pop_table[, "cohort" := "1,000 Genomes"],
    by = "iid",
    all = TRUE
  )

  cohort <- NULL
  pca_dfxy[!cohort %in% "1,000 Genomes", (c("pop", "super_pop", "cohort")) := list(cohort_name, cohort_name, cohort_name)]
  pca_dfxy[,
    (c("pop", "super_pop", "cohort")) := lapply(
      X = .SD,
      FUN = function(x) factor(x, levels = c(cohort_name, sort(setdiff(x, cohort_name))))
    ),
    .SDcols = c("pop", "super_pop", "cohort")
  ]

  pop_centre <- pca_dfxy[
    cohort %in% "1,000 Genomes",
    lapply(.SD, mean),
    .SDcols = sprintf("PC%02d", 1:2),
    by = "pop"
  ]

  PC01 <- PC02 <- dist <- which_closest <- pop <- super_pop <- NULL
  pca_dfxy[,
    "pop_closest" := (function(.x, .y) {
      as.character(
        pop_centre[,
          list("dist" = sqrt((.x - PC01)^2 + (.y - PC02)^2)),
          by = "pop"
        ][,
          "which_closest" := dist == min(dist)
        ][
          (which_closest)
        ][["pop"]]
      )
    })(PC01, PC02),
    by = "iid"
  ]

  pca_dfxy <- merge(
    x = pca_dfxy,
    y = unique(ref_pop_table[, list("pop_closest" = pop, "super_pop_closest" = super_pop)]),
    by = "pop_closest",
    all.x = TRUE
  )

  p <- lapply(
    X = c("pop" = "pop", "super_pop" = "super_pop"),
    FUN = function(ipop) {
      ggplot2::ggplot(
        data = pca_dfxy,
        mapping = ggplot2::aes(x = .data[["PC01"]], y = .data[["PC02"]], colour = .data[[ipop]])
      ) +
        ggplot2::theme_light(base_size = 11) +
      	ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 0.5, na.rm = TRUE) +
      	ggplot2::geom_vline(xintercept = 0, linetype = 2, size = 0.5, na.rm = TRUE) +
        ggforce::geom_mark_ellipse(mapping = ggplot2::aes(fill = .data[[ipop]]), con.cap = 0, alpha = 0.1) +
        ggplot2::geom_point(mapping = ggplot2::aes(shape = .data[[ipop]]), na.rm = TRUE) +
      	ggplot2::scale_colour_viridis_d(na.translate = FALSE, drop = FALSE, end = 0.9) +
        ggplot2::scale_fill_viridis_d(na.translate = FALSE, drop = FALSE, end = 0.9) +
        ggplot2::scale_shape_manual(values = c(3, rep(1, length(unique(pca_dfxy[[ipop]]))))) +
        ggplot2::labs(
          x = pca_contrib["PC01"],
          y = pca_contrib["PC02"],
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

  p_zoom <- lapply(
    X = p,
    FUN = function(.p) {
      .p +
        ggplot2::coord_cartesian(
          xlim = range(pca_dfxy[!cohort %in% "1,000 Genomes"][["PC01"]]),
          ylim = range(pca_dfxy[!cohort %in% "1,000 Genomes"][["PC02"]])
        )
    }
  )

  p_subtitle <- paste0(
    "Principal Component Analysis using ",
    format(nrow(data.table::fread(paste0(input_plink, ".bim"))), big.mark = ",", digits = 0),
    " SNPs<br>",
    "With A) population level and B) super population level"
  )

  p_ethni <- patchwork::wrap_plots(
    patchwork::wrap_plots(p[["pop"]], p_zoom[["pop"]]) +
      patchwork::plot_layout(tag_level = "new", guides = "collect"),
    patchwork::wrap_plots(p[["super_pop"]], p_zoom[["super_pop"]]) +
      patchwork::plot_layout(tag_level = "new", guides = "collect")
  ) +
    patchwork::plot_layout(nrow = 2, ncol = 1) +
    patchwork::plot_annotation(
      tag_levels = c("A", "1"),
      title = "Ethnicty Inference Based On 1,000 Genomes Project Data",
      subtitle = gsub("<br>W", " w", p_subtitle),
      theme = ggplot2::theme(plot.title.position = "plot", plot.subtitle = ggtext::element_markdown())
    )

  p_cohort_base <- ggplot2::ggplot(
    data = pca_dfxy[!cohort %in% "1,000 Genomes"],
    mapping = ggplot2::aes(x = .data[["PC01"]], y = .data[["PC02"]])
  ) +
    ggplot2::theme_light(base_size = 11) +
  	ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 0.5, na.rm = TRUE) +
  	ggplot2::geom_vline(xintercept = 0, linetype = 2, size = 0.5, na.rm = TRUE) +
  	ggplot2::scale_colour_viridis_d(na.translate = FALSE, drop = FALSE, begin = 0.10, end = 0.90) +
    ggplot2::scale_fill_viridis_d(na.translate = FALSE, drop = FALSE, begin = 0.10, end = 0.90) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = ifelse(nrow(unique(pca_dfxy[!cohort %in% "1,000 Genomes", "pop_closest"])) > 5, 2, 1)
      )
    ) +
    ggplot2::theme(plot.title.position = "plot") +
    ggplot2::coord_cartesian(
      xlim = range(pca_dfxy[!cohort %in% "1,000 Genomes"][["PC01"]]),
      ylim = range(pca_dfxy[!cohort %in% "1,000 Genomes"][["PC02"]])
    )

  p_cohort <- patchwork::wrap_plots(
    p_cohort_base +
      ggplot2::geom_point(mapping = ggplot2::aes(colour = .data[["pop_closest"]]), shape = 1, na.rm = TRUE) +
      ggforce::geom_mark_ellipse(
        mapping = ggplot2::aes(
          colour = .data[["pop_closest"]],
          group = .data[["pop_closest"]]
        ),
        fill = "transparent",
        alpha = 0.10
      ) +
      ggplot2::labs(
        x = pca_contrib["PC01"],
        y = pca_contrib["PC02"],
        shape = "Population",
        colour = "Population",
        fill = "Population"
      ),
    p_cohort_base +
      ggplot2::geom_point(mapping = ggplot2::aes(colour = .data[["super_pop_closest"]]), shape = 1, na.rm = TRUE) +
      ggforce::geom_mark_ellipse(
        mapping = ggplot2::aes(
          colour = .data[["super_pop_closest"]],
          group = .data[["super_pop_closest"]]
        ),
        fill = "transparent",
        alpha = 0.10
      ) +
      ggplot2::labs(
        x = pca_contrib["PC01"],
        y = pca_contrib["PC02"],
        shape = "Super Population",
        colour = "Super Population",
        fill = "Super Population"
      )
  ) +
    patchwork::plot_layout(nrow = 2, ncol = 1) +
    patchwork::plot_annotation(
      tag_levels = "A",
      title = "Ethnicty Inference Based On 1,000 Genomes Project Data",
      subtitle = p_subtitle,
      theme = ggplot2::theme(plot.title.position = "plot", plot.subtitle = ggtext::element_markdown())
    )

  message(message_prefix, "Exporting ...")
  ggplot2::ggsave(
    filename = file.path(output_directory, paste0(cohort_name, "_ethnicity_1kg.pdf")),
    plot = p_ethni,
    width = 29.7 - 5,
    height = 21 - 5,
    units = "cm",
    dpi = 120
  )
  ggplot2::ggsave(
    filename = file.path(output_directory, paste0(cohort_name, "_ethnicity.pdf")),
    plot = p_cohort,
    width = 16,
    height = 16,
    units = "cm",
    dpi = 120
  )

  invisible(
    data.table::fwrite(
      x = pca_dfxy[
        !cohort %in% "1,000 Genomes",
        .SD,
        .SDcols = c("iid", sprintf("PC%02d", 1:n_comp), "cohort", "pop_closest", "super_pop_closest")
      ],
      file = file.path(output_directory, paste0(cohort_name, "_ethnicity.csv"))
    )
  )
}
