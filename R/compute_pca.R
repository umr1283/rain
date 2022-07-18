#' Compute the genomic components (and some figures) for ethnicity based on VCF files.
#'
#' @inheritParams estimate_ethnicity
#' @param input_plink A `character`. The path to plink format files (i.e., `.bed`, `.bim` and `.fam` files).
#'
#' @return A `data.frame`.
#'
#' @export
compute_pca <- function(cohort_name, input_plink, output_directory, ref1kg_population, n_comp = 10) {
  PC01 <- PC02 <- dist <- which_closest <- NULL # For global variable warnings
  pop <- super_pop <- pop_closest <- cohort <- NULL # For global variable warnings
  .data <- ggplot2::.data

  message_prefix <- "[rain] "

  message(message_prefix, "Performing PCA ...", appendLF = TRUE)

  fid_iid <- data.table::fread(
    file = paste0(input_plink, ".fam"),
    colClasses = c("V1" = "character", "V2" = "character"),
    header = FALSE
  )

  pca_res <- flashpcaR::flashpca(input_plink, ndim = n_comp)

  pca_contrib <- compute_pve(pca_res[["pve"]])

  pca_dfxy <- compute_vec(pca_res[["vectors"]])

  if (all.equal(fid_iid[[1]], fid_iid[[2]])) {
    pca_dfxy[j = "iid" := fid_iid[[2]]]
  } else {
    pca_dfxy[j = "iid" := paste(fid_iid[[1]], fid_iid[[2]], sep = "/")]
  }

  ref_pop_table <- data.table::fread(
    file = ref1kg_population,
    colClasses = "character",
    header = FALSE,
    col.names = c("iid", "pop", "super_pop", "sex")
  )[j = "cohort" := "1,000 Genomes"]

  pca_dfxy <- merge(x = pca_dfxy, y = ref_pop_table, by = "iid", all = TRUE)

  cohort <- NULL
  pca_dfxy[
    i = !cohort %in% "1,000 Genomes",
    j = (c("pop", "super_pop", "cohort")) := list(cohort_name, cohort_name, cohort_name)
  ][
    j = (c("pop", "super_pop", "cohort")) := lapply(
      X = .SD,
      FUN = function(x) factor(x, levels = c(cohort_name, sort(setdiff(x, cohort_name))))
    ),
    .SDcols = c("pop", "super_pop", "cohort")
  ]

  pop_centre <- pca_dfxy[
    i = cohort %in% "1,000 Genomes",
    j = lapply(.SD, mean),
    .SDcols = sprintf("PC%02d", 1:2),
    by = "pop"
  ]

  pca_dfxy[
    j = "pop_closest" := (function(.x, .y) {
      as.character(
        pop_centre[
          j = list("dist" = sqrt((.x - PC01)^2 + (.y - PC02)^2)),
          by = "pop"
        ][
          j = "which_closest" := dist == min(dist)
        ][(which_closest)][["pop"]]
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

  p1 <- ggplot2::ggplot(
    data = data.table::setnames(
      x = pca_dfxy[
        j = .SD[
          pop %in% .SD[!cohort %in% "1,000 Genomes", unique(pop_closest)] |
            !cohort %in% "1,000 Genomes"
        ]
      ],
      old = names(pca_contrib),
      new = pca_contrib
    )
  ) +
    ggplot2::aes(
      x = .data[[pca_contrib[["PC01"]]]],
      y = .data[[pca_contrib[["PC02"]]]],
      colour = pop, fill = pop, label = pop
    ) +
  	ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 0.5, na.rm = TRUE) +
  	ggplot2::geom_vline(xintercept = 0, linetype = 2, size = 0.5, na.rm = TRUE) +
    ggforce::geom_mark_hull(
      data = ~ .x[cohort %in% "1,000 Genomes"],
      concavity = 2,
      expand = ggplot2::unit(1, "mm"),
      radius = ggplot2::unit(1, "mm"),
      con.cap = ggplot2::unit(0, "mm"),
      con.arrow = ggplot2::arrow(angle = 45, length = ggplot2::unit(2, "mm")),
      label.colour = "grey50",
      con.colour = "grey50",
      label.fontsize = 5,
      label.buffer = ggplot2::unit(2.5, "mm")
    ) +
    ggplot2::geom_point(
      data = ~ .x[!cohort %in% "1,000 Genomes"],
      colour = "black",
      shape = 21
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
  	ggplot2::scale_colour_viridis_d(na.translate = FALSE, drop = FALSE, begin = 0.10, end = 0.90) +
    ggplot2::scale_fill_viridis_d(na.translate = FALSE, drop = FALSE, begin = 0.10, end = 0.90) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position = "none"
    )

  p2 <- p1 + ggplot2::aes(colour = super_pop, fill = super_pop, label = super_pop)

  ragg::agg_png(
    filename = file.path(output_directory, sprintf("%s_ethnicity_1kg.png", cohort_name)),
    width = 16, height = 9, units = "cm", res = 300, scaling = 0.75
  )
    print(
      patchwork::wrap_plots(p1, p2, ncol = 2) +
        patchwork::plot_annotation(
          title = "Ethnicity Inference Based On 1,000 Genomes Project Data",
          subtitle = paste0(
            "Principal Component Analysis using ",
            format(
              x = nrow(data.table::fread(paste0(input_plink, ".bim"))),
              big.mark = ",",
              digits = 1L,
              nsmall = 0L
            ),
            " SNPs, with <b>A</b>) population level and <b>B</b>) super population level"
          ),
          tag_levels = "A",
          theme = ggplot2::theme(
            plot.subtitle = ggtext::element_markdown(face = "italic", size = ggplot2::rel(0.80))
          )
        )
    )
  invisible(grDevices::dev.off())

  invisible(
    data.table::fwrite(
      x = pca_dfxy[
        # i = !cohort %in% "1,000 Genomes",
        j = .SD,
        .SDcols = c(
          "iid", sprintf("PC%02d", 1:n_comp),
          "cohort", "pop_closest", "super_pop_closest", "pop", "super_pop"
        )
      ],
      file = file.path(output_directory, paste0(cohort_name, "_ethnicity_1kg.csv"))
    )
  )
}

#' @keywords internal
compute_pve <- function(pve) {
  pve_dt <- sprintf(
    fmt = "PC%02d (%s %%)",
    seq_along(pve),
    format(pve * 100, digits = 2L, nsmall = 2L, trim = TRUE)
  )
  names(pve_dt) <- sprintf("PC%02d", seq_along(pve))
  pve_dt
}

#' @keywords internal
compute_vec <- function(vec) {
  vec_dt <- data.table::as.data.table(vec, keep.rownames = "iid")
  data.table::setnames(
    x = vec_dt,
    old = setdiff(names(vec_dt), "iid"),
    new = sprintf("PC%02d", as.numeric(gsub("V", "", setdiff(names(vec_dt), "iid"))))
  )
  vec_dt
}
