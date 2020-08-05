#' Compute an analysis report using principal component analysis from flashpca tool.
#'
#' The function can be used in a chunk within a Rmarkdown document/script with results="asis" to render the report.
#'
#' @param data A `vector` or `data.frame`. The numeric data on which the PCA has to be performed.
#' @param design A `data.frame`. Additional variables to be used with factorial planes.
#' @param id_var A `character`. The identifier column used to merge the data.
#' @param technical_vars A `vector(character)`. Variables from design to be used with factorial planes.
#' @param n_comp A `numeric`. The number of principal components to be computed.
#' @param fig_n_comp A `numeric`. The number of principal components to be used for figures.
#' @param outliers_threshold A `numeric`. The threshold to define outliers.
#' @param title_level A `numeric`. The markdown title level, *i.e.*, the number of `#` preceding the section.
#'
#' @return A `data.frame`.
#' @export
#'
#' @import data.table
#' @import ggplot2
#' @import gt
#' @import patchwork
#'
#' @examples
#'
#' if (interactive()) {
#'   pca_report(
#'     data = t(mtcars),
#'     design = as.data.table(mtcars, keep.rownames = "Sample_ID"),
#'     id_var = "Sample_ID",
#'     technical_vars = c("cyl", "gear", "vs"),
#'     n_comp = 5,
#'     fig_n_comp = 5,
#'     outliers_threshold = 3,
#'     title_level = 0
#'   )
#' }
#'
pca_report <- function(
  data,
  design,
  id_var = "Sample_ID",
  technical_vars,
  n_comp = 10,
  fig_n_comp = n_comp,
  outliers_threshold = 1.5,
  title_level = 3
) {
  .data <- ggplot2::.data
  `%>%` <- gt::`%>%`

  message_prefix <- "[rain] "
  message(message_prefix, "PCA started ...", appendLF = TRUE)

  section_prefix <- paste(c("\n", rep("#", title_level)), collapse = "")

  pca_theme <- ggplot2::theme(
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.title = ggtext::element_markdown(),
    plot.subtitle = ggtext::element_markdown(face = "italic", size = ggplot2::rel(0.8)),
    plot.caption = ggtext::element_markdown(face = "italic", size = ggplot2::rel(0.5)),
    axis.text.x = ggtext::element_markdown()
  )

  pca_methylation <- data[rowSums(is.na(data)) == 0, ]
  design_idx <- design[[id_var]] %in% colnames(data)
  pca_phenotypes <- data.table::as.data.table(design[design_idx, ])

  n_comp <- min(c(n_comp, 10, ncol(pca_methylation)))
  fig_n_comp <- min(c(n_comp, 3, ncol(pca_methylation)))

  keep_technical <- names(which(sapply(pca_phenotypes[,
    lapply(.SD, function(x) (data.table::uniqueN(x) > 1 & data.table::uniqueN(x) < length(x))) | is.numeric(x),
    .SDcols = technical_vars
  ], isTRUE)))

  variables_excluded <- setdiff(technical_vars, keep_technical)
  if (length(variables_excluded) != 0) {
    cat(
      "The following variables have been excluded (null variances or confounding with samples):\n",
      paste("+", variables_excluded),
      "\n",
      sep = "\n"
    )
  }

  pca_res <- flashpcaR::flashpca(X = t(pca_methylation), stand = "sd", ndim = n_comp)

  pca_dfxy <- data.table::as.data.table(pca_res[["vectors"]], keep.rownames = id_var)
  data.table::setnames(
    x = pca_dfxy,
    old = setdiff(names(pca_dfxy), id_var),
    new = sprintf("PC%02d", as.numeric(gsub("V", "", setdiff(names(pca_dfxy), id_var))))
  )
  pca_dfxy <- merge(x = pca_dfxy, y = pca_phenotypes, by = id_var)

  p_inertia <- ggplot2::ggplot(
    data = data.table::data.table(
      y = pca_res[["pve"]],
      x = sprintf("PC%02d", seq_along(pca_res[["pve"]]))
    )[1:fig_n_comp],
    mapping = ggplot2::aes(
      x = paste0(
        .data[["x"]],
        "<br><i style='font-size:5pt;'>(",
          format(.data[["y"]] * 100, digits = 2, nsmall = 2),
        " %)</i>"
      ),
      y = .data[["y"]]
    )
  ) +
    ggplot2::geom_col(width = 1, colour = "white", fill = "#21908CFF", na.rm = TRUE) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste(format(x * 100, digits = 2, nsmall = 2), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      x = "Principal Components",
      y = "Contribution"
    ) +
    pca_theme

  if (length(keep_technical) > 0) {
    cat(section_prefix, "## Association Tests\n\n", sep = "")
    asso_dt <- data.table::melt(
      data = pca_dfxy,
      measure.vars = grep("^PC[0-9]+$", names(pca_dfxy), value = TRUE),
      variable.name = "pc",
      value.name = "values"
    )
    asso_dt <- asso_dt[asso_dt[["pc"]] %in% sprintf("PC%02d", 1:n_comp)]
    asso_dt <- asso_dt[,
      data.table::as.data.table(
        stats::anova(
          stats::lm(
            formula = stats::as.formula(paste0("values ~ ", paste(keep_technical, collapse = " + "))),
            data = .SD
          )
        ),
        keep.rownames = "term"
      ),
      by = "pc"
    ]
    p_association <- ggplot2::ggplot(
      data = asso_dt[asso_dt[["term"]] != "Residuals"],
      mapping = ggplot2::aes(x = factor(.data[["pc"]]), y = .data[["term"]], fill = .data[["Pr(>F)"]])
    ) +
      ggplot2::geom_tile(colour = "white", na.rm = TRUE) +
      ggtext::geom_richtext(
        mapping = ggplot2::aes(
          label = gsub(
            pattern = "(.*)e([-+]*)0*(.*)",
            replacement = "\\1<br>&times;<br>10<sup>\\2\\3</sup>",
            x = format(.data[["Pr(>F)"]], digits = 2, nsmall = 2, scientific = TRUE)
          )
        ),
        colour = "white",
        fill = NA,
        label.colour = NA,
        size = 2.5,
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_viridis_c(name = "P-Value", na.value = "grey85", end = 0.75, limits = c(0, 0.1)) +
      ggplot2::theme(panel.grid = ggplot2::element_blank()) +
      ggplot2::scale_x_discrete(
        expand = c(0, 0),
        labels = function(x) {
          paste0(
            x, "<br><i style='font-size:5pt;'>(",
            format(
              x = pca_res[["pve"]][as.numeric(gsub("PC", "", x))] * 100,
              digits = 2,
              nsmall = 2
            ),
            " %)</i>"
          )
        }
      ) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::labs(
        x = "Principal Components",
        y = "Variables",
        title = "Association Tests Between Variables And Principal Components",
        caption = paste0(
          "Contribution computed using ", n_comp," principal components.<br>",
          "Variables are tested against principal components using ANOVA."
        )
      ) +
      pca_theme

    print(p_association)
    cat("\n")

    cat(section_prefix, "## Factorial Planes\n\n", sep = "")
    for (ivar in keep_technical) {
      cat(section_prefix, "### `", ivar, "`\n\n", sep = "")
      p <- patchwork::wrap_plots(
        c(
          apply(
            X = utils::combn(sprintf("PC%02d", 1:fig_n_comp), 2),
            MARGIN = 2,
            FUN = function(x) {
              ggplot2::ggplot(
                data = pca_dfxy[, .SD, .SDcols = c(ivar, x)],
                mapping = ggplot2::aes(x = .data[[x[1]]], y = .data[[x[2]]], colour = .data[[ivar]])
              ) +
                ggplot2::geom_hline(yintercept = 0, linetype = 2, na.rm = TRUE) +
                ggplot2::geom_vline(xintercept = 0, linetype = 2, na.rm = TRUE) +
                ggplot2::geom_point(na.rm = TRUE) +
                ggplot2::stat_ellipse(type = "norm", na.rm = TRUE, show.legend = FALSE) +
                {
                  if (is.numeric(pca_dfxy[[ivar]])) {
                    ggplot2::scale_colour_viridis_c(
                      name = NULL,
                      begin = 0,
                      end = 0.75
                    )
                  } else {
                    ggplot2::scale_colour_viridis_d(
                      name = NULL,
                      begin = if (data.table::uniqueN(pca_dfxy[[ivar]]) == 2) 0.25 else 0,
                      end = 0.75,
                      guide = ggplot2::guide_legend(override.aes = list(size = 4))
                    )
                  }
                }
            }
          ),
          list(p_inertia)
        ),
        guides = "collect"
      ) +
        patchwork::plot_annotation(
          title = paste0("Structure Detection For: '<i>", ivar, "</i>'"),
          caption = paste0("Contribution computed using ", n_comp," principal components."),
          tag_levels = "A",
          theme = pca_theme
        )

      print(p)
      cat("\n")
    }
  }

  cat(section_prefix, "## Outliers Detection\n\n", sep = "")
  pca_dfxy[,
    "dist_centre" := rowSums(sapply(.SD, function(x) as.vector(scale(x))^2)),
    .SDcols = sprintf("PC%02d", 1:fig_n_comp)
  ]
  pca_dfxy[,
    "is_outlier" := lapply(.SD, function(x) {
      factor(
        x = x > (
          stats::quantile(x, 0.75, na.rm = TRUE) + outliers_threshold * stats::IQR(x, na.rm = TRUE)
        ),
        levels = c(FALSE, TRUE),
        labels = c("No", "Yes")
      )
    }),
    .SDcols = "dist_centre"
  ]
  p <- patchwork::wrap_plots(
    c(
      apply(
        X = utils::combn(sprintf("PC%02d", 1:fig_n_comp), 2),
        MARGIN = 2,
        FUN = function(x) {
          ggplot2::ggplot(
            data = pca_dfxy[, .SD, .SDcols = c("is_outlier", x)],
            mapping = ggplot2::aes(
              x = .data[[x[1]]],
              y = .data[[x[2]]],
              colour = .data[["is_outlier"]],
              shape = .data[["is_outlier"]]
            )
          ) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2, na.rm = TRUE) +
            ggplot2::geom_vline(xintercept = 0, linetype = 2, na.rm = TRUE) +
            ggplot2::geom_point(na.rm = TRUE) +
            ggplot2::stat_ellipse(type = "norm", na.rm = TRUE, show.legend = FALSE) +
            ggplot2::scale_colour_viridis_d(
              name = "Outlier",
              begin = 0.25,
              end = 0.75,
              guide = ggplot2::guide_legend(override.aes = list(size = 4))
            ) +
            ggplot2::scale_shape_manual(name = "Outlier", values = c(1, 4))
        }
      ),
      list(p_inertia)
    ),
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = "Outliers Detection In Factorial Planes",
      caption = paste0(
        "Outliers defined for a Euclidean distance from cohort centroid (based on the principal components up to ", fig_n_comp, ")<br>",
        "higher than ", outliers_threshold, " times the interquartile range above the 75<sup>th</sup> percentile.<br>",
        "Contribution computed using ", n_comp," principal components."
      ),
      tag_levels = "A",
      theme = pca_theme
    )
  print(p)

  gt::gt(
    data = pca_dfxy[
      pca_dfxy[["is_outlier"]] == "Yes",
      .SD,
      .SDcols = c(id_var, technical_vars, sprintf("PC%02d", 1:fig_n_comp))
    ],
    auto_align = "center"
  ) %>%
    gt::tab_header(title = "Samples Identified As Possible Outliers") %>%
    gt::fmt_number(columns = sprintf("PC%02d", 1:fig_n_comp), decimals = 2) %>%
    gt::opt_row_striping() %>%
    gt::opt_all_caps() %>%
    print()

  message(message_prefix, "PCA ended.", appendLF = TRUE)

  invisible(pca_dfxy)
}
