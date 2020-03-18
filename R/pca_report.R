#' Compute an analysis report using principal component analysis from flashpca tool.
#'
#' The function can be used in a chunk within a Rmarkdown document/script with results="asis" to render the report.
#'
#' @param data A `vector` or `data.frame`. The data on which the PCA has to be performed.
#' @param design A `data.frame`. Additional variables to be used with factorial planes.
#' @param id_var A `character`. The identifier column used to merge the data.
#' @param technical_vars A `vector(character)`. Variables from design to be used with factorial planes.
#' @param n_comp A `numeric`. The number of principal components to be computed.
#' @param fig_n_comp A `numeric`. The number of principal components to be used for figures.
#' @param outliers_component A `logical`. The principal components to be used to outliers detection.
#' @param outliers_threshold A `numeric`. The threshold to define outliers.
#' @param title_level A `numeric`. The markdown title level, *i.e.*, the number of `#` preceding the section.
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#'
#' pca_report(
#'   data = t(mtcars),
#'   design = tibble::rownames_to_column(mtcars, "Sample_ID"),
#'   id_var = "Sample_ID",
#'   technical_vars = c("cyl", "gear", "vs"),
#'   n_comp = 5,
#'   fig_n_comp = 5,
#'   outliers_component = 1:2,
#'   outliers_threshold = 3,
#'   title_level = 2
#' )
#'
pca_report <- function(
  data,
  design,
  id_var = "Sample_ID",
  technical_vars,
  n_comp = 5,
  fig_n_comp = n_comp,
  outliers_component = 1:2,
  outliers_threshold = 3,
  title_level = 2
) {
  message_prefix <- "[rain] "
  message(message_prefix, "PCA started ...")

  if (!inherits(design, "data.frame")) stop(message_prefix, '"design" must be a "data.frame"!')

  design <- design %>%
    dplyr::mutate_at(dplyr::vars(tidyselect::all_of(id_var)), as.character) %>%
    dplyr::filter(.data[[id_var]] %in% !!colnames(data))

  keep_technical <- design %>%
    dplyr::summarise_at(
      .vars = dplyr::vars(tidyselect::all_of(technical_vars)),
      .funs = ~ dplyr::n_distinct(.x) > 1 & dplyr::n_distinct(.x) < length(.x)
    ) %>%
    dplyr::select_if(~ all(isTRUE(.x)), identity) %>%
    colnames()

  variables_excluded <- setdiff(technical_vars, keep_technical)
  if (length(variables_excluded)!=0) {
    message(message_prefix,
      "The following variables have been excluded (null variances or confounding with samples): ",
      glue::glue_collapse(variables_excluded, sep = ", ", last = " and ")
    )
  }

  pca_res <- flashpcaR::flashpca(
    X = t(as.matrix(data)),
    stand = "sd",
    ndim = n_comp
  )

  pca_dfxy <- pca_res[["projection"]] %>%
    tibble::as_tibble(.name_repair = ~ paste0("PC", seq_len(length(.x)))) %>%
    dplyr::mutate(!!id_var := as.character(colnames(data))) %>%
    dplyr::right_join(y = design, by = id_var)

  cat(paste0("\n", paste(rep("#", title_level), collapse = ""), " PCA inertia contribution {-}\n"))
  p <- tibble::tibble(
    y = (pca_res$values / sum(pca_res$values)),
    x = sprintf("PC%02d", seq_along(pca_res$values)),
    cumsum = cumsum(.data[["y"]])
  ) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
    ggplot2::geom_bar(stat = "identity", width = 1, colour = "white", fill = "#3B528BFF", na.rm = TRUE) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = ggplot2::expand_scale(mult = c(0, 0.05))) +
    ggplot2::labs(y = glue::glue("Inertia (% of {n_comp} PCs)"), x = "Principal Components (PCs)")
  print(p)
  cat("\n")

  if (length(keep_technical)>0) {
    cat(paste0("\n", paste(rep("#", title_level), collapse = ""), " PCA factorial planes {- .tabset}\n"))
    for (ivar in keep_technical) {
      cat(paste0("\n", paste(rep("#", title_level + 1), collapse = ""), " ", ivar, " {-}\n"))
      p <- paste0("PC", 1:fig_n_comp) %>%
        utils::combn(2) %>%
        t() %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        stats::setNames(c("X.PC", "Y.PC")) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(
          data = purrr::map2(
            .x = .data[["X.PC"]],
            .y = .data[["Y.PC"]],
            .f = ~ stats::setNames(pca_dfxy[, c(!!ivar, .x, .y)], c(!!ivar, "X", "Y"))
          )
        ) %>%
        tidyr::unnest(.data[["data"]]) %>%
        dplyr::mutate_at(dplyr::vars(ivar), as.character) %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["X"]], y = .data[["Y"]], colour = .data[[ivar]])) +
        ggplot2::geom_hline(yintercept = 0, na.rm = TRUE) +
        ggplot2::geom_vline(xintercept = 0, na.rm = TRUE) +
        ggplot2::geom_point(shape = 4, size = 2, na.rm = TRUE) +
        ggplot2::stat_ellipse(type = "norm", na.rm = TRUE) +
        ggplot2::scale_colour_viridis_d() +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(!!ggplot2::sym("Y.PC")),
          cols = ggplot2::vars(!!ggplot2::sym("X.PC")),
          scales = "free"
        ) +
        ggplot2::guides(colour = ifelse(dplyr::n_distinct(pca_dfxy[[ivar]]) <= 12, "legend", "none"))+
        ggplot2::labs(colour = ivar)
      print(p)
      cat("\n")
    }

    cat(paste0("\n", paste(rep("#", title_level), collapse = ""), " PCA association {-}\n"))
    p <- pca_dfxy %>%
      tidyr:: pivot_longer(names_to = "PC", values_to = "Values", cols = dplyr::num_range("PC", 1:n_comp)) %>%
      dplyr::filter(.data[["PC"]] %in% paste0("PC", 1:fig_n_comp)) %>%
      dplyr::group_by(.data[["PC"]]) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        lm = purrr::map(.data[["data"]], function(data) {
          stats::lm(
            stats::as.formula(paste0("Values ~ ", paste(keep_technical, collapse = " + "))),
            data = data
          ) %>%
            stats::anova() %>%
            tibble::rownames_to_column(var = "term") %>%
            dplyr::filter(.data[["term"]] != "Residuals") %>%
            dplyr::mutate(term = gsub("factor\\((.*)\\)", "\\1", .data[["term"]]))
        })
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data[["PC"]], .data[["lm"]]) %>%
      tidyr::unnest(.data[["lm"]]) %>%
      ggplot2::ggplot(
        mapping = ggplot2::aes(
          x = factor(.data[["PC"]]),
          y = .data[["term"]],
          fill = .data[["Pr(>F)"]]
        )
      ) +
        ggplot2::geom_tile(colour = "white", na.rm = TRUE) +
        ggplot2::geom_text(
          mapping = ggplot2::aes(label = scales::scientific(.data[["Pr(>F)"]], digits = 2)),
          colour = "white",
          size = 3,
          na.rm = TRUE
        ) +
        ggplot2::scale_fill_viridis_c(name = "P-Value", na.value = "grey85", limits = c(0, 0.1)) +
        ggplot2::theme(panel.grid = ggplot2::element_blank()) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::labs(x = "PCA Components", y = NULL)
    print(p)
    cat("\n")
  }


  if (!is.null(outliers_component)) {
    cat(paste0("\n", paste(rep("#", title_level), collapse = ""), " PCA Outliers {-}\n"))
    pcs <- paste0("PC", outliers_component)
    pca_dfxy <- pca_dfxy %>%
      dplyr::mutate(
        dist_centre = sqrt(as.vector(scale(.data[[pcs[1]]]))^2 + as.vector(scale(.data[[pcs[2]]]))^2),
        high = .data[["dist_centre"]] >=
          (stats::median(.data[["dist_centre"]], na.rm = TRUE) +
            !!outliers_threshold * stats::IQR(.data[["dist_centre"]], na.rm = TRUE)),
        low = .data[["dist_centre"]] <=
          (stats::median(.data[["dist_centre"]], na.rm = TRUE) -
             !!outliers_threshold * stats::IQR(.data[["dist_centre"]], na.rm = TRUE)),
        bad_samples_bool = .data[["high"]] | .data[["low"]],
        dist_centre = NULL,
        high = NULL,
        low = NULL,
        bad_samples = factor(ifelse(.data[["bad_samples_bool"]], "BAD", "GOOD"), levels = c("BAD", "GOOD"))
      )

    ivar <- "bad_samples"
    p <- paste0("PC", 1:fig_n_comp) %>%
      utils::combn(2) %>%
      t() %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      tibble::as_tibble() %>%
      stats::setNames(c("X.PC", "Y.PC")) %>%
      dplyr::mutate(
        data = purrr::map2(
          .x = .data[["X.PC"]],
          .y = .data[["Y.PC"]],
          .f = ~ stats::setNames(pca_dfxy[, c(!!ivar, .x, .y)], c(!!ivar, "X", "Y"))
        )
      ) %>%
      tidyr::unnest(.data[["data"]]) %>%
      dplyr::mutate_at(dplyr::vars(tidyselect::all_of(ivar)), as.character) %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data[["X"]], y = .data[["Y"]], colour = .data[[ivar]])) +
      ggplot2::geom_hline(yintercept = 0, na.rm = TRUE) +
      ggplot2::geom_vline(xintercept = 0, na.rm = TRUE) +
      ggplot2::geom_point(shape = 4, size = 2, na.rm = TRUE) +
      ggplot2::stat_ellipse(type = "norm", na.rm = TRUE) +
      ggplot2::scale_colour_viridis_d() +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::facet_grid(
        rows = ggplot2::vars(!!ggplot2::sym("Y.PC")),
        cols = ggplot2::vars(!!ggplot2::sym("X.PC")),
        scales = "free"
      ) +
      ggplot2::guides(colour = ifelse(dplyr::n_distinct(pca_dfxy[[ivar]]) <= 12, "legend", "none")) +
      ggplot2::labs(colour = ivar)
    print(p)
    cat("\n")

  }

  message(message_prefix, "PCA ended.")

  invisible(pca_dfxy)
}
