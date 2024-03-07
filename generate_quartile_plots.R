# Make Quartile Plots

library(tidyverse)
library(creedenzymatic)

quartile_figure <- function(df, grouping = "KinaseFamily") {
  df |>
    dplyr::select(hgnc_symbol, one_of(grouping), Qrt, Method) |>
    tidyr::pivot_wider(
      names_from = Method,
      values_from = Qrt,
      values_fn = unique
    ) |>
    tidyr::pivot_longer(where(is.numeric), names_to = "Method", values_to = "Qrt") |>
    dplyr::mutate(
      present = ifelse(is.na(Qrt), "No", "Yes"),
      Qrt = ifelse(present == "No", 2, Qrt),
      present = as.factor(present),
      Qrt = as.factor(Qrt),
      Method = as.factor(Method)
    ) |>
    ggplot2::ggplot(ggplot2::aes(hgnc_symbol, Method)) +
    ggplot2::geom_point(ggplot2::aes(size = Qrt, shape = present)) +
    ggplot2::scale_size_manual(values = c(
      `4` = 4,
      `3` = 3,
      `2` = 2,
      `1` = 1
    )) +
    ggplot2::theme_bw() +
    { # nolint: brace_linter.
      if (grouping == "subfamily") {
        ggplot2::facet_grid(. ~ subfamily, scales = "free", space = "free")
      } else if (grouping == "group") {
        ggplot2::facet_grid(. ~ group, scales = "free", space = "free")
      } else {
        ggplot2::facet_grid(. ~ KinaseFamily,
          scales = "free",
          space = "free"
        )
      }
    } +
    ggplot2::scale_shape_manual(values = c(Yes = 19, No = 1)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 30,
      size = 7.5,
      vjust = 0.7
    ), axis.ticks = ggplot2::element_blank(), legend.position = "bottom") +
    ggplot2::labs(x = "", y = "") +
    ggplot2::guides(shape = "none")
}

generate_quartile_plot <- function(datafile) {
  creeden_data <-
    readr::read_csv(file.path("results", datafile), show_col_types = FALSE)

  sig_kinases <- creeden_data |>
    extract_top_kinases(min_qrt = 4, min_counts = 2) |>
    c("AATK", "LTK", "INSRR", "TNK1", "EPHA6")

  creeden_data |>
    dplyr::filter(hgnc_symbol %in% sig_kinases) |>
    quartile_figure()
}

creedenzymatic_files <- list.files("results", "creedenzymatic") |>
  set_names(~ str_remove(.x, "_.*")) |>
  map(generate_quartile_plot) |>
  imap(~ ggsave(
    str_glue("{.y}-creedenzymatic.svg"),
    path = "figures",
    plot = .x,
    width = 10,
    height = 3
  ))
