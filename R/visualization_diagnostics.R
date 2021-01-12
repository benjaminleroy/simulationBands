#' compare all densities visually (true, fit, conformal_pdf)
#'
#' @param stepwise_scaling new scaling for each level set
#' @param individual_df_rich data for individual x's potential y values and
#' true, fit values
#' @param cumulative_comparisons_df cumlative information for level sets
#' @param scale boolean if we should normalize the solution
#'
#' @return visualization plot
vis_all_densities <- function(stepwise_scaling, individual_df_rich,
                              cumulative_comparisons_df,
                              scale = T){
  cumulative_comparisons_df$multiplier <- stepwise_scaling

  all_df2 <- individual_df_rich %>%
    dplyr::left_join(cumulative_comparisons_df,
              by = c("g_id", "g_id2")) %>%
    dplyr::mutate(conformal_pdf_density = .data$fit_density * .data$multiplier)

  if (scale){
    # delta_y
    delta_y <- unique(diff(individual_df_rich$yy))
    assertthat::assert_that(length(delta_y) == 1 |
                              all(abs(diff(delta_y) < 1e-10)),
                            msg = paste("difference between y values not the",
                                        "same - please correct"))
    if (length(delta_y) > 1){
      delta_y <- delta_y[1]
    }

    total <- sum(all_df2$conformal_pdf_density)*delta_y
    all_df2 <- all_df2 %>%
      dplyr::mutate(conformal_pdf_density = .data$conformal_pdf_density / total)
  }

  scaled_string <- c(""," (update correctly scaled)")[scale + 1]

  all_df2_long <- all_df2 %>%
    dplyr::select(.data$yy,
                  .data$true_density,
                  .data$fit_density,
                  .data$conformal_pdf_density) %>%
    dplyr::rename(conformal_density = "conformal_pdf_density",
           fit_density = "fit_density") %>%
    tidyr::pivot_longer(cols = c(.data$true_density,
                                 .data$fit_density,
                                 .data$conformal_density),
                        names_to = "name",
                        values_to = "value"
    ) %>%
    dplyr::mutate(name = stringr::str_replace(.data$name, "_", " "))

  vis2 <-  all_df2_long %>% ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(x = .data$yy,
                   y = .data$value,
                   color = .data$name), size = 1) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$yy,
                   y = .data$value,
                   color = .data$name,
                   size = .data$name)) +
    ggplot2::scale_size_manual
  (values = c("true density" = 0, "fit density" = 0, "conformal density" = 2)) +
    ggplot2::scale_color_manual(
      values = c("true density" = "black", "fit density" = "red",
                 "conformal density" = "purple")) +
    ggplot2::labs(x = "range of y values",
                  y = paste0("probability distribution", scaled_string),
                  color = "density type",
                  size = "density type")

  return(vis2)

}

#' visualize scaling values
#'
#' @param stepwise_scaling scaling values (for each level)
#' @param x x axis values, e.g. the level set values, proportion of mass
#' values, etc. The default is just an index.
#'
#' @return \code{ggplot} visual
#' @export
vis_stepwise_scaling_only <- function(stepwise_scaling,
                                      x = 1:length(stepwise_scaling)){
  vis1 <- data.frame(multiplier = stepwise_scaling,
                     x = x) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$multiplier)) +
    ggplot2::geom_line(ggplot2::aes(x = .data$x, y = .data$multiplier)) +
    ggplot2::labs(x = "threshold for cde",
         y = "scaling")
  return(vis1)
}
