#' plot_pah
#'
#' Create a bar plot of PAH concentrations. The function allows you to group the data by one
#' or more variables, and decide which grouping variables will determine order or color.
#'
#' @export
#' @param pah_dat dataframe with PAH concentrations
#' @param conc_column column that contains PAH concentrations
#' @param sample_id_column column name that contains unique sample id
#' @param compound_column column name that contains PAH compound names. This can also include other chemicals, or sums of chemicals.
#' @param compound_plot a vector of strings identifying which compounds from compound_column to include in the plot. If more than one compound is given, the plot will be faceted by compound.
#' @param color_column a column with group variable by which to color code bars
#' @param group_column a vector of one or more columns with grouping variables that will be used to order the bars. Groups should be listed from highest to lowest order.
#' @import ggplot2
#' @import dplyr
#' @examples
#'
plot_pah <- function(pah_dat, conc_column = "Value", sample_id_column = "Sample",
                     compound_column = "Parameter", compound_plot = "Total PAH",
                     color_column = NA, group_column = NA, conc_units = "ppb") {

  pah_dat_temp <- filter(pah_dat, compound_column %in% compound_plot)
  barchart_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1),
                          legend.position = "bottom", panel.spacing=unit(5,"pt"), strip.background = element_blank(),
                          strip.text = element_text(margin = margin(4,10,4,10, unit="pt")))
  ybreaks <- c(round(min(log10(pah_dat_temp[[conc_column]])), 0), round(max(log10(pah_dat_temp[[conc_column]])), 0))
  yticklabs <- 10^(ybreaks[1]:ybreaks[2])

  if (anyNA(group_column)) {
    # plotting if there is no grouping column given
    p <- ggplot(pah_dat_temp, aes(x = reorder(sample_id_column, Value), y = Value)) +
      geom_bar(stat="identity", position="identity", colour="black") +
      barchart_theme() +
      labs(x = "", y = paste0("Concentration (", conc_units, ")"))
    if (!is.na(color_column)) {
      # add color if color column is given
    p <- p + aes(fill = color_column)
    }

    if (length(compound_column) > 1) {
      # facet by compound if more than one compound ID is given
      p <- p + facet_wrap(~compound_column, ncol = ifelse(length(compound_plot) > 3, 2, 1))
    }
  } else {
    # now handle groups if groups are given
    # first group data, calculate mean by group to determine group order
    temp.order <- pah_dat_temp %>%
      group_by(group_column[1]) %>%
      summearize(mean = mean(conc_column)) %>%
      arrange(mean)

    pah_dat_temp[group_column[1]] <- factor(pah_dat_temp[[group_column[1]]],
                                                         levels = temp.order[[group_column[1]]])
    p <- ggplot(pah_dat_temp, aes(x = sample_id_column, y = Value)) +
      geom_bar(width = 0.8, position = position_dodge(width = 2), colour="black" ) +
      # need to test this: does this work of PARAM_SYNYNOYM is one variable?
      facet_grid(PARAM_SYNONYM~State, space = "free_x", scales = 'free_x') +
      barchart_theme() +
      labs(x = "", y = paste0("Concentration (", conc_units, ")")) +
      scale_y_log10(breaks = seq(round(min(pah_dat_temp[[conc_column]]), 0)),
                                 labels = c(10, 100,1000, 10000,100000))

      if (!is.na(color_column)) {
        p <- p + aes(fill = color_column)
      }

  }
}
