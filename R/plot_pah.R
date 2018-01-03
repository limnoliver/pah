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
#' @importFrom ggplot2
#' @importFrom dplyr order_by
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom forcats fct
#' @examples
#'
plot_pah <- function(pah_dat, conc_column = "Value", sample_id_column = "Sample",
                     compound_column = "Parameter", compound_plot = "Total PAH",
                     color_column = NA, group_column = NA) {
  pah_dat_temp <- filter(pah_dat, compound_column %in% compound_plot)
  if (anyNA(group_column)) {
    p <- ggplot(pah_dat_temp, aes(x = reorder(sample_id_column, Value), y = Value)) +
      geom_bar(stat="identity", position="identity", colour="black")
    if (!is.na(color_column)) {
    p <- p + aes(fill = color_column)
    }

    if (length(compound_column) > 1) {
      p <- p + facet_wrap(~compound_column, ncol = ifelse(length(compound_plot) > 3, 2, 1))
    }
  } else {


  }
}
