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
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @examples
#'
plot_pah <- function(pah_dat, conc_column = "Value", sample_id_column = "Sample",
                     compound_column = "Parameter",
                     color_column = NA, group_column = NA) {
  if (is.na(color_column) & is.na(group_column)) {
    p <- ggplot(pah_dat, aes(x = conc_column, y = sample_id_column)) +
      geom_bar(stat="identity", position="identity", colour="black")

  }
}
