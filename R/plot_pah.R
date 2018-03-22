#' plot_pah
#'
#' Create a bar plot of PAH concentrations. The function allows you to group the data by one
#' or more variables, and decide which grouping variables will determine order or color.
#'
#' @export
#' @param pah_dat dataframe with PAH concentrations. It is recommended that entire samples with censored observations
#' (e.g., below detection limit) are removed from this data frame prior to calculating the difference between source and samples.
#' This function assumes these samples have already been removed.
#' @param conc_column column that contains PAH concentrations
#' @param sample_id_column column name that contains unique sample id
#' @param compound_column column name that contains PAH compound names. This can also include other chemicals, or sums of chemicals.
#' @param compound_plot a vector of strings identifying which compounds from compound_column to include in the plot. If more than one compound is given, the plot will be faceted by compound.
#' @param color_column a column with group variable by which to color code bars
#' @param group_column a column by which to group and order the bars.
#' @param order_column a column by which to order bars within groups. If left NA, this will default
#' to the conc_column so bars are ordered low to high within groups or across all
#' bars if no groups are given.
#' @param conc_units character string, units for PAH concentration that will be used
#' in labels on plots
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom rlang sym
#' @examples
#'
plot_pah <- function(pah_dat, conc_column = 'Value', sample_id_column = 'Sample',
                     compound_column = 'Parameter', compound_plot = "Total PAH",
                     color_column = NA, group_column = NA, order_column = NA, conc_units = "ppb") {

  quo_compound_column <- sym(compound_column)
  quo_conc_column <- sym(conc_column)
  quo_sample_id_column <- sym(sample_id_column)
  quo_group_column <- sym(group_column)

  order_column <- ifelse(is.na(order_column), conc_column, order_column)
  quo_order_column <- sym(order_column)

  pah_dat_temp <- filter(pah_dat, (!!quo_compound_column) %in% compound_plot)

  barchart_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                          legend.position = "top", panel.spacing=unit(5,"pt"), strip.background = element_blank(),
                          strip.text = element_text(margin = margin(4,10,4,10, unit="pt")))

  yvals <- c(round(min(log10(pah_dat_temp[[conc_column]])), 0), round(max(log10(pah_dat_temp[[conc_column]])), 0))
  ybreaks <- yvals[1]:yvals[2]
  yticklabs <- 10^(yvals[1]:yvals[2])

  if (anyNA(group_column)) {
    # plotting if there is no grouping column given
    p <- ggplot(pah_dat_temp, aes_string(x = paste0('reorder(',sample_id_column,',', order_column,')'), y = conc_column)) +
      geom_bar(stat="identity", position="identity", colour="black") +
      barchart_theme +
      labs(x = "", y = paste0("Concentration (", conc_units, ")"))
    if (!is.na(color_column)) {
      # add color if color column is given
    p <- p + aes_string(fill = color_column)
    }

    if (length(compound_column) > 1) {
      # facet by compound if more than one compound ID is given
      p <- p + facet_wrap(as.formula(paste('~',compound_column)), ncol = ifelse(length(compound_plot) > 3, 2, 1))
    }
  } else {
    # now handle groups if groups are given
    # first group data, calculate mean by group to determine group order
    temp.order <- pah_dat_temp %>%
      group_by(!!quo_group_column) %>%
      summarize(max = max(!!quo_conc_column)) %>%
      arrange(max)

    pah_dat_temp[group_column] <- factor(pah_dat_temp[[group_column]],
                                                         levels = temp.order[[group_column]])

    pah_dat_temp <- pah_dat_temp %>%
      ungroup() %>%
      mutate(thresholdPEC = ifelse((!!quo_compound_column) == "Priority Pollutant PAH", 22800, NA),
             thresholdTEC = ifelse((!!quo_compound_column) == "Priority Pollutant PAH", 1610, NA))

    dummy.dat <- pah_dat_temp %>%
      select(!!quo_group_column, !!quo_compound_column, thresholdPEC, thresholdTEC) %>%
      distinct()

    p <- ggplot(pah_dat_temp, aes_string(x = paste0('reorder(',sample_id_column,',',order_column,")"), y = conc_column)) +
      geom_bar(stat = 'identity', width = 0.8, position = position_dodge(width = 2), colour="black" ) +
      # need to test this: does this work of PARAM_SYNYNOYM is one variable?
      geom_hline(data = dummy.dat, aes(yintercept = thresholdPEC)) +
      geom_hline(data = dummy.dat, aes(yintercept = thresholdTEC)) +
      #geom_text(aes(y = 50000, label = c(rep(NA, 4), "TEC", rep(NA, 66)))) +
      facet_grid(as.formula(paste(".~", group_column)),
                 space = "free_x", scales = 'free_x', switch = 'x') +
      barchart_theme +
      labs(x = "", y = paste0("Concentration (", conc_units, ")")) +
      scale_y_continuous(breaks = 10^ybreaks, labels = yticklabs, trans = 'log10') +
      theme(strip.text.x = element_text(angle = 90, hjust = 1),
            panel.spacing = unit(0.2, "lines"), strip.placement = 'outside')


      if (!is.na(color_column)) {
        p <- p + aes_string(fill = color_column)
      }

  }
  return(p)
}
