#' pah_mw_parent
#'
#' Summarizes compound characteristics, including low versus high molecular weight and
#' parent versus alkylated, as one line of evidence for potential sources. Samples that
#' contain high molecular weight and parent compounds indicate a pyrogenic rather than petrogenic
#' PAH source.
#'
#' @export
#' @param compound_info The output dataframe from `get_compound_info`, which contains sample concentrations
#' as well as compound-specific information, including low versus high molecular weight
#' (column "molwt_highlow") and parent versus alkylated (column "parentAlkyl").
#' @param sample_column string, column that contains unique sample identifier
#' @param conc_column string, column that contains sample concentrations
#' @param statistic string, either "sum" or "average". If statistic = "sum", compounds are summed by sample and
#' category (low vs. high, parent vs. alkyl for each sample). Because different analyses may measure a different number
#' of chemicals or focus on different types of chemicals, these numbers are hard to interpret both within
#' and across studies. The raw sums are reported ("sum") along with the number of compounds in each category.
#' Alternatively, the average concentration of compounds in each category ("average") is a way to standardize
#' the number of compounds across categories. If plot = FALSE, both statistics  are calculated and exported
#' in the first table in the list ("all_dat"), but only the chosen statistic is summarized and exported in
#' the second table that summarizes HMW:LMW and parent:alkyl for each sample ("summarized_dat"). If plot = TRUE,
#' only the statistic of choice is plotted.
#' @param plot logical, whether the data should be exported as two dataframes within a list (FALSE) or
#' plotted as a boxplot (TRUE).
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang sym
#' @examples

pah_mw_parent <- function(compound_info, sample_column, conc_column, statistic = "average", plot = TRUE) {
  # make column names dplyr-ready
  quo_conc_column <- sym(conc_column)
  quo_sample_column <- sym(sample_column)

  mw <- compound_info %>%
    filter(!is.na(molwt_highlow)) %>%
    group_by(!!quo_sample_column, molwt_highlow) %>%
    summarize(totals = sum(!!quo_conc_column), means = mean(!!quo_conc_column), counts = n()) %>%
    rename(variable = molwt_highlow)

  parent <- compound_info %>%
    filter(!is.na(parentAlkyl)) %>%
    group_by(!!quo_sample_column, parentAlkyl) %>%
    summarize(totals = sum(!!quo_conc_column), means = mean(!!quo_conc_column), counts = n()) %>%
    rename(variable = parentAlkyl)

  all <- bind_rows(mw, parent)

  all$variable <- factor(all$variable, levels = c('parent', 'alkyl', 'LMW', 'HMW'))

  all.text <- ungroup(all) %>%
    select(variable, counts) %>%
    distinct() %>%
    arrange(variable)

  my_labels <- paste0(all.text$variable, " \nn = ", all.text$counts)

  p <- ggplot(all, aes_string(x = 'variable', y = ifelse(statistic == "sum", 'totals', 'means'))) +
    scale_y_continuous(trans = 'log10') +
    geom_boxplot() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Compound Type", y = ifelse(statistic == "sum", "Sum Sample Concentration (ppb)", "Avg. Sample Concentration (ppb)")) +
    scale_x_discrete(labels = my_labels)

  # create sample-specific parent:alkyl and hmw:lmw variables
  all.wide <- select(all, -counts) %>%
    select_(ifelse(statistic == "sum", "-means", "-totals")) %>%
    spread(key = variable, value = ifelse(statistic == "sum", totals, means)) %>%
    mutate(parent_alkyl = parent/alkyl, HMW_LMW = HMW/LMW)

  if (plot == T) {
    out <- p
  } else {
    out <- list(all, all.wide)
    names(out) <- c("all_dat", "summarized_dat")
  }
  return(out)
}
