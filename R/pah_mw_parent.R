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
#' (column "molwt_highlow").
#' @param parent_compounds A vector of all parent compounds in analysis that have associated alkylated compounds. If two alkylated compounds
#' can not be distinguished from each other (e.g., have the same molecular weight), and therefore the parents cannot be distinguished,
#' all possible parents should be listed. These names should reflect the exact names that are in your compound_info table.
#' @param alkylated_compounds A vector of all alkylated compounds in the analysis that are associated with the listed parents. These
#' names should reflect the exact names that are in your compound_info table.
#' @param sample_column string, column that contains unique sample identifier
#' @param conc_column string, column that contains sample concentrations
#' @param compound_column string, column that contains compound names.
#' @param statistic string, either "sum" or "average" which will apply to molecular weight ratios.
#' If statistic = "sum", compounds in each category (low vs. high molecular weight) are summed by sample.
#' Because different analyses may measure a different number
#' of chemicals or focus on different types of chemicals, these numbers are hard to interpret both within
#' and across studies. For example, the analysis may measure many more high versus low molecular weight compounds,
#' skewing this ratio. Alternatively, the average concentration of compounds in each category ("average")
#' is a way to standardize the number of compounds in each category. If plot = FALSE, both statistics  are calculated
#' and exported in the first table in the list ("all_dat"), but only the chosen statistic is summarized and exported in
#' the second table that summarizes HMW:LMW for each sample ("summarized_dat"). If plot = TRUE,
#' only the statistic of choice is plotted. The sum is always used for the parent:alkyl calculation, assuming
#' the user has provided only parents that have corresponding alkyls, and vice versa.
#' @param plot logical, whether the data should be exported as two dataframes within a list (FALSE) or
#' plotted as a boxplot (TRUE).
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang sym
#' @examples

pah_mw_parent <- function(compound_info, sample_column, conc_column, compound_column, statistic = "average", plot = TRUE) {

  # test to be sure column and compound names exist
  if (!sample_column %in% names(compound_info)) {
    stop(paste(sample_column, "column is not in compound_info dataframe."))
  }

  if (!conc_column %in% names(compound_info)) {
    stop(paste(conc_column, "column is not in compound_info dataframe."))
  }

  if (!compound_column %in% names(compound_info)) {
    stop(paste(compound_column, "column is not in compound_info dataframe."))
  }

  if (!all(parent_compounds %in% unique(compound_info[[compound_column]]))) {
    stop(paste0("Not all parent compounds provided match compound names in column", compound_column, "."))
  }

  if (!all(alkylated_compounds %in% unique(compound_info[[compound_column]]))) {
    stop(paste0("Not all alkylated compounds provided match compound names in column", compound_column, "."))
  }


  # make column names dplyr-ready
  quo_conc_column <- sym(conc_column)
  quo_sample_column <- sym(sample_column)
  quo_compound_column <- sym(compound_column)

  mw <- compound_info %>%
    filter(!is.na(molwt_highlow)) %>%
    group_by(!!quo_sample_column, molwt_highlow) %>%
    summarize(totals = sum(!!quo_conc_column), means = mean(!!quo_conc_column), counts = n()) %>%
    rename(variable = molwt_highlow)

  compound_info <- compound_info %>%
    mutate(parentAlkyl = case_when(
      !!quo_compound_column %in% parent_compounds ~ 'parent',
      !!quo_compound_column %in% alkylated_compounds ~ 'alkyl'))

  parent <- compound_info %>%
    filter(!is.na(parentAlkyl)) %>%
    group_by(!!quo_sample_column, parentAlkyl) %>%
    summarize(totals = sum(!!quo_conc_column), means = sum(!!quo_conc_column), counts = n()) %>%
    rename(variable = parentAlkyl)

  all <- bind_rows(mw, parent)

  all$variable <- factor(all$variable, levels = c('parent', 'alkyl', 'LMW', 'HMW'))

  all.text <- ungroup(all) %>%
    select(variable, counts) %>%
    distinct() %>%
    arrange(variable)

  my_labels <- paste0(all.text$variable, " \nn = ", all.text$counts)

  p <- ggplot(all, aes_string(x = 'variable', y = ifelse(statistic == "sum", 'totals', 'means'))) +
    scale_y_log10() +
    geom_boxplot() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Compound Type", y = "Sample Concentration (ppb)") +
    scale_x_discrete(labels = my_labels) +
    annotation_logticks(sides = 'l')

  # create sample-specific parent:alkyl and hmw:lmw variables
  all.wide <- select(all, -counts) %>%
    select_(ifelse(statistic == "sum", "-means", "-totals")) %>%
    spread(key = variable, value = ifelse(statistic == "sum", totals, means)) %>%
    mutate(parent_alkyl = parent/alkyl, HMW_LMW = HMW/LMW)

  all$means[all$variable %in% c('alkyl', 'parent')] <- NA

  if (plot == T) {
    out <- p
  } else {
    out <- list(all, all.wide)
    names(out) <- c("all_dat", "summarized_dat")
  }
  return(out)
}
