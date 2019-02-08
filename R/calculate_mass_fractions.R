#' calc_mass_fractions
#'
#' A mass balance approach to ruling out sources of contamination. This function uses published
#' source PAH concentrations, and calculates the mass fraction (as a percent) of the source that
#' would be required to account for the reported concentration in the samples. For sources with
#' multiple observations, the mean value is used. Mass fraction can be calculated for each sample,
#' or for summary statistics across all samples (min, quartiles, mean, max).
#'
#' @export
#' @param compound_info The output dataframe from `get_compound_info`, which contains sample
#' concentrations as well as compound-specific information, including whether the compound is one of
#' 16 EPA priority compounds and compound-specific toxicity.
#' @param sample_column string, column that contains unique sample identifier
#' @param conc_column string, column that contains sample concentrations
#' @param conc_unit string, the units of PAH concentrations,
#' either "ppb" (ug/kg) or "ppm" (mg/kg).
#' @param compound_column string, column that contains compound names.
#' @param calc_type how to calculate mass fractions, either for each individual sample ('by_sample'),
#' or by summary statistics across all samples 'summary'. Summary calculates mass fractions for all quartiles,
#' minimum, mean, and maximum of sample concentrations.
#' @param plot logical, whether 'by_sample' should be summarized as a tile plot rather than table.
#' @param sample_order string, how the samples should be ordered, either by 'pah_conc', which is
#' the sum of the EPA 16 priority compounds, or 'norm_pah_conc' which is the TOC-normalized PAH 16 concentration.
#' Sources are considered "unlikely" when the percent of source to sample is greater than the percent TOC
#' in the sample, given that PAHs are limited to the organic fraction. Ordering by TOC-normalized PAH concentration
#' gives a smoother look to the figure, but is less intuitive in terms of sample ordering.
#' @return If calc_type is "summary", each row represents a source, and source mean concentrations,
#' number of PAHs used, and references are reported alongside percent mass fractions calculated for all quartiles,
#' minimum, and maximum of all sample concentrations. If calc_type is 'by_sample', a data frame of
#' n samples x j sources is given, where each cell represents the mass fraction for that sample-source combination.
#' Sample IDs are given in a column, and columns are named by source ID.
#' @importFrom rlang sym
#' @importFrom tidyr gather
#' @import dplyr
#' @import ggplot2
#' @examples

calc_mass_fractions <- function(compound_info, sample_column, conc_column, compound_column,
                                conc_unit = "ppb", calc_type = 'summary', plot = FALSE,
                                sample_order = 'norm_pah_conc') {
  # make column names dplyr-ready
  quo_sample_column <- sym(sample_column)
  quo_conc_column <- sym(conc_column)
  quo_compound_column <- sym(compound_column)

  dat <- filter(compound_info, EPApriority16 %in% TRUE) %>%
    group_by(!!quo_sample_column) %>%
    summarize(sum_EPA16 = sum(!!quo_conc_column))

  toc <- filter(compound_info, !!quo_compound_column %in% 'TOC') %>%
    select(!!quo_sample_column, TOC = !!quo_conc_column)

  if (conc_unit == 'ppb') {
    dat <- mutate(dat, concentration = (sum_EPA16/1000))
  } else {
    dat <- mutate(dat, concentration = sum_EPA16)
  }

  dat <- left_join(dat, toc) %>%
    mutate(norm_concentration = concentration/(TOC/100))

  toc <- rename(toc, id = !!quo_sample_column)

  source_mean <- group_by(source_conc, source) %>%
    summarize(concentration = mean(conc_mgkg),
              PAHs_used = mean(PAHs_used),
              reference = paste0(reference, collapse = ", ")) %>%
    mutate(units = "mg/kg")

  # function to calculate mass fraction as a percent, round to two digits, and then
  # format to say ">100" for all values over 100.
  mass_fraction <- function(x) {
    temp <- round((x/source_mean$concentration)*100, 2)
    temp <- ifelse(temp <= 100, temp, ">100")
  }

  if (calc_type == 'summary') {
    samp_summary <- as.numeric(summary(dat$concentration))

    source_mean$mass_frac_min <- mass_fraction(samp_summary[[1]])
    source_mean$mass_frac_Q1 <- mass_fraction(samp_summary[[2]])
    source_mean$mass_frac_median <- mass_fraction(samp_summary[[3]])
    source_mean$mass_frac_mean <- mass_fraction(samp_summary[[4]])
    source_mean$mass_frac_Q3 <- mass_fraction(samp_summary[[5]])
    source_mean$mass_frac_max <- mass_fraction(samp_summary[[6]])

    source_mean$concentration <- round(source_mean$concentration, 1)
    source_mean$PAHs_used <- round(source_mean$PAHs_used, 1)

    out <- source_mean
  }

  if (calc_type == "by_sample") {
    # create empty data frame
    frac_by_site <- as.data.frame(matrix(ncol = 18, nrow = nrow(dat)))

    # name the columns by source names, create a column for sample IDs
    names(frac_by_site) <- source_mean$source

    source_concentrations <- source_mean$concentration
    sample_concentrations <- dat$concentration

    for (i in 1:nrow(frac_by_site)) {
      frac_by_site[i,] <- round((sample_concentrations[i]/source_concentrations)*100, 2)
    }

    if (plot == FALSE) {

      convert_100 <- function(x){
        dat <- ifelse(x <= 100, x, ">100")
        return(dat)
      }

      # convert vals > 100 to ">100"
      frac_by_site <- as.data.frame(apply(frac_by_site, 2, convert_100))

      # add sample ID column and put as first column
      frac_by_site <- mutate(frac_by_site, id = dat[[sample_column]]) %>%
        select(id, 2:ncol(frac_by_site))

      names(frac_by_site)[1] <- sample_column

      out <- frac_by_site

    } else {
      frac_for_plot <- mutate(frac_by_site, id = dat[[sample_column]]) %>%
        gather(key = "source", value = "mass_fraction", -id)

      frac_for_plot <- left_join(frac_for_plot, toc) %>%
        mutate(category = case_when(TOC < mass_fraction ~ 'impossible (> %TOC)',
                                    TOC >= mass_fraction & TOC*0.5 < mass_fraction ~ 'unlikely (> 0.5 x %TOC)',
                                    TOC *0.5 >= mass_fraction ~ 'possible (< 0.5 x %TOC)'))

      if (sample_order == 'pah_conc') {
        plot_sample_order <- arrange(dat, concentration)
        frac_for_plot$id <- factor(frac_for_plot$id, levels = plot_sample_order[[sample_column]])
      } else if (sample_order == 'norm_pah_conc') {
        plot_sample_order <- arrange(dat, norm_concentration)
        frac_for_plot$id <- factor(frac_for_plot$id, levels = plot_sample_order[[sample_column]])
      }


      source_order <- arrange(source_mean, concentration)
      frac_for_plot$source <- factor(frac_for_plot$source, levels = source_order[['source']])

      # drop any NA values where TOC is missing
      frac_for_plot <- filter(frac_for_plot, !is.na(category))

      # plot
      p <- ggplot(frac_for_plot, aes(y = source, x = id)) +
        geom_tile(aes(fill = category), color = 'white') +
        scale_fill_manual(name = 'Mass Fraction', values = c("possible (< 0.5 x %TOC)" = '#1a9850',
                                                             "unlikely (> 0.5 x %TOC)" = '#fdae61',
                                                             "impossible (> %TOC)" = '#d73027')) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = "", y = "") +
        geom_vline(xintercept = round(nrow(dat)/4, 0)*1:3, size = 1)

      out <- p
    }
  }

  return(out)
}
