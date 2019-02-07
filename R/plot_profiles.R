#' plot_profiles
#'
#' Creates common figures from the calculated distance between sample and source profiles.
#' Possible plots include a boxplot of the sum of the distance metric across all sources and
#' proportional concentration profile figures comparing an individual sample to a source.
#'
#' @export
#' @param profile_dat The output of `pah_profiler` which includes the proportional concentration
#' profiles from the samples and sources, as well the sum of the distance metric for each
#' sample/source combination. These are combined in a list.
#' @param plot_type The desired plot type. This can either be a boxplot ('boxplot') of all distance measures
#' by source, or various profile ('profile') plots.
#' @param sources_plot If `plot_type` = `profile`, include a vector of abbreviated sources that you want to include
#' in a profile plot to compare to samples. Set equal to "all" to include paneled plots of all source profiles.
#' @param sample_column Column that contains unique sample IDs.
#' @param samples_plot If `plot_type` = `profile`, a unique sample ID to use to plot against source profiles.
#' Can either be a single unique ID or "all" to indicate taking the mean and standard deviation of all samples.
#' @param source_abbreviation logical, whether source abbreviations should be used when plot_type = 'distance_boxplot'.
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang sym
#' @examples

plot_profiles <- function(profile_dat, plot_type = 'boxplot', sources_plot = NA, samples_plot = 'all',
                          sample_column = 'sample_id', source_abbreviation = FALSE) {
  quo_sample_column <- sym(sample_column)
  if (plot_type == 'boxplot') {

    if (source_abbreviation == FALSE) {
      distance <- left_join(profile_dat,
                            select(pah::sources, source_abbrev, source_short_no_ref),
                            by = c('source' = 'source_abbrev')) %>%
        select(-source) %>%
        rename(source = source_short_no_ref)
    }

      sum_chi2 <- profile_dat[[2]]

      order.vals <- sum_chi2 %>%
        group_by(source) %>%
        summarize(med = median(sum_chi2)) %>%
        arrange(med)

      sum_chi2$source <- factor(sum_chi2$source, levels = order.vals$source)

    # create boxplot
    p <- ggplot(sum_chi2, aes(x = source, y = sum_chi2)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 0))) +
      labs(x = '', y = 'Sum Chi2')

  } else if (plot_type == 'profile') {
    pro_dat <- profile_dat[[1]]

    if (samples_plot == 'all') {
      # calculate sample means and sds if samples == 'all'
      sample_pro_dat <- group_by(pro_dat, Compound) %>%
        summarize(mean_prop_conc = mean(prop_conc),
                  sd_prop_conc = sd(prop_conc)) %>%
        rename(sample_prop_conc = mean_prop_conc)
    } else {
      sample_pro_dat <- filter(pro_dat, (!!quo_sample_column) %in% samples_plot) %>%
        rename(sample_prop_conc = prop_conc)
    }
    sources_dat <- select(pro_dat, Compound, source, source_prop_conc) %>%
      rename(prop_conc = source_prop_conc, profile = source) %>%
      distinct()

    profiles_all <- left_join(sources_dat, sample_pro_dat) %>%
      left_join(distinct(pro_dat[,c('Compound', 'molwt')]))

    if (!('all' %in% sources_plot)){
      profiles_all <- filter(profiles_all, profile %in% sources_plot)
    }
    if (samples_plot == 'all') {
      p <- ggplot(profiles_all, aes(x = reorder(Compound, molwt))) +
        geom_point(alpha = 0.7, aes(y = prop_conc)) +
        geom_line(aes(y = prop_conc, group = profile)) +
        geom_point(aes(x = reorder(Compound, molwt), y = sample_prop_conc), color = '#e41a1c', alpha = 0.8) +
        geom_errorbar(aes(ymin = sample_prop_conc-sd_prop_conc, ymax = sample_prop_conc+sd_prop_conc),
                      color = '#e41a1c', alpha = 0.7) +
        geom_line(aes(y = sample_prop_conc, group = 1), color = '#e41a1c') +
        facet_wrap(~profile, scales = 'free_y') +
        theme_bw() +
        labs(x = "", y = "Proportional Concentration") +
        theme(axis.text.x = element_text(size = rel(1.3), angle = 45, hjust = 1, margin = margin(t = 0)),
              axis.text.y = element_text(size = rel(1.3)),
              axis.title.y = element_text(size = rel(1.4)),
              strip.text = element_text(size = rel(1.4)),
              panel.grid.minor.y = element_blank())
    } else {
      p <- ggplot(profiles_all, aes(x = reorder(Compound, molwt))) +
        geom_point(alpha = 0.7, aes(y = prop_conc)) +
        geom_line(aes(y = prop_conc, group = profile)) +
        geom_point(aes(x = reorder(Compound, molwt), y = sample_prop_conc),
                   color = '#e41a1c', alpha = 0.8) +
        geom_line(aes(y = sample_prop_conc, group = 1), color = '#e41a1c') +
        facet_wrap(~profile, scales = 'free_y') +
        theme_bw() +
        labs(x = "", y = "Proportional Concentration") +
        theme(axis.text.x = element_text(size = rel(1.3), angle = 45, margin = margin(t = 0), hjust = 1),
              axis.text.y = element_text(size = rel(1.3)),
              axis.title.y = element_text(size = rel(1.4)),
              strip.text = element_text(size = rel(1.4)),
              panel.grid.minor.y = element_blank())
    }


  } else {
    warning('plot_type does not equal "boxplot" or "profile"')
  }
  return(p)
}
