#' plot_ratio_dist
#'
#' This function uses  the output of `calc_ratio_dist` to create summary visualizations, either summarized
#' by the samples or sources.
#'
#' @export
#' @param ratio_dist_dat  The output from `calc_ratio_dist`, a list of three dataframes.
#' @param plot_type Options for different summary plots. Options include "sample", which color codes the
#' sources and shows top source for each sample, "source-mean" which creates a 3-paneled plot of
#' the calculated mean Euclidean distances between samples and sources for each double ratio plot and color-
#' codes each bar by the mean rank of that source across all samples, and "source-top" which creates a
#' bar chart of the proportion of source-sample comparisons where each source was the closest source
#' to a sample.
#' @param percent_cutoff Numeric, from 0 to 100. A cutoff for plot_type = "sample", below which sources will
#' be pooled into category "other". The percent refers to the percent of source-sample comparisons where the
#' source was ranked as the top source by Euclidean distance in a double ratio plot space. This allows the
#' user to bin "unimportant" sources into "other" and reduces the color complexity of the plot.
#' @param sample_order A vector of sample IDs (that are equivalent to the IDs stored in ratio_dist_dat) in
#' the order they are to appear when plot_type = "sample". For example, samples could be sorted from low to
#' high PAH concentration to look for patterns in likely sources across concentration gradients.
#' @return One of three plots. "sample" returns a tiled plot showing the top source by each sample (y-axis)
#' in all three double ratio comparisons (x-axis). "source-mean" creates a bar chart with source on the x-axis
#' and mean Euclidean distance across all samples on the y-axis. The plot is organized into three panels, with
#' each panel representing a different double ratio comparison. "source-top" creates a bar chart with source
#' on the x-axis and the proportion of comparisons each source was considered the top or closest source for a
#' given sample. A proportion rather than a raw count is given because not all sources were in each double
#' ratio comparison. The bars are colored by the number of source to sample distances that were calculated for
#' that source.
#' @importFrom dplyr filter
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra grid.arrange
#' @examples

plot_ratios <- function(ratio_dist_dat, plot_type = "source-mean", percent_cutoff = 5, sample_order = NA) {
  if (plot_type == 'source-mean') {

    # convert data from wide to long
    mean.sources <- ratio_dist_dat$source

    long_diff <- select(mean.sources, source:BBC_IIB_meandiff) %>%
      gather(double_ratio, double_ratio_meandiff, -source) %>%
      mutate(double_ratio = gsub('_meandiff', '', double_ratio))

    long_rank <- select(mean.sources, source, FFP_AAP_meanrank:BBC_IIB_meanrank) %>%
      gather(double_ratio, double_ratio_meanrank, -source) %>%
      mutate(double_ratio = gsub('_meanrank', '', double_ratio))

    mean.sources.long <- left_join(long_diff, long_rank, by = c('source', 'double_ratio'))

    # order the sources by the first double ratio comparison
    order.sources <- filter(mean.sources.long, double_ratio == "FFP_AAP") %>%
      arrange(double_ratio_meandiff)

    # modify double ratio comparisons to include full ratio description
    mean.sources.long$source <- as.factor(mean.sources.long$source)
    mean.sources.long$source <- factor(mean.sources.long$source, levels = order.sources[[1]])
    mean.sources.long$double_ratio <- as.factor(mean.sources.long$double_ratio)
    levels(mean.sources.long$double_ratio) <- c('BaA/(BaA+Ch) : IndPy/(IndPy+BghiP)',
                                                'FluA/(FluA+Pyr) : Anth/(Anth+Phen)',
                                                'FluA/(FluA+Pyr) : IndPy/(IndPy+BghiP)')
    mean.sources.long$double_ratio <- factor(mean.sources.long$double_ratio,
                                             levels(mean.sources.long$double_ratio)[c(2,3,1)])

    # create barchart that plots mean distance by source, color codes by rank
    p <- ggplot(mean.sources.long, aes(x = source, y = double_ratio_meandiff)) +
      geom_bar(stat = 'identity', position = 'dodge', aes(fill = double_ratio_meanrank), color = 'black', size = 0.4) +
      scale_fill_gradient2(low = '#67001f', mid = '#f7f7f7', high = '#053061',
                           midpoint = median(mean.sources.long$double_ratio_meanrank, na.rm = T),
                           name = "Mean Rank") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), strip.background = element_blank()) +
      facet_wrap(~double_ratio, nrow = 3) +
      labs(x = 'Source', y = 'Mean distance between source and samples')

    return(p)

  } else if (plot_type == "source-top") {
    mean.sources <- ratio_dist_dat$source

    p <- ggplot(mean.sources, aes(x = reorder(source, n_prop), y = n_prop)) +
      geom_bar(stat = 'identity', aes(fill = factor(n_poss))) +
      scale_fill_brewer(palette = 'Dark2', type = 'qual', name = 'Number of Comparisons') +
      labs(x = 'Sources', y = 'Percent Times Top Source in\nDouble Ratio Comparisons') +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = c(0,1), legend.justification = c(0,1),
            legend.box.background = element_rect(color = 'black'),
            legend.box.margin = margin(rep(1,4))) +
      guides(fill = guide_legend(nrow = 1))

    return(p)

  } else if (plot_type == "sample") {

    # find any sources that have were the top source
    # at least 5% of the time, classify all other sources as "other"
    mean.sources <- ratio_dist_dat$source
    top_sources <- mean.sources$source[mean.sources$percent_top >= percent_cutoff]
    top_sources <- na.omit(top_sources)

    ##################################
    # plot top sources by sample
    # use this code only for plotting
    by.samples <- ratio_dist_dat$sample
    by.samples.long <- by.samples %>%
      gather(key = comparison, value = top_var, -sample) %>%
      mutate(top_var_cat = ifelse(top_var %in% top_sources, top_var, 'other'))

    if (!is.na(sample_order)) {
      if (length(sample_order) != nrow(by.samples)) {
        warning("sample_order length does not match the number of unique sample IDs")
      }
      if (!all(by.samples$sample %in% sample_order)) {
        warning("some sample IDs missing from sample_order")
      }

      by.samples.long$sample <- as.factor(by.samples.long$sample)
      by.samples.long$sample <- factor(by.samples.long$sample, levels = sample_order)
    }

    p <- ggplot(by.samples.long, aes(x = comparison, y = sample)) +
      geom_tile(aes(fill = top_var_cat), color = 'white') +
      scale_fill_manual(values = c('#8c510a', '#bf812d', '#c51b7d',
                                   '#01665e', '#35978f', 'gray', '#b2182b'),
                        name = "Source") +
      labs(x = "", y = "")
  }
}
