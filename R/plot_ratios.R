#' plot_ratios
#'
#' This function uses takes the output of `calc_ratios` and creates double-ratio plots. This allows
#' comparison between samples and source ratios, and calculate which source ratio is closest (in Euclidean
#' distance) to the sample.
#'
#' @export
#' @param ratio_dat A dataframe which is the output from `calc_ratios`
#' @param source_drop A vector of source abbreviations to drop from plotting. This allows the user
#' to drop source outliers, or sources that are not near any samples.
#' @return A 1x3 plot of relevant diagnostic double ratio plots. Note that this function defaults to drop
#' the source "COA2" because it is an outlier relative to other sources. Use 'source_drop = NA' to include
#' all sources.
#' @importFrom dplyr filter
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom gridExtra grid.arrange
#' @examples

plot_ratios <- function(ratio_dat, source_drop = "COA2") {
  # drop outliers specificed by user
  ratio_dat <- filter(ratio_dat, !(sample_id %in% source_drop))

  p1 <- ggplot(ratio_dat, aes(y = Anth_AnthPhen, x = Flua_FluaPyr)) +
    geom_point(aes(color = sample_type), alpha = 0.6, shape = 16, size = 2) +
    scale_color_manual(values = c('gray', 'red'), name = "Sample Type") +
    #scale_fill_manual(values=c( "purple", "pink", "red"),   # for coloring pts according to conc
    #                     breaks=c("high","medium","low"),
    #                     labels=c(">100mg/kg", "10.1-100mg/kg","< 10mg/kg")) +
    geom_text_repel(data=subset(ratio_dat, sample_type == 'source'),
                    aes(y = Anth_AnthPhen, x = Flua_FluaPyr, label=sample_id)) +
    theme_bw()+
    theme(legend.position=c(0,1), legend.justification=c(0,1),
          legend.title=element_text(size=rel(.9)),
          legend.text=element_text(size=rel(.8)),
          legend.box.background = element_rect(colour = "black"),
          axis.text.x = element_text(size=rel(1)),
          axis.text.y = element_text(size=rel(1)),
          axis.title.y = element_text(size=rel(1.1)),
          axis.title.x = element_text(size=rel(1.1)),
          panel.grid = element_blank()) +
    labs(x = "FluA / FluA + Pyr", y = "Anth / Anth + Phen")

  p2 <- ggplot(ratio_dat, aes(y = Indpy_IndpyBghip, x = Flua_FluaPyr)) +
    geom_point(aes(color = sample_type), alpha = 0.6, shape = 16, size = 2) +
    scale_color_manual(values = c('gray', 'red'), guide = F) +
    #scale_fill_manual(values=c( "purple", "pink", "red"),   # for coloring pts according to conc
    #                     breaks=c("high","medium","low"),
    #                     labels=c(">100mg/kg", "10.1-100mg/kg","< 10mg/kg")) +
    geom_text_repel(data=subset(ratio_dat, sample_type == 'source'),
                    aes(y = Indpy_IndpyBghip, x = Flua_FluaPyr, label=sample_id)) +
    theme_bw()+
    theme(axis.text.x = element_text(size=rel(1)),
          axis.text.y = element_text(size=rel(1)),
          axis.title.y = element_text(size=rel(1.1)),
          axis.title.x = element_text(size=rel(1.1)),
          panel.grid = element_blank()) +
    labs(x = "FluA / FluA + Pyr", y = "IndPy / IndPy + BghiP")

  p3 <- ggplot(ratio_dat, aes(y = Indpy_IndpyBghip, x = Baa_BaaCh)) +
    geom_point(aes(color = sample_type), alpha = 0.6, shape = 16, size = 2) +
    scale_color_manual(values = c('gray', 'red'), name = "Sample Type", guide = F) +
    #scale_fill_manual(values=c( "purple", "pink", "red"),   # for coloring pts according to conc
    #                     breaks=c("high","medium","low"),
    #                     labels=c(">100mg/kg", "10.1-100mg/kg","< 10mg/kg")) +
    geom_text_repel(data=subset(ratio_dat, sample_type == 'source'),
                    aes(y = Indpy_IndpyBghip, x = Baa_BaaCh, label=sample_id)) +
    theme_bw()+
    theme(axis.text.x = element_text(size=rel(1)),
          axis.text.y = element_text(size=rel(1)),
          axis.title.y = element_text(size=rel(1.1)),
          axis.title.x = element_text(size=rel(1.1)),
          panel.grid = element_blank()) +
    labs(x = "BaA / BaA + Ch", y = "IndPy / IndPy + BghiP")

  p4 <- grid.arrange(p1, p2, p3, nrow = 1)

  return(p4)
}
