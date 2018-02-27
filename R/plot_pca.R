#' plot_pca
#'
#' Creates various plots from the PCA components as well as the calculated PCA Euclidean distances.
#'
#' @export
#' @param pca_dat The output from 'pah_pca' which includes the component data as well as the calculated Euclidean
#' distances data.
#' @param plot_type Which plot to create. Options include 'distance_boxplot' which summarizes the Euclidean
#' distances by source, 'pca_components' which shows where the sources and samples lie in all combinations
#' of the chosen components space.
#' @param sources a dataframe of source profiles. The default is to use the built-in `source_profiles` table,
#' but users can provide their own table. This is useful if the user has a source profile to add to the built-in table.
#' @return If plot_type = 'distance_boxplot', a single boxplot is returned with source IDs on the x axis and
#' Euclidean distances on the y axis. If plot_type = 'pca_components', a scatterplot of
#' all possible combinations of chosen components are included in a panel matrix, with samples as black dots
#' and sources as red dots that are labeled with the source abbreviation. To see the full names of each source
#' abbrevation, see table source_ratios.
#' @import ggplot2
#' @import dplyr
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @examples

plot_pca <- function(pca_dat, plot_type = "distance_boxplot") {
  if (plot_type == "distance_boxplot") {
    distance <- pca_dat$pca_distance

    median_diff <- distance %>%
      group_by(source) %>%
      summarize(median = median(euc_dist)) %>%
      arrange(median)

    distance$source <- factor(distance$source, levels = median_diff$source)

    p <- ggplot(distance, aes(x = source, y = euc_dist)) +
      geom_boxplot() +
      theme_bw() +
      labs(y = "Euclidean distance\n(zero = identical to sample)", x = "") +
      theme(panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    return(p)

  } else if (plot_type == "pca_components") {
    comp_dat <- pca_dat$pca_dat

    plot_num <- 1
    plot_list <- list()
    my_labels <- comp_dat$sample_id
    my_labels[comp_dat$type == 'sample'] <- ""

    n_pca <- ncol(comp_dat) - 2

    for(xcol in 1:(n_pca-1)){
      for(ycol in 2:n_pca){
        if (xcol >= ycol) {next}

        temp_dat <- data.frame(x = comp_dat[,xcol],
                               y = comp_dat[,ycol],
                               type = comp_dat$type)

        p <- ggplot(data = temp_dat, aes(x = x, y = y)) +
          geom_point(aes_string(color = 'type'), alpha = 0.5, show.legend = F) +
          geom_text_repel(data = temp_dat,
                          aes(x = x, y = y, label = my_labels)) +
          scale_color_manual(values = c('black', 'red')) +
          labs(x = paste('Component', xcol), y = paste('Component', ycol)) +
          theme_classic()

        plot_list[[plot_num]] <- p
        plot_num <- plot_num + 1
      }
    }

    if (n_pca < 4) {
      p_final <- cowplot::plot_grid(plotlist = plot_list, nrow = 1)
    } else if(n_pca == 4) {
      p_final <- cowplot::plot_grid(plotlist = plot_list, nrow = 3)
    } else if (n_pca > 4) {
      p_final <- cowplot::plot_grid(plotlist = plot_list, nrow = 4)
    }

    return(p_final)

  } else {
    warning('That plot_type does not exist.')
  }
}
