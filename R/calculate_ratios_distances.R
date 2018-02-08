#' calc_ratio_distance
#'
#' Calculates the Euclidean distances between samples and sources from
#' the double ratio plots.
#'
#' @export
#' @param ratio_dat The output of calc_ratios, which includes the four ratios of interest for both sources
#' and samples.
#' @import dplyr
#' @return A list of three data frames, which includes 1) distances for each
#' double ratio plot, source, and sample combination, along with the ranking of each source by sample (1 = closest) for each
#' double ratio plot ("raw"), 2) distances summarized by source, which contains the mean
#' distance, rank, and the number of top ranks calculated for each source and double ratio plot combination ("source"),
#' 3) distances summarized by sample, which reports the top-ranked source for each sample and double
#' ratio combination ("sample"). Abbreviations for the ratios include: Anth_AnthPhen or AAP (ratio of anthracene
#' to anthracene plus phenanthrene), Flua_FluaPyr or FFP (ratio of fluoranthene to fluoranthene plus pyrene),
#' Baa_BaaCh or BBC (ratio of benzo(a)anthracene to benzo(a)anthracene plus chrysene), Indpy_IndpyBghip or IIB (ratio
#' of indeno(1,2,3-cd)pyrene to indeno(1,2,3-cd)pyrene plus benzo(g,h,i)perylene). Each double ratio comparison
#' is distinguished by the shorter of the abbreviations separated by an underscore (e.g., FFP_AAP is one double ratio
#' comparison).
#' @examples

calc_ratio_dist <- function(ratio_dat) {

  structure1 <- select(ratio_dat, sample_id, sample_type) %>%
    filter(sample_type == 'sample') %>%
    select(sample_id)

  structure2 <- select(ratio_dat, sample_id, sample_type) %>%
    filter(sample_type == 'source') %>%
    select(sample_id)

  structure <- expand.grid(structure1[[1]], structure2[[1]])

  ratio_dat_sub <- select(ratio_dat, -sample_type)

  comp <- rename(structure, sample = Var1, source = Var2) %>%
    left_join(ratio_dat_sub, by = c('sample' = 'sample_id')) %>%
    rename(Anth_AnthPhen_sam = Anth_AnthPhen,
           Flua_FluaPyr_sam = Flua_FluaPyr,
           Baa_BaaCh_sam = Baa_BaaCh,
           Indpy_IndpyBghip_sam = Indpy_IndpyBghip) %>%
    left_join(ratio_dat_sub, by = c('source' = 'sample_id')) %>%
    rename(Anth_AnthPhen_src = Anth_AnthPhen,
           Flua_FluaPyr_src = Flua_FluaPyr,
           Baa_BaaCh_src = Baa_BaaCh,
           Indpy_IndpyBghip_src = Indpy_IndpyBghip)

  comp <- comp %>%
    mutate(diff1 = sqrt((Anth_AnthPhen_sam - Anth_AnthPhen_src)^2 +
                          (Flua_FluaPyr_sam - Flua_FluaPyr_src)^2),
           diff2 = sqrt((Indpy_IndpyBghip_sam - Indpy_IndpyBghip_src)^2 +
                          (Flua_FluaPyr_sam - Flua_FluaPyr_src)^2),
           diff3 = sqrt((Indpy_IndpyBghip_sam - Indpy_IndpyBghip_src)^2 +
                          (Baa_BaaCh_sam - Baa_BaaCh_src)^2))

  ranks <- comp %>%
    group_by(sample) %>%
    mutate(rank1 = row_number(diff1),
           rank2 = row_number(diff2),
           rank3 = row_number(diff3))

  # find how many times each source is the top source by sample
  # then standardize that by the number of comparisons were made
  # between that source and samples (some sources do not have
  # all ratios, so sum of top ranks is hard to compare)

  by.sources.rank.counts <- ranks %>%
    group_by(source) %>%
    summarize(top_diff1 = ifelse(anyNA(rank1), NA, length(which(rank1 %in% 1))),
              top_diff2 = ifelse(anyNA(rank2), NA, length(which(rank2 %in% 1))),
              top_diff3 = ifelse(anyNA(rank3), NA, length(which(rank3 %in% 1))))

  by.sources.rank.counts$n_poss <- (3-rowSums(is.na(by.sources.rank.counts)))*length(unique(ranks$sample))
  by.sources.rank.counts$n_prop <- (rowSums(by.sources.rank.counts[,c('top_diff1', 'top_diff2', 'top_diff3')], na.rm = T)/by.sources.rank.counts$n_poss)*100
  by.sources.rank.counts <- filter(by.sources.rank.counts, !is.na(n_prop))


  # calculate the mean the distances and ranks for each source
  mean.sources <- ranks %>%
    select(source, diff1, diff2, diff3, rank1, rank2, rank3) %>%
    group_by(source) %>%
    summarize_all(mean) %>%
    ungroup() %>%
    mutate_if(is.numeric, funs(round(., digits = 2))) %>%
    rename(FFP_AAP_meandiff = diff1, FFP_IIB_meandiff = diff2, BBC_IIB_meandiff = diff3,
           FFP_AAP_meanrank = rank1, FFP_IIB_meanrank = rank2, BBC_IIB_meanrank = rank3)

  # merge back in with count of top source
  mean.sources <- left_join(mean.sources, by.sources.rank.counts, by = 'source') %>%
    rename(FFP_AAP_n_toprank = top_diff1, FFP_IIB_n_toprank = top_diff2, BBC_IIB_n_toprank = top_diff3)

  by.samples <- ranks %>%
    group_by(sample) %>%
    summarize(top_source1 = source[which.min(diff1)],
              top_source2 = source[which.min(diff2)],
              top_source3 = source[which.min(diff3)])

  top.sources <- by.sources$source[by.sources$n_prop > 5]
  top.sources <- na.omit(top_sources)

  return_dat <- list(ranks, mean.sources, top.sources)
  names(return_dat) <- c('raw', 'source', 'sample')

  return(return_dat)

}
