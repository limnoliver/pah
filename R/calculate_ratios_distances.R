#' calc_ratio_distance
#'
#' Calculates the Euclidean distances between samples and sources from
#' the double ratio plots.
#'
#' @export
#' @param ratio_dat The output of calc_ratios, which includes the four ratios of interest for both sources
#' and samples.
#' @param summary_type Whether the distances should be exported as the raw distances and ranks for each
#' double ratio plot, source, and sample combination ("raw"),
#' or summarized by "sample" or by "source". If the data are summarized by source,
#' the generated table contains a source in each row, and the mean distance and rank between that source
#' and all samples. This will help discern which sources are closest to the samples overall (or on average). If the distances are summarized
#' If the data are summarized by sample, the generated table contains a sample in each row, and the top ranked
#' source for each double ratio comparison. This will help discern the likely source for each sample.
#' @import dplyr
#' @return A data frame that contains a column `sample_id` which contains the original sample identifiers,
#' as well as the unique identifiers of each source contained in `source_profiles`. The source and sample
#' ratios are defined in the `sample_type` column which is defined as either `source` or `sample`. Additionally,
#' four columns of ratios are included: Anth_AnthPhen (ratio of anthracene to anthracene plus phenanthrene),
#' Flua_FluaPyr (ratio of fluoranthene to fluoranthene plus pyrene), Baa_BaaCh (ratio of benzo(a)anthracene
#' to benzo(a)anthracene plus chrysene), Indpy_IndpyBghip (ratio of indeno(1,2,3-cd)pyrene to
#' indeno(1,2,3-cd)pyrene plus benzo(g,h,i)perylene).
#' @examples

calc_ratio_dist <- function(ratio_dat, summary_type = 'sample') {

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
           rank3 = row_number(diff3)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(sumranks = sum(rank1, rank2, rank3, na.rm = T))

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

  if (summary_type == "sources") {

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

  return(mean.sources)
  } else if (summary_type == "samples") {

    by.samples <- ranks %>%
      group_by(sample) %>%
      summarize(top_source1 = source[which.min(diff1)],
                top_source2 = source[which.min(diff2)],
                top_source3 = source[which.min(diff3)])

    top_sources <- by.sources$source[by.sources$n_prop > 5]
    top_sources <- na.omit(top_sources)

    return(top_sources)
  } else if (summary_type == "raw") {
    return(ranks)
  } else {
    warning("summary_type must be set to 'raw', 'sources', or 'samples'.")
  }
}
