#' calc_ratios
#'
#' Calculates diagnostic ratios and compares them to known source ratios. Double ratio plots are then
#' used to calculate the distance between source and sample. All double ratio Euclidean distances between source
#' and sample are averaged by source to calculate a mean distance between sample and source.
#'
#' @export
#' @param pah_dat dataframe with PAH concentrations, where there is one column of compounds, and each sample
#' is contained in a column. If individual samples are to be averaged before computation, the user should supply
#' the average values. The dataframe must contained a column named 'casrn' which contains the casrn numbers
#' for each compound. Samples below detection limit should be coded as "0" in the concentration column. The function
#' removes any samples that report a 0 for any of the compounds required to calculate ratios.
#' @param sample_column name of the column that contains your unique sample identifier
#' @param conc_column name of the column that contains the results (concentrations)
#' @return A data frame that contains a column `sample_id` which contains the original sample identifiers,
#' as well as the unique identifiers of each source contained in `source_profiles`. The source and sample
#' ratios are defined in the `sample_type` column which is defined as either `source` or `sample`. Additionally,
#' four columns of ratios are included: Anth_AnthPhen (ratio of anthracene to anthracene plus phenanthrene),
#' Flua_FluaPyr (ratio of fluoranthene to fluoranthene plus pyrene), Baa_BaaCh (ratio of benzo(a)anthracene
#' to benzo(a)anthracene plus chrysene), Indpy_IndpyBghip (ratio of indeno(1,2,3-cd)pyrene to
#' indeno(1,2,3-cd)pyrene plus benzo(g,h,i)perylene).
#' @import dplyr
#' @importFrom tidyr spread
#' @importFrom rlang sym
#' @examples

calc_ratios <- function(pah_dat, sample_column = 'sample_id', conc_column) {

  # make column names dplyr-ready
  quo_sample_column <- sym(sample_column)
  quo_conc_column <- sym(conc_column)

  # remove 0 (BDL) values for compounds that we will use
  # this includes pcodes c(63208, 63180, 64111, 63610, 64118, 64113)
  # or casrn
  casrn.keep <- c('206-44-0', '120-12-7', '205-99-2', '56-55-3', '191-24-2')
  pah_temp <- filter(pah_dat, casrn %in% casrn.keep)

  # get a list of sites to remove - that is, those that have conc = 0
  sites.drop <- unique(pah_temp[, sample_column][pah_temp[, conc_column] == 0])
  pah_temp <- filter(pah_temp, !((!!quo_sample_column) %in% sites.drop)) %>%
    select(!!quo_sample_column, casrn, !!quo_conc_column)

  pah_dat_wide <- pah_temp %>%
    spread(key = casrn, value = !!quo_conc_column)

  pah_dat_wide <- pah_dat_wide %>%
    mutate(Anth_AnthPhen = `120-12-7`/(`120-12-7` + `85-01-8`),
           Flua_FluaPyr = `206-44-0`/ (`206-44-0` + `129-00-0`),
           Indpy_IndpyBghip = `193-39-5`/(`193-39-5` + `191-24-2`),
           Baa_BaaCh = `56-55-3`/(`56-55-3` + `218-01-9`)) %>%
    select(!!quo_sample_column, Anth_AnthPhen:Indpy_IndpyBghip) %>%
    mutate(sample_type = 'sample') %>%
    rename(sample_id = !!quo_sample_column)

  # Combine sample ratios with source ratios
  ratios <- rename(source_ratios, sample_id = abbrev) %>%
    select(sample_id, Anth_AnthPhen:Indpy_IndpyBghip) %>%
    mutate(sample_type = 'source') %>%
    bind_rows(pah_dat_wide)

  return(ratios)
}
