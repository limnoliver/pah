#' pah_profiler
#'
#' Creates profiles of PAH compounds in the sample and compares them to source profiles. The
#' difference between the sample and the source is calculated using a chi-squared statistic. To see
#' which 12 compounds are used from the source profiles, see the table `source_profiles`.
#'
#' @export
#' @param pah_dat dataframe with PAH concentrations, where there is one column of compounds, and each sample
#' is contained in a column. If individual samples are to be averaged before computation, the user should supply
#' the average values.
#' @param compound_column column name which will be used to merge with source profiles. This can be a USGS pcode ('pcode'),
#' CAS registration number ('casrn'), or compound name ('Compound').
#' @param sources a dataframe of source profiles. The default is to use the built-in `source_profiles` table,
#' but users can provide their own table. This is useful if the user has a source profile to add to the built-in table.
#' @return Returns two data frames. The first is a long dataframe where observations are repeated for each sample/compound combination. Adds to the input
#' dataframe all columns in `source_profiles`, as well as the chi-squared difference between source and sample, and the proportional
#' concentrations. The second dataframe is a long data frame where observations are repeated for
#' each sample/source combination, and the chi2 value is summed across all compounds to create a distance
#' metric between each sample and source profile.
#' @import dplyr
#' @importFrom rlang sym
#' @examples

pah_profiler <- function(pah_dat, compound_column = 'casrn', sample_column = 'sample_id',
                         conc_column = 'RESULT', sources = source_profiles) {
  # make column names dplyr-ready
  quo_compound_column <- sym(compound_column)
  quo_conc_column <- sym(conc_column)
  quo_sample_column <- sym(sample_column)

  # pull out all 12 compounds
  profile_compounds <- select(sources, !!quo_compound_column)

  # filter user samples to include only those in the source profiles
  samp.prof <- filter(pah_dat, (!!quo_compound_column) %in% profile_compounds[[1]])

  # group by sample - calculate sum total and whether or not there are any samples = 0
  samp.prof.bysample <- group_by(samp.prof, !!quo_sample_column) %>%
    summarize(total_pah = sum(!!quo_conc_column))

  samp.prof <- left_join(samp.prof, samp.prof.bysample) %>%
    mutate(prop_conc = (!!quo_conc_column)/total_pah)

  # merge in source compound info
  all.profs <- full_join(samp.prof, sources, by = compound_column) %>%
    select(-RESULT) %>%
    gather(key = source, value = source_prop_conc, -(1:8), -Compound, -Abbreviation, -pcode, -molwt)

  # calculate the chi squared difference
  all.profs <- mutate(all.profs, chi2 = (abs(prop_conc - source_prop_conc)^2)/((prop_conc + source_prop_conc)/2))

  # sum over the sources by sample
  all.summary <- group_by(all.profs, sample_id, source) %>%
    summarize(sum_chi2 = sum(chi2))

  out <- list(all.profs, all.summary)
  names(out) <- c('profiles', 'sum_chi2')
  return(out)
}
