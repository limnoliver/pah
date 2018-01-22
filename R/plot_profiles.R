#' plot_profiles
#'
#' Creates common figures from the calculates sample and source profiles. Possible plots include
#' a boxplot of the sum of the distance metric across all sources, proportional concentration
#' profile figures comparing an individual sample to a source.
#'
#' @export
#' @param profile_dat The output of `pah_profiler` which includes the proportional concentration
#' profiles from the samples and sources, as well the sum of the distance metric for each
#' sample/source combination.
#' @param plot_type The desired plot type. This can either be a boxplot ('boxplot') of all distance measures
#' by source, or various profile ('profile') plots.
#' @param source If `plot_type` = `profile`, include a vector of sources that you want to include
#' in a profile plot to compare to samples.
#' @return Returns two data frames. The first is a long dataframe where observations are repeated for each sample/compound combination. Adds to the input
#' dataframe all columns in `source_profiles`, as well as the chi-squared difference between source and sample, and the proportional
#' concentrations. The second dataframe is a long data frame where observations are repeated for
#' each sample/source combination, and the chi2 value is summed across all compounds to create a distance
#' metric between each sample and source profile.
#' @import dplyr
#' @importFrom rlang sym
#' @examples

plot_profiles <- function(profile_dat, compound_column = 'casrn', sample_column = 'sample_id',
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
