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
#' @param sample_column column name of unique sample identifier
#' @param conc_column column name of reported concentrations
#' @param source_profs a dataframe of source profiles. The default is to use the built-in `source_profiles` table,
#' but users can provide their own table. This is useful if the user has a source profile to add to the built-in table.
#' @param include_creosote Logical, whether to include the source profiles for creosote (n = 2). The source profiles
#' for creosote only include 11 compounds, and the missing 12th compound will be dropped in all sample-source
#' comparisons if include_creosote = T. It is recommended that users first conduct analysis with creosote to determine
#' if creosote is an important source. If not, continue analysis without the creosote sources using the 12-compound profiles.
#' @return Returns two data frames. The first (profiles) is a long dataframe where observations are repeated for
#' each sample/compound/source combination, and reports the proportional concentration of that unique compound/sample combination,
#' and chi-squared distance between the source and sample. Additionally, the function adds
#' all columns in `source_profiles`.  The second dataframe is a long data frame where observations are repeated for
#' each unique sample/source combination, and the chi2 value is summed across all compounds to create a distance
#' metric between each sample and source profile.
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom rlang sym
#' @examples

pah_profiler <- function(pah_dat, compound_column = 'casrn', sample_column,
                         conc_column, source_profs = source_profiles, include_creosote = T) {
  # make column names dplyr-ready
  quo_compound_column <- sym(compound_column)
  quo_conc_column <- sym(conc_column)
  quo_sample_column <- sym(sample_column)

  # pull out all 12 compounds
  # filter to 11 compounds if using creosote

  if (include_creosote == T) {

    profile_compounds <- filter(source_profs, Abbreviation != 'BeP') %>%
      select(!!quo_compound_column)

    # because other proportions were calculated based on 12 compounds,
    # need to adjust for the proportion that is made up of BeP (the
    # compound being dropped)
    BeP_prop <- filter(source_profs, Abbreviation == 'BeP') %>%
      select(PPLT:CTD7)

    BeP_prop <- as.numeric(BeP_prop[1,])

    prop_fixed <- filter(source_profs, Abbreviation != 'BeP') %>%
      select(PPLT:CTD7)

    prop_fixed <- data.frame(t(round(t(prop_fixed)/(1-BeP_prop), 3)))

    prop_rest <- filter(source_profs, Abbreviation != 'BeP') %>%
      select(-(PPLT:CTD7)) %>%
      bind_cols(prop_fixed) %>%
      select(Compound, Abbreviation, pcode, PPLT:CTD7, CRE2, CRE4, casrn, molwt)

    source_profs <- prop_rest
  } else {
    profile_compounds <- select(source_profs, !!quo_compound_column)
    source_profs <- select(source_profs, -CRE2, -CRE4)
  }

  # filter user samples to include only those in the source profiles
  # only include necessary columns
  samp.prof <- filter(pah_dat, (!!quo_compound_column) %in% profile_compounds[[1]]) %>%
    select(!!quo_sample_column, !!quo_compound_column, !!quo_conc_column)


  # group by sample - calculate sum total and whether or not there are any samples = 0
  samp.prof.bysample <- group_by(samp.prof, !!quo_sample_column) %>%
    summarize(total_pah = sum(!!quo_conc_column))

  samp.prof <- left_join(samp.prof, samp.prof.bysample) %>%
    mutate(prop_conc = (!!quo_conc_column)/total_pah)

  # merge in source compound info
  all.profs <- full_join(samp.prof, source_profs, by = compound_column) %>%
    select(-!!quo_conc_column) %>%
    gather(key = source, value = source_prop_conc, -!!quo_sample_column, -casrn, -total_pah, -prop_conc, -Compound, -Abbreviation, -pcode, -molwt)

  # calculate the chi squared difference
  all.profs <- mutate(all.profs, chi2 = (abs(prop_conc - source_prop_conc)^2)/((prop_conc + source_prop_conc)/2))

  # sum over the sources by sample
  all.summary <- group_by(all.profs, !!quo_sample_column, source) %>%
    summarize(sum_chi2 = sum(chi2))

  out <- list(all.profs, all.summary)
  names(out) <- c('profiles', 'sum_chi2')
  return(out)
}
