#' plot_profiles
#'
#' Creates common figures from the calculated distance between sample and source profiles.
#' Possible plots include a boxplot of the sum of the distance metric across all sources and
#' proportional concentration profile figures comparing an individual sample to a source.
#'
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

}
