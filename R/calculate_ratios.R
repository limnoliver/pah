#' calc_ratios
#'
#' Calculates diagnostic ratios and compares them to known source ratios. Double ratio plots are then
#' used to calculate the distance between source and sample. All double ratio Euclidean distances between source
#' and sample are averaged by source to calculate a mean distance between sample and source.
#'
#' @export
#' @param pah_dat dataframe with PAH concentrations, where there is one column of compounds, and each sample
#' is contained in a column. If individual samples are to be averaged before computation, the user should supply
#' the average values.
#' @return
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom rlang sym
#' @examples

pah_profiler <- function(pah_dat, compound_column = 'casrn', sample_column = 'sample_id',
                         conc_column = 'RESULT', sources = source_profiles) {


  # make column names dplyr-ready
  quo_compound_column <- sym(compound_column)

  # remove 0 (BDL) values for compounds that we will use
  # this includes pcodes 63208, 63227, 63180, 63224, 64111, 63610, 64118, 64114, 64115, 64113

  # create bins based on PAH16 concentration (do this in workflow?)
  # remove source COA2 - it creates too much dead space -- look at graphs before doing this
  # compute ratios -- in Austin's version, each compound has its own column
  ratios$Fl_Py <- ratios$p63208 / ratios$p63227  # flouranthene / pyrene
  ratios$BbF_BkF <- ratios$p64111 / ratios$p64114  # benzo[b]fluoranthene / benzo[k]fluoranthene
  ratios$An_178 <- ratios$p63180 / (ratios$p63180 + ratios$p64111) # anthracene / anthracene + phenanthrene (MW178)
  ratios$Fl_FlPy <- ratios$p63208 / (ratios$p63208 + ratios$p63227) # flouranthene / flouranthene + pyrene (MW202)
  ratios$BaA_BaACh <- ratios$p63610 / (ratios$p63610 + ratios$p64115) # benz[a]anthracene / benz[a]anthracene + chrysene (MW228)
  ratios$IP_IPBghiP <- ratios$p64118 / (ratios$p64118 + ratios$p64113)  # indeno[1,2,3-cd]pyrene / indeno[1,2,3-cd]pyrene+ benzo[g,h,i]perylene (MW276)



}
