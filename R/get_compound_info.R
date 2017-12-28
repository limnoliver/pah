#' get_compound_info
#'
#' Merges user PAH data with compound-specific metadata, including whether the compound is: one of
#' the 16 EPA priority compounds, low or high molecular weight, a parent or alkylated PAH, and whether
#' the compound is in the source profiles.
#'
#' @export
#' @param pah_dat data frame of user PAH data, either retrieved using the get_pah data, or
#' the user's own PAH concentration data. Data frame must include a PAH identifier column, which can be
#' either the USGS parameter code, casrn, or compound name.
#' @param merge_type whether the merge is being done by USGS parameter code ('pcode'), Chemical Abstracts
#' Service Registry Number ('casrn'), or compound name ('name').
#' @param merge_col column name that contains the PAH identifier by which to merge
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @examples
get_compound_info <- function(pah_dat, ){

}
