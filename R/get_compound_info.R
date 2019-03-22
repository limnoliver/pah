#' get_compound_info
#'
#' Merges user PAH data with compound-specific metadata, including whether the compound is: one of
#' the 16 EPA priority compounds, low or high molecular weight, a parent or alkylated PAH, and whether
#' the compound is in the source profiles.
#'
#' @param pah_dat data frame of user PAH data, either retrieved using the get_pah data, or
#' the user's own PAH concentration data. Data frame must include a PAH identifier column, which can be
#' either the USGS parameter code, casrn, or compound name. Not all compounds have casrn numbers. This function
#' requires that pah_dat be in long format; that is, there is a parameter column and site column, where site
#' or sample id is repeated for each compound.
#'
#' @param merge_type whether the merge should being done by USGS parameter codes ('pcode'), Chemical Abstracts
#' Service Registry Number ('casrn'), or compound name ('name'). Note that alkylated PAHs
#' (e.g., C2-alkylated naphthalene) do not have
#' casrn numbers, and metadata on these compounds will not be merged if merge_type = "casrn".
#' @param merge_col column name that contains the PAH identifier by which to merge
#' @return A data frame with all of the columns from the user-supplied pah_dat, as well as columns from the
#' built in data frame pah_compounds. These included...
#' @export
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @examples

get_compound_info <- function(pah_dat, merge_type = "name", merge_col = "parameter_nm"){
  temp_dat <- pah_dat
  if (!merge_type %in% c("name", "pcode", "casrn")) {
    stop('merge_type must be set to either "name", "pcode", or "casrn"')
  }
  if (!merge_col %in% names(pah_dat)) {
    stop(paste("Column", merge_col, "not found in pah_dat."))
  }

  if (merge_type == "name") {
    pah_compounds$Parameter_stripped <- tolower(gsub(pattern = "[[:punct:]]|[[:blank:]]|alkylated|s", replacement = "", pah_compounds$Parameter))
    temp_dat$Parameter_stripped <- tolower(gsub(pattern = "[[:punct:]]|[[:blank:]]|alkylated|s", replacement = "", pah_dat[[merge_col]]))
    temp_dat <- left_join(temp_dat, pah_compounds, by = "Parameter_stripped") %>%
      select(-Parameter_stripped)

    if (any(is.na(temp_dat$EPApriority16))) {
      warning(paste0("Not all compound names in column ", merge_col, " matched a compound name in the built-in data frame pah_compounds, so function did not perform a complete merge. Rows with compound names that did not match are included and were given NA values for the new columns."))
    }
  }
  if (merge_type == "pcode") {
    if (is.numeric(temp_dat[[merge_col]])) {
      warning("Pcodes provided are numeric, converted to character string to match data type in pah_compounds")
      temp_dat[merge_col] <- as.character(temp_dat[[merge_col]])
    }
    temp_dat <- left_join(temp_dat, pah_compounds, by = setNames('pcode', merge_col))

    if (any(is.na(temp_dat$EPApriority16))) {
      warning(paste0("Not all codes in column ", merge_col, " matched a code in the built-in data frame pah_compounds, so function did not perform a complete merge. Rows with compound names that did not match are included and were given NA values for the new columns."))
    }
  }
  if (merge_type == "casrn") {

    temp_dat <- left_join(temp_dat, pah_compounds$pcode, by = setNames("casrn" = merge_col))

    if (any(is.na(temp_dat$EPApriority16))) {
      warning(paste0("Not all codes in column ", merge_col, " matched a code in the built-in data frame pah_compounds, so function did not perform a complete merge. Rows with codes that did not match are included and were given NA values for the new columns."))
    }
  }
  return(temp_dat)

}
