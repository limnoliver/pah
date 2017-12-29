#' Default PAH pcodes
#'
#' A vector of pcodes related to PAH analysis identified in Baldwin et al. (2017).
#'
#' @examples
#' length(pah_pcodes)
"pah_pcodes"

#' PAH compound metadata
#'
#' A table of PAH compound names, abbreviations, pcodes, and supporting information such as
#' molecular weights and whether the compound is a parent or alkylated compound.
#'
#' @format A data frame with 84 observations of 18 variables:
#' \describe{
#'   \item{pcode}{USGS parameter code}
#'   \item{schedule}{USGS lab schedule}
#'   \item{Parameter}{Compound name}
#'   \item{Abbreviation}{Compound short name}
#'   \item{EPApriority16}{logical, is this compound one of the 16 EPA priority compounds?}
#'   }
#' @examples
#' head(pah_compounds$Parameter)
"pah_compounds"


