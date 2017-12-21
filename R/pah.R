#' Default PAH pcodes
#'
#' A vector of pcodes related to PAH analysis identified in Baldwin et al. (2017).
#'
#' @name pah_pcodes
#' @rdname sampleData
#' @docType data
#' @examples
#' length(pah_pcodes)
"pah_pcodes"

#' A table of PAH compound names, abbreviations, pcodes, and supporting information such as
#' molecular weights and whether the compound is a parent or alkylated compound.
#'
#' @name pah_compounds
#' @rdname sampleData
#' @docType data
#' @format A data frame with 85 observations of 20 variables:
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


