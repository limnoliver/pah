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

#' PAH source profiles
#'
#' A table of various PAH sources and the associated 12-PAH compound profiles. These are the source
#' profiles used in Baldwin et al. (2017).
#'
#' @format A data frame with 12 observations of 24 variables:
#' \describe{
#'   \item{Compound}{USGS parameter code}
#'   \item{Abbreviation}{Compound abbreviation used for plotting}
#'   \item{pcode}{USGS parameter code}
#'   \item{Power_plant_emissions}{Power plant emissions profile, from ...}
#'   \item{Residential_heating}{Residential heating profile, from...}
#'   \item{Coal_average}{The average coal profile, from ...}
#'   \item{Coke_oven_emissions}{Profile from coking emissions, from ...}
#'   \item{Diesel_vehicle}{PAH profile from diesel vehicles, from ...}
#'   \item{Gasoline_vehicle}{PAH profile from gasoline vehicles, from ...}
#'   \item{Traffic_tunnel_air}{PAH profiles of traffic tunnel air, from ...}
#'   \item{Vehicle_traffic_avg}{The average PAH profile from vehicle traffic emissions, from ...}
#'   \item{casrn}{The chemical identification number.}
#'   }
#' @examples
#' head(source_profiles)
"source_profiles"


