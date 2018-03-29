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

#' PAH source ratios
#'
#' A table of various PAH sources and associated PAH ratios. These source ratios come from various sources,
#' and the citation is included in the table. These ratios can be used to calculate distance in a
#' double ratio plot.
#'
#' @format A data frame with 36 observations of 9 variables:
#' \describe{
#'   \item{sourceCategory}{Type of material sampled.}
#'   \item{sourceNo}{Unique identifier for the source material.}
#'   \item{sourceID}{Source descriptor.}
#'   \item{abbrev}{Source abbreviation (for plotting).}
#'   \item{Reference}{Published paper which contains the source ratios.}
#'   \item{Anth_AnthPhen}{The ratio of anthracene to anthracene plus phenanthrene.}
#'   \item{Flua_FluaPyr}{The ratio of fluoranthene to fluoranthene plus pyrene.}
#'   \item{Baa_BaaCh}{The ratio of benzo(a)anthracene to benzo(a)anthracene plus chrysene.}
#'   \item{Indpy_IndpyBghip}{The ratio of indeno(1,2,3-cd)pyrene to indeno(1,2,3-cd)pyrene plus benzo(g,h,i)perylene.}
#'   }
#' @examples
#' head(source_ratios)
"source_ratios"

#' PAH source concentrations
#'
#' A table of various PAH sources and the concentration of PAHs measured in each source. Most source
#' PAH concentrations are the EPA 16 priority compounds, but some concentrations are a sum of more than
#' 16 compounds. These can be used to rule out various sources by comparing to sample concentrations; e.g.,
#' if source concentrations are less than sample concentrations, that source by definition cannot be the
#' primary source of contamination. These concentrations have been extracted from a number of published
#' papers.
#'
#' @format A data frame with 36 observations of 9 variables:
#' \describe{
#'   \item{source}{Brief description of the source material.}
#'   \item{reference}{Citation for paper from which concentration was extracted.}
#'   \item{units}{Unit if source concentration.}
#'   \item{conc_mgkg}{Measured concentration.}
#'   \item{PAHs_used}{The number of PAH compounds that were summed to calculate concentration.}
#'   }
#' @examples
#' head(source_conc)
"source_conc"

