#' get_pah
#'
#' Retrieve PAH data from NWIS. If you already have your PAH data in-hand, there is no
#' need to run this function.
#'
#' @export
#' @param site vector of character strings representing USGS site IDs
#' @param pcodes vector of pcodes to call from NWIS, default set to 84 pcodes that specify PAH-related pcodes
#' @param startDate start date from which to pull data, in YYY-MM-DD format
#'
#' @importFrom dataRetrieval readNWISqw
#' @examples

getPAH <- function(sites, pcodes = pah_pcodes, startDate, endDate) {

  pahs <- readNWISqw(siteNumbers = sites, parameterCd = pcodes, startDate, endDate)

}

