#' summary_dv
#'
#' Get max and min of discharge for day of the year.
#'
#' @export
#' @param site character USGS site ID
#' @importFrom dataRetrieval readNWISdv
#' @importFrom dataRetrieval renameNWISColumns
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @examples
#' site <- "08279500"
#' sum_dv <- summary_dv(site)
summary_dv <- function(site){
  Date <- doy <- Flow <- ".dplyr.var"

  dv_data <- readNWISdv(siteNumbers=site,
                        parameterCd = "00060", startDate = "", endDate = "")

  dv_summ <- renameNWISColumns(dv_data)
  dv_summ <- mutate(dv_summ, doy = as.numeric(strftime(Date, format = "%j")))
  dv_summ <- group_by(dv_summ, doy)
  dv_summ <- summarise(dv_summ,
                       max = max(Flow, na.rm = TRUE),
                       min = min(Flow, na.rm = TRUE))
  return(dv_summ)

}
