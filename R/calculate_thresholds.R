#' calc_tox_thresholds
#'
#' Summarizes PAH concentrations relative to toxicity thresholds including
#' 1) the concentration above which adverse outcomes to aquatic biota are expected
#' (probable effect concentration or PEC), 2) the concentration below which adverse
#' outcomes to aquatic biota are unlikely (threshold effect concentration or TEC),
#' and 3) a threshold that accounts for the bioavailability of PAHs based on the
#' organic carbon content of the sediment (sum equilibrium-partitioning sediment
#' benchmark toxicity unit or ESBTU).
#'
#' @export
#' @param compound_info The output dataframe from `get_compound_info`, which contains sample
#' concentrations as well as compound-specific information, including whether the compound is one of
#' 16 EPA priority compounds and compound-specific toxicity.
#' @param sample_column string, column that contains unique sample identifier
#' @param conc_column string, column that contains sample concentrations
#' @param compound_column string, column that contains compound names. For this function to work properly,
#' there must be a compound named "TOC" for total organic carbon (for TOC normalization in ESBTU calculation),
#' with units of percent TOC.
#' @param conc_unit string, the units of PAH concentrations,
#' either "ppb" (ug/kg) or "ppm" (mg/kg).
#' @return a list of two dataframes. 'results_bysite' reports sample-specific thresholds, while
#' 'results_summary' reports summary statistics of thresholds across all samples.
#' @importFrom rlang sym
#' @import dplyr
#' @examples

calc_tox_thresholds <- function(compound_info, sample_column, conc_column, compound_column, conc_unit = "ppb") {
  # make column names dplyr-ready
  quo_sample_column <- sym(sample_column)
  quo_conc_column <- sym(conc_column)
  quo_compound_column <- sym(compound_column)

  if (conc_unit != "ppb" & conc_unit != "ppm") {
    warning('conc_unit must be set to either ppb or ppm')
  }

  toc_test <- which(compound_info[,compound_column] %in% "TOC")
  tec <- ifelse(conc_unit == 'ppb', 1610, 1.610)
  pec <- ifelse(conc_unit == 'ppb', 22800, 22.8)

  pec_tec <- filter(compound_info, which(EPApriority16 == TRUE)) %>%
    group_by(!!quo_sample_column) %>%
    summarize(sum_EPA16 = sum(!!quo_conc_column)) %>%
    mutate(tec_ratio = sum_EPA16/tec,
           pec_ratio = sum_EPA16/pec) %>%
    select(!!quo_sample_column, sum_EPA16, tec_ratio, pec_ratio)

  toc_dat <- filter(compound_info, !!quo_compound_column == "TOC") %>%
    mutate(f_TOC = (!!quo_conc_column)/100) %>%
    select(!!quo_sample_column, f_TOC)

  esbtu_dat <- filter(compound_info, !is.na(coc_pah_fcv)) %>%
    left_join(toc_dat, by = !!quo_sample_column) %>%
    mutate(conc_ug_g = ifelse(conc_unit == 'ppb', round((!!quo_conc_column)/f_TOC, 0)/1000, round((!!quo_conc_column)/f_TOC, 0))) %>%
    mutate(esbtu = conc_ug_g/coc_pah_fcv) %>%
    group_by(!!quo_sample_column) %>%
    summarize(n_esbtu = n(),
              sum_esbtu = sum(esbtu, na.rm = T))

  site_results <- left_join(pec_tec, esbtu_dat, by = !!quo_sample_column)

  site_results_summary <- data.frame(unique_id = c("TEC", 'PEC', 'ESBTU'),
                                     mean_EPApriority16_conc = rep(mean(site_results$sum_EPA16), 3),
                                     n_sites = rep(nrow(site_results), 3),
                                     mean_ratio = c(mean(site_results$tec_ratio), mean(site_results$pec_ratio), mean(site_results$sum_esbtu)),
                                     median_ratio = c(median(site_results$tec_ratio), median(site_results$pec_ratio), median(site_results$sum_esbtu)),
                                     sd_ratio = c(sd(site_results$tec_ratio), sd(site_results$pec_ratio), sd(site_results$sum_esbtu)),
                                     min_ratio = c(min(site_results$tec_ratio), min(site_results$pec_ratio), min(site_results$sum_esbtu)),
                                     max_ratio = c(max(site_results$tec_ratio), max(site_results$pec_ratio), max(site_results$sum_esbtu)))

  out <- list(site_results, site_results_summary)
  names(out) <- c('results_bysample', 'results_summary')

  return(out)
}
