#' pah_pca
#'
#' Conducts principle components analysis on the source and sample PAH profiles. This reduces the profiles
#' to manageable dimensions. All important PCA axes are then used to calculate the Euclidean distance between
#' source and samples to determine likely sources.
#'
#' @export
#' @param profiles The output of pah_profiles, which contains PAH profile data for both source and samples.
#' The reported values are proportional concentrations for each compound.
#' @param perc_cutoff The threshold for keeping PCA components in subsequent analyses; all
#' individual components that explain greater than or equal to this cutoff of the total
#' variance in the dataset will be kept in the analysis.
#' @param sample_column column name of unique sample identifier
#' @return A list of three data frames. "pca_dat" contains all PCA components that met the perc_cutoff,
#' as well as a sample_id column and type which identifies whether the values correspond to a sample
#' or source profile. "pca_summary" prints the standard deviation, proportion of variance, and cumulative
#' proportion of variance explained for all components of the PCA. "pca_distance" is a long data frame of all
#' source-sample compindations, and the euclidean distance between the source and sample in the PCA
#' components space.
#' @import dplyr
#' @importFrom tidyr spread drop_na
#' @importFrom stats prcomp
#' @examples

pah_pca <- function(profiles, perc_cutoff = 10, sample_column) {

  quo_sample_column <- sym(sample_column)

  # subset list that comesout of pah_profiles
  all_profiles <- profiles$profiles

  # get sample profiles transposed and ready for pca
  sample_profiles_t <- select(all_profiles, !!quo_sample_column, Abbreviation, prop_conc) %>%
    distinct() %>%
    spread(key = Abbreviation, value = prop_conc)

  sample_profiles_t <- as.data.frame(sample_profiles_t)
  #row.names(sample_profiles_t) <- sample_profiles_t$sample_id
  #sample_profiles_t <- select(sample_profiles_t, -!!quo_sample_column)

  sample_profiles_t$type <- 'sample'

  # get source profiles tansposed and ready for pca
  source_profiles_t <- select(all_profiles, source, Abbreviation, source_prop_conc) %>%
    distinct() %>%
    rename(prop_conc = source_prop_conc)

  names(source_profiles_t)[names(source_profiles_t) %in% 'source'] <- sample_column

  source_profiles_t <- spread(source_profiles_t, key = Abbreviation, value = prop_conc)

  #source_profiles_t <- as.data.frame(source_profiles_t)
  #row.names(source_profiles_t) <- source_profiles_t[, sample_column]
  #source_profiles_t <- select(source_profiles_t, -!!quo_sample_column)

  source_profiles_t$type <- 'source'

  # put together
  profiles_t <- rbind(source_profiles_t, sample_profiles_t)
  #profiles_t[,sample_column] <- row.names(profiles_t)

  # drop any rows (samples) with 0 values (BDL)
  profiles_t <- filter_all(profiles_t, all_vars(. != 0))

  sample_id <- profiles_t[, sample_column]
  #sample_id <- sample_id[[1]]
  type <- profiles_t$type
  #
  profiles_t <- select(profiles_t, -type, -!!quo_sample_column)
  profiles_t <- data.frame(profiles_t)
  row.names(profiles_t) <- sample_id

  # run pca
  pca <- prcomp(profiles_t, center = T, scale = T)
  pca_out <- summary(pca)$importance

  # determine what dimensions the plotting space should be in
  n_pca <- length(which(pca_out[2,] >= (perc_cutoff/100)))
  if (n_pca == 1) {
    n_pca <- 2
  } else {
    n_plots <- factorial(n_pca)/(factorial(2)*factorial(n_pca - 2))
  }

  cum_variance <- round(pca_out[3,n_pca], 2)*100

  # get data in format to be used with ggplot
  pca_plot_dat <- data.frame(pca$x[,1:n_pca])
  pca_plot_dat[, sample_column] <- row.names(pca_plot_dat)
  pca_plot_dat$type <- type

  # now calculate euclidean distances
  structure1 <- select(pca_plot_dat, !!quo_sample_column, type) %>%
    filter(type == 'sample') %>%
    select(!!quo_sample_column)

  structure2 <- select(pca_plot_dat, !!quo_sample_column, type) %>%
    filter(type == 'source') %>%
    select(!!quo_sample_column)

  structure <- expand.grid(structure1[[1]], structure2[[1]])

  pca_plot_dat_sub <- select(pca_plot_dat, -type)

  # left join by sample ids and rename columns
  comp <- rename(structure, sample = Var1, source = Var2) %>%
    left_join(pca_plot_dat_sub, by = c('sample' = sample_column))

  sample_colnames <- paste0('PC', 1:n_pca, '_sam')
  names(comp)[3:(2+n_pca)] <- sample_colnames

  #left join by source ids and rename colums
  comp <- comp %>%
    left_join(pca_plot_dat_sub, by = c('source' = sample_column))

  source_colnames <- paste0('PC', 1:n_pca, '_src')
  names(comp)[(3+n_pca):ncol(comp)] <- source_colnames

  # find squared differences between source and sample of
  # each component to reduce dimensions
  pca_names <- names(pca_plot_dat)[1:n_pca]
  comp_diff <- data.frame(sample = comp$sample,
                          source = comp$source)
  for (i in 1:n_pca) {
    temp_col <- grep(pca_names[i], names(comp))
    temp_diff <- (comp[,temp_col[1]] - comp[,temp_col[2]])^2
    comp_diff[,2+i] <- temp_diff
  }

  # sum squared differences across all pca dimensions,
  # take square root to get euclidean distance
  comp_diff <- comp_diff %>%
    mutate(euc_dist = sqrt(rowSums(.[3:(2+n_pca)]))) %>%
    select(sample, source, euc_dist)

  # outputs
  sum_desc <- paste0("PCA components 1 through ", n_pca, " each explained more than the cutoff of ",
                     perc_cutoff, "% of the total variance in the data and together explained ",
                     cum_variance, "% of the total variance.")
  message(sum_desc)
  out <- list(pca_plot_dat, pca_out, comp_diff)
  names(out) <- c('pca_dat', 'pca_summary', 'pca_distance')
  return(out)
}
