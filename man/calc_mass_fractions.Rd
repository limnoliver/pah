% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_mass_fractions.R
\name{calc_mass_fractions}
\alias{calc_mass_fractions}
\title{calc_mass_fractions}
\usage{
calc_mass_fractions(compound_info, sample_column, conc_column,
  compound_column, conc_unit = "ppb", calc_type = "summary",
  plot = FALSE, sample_order = "norm_pah_conc")
}
\arguments{
\item{compound_info}{The output dataframe from `get_compound_info`, which contains sample
concentrations as well as compound-specific information, including whether the compound is one of
16 EPA priority compounds and compound-specific toxicity.}

\item{sample_column}{string, column that contains unique sample identifier}

\item{conc_column}{string, column that contains sample concentrations}

\item{compound_column}{string, column that contains compound names.}

\item{conc_unit}{string, the units of PAH concentrations,
either "ppb" (ug/kg) or "ppm" (mg/kg).}

\item{calc_type}{how to calculate mass fractions, either for each individual sample ('by_sample'),
or by summary statistics across all samples 'summary'. Summary calculates mass fractions for all quartiles,
minimum, mean, and maximum of sample concentrations.}

\item{plot}{logical, whether 'by_sample' should be summarized as a tile plot rather than table.}

\item{sample_order}{string, how the samples should be ordered, either by 'pah_conc', which is
the sum of the EPA 16 priority compounds, or 'norm_pah_conc' which is the TOC-normalized PAH 16 concentration.
Sources are considered "unlikely" when the percent of source to sample is greater than the percent TOC
in the sample, given that PAHs are limited to the organic fraction. Ordering by TOC-normalized PAH concentration
gives a smoother look to the figure, but is less intuitive in terms of sample ordering.}
}
\value{
If calc_type is "summary", each row represents a source, and source mean concentrations,
number of PAHs used, and references are reported alongside percent mass fractions calculated for all quartiles,
minimum, and maximum of all sample concentrations. If calc_type is 'by_sample', a data frame of
n samples x j sources is given, where each cell represents the mass fraction for that sample-source combination.
Sample IDs are given in a column, and columns are named by source ID.
}
\description{
A mass balance approach to ruling out sources of contamination. This function uses published
source PAH concentrations, and calculates the mass fraction (as a percent) of the source that
would be required to account for the reported concentration in the samples. For sources with
multiple observations, the mean value is used. Mass fraction can be calculated for each sample,
or for summary statistics across all samples (min, quartiles, mean, max).
}
