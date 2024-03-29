% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization_diagnostics.R
\name{vis_all_densities}
\alias{vis_all_densities}
\title{compare all densities visually (true, fit, conformal_pdf)}
\usage{
vis_all_densities(
  stepwise_scaling,
  individual_df_rich,
  cumulative_comparisons_df,
  scale = T
)
}
\arguments{
\item{stepwise_scaling}{new scaling for each level set}

\item{individual_df_rich}{data for individual x's potential y values and
true, fit values}

\item{cumulative_comparisons_df}{cumlative information for level sets}

\item{scale}{boolean if we should normalize the solution}
}
\value{
visualization plot
}
\description{
compare all densities visually (true, fit, conformal_pdf)
}
