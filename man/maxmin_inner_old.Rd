% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate_funtional.R
\name{maxmin_inner_old}
\alias{maxmin_inner_old}
\title{minimum distance to contain all points of first set with union of balls of
second set}
\usage{
maxmin_inner_old(df_row, df_col)
}
\arguments{
\item{df_row}{data.frame, rows are observations, all columns are used in
distance calculation.}

\item{df_col}{data.frame, rows are observations, all columns are used in
distance calculation.}
}
\value{
single minimum distance scalar
}
\description{
Inner function, slow version. For each row in \code{df_row}, the minimum
distance to a point in \code{df_col} is calculated and then the maximum of
this vector is returned.
}
