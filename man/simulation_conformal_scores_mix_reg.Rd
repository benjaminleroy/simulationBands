% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flexmix_poor_fit_Y.R
\name{simulation_conformal_scores_mix_reg}
\alias{simulation_conformal_scores_mix_reg}
\title{calculate some conformal scores related to simulation based CDE}
\usage{
simulation_conformal_scores_mix_reg(
  x,
  y,
  model_params,
  sim_num = 100,
  h = 0.4,
  grid = seq(-15, 15, by = 0.01)
)
}
\arguments{
\item{x}{true x}

\item{y}{true y}

\item{model_params}{list of model parameters of with names as follows
\itemize{
\item{\code{beta}:}{ a vector, length 2}
\item{\code{intercept}:}{ a vector, length 2}
\item{\code{sd}:}{ a vector, length 2}
\item{\code{proportions}:}{ a vector, length 2}
}
that represent a conditional distribution defined by
\deqn{y \sim p_1 M_1(x) + p_2 M_2(x)}{y ~ p1 M1(x) + p2 M2(x)}
where \eqn{(M_1,M_2)(x) \sim right\_intercept + N(right\_beta \cdot x,
right\_sd)}{(M1,M2)(x) ~ right_intercept + N(right_beta*x, right_sd)} and
\eqn{(p_1,p_2) = right]_proportions}{(p1,p2) = right_proportions}}

\item{sim_num}{number of simulations to base the kernel density estimate of
the cde}

\item{h}{kernel density estimate bandwidth parameter}

\item{grid}{discretizes grid that y should range over (needs to be equally
spaced)}
}
\value{
a data.frame with 3 columns (and 1 row)
\itemize{
\item{\code{non_conformal_mass_full_sim_based}:}{Estimated amount of mass
with CDE above \code{y} value using just the simulations.}
\item{\code{empirical_cde}:}{Estimated CDE value.}
\item{\code{non_conformal_mass_ecde_and_sim_based}:}{Estimated amount of
mass with CDE above \code{y} value (using \code{grid})}
}
}
\description{
#' this wrapper assumes a model has been fit relative to the 2 mode mixture
of linear regression (similar to \code{flexmix}) for the underlying mass
non-conformal mass scores
}
