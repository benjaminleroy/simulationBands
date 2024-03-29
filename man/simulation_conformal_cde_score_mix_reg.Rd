% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flexmix_poor_fit_Y.R
\name{simulation_conformal_cde_score_mix_reg}
\alias{simulation_conformal_cde_score_mix_reg}
\title{simulation based conformal cde scores}
\usage{
simulation_conformal_cde_score_mix_reg(x, y, model_params, sim_num = 100)
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
where \eqn{(M_1,M_2)(x) \sim right\_intercept + N(right\_beta \cdot x, right\_sd)}{(M1,M2)(x) ~ right_intercept + N(right_beta*x, right_sd)} and
\eqn{(p_1,p_2) = right]_proportions}{(p1,p2) = right_proportions}}

\item{sim_num}{number of simulations to base the kernel density estimate of
the cde}
}
\value{
conformal score for (x,y) pair under specified model
}
\description{
this wrapper assumes a model has been fit relative to the 2 mode mixture
of linear regression (similar to \code{flexmix}) for the underlying mass
non-conformal mass scores
}
