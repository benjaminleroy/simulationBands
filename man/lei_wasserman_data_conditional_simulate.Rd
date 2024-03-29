% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_lei_wasserman.R
\name{lei_wasserman_data_conditional_simulate}
\alias{lei_wasserman_data_conditional_simulate}
\title{Generate potential Y|X=x from example in Lei and Wasseerman (2014) paper}
\usage{
lei_wasserman_data_conditional_simulate(
  x,
  n = 200,
  sigma_function = inner_sigma,
  verbose = T
)
}
\arguments{
\item{x}{x value (can be a vector)}

\item{n}{number of simulations per x value}

\item{sigma_function}{a function for the sigma function. The default is
described above, but using \code{function(x){0}} would return the means of
the functions}

\item{verbose}{boolean, if we should enumerate the how many \code{x} values
it has created simulations for.}
}
\value{
data.frame for (x,y) pairs of \code{n} simulations from the
generative process

list of vectors of simulations (ordered by the x values)
}
\description{
Generate data from example in Lei and Wasseerman (2014) paper
}
\details{
The generative model is defined as
\deqn{(Y|X=x) ~ .5 N(f(x)-g(x), sigma^2(x)) + .5 N(f(x)+g(x),
     sigma^2(x))}{(Y|X=x) \sim .5 N(f(x)-g(x),
     \sigma^2(x)) + .5 N(f(x)+g(x),\sigma^2(x))}
where
\enumerate{
\item \eqn{f(x) = (x-1)^2(x+1)},
\item \eqn{g(x) = 2 sqrt(x+.5) II(x >= -.5)}{
    g(x) = 2\sqrt(x+.5)\mathbb{I}(x\geq -.5)}, and
\item \eqn{sigma(x) = 1/4 + |x|.}{\sigma(x) = 1/4 + |x|.}
}
}
