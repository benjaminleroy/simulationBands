
#' f function for data generation for example in Lei and Wasserman (2014) paper
#'
#' \eqn{f(x) = (x-1)^2(x+1)}
#'
#' @param x x value (can be a vector)
#'
#' @return f(x)
inner_f <- function(x){
  return((x-1)^2*(x+1))
}

#' g function for data generation for example in Lei and Wasserman (2014) paper
#'
#' \eqn{g(x) = 2 sqrt(x+.5) II(x >= -.5)}{
#'      g(x) = 2\sqrt(x+.5)\mathbb{I}(x\geq -.5)}
#'
#' @param x x value (can be a vector)
#'
#' @return g(x)
inner_g <- function(x){
  value <- 2*sqrt(ifelse(x >= -.5,
                         x+.5, 0))
  return(value)
}

#' inner mixture creation function for \eqn{Y|X = x} for data generation for
#' example in Lei and Wasseerman (2014) paper
#'
#' \deqn{f(x) + [-1,1](\code{group} == 1) g(x)}
#'
#' where
#' 1. \eqn{f(x) = (x-1)^2(x+1)},
#' 2. \eqn{
#'     g(x) = 2 sqrt(x+.5) II(x >= -.5)}{g(x) = 2\sqrt(x+.5)\mathbb{I}(x\geq -.5)},
#'
#' @param x x value (can be a vector)
#' @param group group value (1,2) which mixture the y is from  (can be a vector)
#'
#' @return mixture value f(x) - g(x) or f(x) + g(x)
inner_mixture <- function(x, group){
  value <- ifelse(group == 1,
                  inner_f(x) - inner_g(x),
                  inner_f(x) + inner_g(x))
}

#' sigma function for data generation for example in Lei and Wasserman (2014)
#' paper
#'
#' \deqn{sigma(x) = 1/4 + |x|.}{\sigma(x) = 1/4 + |x|.}
#'
#' @param x x value (can be a vector)
#'
#' @return sigma(x)
#' @export
inner_sigma <- function(x){
  return(1/4 + abs(x))
}


#' Generate data from example in Lei and Wasseerman (2014) paper
#'
#' The generative model is defined as
#' \deqn{X ~ Unif(-1.5, 1.5)}{X \sim Unif(-1.5, 1.5)}
#' \deqn{(Y|X=x) ~ .5 N(f(x)-g(x), sigma^2(x)) + .5 N(f(x)+g(x),
#'       sigma^2(x))}{(Y|X=x) \sim .5 N(f(x)-g(x),
#'        \sigma^2(x)) + .5 N(f(x)+g(x), \sigma^2(x))}
#' where
#' 1. \eqn{f(x) = (x-1)^2(x+1)},
#' 2. \eqn{g(x) = 2 sqrt(x+.5) II(x >= -.5)}{
#' g(x) = 2\sqrt(x+.5)\mathbb{I}(x\geq -.5)}, and
#' 3. \eqn{sigma(x) = 1/4 + |x|.}{\sigma(x) = 1/4 + |x|.}

#' @param n number of observations
#' @param sigma_function a function for the sigma function. The default is
#' described above, but using \code{function(x){0.0001}} would return the means
#' of the functions
#'
#' @return data.frame for (x,y) pairs of \code{n} simulations from the
#' generative process
#' @export
#'
lei_wasserman_data <- function(n = 1000, sigma_function = inner_sigma){
  x <- runif(n = n, min = -1.5, 1.5)
  group_id <- sample(c(1,2), replace = T, size = n)

  mu <- inner_mixture(x, group_id)
  sigma <- sigma_function(x)
  y <- rnorm(n = n,
             mean = mu, sd = sigma )

  return(data.frame(x = x,
                    y = y))

}

#' Generate potential Y|X=x from example in Lei and Wasseerman (2014) paper
#'
#' Generate data from example in Lei and Wasseerman (2014) paper
#'
#' The generative model is defined as
#' \deqn{(Y|X=x) ~ .5 N(f(x)-g(x), sigma^2(x)) + .5 N(f(x)+g(x),
#'      sigma^2(x))}{(Y|X=x) \sim .5 N(f(x)-g(x),
#'      \sigma^2(x)) + .5 N(f(x)+g(x),\sigma^2(x))}
#' where
#' 1. \eqn{f(x) = (x-1)^2(x+1)},
#' 2. \eqn{g(x) = 2 sqrt(x+.5) II(x >= -.5)}{
#'     g(x) = 2\sqrt(x+.5)\mathbb{I}(x\geq -.5)}, and
#' 3. \eqn{sigma(x) = 1/4 + |x|.}{\sigma(x) = 1/4 + |x|.}
#'
#' @param x x value (can be a vector)
#' @param n number of simulations per x value
#' @param sigma_function a function for the sigma function. The default is
#' described above, but using \code{function(x){0}} would return the means of
#' the functions
#' @param verbose boolean, if we should enumerate the how many \code{x} values
#' it has created simulations for.
#'
#' @return data.frame for (x,y) pairs of \code{n} simulations from the
#' generative process
#' @export
#'
#' @return list of vectors of simulations (ordered by the x values)
#' @export
lei_wasserman_data_conditional_simulate <- function(x, n = 200,
                                                    sigma_function = inner_sigma,
                                                    verbose = T){

  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Simulating [:bar] :percent eta: :eta",
      total = length(x), clear = FALSE, width = 38)
  }

  sim_list <- list()
  id <- 1
  for (x_value in x){

    xx <- rep(x_value, n)
    group_id <- sample(c(1,2), replace = T, size = n)

    mu <- inner_mixture(xx, group_id)
    sigma <- sigma_function(xx)
    yy <- rnorm(n = n,
                mean = mu, sd = sigma )

    sim_list[[id]] <- data.frame(x = xx,
                                 sim = yy)
    if (verbose) {
      pb$tick()
    }
    id <- id + 1
  }

  return(sim_list)

}

#' 0- sigma function for data generation for example in Lei and Wasserman (2014)
#' paper. To help just return the mean lines
#'
#' \deqn{sigma(x) = .0001}{\sigma(x) = .0001}
#'
#' @param x x value (can be a vector)
#'
#' @return sigma_zero(x)
inner_sigma_zero <- function(x){
  return(0.0001)
}


#' Create the CDE Y|X=x for data example from in Lei and Wasseerman (2014) paper
#'
#' The generative model is defined as
#' \deqn{(Y|X=x) ~ .5 N(f(x)-g(x), sigma^2(x)) + .5 N(f(x)+g(x),
#'      sigma^2(x))}{(Y|X=x) \sim .5 N(f(x)-g(x),
#'      \sigma^2(x)) + .5 N(f(x)+g(x),\sigma^2(x))}
#' where
#' 1. \eqn{f(x) = (x-1)^2(x+1)},
#' 2. \eqn{g(x) = 2 sqrt(x+.5) II(x >= -.5)}{
#'     g(x) = 2\sqrt(x+.5)\mathbb{I}(x\geq -.5)}, and
#' 3. \eqn{sigma(x) = 1/4 + |x|.}{\sigma(x) = 1/4 + |x|.}
#'
#' @param x x value (a scalar)
#'
#' @return function to calculate the CDE for any y value
#' @export
cde_lei_wassserman <- function(x) {
  mu <- inner_mixture(rep(x,2), c(1,2))
  sigma <- inner_sigma(x)

  marginal_cde <- function(y) {
    1/2 * dnorm(y, mean = mu[1], sd = sigma) +
      1/2 * dnorm(y, mean = mu[2], sd = sigma)
  }
  return(marginal_cde)
}


#' Caculate CDE Y|X=x from data example from in Lei and Wasseerman (2014) paper
#' on a grid of y values
#'
#' The generative model is defined as
#' \deqn{(Y|X=x) ~ .5 N(f(x)-g(x), sigma^2(x)) + .5 N(f(x)+g(x),
#'      sigma^2(x))}{(Y|X=x) \sim .5 N(f(x)-g(x),
#'      \sigma^2(x)) + .5 N(f(x)+g(x),\sigma^2(x))}
#' where
#' 1. \eqn{f(x) = (x-1)^2(x+1)},
#' 2. \eqn{g(x) = 2 sqrt(x+.5) II(x >= -.5)}{
#'     g(x) = 2\sqrt(x+.5)\mathbb{I}(x\geq -.5)}, and
#' 3. \eqn{sigma(x) = 1/4 + |x|.}{\sigma(x) = 1/4 + |x|.}
#'
#' @param x x value (can be a vector)
#' @param grid vector of y values to be evaluated (must be equally spaced)
#' @param verbose boolean, if we should enumerate the how many \code{x} values
#' it has examined values for.
#'
#' @return data frame with columns
#' \itemize{
#'   \item{\code{x}: }{the associated x values}
#'   \item{\code{y}: }{the associated y values}
#'   \item{\code{cde}: }{CDE for Y=\code{y}|X=\code{x} as described above.}
#'   \item{\code{mass}: }{Approximation of P(z|x and f(z|x) >= f(y|x))}
#'}
#' @export
cde_and_mass_lei_wasserman_discrete <- function(x,
                                                grid = seq(-8,8,
                                                              length.out = 100),
                                                verbose = T){


  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "Simulating [:bar] :percent eta: :eta",
      total = length(x), clear = FALSE, width = 38)
  }

  # parameter creation
  delta_y <- unique(diff(grid))
  assertthat::assert_that(length(delta_y) == 1 |
                            all(abs(diff(delta_y) < 1e-10)),
                                msg = paste("difference between grid values not the",
                                            "same - please correct"))

  if (length(delta_y) > 1){
    delta_y <- delta_y[1]
  }


  # analysis
  df_out <- data.frame()
  for (idx in 1:length(x)) {
    inner_cde_values <- cde_lei_wassserman(x[idx])(grid)

    inner_mass_values <- inner_discrete_mass_cde(inner_cde_values, delta_y)

    df_out <- rbind(df_out,
                    data.frame(x = x[idx],
                               y = grid,
                               cde = inner_cde_values,
                               mass = inner_mass_values))
    if (verbose) {
      pb$tick()
    }
  }
  return(df_out)
}


#' inner function for estimation of probability of mass with CDE above or equal
#' to value
#'
#' @param cde_vec vector of CDE values for Y=y|X=x
#' @param delta_y space between each y value
#'
#' @return mass vector
inner_discrete_mass_cde <- function(cde_vec, delta_y = 16/99){
  ordering <- order(cde_vec,decreasing = F)
  mass_estimate <- cumsum(sort(cde_vec)) * delta_y
  mass_from_cde <- rep(NA, length(cde_vec))
  mass_from_cde[ordering] <- mass_estimate
  return(mass_from_cde)
}


#' calculate true CDE conformal score and approximate mass non-conformal score
#' under Lei & Wasserman 2014 perfect data fit.
#'
#' @param x x value
#' @param y observed y value
#' @param grid good range of potential y values that take up most of mass of
#' CDE, should also be equally spaced.
#'
#' @return data frame with
#' \itemize{
#'   \item{\code{x}: }
#'   \item{\code{y}: }
#'   \item{\code{cde}: }{CDE of true y value}
#'   \item{\code{mass}: }{approximation of P(z|x and f(z|x) >= f(y|x))}
#' }
#' @export
conformal_scores_lei_wasserman <- function(x, y,
                                           grid = seq(-8,8,
                                                      length.out = 100)){

  # parameter creation
  delta_y <- unique(diff(grid))
  assertthat::assert_that(length(delta_y) == 1 |
                            all(abs(diff(delta_y) < 1e-10)),
                          msg = paste("difference between grid values not the",
                                      "same - please correct"))

  if (length(delta_y) > 1){
    delta_y <- delta_y[1]
  }

  y_cde <- cde_lei_wassserman(x)(y)
  y_range_cde <- cde_lei_wassserman(x)(grid)

  df_out <-  data.frame(x = x, y = y,
                        cde = y_cde,
                        mass = sum(y_range_cde[y_range_cde >= y_cde])*delta_y)

  return(df_out)
}

