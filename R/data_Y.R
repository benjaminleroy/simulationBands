#' generate 2d "Y" distribution
#'
#' data is a mixture of 1/3 from y ~ norm(mu = 0, sd = 1),
#' x ~ unif(min = -6,max = -1)
#'
#' and 2/3 from x ~ unif(min = -1, max = 10), y ~ p1 M1(x) + p2 M2(x)
#' where (M1,M2)(x) ~ right_intercept + N(right_beta*x, right_sd) and
#' (p1,p2) = right_proportions
#'
#' @param n number of points (2/3 right side, 1/3 left side)
#' @param right_beta vector length 2, see above
#' @param right_intercept vector length 2, see above
#' @param right_sd vector length 2, see above
#' @param right_proportions vector length 2, see above
#'
#' @return data frame (x,y) sampled from the above distribution
#' @export
#'
#' @examples
#' library(ggplot2)
#' data_all <- generate_split_mixture(n = 1500)
#' ggplot(data_all) +
#'   geom_point(aes(x = x, y = y), alpha = .3)
generate_split_mixture <- function(n = 1500,
                                   right_beta = c(-1,1),
                                   right_intercept = c(0,0),
                                   right_sd = c(1,1),
                                   right_proportions = c(.5,.5)){
  right_X <- stats::runif(n = floor(2/3*n), min = -1, max = 10)
  right_y <- gen_mix_reg(x = right_X,
                         beta = right_beta,
                         intercept = right_intercept,
                         sd = right_sd,
                         proportions = right_proportions)

  left_X <- stats::runif(n = n - floor(2/3*n), min = -6,-1)
  left_y <- stats::rnorm(n = n - floor(2/3*n))

  data_all <- data.frame(x = c(right_X, left_X),
                         y = c(right_y, left_y))

  return(data_all)
}

#' generate 1d "v" distribution from x values
#'
#' Data is generated from y ~ p1 M1(x) + p2 M2(x)
#' where (M1,M2)(x) ~ right_intercept + N(right_beta*x, right_sd) and
#' (p1,p2) = right_proportions
#'
#' @param x x values
#' @param beta vector length 2, see above
#' @param intercept vector length 2, see above
#' @param sd vector length 2, see above
#' @param proportions vector length 2, see above
#'
#' @return data frame (x,y) sampled from the above distribution
#' @export
gen_mix_reg <- function(x, beta = c(-1,1),
                        intercept = c(0,0),
                        sd = c(1,1),
                        proportions = c(.5,.5)) {
  assertthat::assert_that(all(length(beta) == length(intercept),
                              length(beta) == length(sd),
                              length(beta) == length(proportions)),
                          msg = "parameters should be of the correct length")

  n <- length(beta)

  group_id <- sample(1:n, size = length(x), prob = proportions, replace = T)

  beta_x <- beta[group_id]
  intercept_x <- intercept[group_id]
  sd_x <- sd[group_id]

  return(intercept_x + beta_x * x + stats::rnorm(n = length(x), sd = sd_x))
}


#' obtain conditional pdf value for y values for "v" distribution
#'
#' conditional pdf is defined under y ~ p1 M1(x) + p2 M2(x)
#' where (M1,M2)(x) ~ right_intercept + N(right_beta*x, right_sd) and
#' (p1,p2) = right_proportions
#'
#' @param y can be a vector
#' @param x_val a single value
#' @param beta vector length 2, see above
#' @param intercept vector length 2, see above
#' @param sd vector length 2, see above
#' @param proportions vector length 2, see above
#'
#' @return conditional pdf values for y vector
multimode_density_f <- function(y, x_val, beta = c(-1,1),
                                intercept = c(0,0),
                                sd = c(1,1),
                                proportions = c(.5,.5)){
  means <- beta * x_val + intercept

  prob <- t(stats::dnorm(t(matrix(rep(y, each = 2), ncol = 2, byrow = T)),
                  sd = sd, mean =  means)) %*%
    matrix(proportions, nrow = 2)
  return(prob)
}

#' obtain conditional pdf value for y values for "Y" distribution
#'
#' conditional pdf is defined as following:
#'
#' If x < -1 we have y ~ norm(mu = 0, sd = 1), else
#' y ~ p1 M1(x) + p2 M2(x)
#' where (M1,M2)(x) ~ right_intercept + N(right_beta*x, right_sd) and
#' (p1,p2) = right_proportions
#'
#' @param y can be a vector
#' @param x_val a single value
#' @param right_beta vector length 2, see above
#' @param right_intercept vector length 2, see above
#' @param right_sd vector length 2, see above
#' @param right_proportions vector length 2, see above
#'
#' @return conditional pdf values for y vector
multimodel_density_split_mixture <- function(y, x_val,
                                             right_beta = c(-1,1),
                                             right_intercept = c(0,0),
                                             right_sd = c(1,1),
                                             right_proportions = c(.5,.5)){
  if (x_val >= -1){
    prob <- multimode_density_f(y, x_val,
                                beta = right_beta,
                                intercept = right_intercept,
                                sd = right_sd,
                                proportions = right_proportions)
  } else {
    prob <- stats::dnorm(y, mean = 0, sd = 1)
  }

  return(prob)
}



