#' simulation based non-conformal mass scores
#'
#' this wrapper assumes a model has been fit relative to the 2 mode mixture
#' of linear regression (similar to \code{flexmix}) for the underlying mass
#' non-conformal mass scores
#'
#' @param x true x
#' @param y true y
#' @param model_params list of model parameters of with names as follows
#' \itemize{
#'  \item{\code{beta}:}{ a vector, length 2}
#'  \item{\code{intercept}:}{ a vector, length 2}
#'  \item{\code{sd}:}{ a vector, length 2}
#'  \item{\code{proportions}:}{ a vector, length 2}
#' }
#' that represent a conditional distribution defined by
#' \deqn{y \sim p_1 M_1(x) + p_2 M_2(x)}{y ~ p1 M1(x) + p2 M2(x)}
#' where \eqn{(M_1,M_2)(x) \sim right\_intercept + N(right\_beta \cdot x, right\_sd)}{(M1,M2)(x) ~ right_intercept + N(right_beta*x, right_sd)} and
#' \eqn{(p_1,p_2) = right]_proportions}{(p1,p2) = right_proportions}
#' @param sim_num number of simulations to base the kernel density estimate of
#' the cde
#'
#' @return nonconformal score for (x,y) pair under specified model
#' @export
#'
simulation_nonconformal_mass_score_mix_reg <- function(x,y, model_params, sim_num = 100){
  y_values <- gen_mix_reg(x = rep(x, sim_num), beta = model_params$beta,
                          intercept = model_params$intercept,
                          sd = model_params$sd,
                          proportions = model_params$prop)

  density_estimate <- ks::kde(y_values, h = .1)
  prob <- predict(density_estimate, x = y_values)
  prob_x <- predict(density_estimate, x = y)

  return(mean(prob <= prob_x))
}


#' simulation based conformal cde scores
#'
#' this wrapper assumes a model has been fit relative to the 2 mode mixture
#' of linear regression (similar to \code{flexmix}) for the underlying mass
#' non-conformal mass scores
#'
#' @param x true x
#' @param y true y
#' @param model_params list of model parameters of with names as follows
#' \itemize{
#'  \item{\code{beta}:}{ a vector, length 2}
#'  \item{\code{intercept}:}{ a vector, length 2}
#'  \item{\code{sd}:}{ a vector, length 2}
#'  \item{\code{proportions}:}{ a vector, length 2}
#' }
#' that represent a conditional distribution defined by
#' \deqn{y \sim p_1 M_1(x) + p_2 M_2(x)}{y ~ p1 M1(x) + p2 M2(x)}
#' where \eqn{(M_1,M_2)(x) \sim right\_intercept + N(right\_beta \cdot x, right\_sd)}{(M1,M2)(x) ~ right_intercept + N(right_beta*x, right_sd)} and
#' \eqn{(p_1,p_2) = right]_proportions}{(p1,p2) = right_proportions}
#' @param sim_num number of simulations to base the kernel density estimate of
#' the cde
#'
#' @return conformal score for (x,y) pair under specified model
#' @export
#'
simulation_conformal_cde_score_mix_reg <- function(x,y, model_params,
                                                   sim_num = 100){
  y_values <- gen_mix_reg(x = rep(x, sim_num), beta = model_params$beta,
                          intercept = model_params$intercept,
                          sd = model_params$sd,
                          proportions = model_params$prop)

  density_estimate <- ks::kde(y_values, h = .1)
  prob_x <- predict(density_estimate, x = y)

  return(prob_x)
}



#' calculate some conformal scores related to simulation based CDE
#'
#' #' this wrapper assumes a model has been fit relative to the 2 mode mixture
#' of linear regression (similar to \code{flexmix}) for the underlying mass
#' non-conformal mass scores
#'
#'
#' @param x true x
#' @param y true y
#' @param model_params list of model parameters of with names as follows
#' \itemize{
#'  \item{\code{beta}:}{ a vector, length 2}
#'  \item{\code{intercept}:}{ a vector, length 2}
#'  \item{\code{sd}:}{ a vector, length 2}
#'  \item{\code{proportions}:}{ a vector, length 2}
#' }
#' that represent a conditional distribution defined by
#' \deqn{y \sim p_1 M_1(x) + p_2 M_2(x)}{y ~ p1 M1(x) + p2 M2(x)}
#' where \eqn{(M_1,M_2)(x) \sim right\_intercept + N(right\_beta \cdot x,
#' right\_sd)}{(M1,M2)(x) ~ right_intercept + N(right_beta*x, right_sd)} and
#' \eqn{(p_1,p_2) = right]_proportions}{(p1,p2) = right_proportions}
#' @param sim_num number of simulations to base the kernel density estimate of
#' the cde
#' @param h kernel density estimate bandwidth parameter
#' @param grid discretizes grid that y should range over (needs to be equally
#' spaced)
#'
#' @return a data.frame with 3 columns (and 1 row)
#' \itemize{
#'  \item{\code{non_conformal_mass_full_sim_based}:}{Estimated amount of mass
#'  with CDE above \code{y} value using just the simulations.}
#'  \item{\code{empirical_cde}:}{Estimated CDE value.}
#'  \item{\code{non_conformal_mass_ecde_and_sim_based}:}{Estimated amount of
#'  mass with CDE above \code{y} value (using \code{grid})}
#' }
#' @export
simulation_conformal_scores_mix_reg <- function(x,y, model_params, sim_num = 100,
                                                h = .4,
                                                grid = seq(-15, 15, by = .01)){



  y_values <- gen_mix_reg(x = rep(x, sim_num), beta = model_params$beta,
                          intercept = model_params$intercept,
                          sd = model_params$sd,
                          proportions = model_params$prop)

  density_estimate <- ks::kde(y_values, h = h)
  prob <- predict(density_estimate, x = y_values)
  prob_x <- predict(density_estimate, x = y)
  prob_xx <- predict(density_estimate, x = grid)

  # delta_xx creation
  delta_xx <- unique(diff(grid))
  assertthat::assert_that(length(delta_xx) == 1 |
                            all(abs(diff(delta_xx) < 1e-10)),
                          msg = paste("difference between grid values not the",
                                      "same - please correct"))

  if (length(delta_xx) > 1){
    delta_xx <- delta_xx[1]
  }

  df_out <- data.frame(non_conformal_mass_full_sim_based = mean(prob >= prob_x),
                       empirical_cde = prob_x,
                       non_conformal_mass_ecde_and_sim_based =
                         sum(prob_xx[prob_xx >= prob_x])*delta_xx)
  return(df_out)
}

