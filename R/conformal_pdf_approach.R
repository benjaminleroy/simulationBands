
#' Create the "S" matrix and "b" vector
#'
#' An internal function that create the "S" matrix and "b" vector for the
#' "conformal pdf".
#'
#' @details
#' The S matrix is a the upper triangle, where the non-zero values have column
#' information that is \eqn{P(y: f(y|x) \in [lambda_{i}, lambda_{i+1})} where
#' \eqn{\lambda_{i}} are the level set ups associated with conformal scores
#' values.
#'
#' The b vector is just (i-1)/n
#'
#' @param n number of observations in the calibration set
#' @param prob_cde \eqn{P(y: f(y|x) \in [lambda_{i}, lambda_{i+1})} vector. See
#' details for more information.
#' @param prime boolean if we should return S' and b'
#'
#' @return list of S matrix and b vector.
#' @export
create_S_b_matrix <- function(n, prob_cde, prime = FALSE){
  # Equalities
  # 1. conformal constraints

  m1.2 <- diag(prob_cde)

  if (prime){
    return(list(S = m1.2, b = rep(1/(n+1), n+1)))
  }

  m1.1 <- matrix(0, nrow = n + 1, ncol = n + 1)
  m1.1[lower.tri(m1.1, diag = T)] <- 1
  m1.1 <- m1.1[, ncol(m1.1):1]

  m1.3 <- m1.1 %*% m1.2

  lhs1 <- m1.3

  rhs1 <- 1 - (n:0)/(n+1)

  return(list(S = lhs1, b = rhs1))
}


#' stepwise transform of CDE to become comformal
#'
#' Solving for the x in the following general equation
#' \deqn{\min_{x, z} 1/2(Sx-b)^T(Sx-b) +
#'       \lambda (\alpha z^Tz + (1-\alpha) 1^Tz)}
#' subject to some of the following equations
#' \tabular{rl}{
#' \eqn{A_1}: & \eqn{z_i \geq x_i - x_{i+1}}  \\
#' \eqn{A_2}: & \eqn{z_i \geq - x_i + x_{i+1}}  \\
#' \eqn{A_3}: & \eqn{x_i \geq 0}  \\
#' \eqn{A_4}: & \eqn{x_i - x_{i+1} \geq 0}  \\
#' }
#'
#' With \eqn{A_1, A_2, A_3} constraints for \code{monotonically_increasing} is
#' \code{FALSE} and \eqn{A_1, A_3, A_4} if \code{monotonically_increasing} is
#' \code{TRUE}.
#'
#' @param n number of conformal values (length of \code{prob_cde} - 1)
#' @param prob_cde vector of mass of probability (defined by fitted density)
#' between each conformal based level grouping
#' @param lambda optimization constraint (weight for smoothing)
#' @param alpha elasticnet penalty (recommend default at \code{0})
#' @param monotonically_increasing constraint to require scaling to increase
#' (recommend default at \code{FALSE})
#' @param delta stabilizing constant (added to diagonal if needed)
#' @param scaling_constraint boolean if we add a constraint to make sure that
#' the final density has mass 1 (currently set at \code{FALSE} - but not for
#' a good reason).
#' @param prime boolean if we should return S' and b'
#'
#' @return solution of above equation
#' @export
stepwise_conformal_cde_update <- function(n, prob_cde,
                                          lambda = -1,
                                          alpha = 0,
                                          monotonically_increasing = F,
                                          delta = 1e-12,
                                          scaling_constraint = F,
                                          prime = FALSE){

  assertthat::assert_that(length(prob_cde) == n + 1,
                          msg = "size correct")
  assertthat::assert_that(lambda == -1 | lambda > 0,
                          msg = "lambda correct")
  if (lambda != -1){
    assertthat::assert_that(alpha >=0 & alpha <= 1,
                            msg = "alpha correct")
  }

  # Pre processing: S, b calculation --------------
  eq_mat <- create_S_b_matrix(n, prob_cde, prime = prime)
  S <- eq_mat[[1]]
  b <- eq_mat[[2]]

  eq_mat_full <- create_S_b_matrix(n, prob_cde, prime = F)
  S_full <- eq_mat_full[[1]]
  b_full <- eq_mat_full[[2]]

  if (lambda == -1){
    solution = solve(S) %*% b
    return(solution)
  }

  if (alpha*lambda >= delta){
    delta <- alpha*lambda
  }

  # loss function ------------------------
  ## c(x) = x^T D_mat x - dvec^T x
  D_mat.1 <- 1/2 * t(S) %*% S
  D_mat <- rbind(cbind(D_mat.1, matrix(0, nrow = nrow(D_mat.1),
                                       ncol = nrow(D_mat.1) - 1)),
                 cbind(matrix(0, nrow = nrow(D_mat.1) - 1,
                              ncol = nrow(D_mat.1)),
                       diag(rep(delta, nrow(D_mat.1) - 1))))

  dvec <- c(t(S) %*% b, rep(-lambda*(1-alpha), nrow(D_mat.1) - 1))

  # inequality constraints ------------------------

  ## z_i >= a_i - a_{i+1}
  Amat1.1 <- diag(1, nrow = nrow(D_mat.1), ncol = nrow(D_mat.1))
  Amat1.2 <- pracma::Diag(rep(-1, nrow(D_mat.1)-1),k = +1)
  Amat1.3 <- Amat1.1 + Amat1.2
  Amat1.4 <- Amat1.3[-nrow(Amat1.3),]
  Amat1 <- cbind(Amat1.4, diag(rep(1, nrow(D_mat.1) - 1)))

  ## z_i >= -a_i + a_{i+1}
  Amat2 <- cbind(-Amat1.4, diag(rep(1, nrow(D_mat.1) - 1)))

  ## x >= 0
  Amat3 <- cbind(diag(rep(1, nrow(D_mat.1))), matrix(0, nrow = nrow(D_mat.1),
                                                     ncol = nrow(D_mat.1) - 1))

  if (monotonically_increasing) {
    ## a_{i+1} - a_i >= 0 & drop z_i >= -a_i + a_{i+1}
    Amat4 <- cbind(-Amat1.4, diag(rep(0, nrow(D_mat.1) - 1)))
    Amat <- rbind(Amat1, Amat3, Amat4)
  } else {
    Amat <- rbind(Amat1, Amat2, Amat3)
  }

  b_0_vec <- rep(0, nrow(Amat))

  # "equality" constraint ------------------------
  ## to get correct scaling
  ## S[final_row]^T x = 1
  if (scaling_constraint){
    Amat5 <- cbind(rbind(S_full[nrow(S_full), ],
                         -S_full[nrow(S_full), ]), matrix(0, nrow = 2, ncol = nrow(D_mat.1) - 1))
    b_0_5 <- c(1,-1)

    Amat <- rbind(Amat, Amat5)
    b_0_vec <- c(b_0_vec, b_0_5)
  }


  # solving qp ------------------------
  qp_solve <- quadprog::solve.QP(Dmat = D_mat, dvec = dvec,
                                 Amat = t(Amat), bvec = b_0_vec, factorize = F)

  solution <- qp_solve$solution[1:ncol(S)]

  return(solution)
}


#' calculate conformal scores and split into level groupings
#'
#' @param individual_df data.frame that contains information from a single (x)
#' with the following columns:
#' \itemize{
#'  \item{\code{yy}:}{ range of y value that conformal inference will explore}
#'  \item{\code{fit_density}:}{ fit density estimates for \code{yy}}
#' }
#' @param split_conformal_info_vector a vector of conformal / non-conformal
#' scores from a calibration set
#' @param conformal_score string to indicate if we will be doing use the CDE
#' conformal score (\code{"cde"}) or mass non-conformal score (\code{"mass"})
#'
#' @details
#' For \code{split_conformal_info_vector} need to feed in acceptable conformal
#' or non-conformal scores depending upon \code{conformal_score} string value.
#'
#' In the output, the "cuts" are assoicated with
#' \eqn{y:\hat{f}(y|x) > \lambda_i}
#'
#' @return \code{individual_df} with additional columns, specifically:
#' \itemize{
#'   \item{\code{g_id}:}{ string defining which conformal "cut" of the conformal
#' score values the fit density lives in}
#'   \item{\code{g_id2}:}{ integer factor representation of \code{g_id}}
#'   \item{\code{g_id3}:}{ \code{g_id2}/\eqn{max}(\code{g_id2})}
#'   \item{\code{cut_off_lower}:}{ lower threshold from \code{g_id} as scalar}
#'   \item{\code{cut_off_upper}:}{ upper threshold from \code{g_id} as scalar}
#'}
#' @export
#'
#' @importFrom rlang .data
conformal_pdf_breaks <- function(individual_df, split_conformal_info_vector,
                                 conformal_score = c("cde", "mass")){
  if (length(conformal_score) >  1){
    conformal_score <- "cde"
  }

  assertthat::assert_that(conformal_score %in% c("cde", "mass"),
                          msg = "please select an acceptable conformal score")

  if (conformal_score == "cde"){
    inner_levels_density <- sort(split_conformal_info_vector)

    individual_df_all <-individual_df %>%
      dplyr::mutate(g_id = cut(.data$fit_density,
                        breaks = c(-Inf, inner_levels_density, Inf)),
             g_id2 = cut(.data$fit_density,
                         breaks = c(-Inf, inner_levels_density, Inf),
                         labels = F),
             g_id3 = .data$g_id2/max(.data$g_id2),
             cut_off_lower = cut_to_numeric(.data$g_id, .lower =  TRUE),
             cut_off_upper = cut_to_numeric(.data$g_id, .lower = FALSE))
  } else {# conformal_score == "mass"
    mass_inner_levels <- sort(split_conformal_info_vector)

    delta_y <- unique(diff(individual_df$yy))

    # parameter clean up
    assertthat::assert_that(length(delta_y) == 1 |
                              all(abs(diff(delta_y) < 1e-10)),
                            msg = paste("difference between y values not the",
                                        "same - please correct"))

    if (length(delta_y) > 1){
      delta_y <- delta_y[1]
    }

    individual_df_all <- individual_df %>%
      dplyr::arrange(.data$fit_density) %>%
      dplyr::mutate(proportion_above = 1 - cumsum(.data$fit_density)*delta_y)

    individual_df_all <- individual_df_all %>%
      dplyr::mutate(g_id = cut(.data$proportion_above,
                        breaks = c(-Inf, mass_inner_levels, Inf)),
             g_id2 = cut(.data$proportion_above,
                         breaks = c(-Inf, mass_inner_levels, Inf),
                         labels = F),
             g_id3 = .data$g_id2/max(.data$g_id2),
             cut_off_lower = cut_to_numeric(.data$g_id, .lower =  TRUE),
             cut_off_upper = cut_to_numeric(.data$g_id, .lower = FALSE))
  }

  return(individual_df_all)
}


#' calculate cumulative comparisons
#'
#' @param individual_df_rich  data.frame that contains information from a single (x)
#' with the following columns:
#' \itemize{
#'  \item{\code{yy}:}{ range of y value that conformal inference will explore}
#'  \item{\code{fit_density}:}{ fit density estimates for \code{yy}}
#' }
#' and also contains conformal grouping structure (identified in
#' \code{group_columns}) and that relative to conformal level groupings
#' @param group_columns tidyverse style column groupings (all columns that
#' identify information about the conformal level groupings)
#' @param delta_y difference betwen \code{yy} values (assumed constant), the
#' default is to calculate them internally.
#'
#' @return \code{cumulative_df} contains cumulative mass information for the
#' fit density and true density (currently assumed known...)
#' @export
cumulative_comparisons <- function(individual_df_rich,
                                   group_columns = c("g_id", "g_id2", "g_id3",
                                                     "cut_off_lower",
                                                     "cut_off_upper"),
                                   delta_y = NULL){
  # parameter clean up
  if (is.null(delta_y)){
    delta_y <- unique(diff(individual_df_rich$yy))
    assertthat::assert_that(length(delta_y) == 1 |
                            all(abs(diff(delta_y) < 1e-10)),
                            msg = paste("difference between y values not the",
                                        "same - please correct"))
    if (length(delta_y) > 1){
      delta_y <- delta_y[1]
    }
  }
  #quos
  group_columns_q <- dplyr::enquos(group_columns)
  group_columns <- unname(tidyselect::vars_select(
    dplyr::tbl_vars(individual_df_rich),
    !!!group_columns_q))



  cumulative_df <- individual_df_rich  %>%
    dplyr::group_by_at(dplyr::vars(group_columns)) %>%
    dplyr::summarize(fit_prop = sum(.data$fit_density)*delta_y,
                     true_prop = sum(.data$true_density)*delta_y) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cum_fit_prop = cumsum(.data$fit_prop), # fit prop
                  fit_prop_lag = dplyr::lag(.data$fit_prop, default = 0),
                  cum_fit_prop_lag = cumsum(.data$fit_prop_lag)) %>%
    dplyr::mutate(fit_prob_above_lower = 1-.data$cum_fit_prop_lag) %>%
    dplyr::mutate(cum_true_prop = cumsum(.data$true_prop), # true prop
                  true_prop_lag = dplyr::lag(.data$true_prop, default = 0),
                  cum_true_prop_lag = cumsum(.data$true_prop_lag)) %>%
    dplyr::mutate(true_prob_above_lower = 1- .data$cum_true_prop_lag)

  return(cumulative_df)
}
