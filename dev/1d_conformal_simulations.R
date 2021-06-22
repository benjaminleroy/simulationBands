
# document hyper parameters

confidence_level <- .6
confidence_level_string <- "60%"
confidence_level_image_string <- "60percent"
n_simulations <- 1000
delta_prop <- .8
delta_prop_string <- "80%"
delta_prop_image_string <- "80percent"

n_sims_containment <- 300
verbose <- T
save_fig <- F
save_image <- save_fig

rerun <- F # to rerun all analysis even if already saved...


# Libraries -----------------

library(tidyverse)
library(devtools)
library(ggrepel)
library(gridExtra)
library(tikzDevice)


load_all("../")
load_all("../../EpiCompare/") # tidy_dist_mat


# Calibration and Test set -------------------

set.seed(1)
df <- tibble(x = runif(1000, min = -1.5, max = 1.5))
df$y <- sapply(df$x, function(x) {simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})

calibration_set <- tibble(x = runif(1000, min = -1.5, max = 1.5))
calibration_set$y <- sapply(calibration_set$x, function(x) {simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})


test_set <- tibble(x = runif(1000, min = -1.5, max = 1.5))
test_set$y <- sapply(test_set$x, function(x) {simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})

test_set_discrete <- tibble(x = seq(-1.5,1.5, by = .01))
test_set_discrete$y <- sapply(test_set_discrete$x, function(x) {simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})


vis_raw <- ggplot(test_set) +
  geom_point(aes(x = x, y = y)) +
  theme_minimal()
vis_raw

if (save_fig){
  ggplot2::ggsave(vis_raw, filename = "quick_images/simulation_distribution_1d.png")

  tikz(file = "quick_images/simulation_distribution_1d.tex")
  print(vis_raw)
  dev.off()
}


# Functions ------------------


#' #' Simulation based conformal score processor
#' #'
#' #' This is a basic one (below is a more complex one...)
#' #'
#' #' @param truth_df data frame of a *single* new observation
#' #' @param simulations_grouped_df grouped data frame of simulated points
#' #' @param data_column_names character vector of column names that define
#' #' locational information of points
#' #'
#' #' @return list of information including the conformal score of the observation
#' #' and internal information about the creation of the prediction regions
#' #'
#' #' @export
#' simulation_based_conformal_1d <- function(truth_df, simulations_grouped_df,
#'                                           data_column_names = c("y"),
#'                                           delta_prop = .8){
#'
#'
#'   assertthat::assert_that(inherits(simulations_grouped_df, "grouped_df"))
#'   group_names <- names(group_keys(simulations_grouped_df))
#'
#'   truth_df_inner <- truth_df %>%
#'     dplyr::select(dplyr::one_of(data_column_names))
#'
#'   simulations_group_df_inner <- simulations_grouped_df %>%
#'     dplyr::select(dplyr::one_of(c(group_names, data_column_names)))
#'
#'   group_info <- simulations_group_df_inner %>% group_keys()
#'
#'   dist_mat <- simulations_group_df_inner %>% group_split() %>%
#'     do.call(rbind, .) %>% # match group_info ordering
#'     select(one_of(data_column_names)) %>%
#'     dist() %>% as.matrix()
#'
#'   tdm_sims <- EpiCompare::tidy_dist_mat(dist_mat, group_info, group_info)
#'
#'   # sigma selection
#'
#'   sigma_size <- c("20%" = .2, "25%" = .25, "30%" = .3,
#'                   "35%" = .35, "40%" = .4, "45%" = .45)
#'
#'   percentage <- names(sigma_size)[stats::quantile(as.matrix(tdm_sims), sigma_size) > 0][1]
#'
#'
#'   # rank_df
#'   pseudo_density_df <- EpiCompare::distance_psuedo_density_function(
#'     tdm_sims,
#'     sigma = percentage, df_out = T) %>%
#'     mutate(ranking = rank(psuedo_density,ties.method = "min")) #spelling error... :(
#'
#'   assertthat::assert_that(all(!is.na(pseudo_density_df$psuedo_density)),
#'                           msg = paste("internal error in",
#'                                       "distance_psuedo_density_function",
#'                                       "function's sigma selection."))
#'
#'   mm_df <- maxmin_distance_vector(truth_df = truth_df_inner,
#'                                    simulations_grouped_df = simulations_group_df_inner,
#'                                    data_column_names = c("y"),
#'                                    .all_info = F)
#'
#'   proportion_points_not_included <- 1 - delta_prop
#'
#'   top_points <- simulations_group_df_inner %>%
#'     left_join(pseudo_density_df, by = group_names) %>%
#'     mutate(keep = ranking > ceiling(proportion_points_not_included*nrow(simulations_group_df_inner))) %>%
#'     ungroup() %>% filter(keep) %>%
#'     select(one_of(data_column_names))
#'
#'
#'   mm_delta <- get_delta_nn(top_points)
#'
#'   containment_df <- pseudo_density_df %>%
#'     left_join(mm_df, by = group_names) %>%
#'     mutate(delta_close = maxmin_dist < mm_delta)
#'
#'
#'   conformal_score <- max(c(containment_df$ranking[containment_df$delta_close],
#'                            0))
#'
#'   return(list(conformal_score = conformal_score, containment_df = containment_df,
#'               mm_delta = mm_delta,
#'               truth_df_inner = truth_df_inner,
#'               simulations_group_df_inner = simulations_group_df_inner,
#'               parameters = c("mm_delta_prop" = proportion_points_not_included,
#'                              "sigma_percentage" = percentage)))
#' }
#'

#' Mode clustering for euclidean data (distance based)
#'
#' Walks up a guassian kernel density estimate to find modes
#'
#' Similar to more complex code developed by Yen-Chi Chen and Chris Genovese
#'
#' @param X_mat matrix of data points (in rows), that define the kernel density
#' @param G_mat matrix of data points (in rows) that need to be "walked" to
#' their modes
#' @param sigma density scalar / sigma value
#' @param maxT int, maximum number of iterations
#' @param eps float, difference between 2 sequential points for which to stop
#' walking toward the mode for a give point's walk
#' @param verbose boolean, if this progression should be verbose
#' @param list_out if we should return every step of the process (for fun
#' visualization purposes and visual checks only)
#'
#' @return list of :
#'  - matrix of \code{G_mat} points walked up all the way
#'  - if \code{list_out} is \code{TRUE} then list of matrices for each step,
#'  else \code{NULL}
#' @export
psuedo_density_mode_cluster_1d <- function(X_mat, G_mat = X_mat,
                            sigma,
                            maxT = 30,
                            eps = 1e-05,
                            verbose = TRUE,
                            list_out = FALSE){
  n_X <- nrow(X_mat)
  n_G <- nrow(G_mat)
  G_mat_out <- G_mat # just a copy

  if (list_out){
    G_mat_list <- list()
    G_mat_list[[1]] <- G_mat_out
  }


  error <- rep(1e08, n_G) # initial error = massive error
  max_error <- 1e08
  t <- 0

  if (verbose){
    pb <- progress::progress_bar$new(
      format = "mode clustering [:bar] :percent eta: :eta",
      total = maxT, clear = FALSE, width = 60)
  }


  while (max_error > eps && t < maxT){
    for (g_idx in 1:n_G) {

      current_vs_X_dist_vals <- apply(
                                  sweep(X_mat, 2, G_mat_out[g_idx,], "-")^2,
                                  1,
                                  function(row_vals) sqrt(sum(row_vals)))

      pseudo_density <- exp(-current_vs_X_dist_vals^2/sigma)

      denominator <- sum(pseudo_density)

      if (denominator == 0) {
        stop("Error: no distance between a observation in G_mat and all observations in X_mat > 0.");
      }

      scaled_props <- pseudo_density / denominator

      new_point <- apply(
                    sweep(X_mat, 1, scaled_props, "*"),
                    2, sum)

      # get error and update
      error[g_idx] <- sqrt(sum((new_point-G_mat_out[g_idx,])^2))
      G_mat_out[g_idx, ] <- new_point

    }
    if (list_out){
      G_mat_list[[t+2]] <- G_mat_out
    }

    t <- t + 1
    max_error <- max(error)
    if (verbose){
      pb$tick()
    }
  }
  if (!list_out){
    G_mat_list <- list()
    }

  return(list(G_mat_out, G_mat_list))
}


testthat::test_that("test psuedo_density_mode_cluster_1d", {
  set.seed(1)
  data <- matrix(c(rnorm(mean = 2,n = 100),
            rnorm(mean = -2, n = 100)),
            byrow = T, nrow = 100)

  sigma <- quantile(as.matrix(dist(data)), .2)

  out <- psuedo_density_mode_cluster_1d(data, sigma = sigma,
                                        maxT = 60, eps = 1e-05,
                                 verbose = FALSE, list_out = TRUE)

  df_out <- do.call(rbind, lapply(1:length(out[[2]]),
                        function(idx) data.frame(out[[2]][idx]) %>%
                          mutate(idx = idx) %>%
                          rownames_to_column()))

  if (FALSE){
    plotly::ggplotly(df_out %>% ggplot() +
                       geom_point(aes(x=X1, y =X2, frame = idx)) +
                       geom_path(aes(x=X1, y=X2, group = rowname))
    )
  }


  dist_mat <- as.matrix(dist(out[[1]]))

  adjmatrix <- dist_mat <= .01
  ig <- igraph::graph_from_adjacency_matrix(adjmatrix, mode = "undirected")
  groupings <- igraph::components(ig, mode = "strong")$membership

  # 2 completely disjoint groups 1:50 and 51:100
  testthat::expect_equal(length(unique(groupings[1:50])), 1)
  testthat::expect_equal(length(unique(groupings[51:100])), 1)
  testthat::expect_equal(length(unique(groupings)), 2)

})

#' Find mode clusters of euclidean objects based on distance - densities
#'
#' @param df_info data frame with rows of observations to mode cluster
#' @param position character vector of columns that define the location of the
#' points
#' @param naming_info a data frame with information about the data points (same
#' number of rows as \code{df_info})
#' @param sigma  density scalar / sigma value
#' @param maxT int, maximum number of iterations
#' @param eps float, difference between 2 sequential points for which to stop
#' walking toward the mode for a give point's walk
#' @param diff_eps float, if the final step of each of the points is within
#' this distance from each-other they will be grouped together.
#' @param verbose boolean, if this progression should be verbose
#'
#' @return data frame with \code{naming_info} and a \code{grouping} column which
#' indicates which mode group each observation is in.
#' @export
mode_clustering_1d <- function(df_info, position, naming_info = NULL, sigma, maxT = 50,
                              eps = 1e-05, diff_eps = 1e-05, verbose = TRUE){
  mat_info_raw <- df_info[,position] %>% as.matrix()
  if (!is.null(naming_info)){
    g_names = df_info[, naming_info, drop = F]
  }

  out <- psuedo_density_mode_cluster_1d(mat_info_raw, sigma = sigma,
                                        maxT = maxT, eps = eps,
                                        verbose = verbose, list_out = FALSE)

  dist_mat <- as.matrix(dist(out[[1]]))

  adjmatrix <- dist_mat <= diff_eps
  ig <- igraph::graph_from_adjacency_matrix(adjmatrix, mode = "undirected")

  groupings <- igraph::components(ig, mode = "strong")$membership

  if (!is.null(naming_info)){
    if (inherits(g_names, "data.frame")){
      return(cbind(g_names, groupings = groupings))
    } else {
      return(data.frame(g_names, groupings = groupings))
    }
  } else {
    return(data.frame(groupings))
  }


}

testthat::test_that("test mode_clustering_1d", {
  set.seed(1)
  data <- matrix(c(rnorm(mean = 2,n = 100),
                   rnorm(mean = -2, n = 100)),
                 byrow = T, nrow = 100)

  sigma <- quantile(as.matrix(dist(data)), .2)

  data_df <- data.frame(data) %>% rownames_to_column()
  names(data_df)[2:3] <- c("x","y")

  out <- mode_clustering_1d(data_df, position = 2:3, naming_info = 1,
                     sigma = sigma, verbose = F)

  testthat::expect_equal(out[,1], data_df$rowname)

  # 2 completely disjoint groups 1:50 and 51:100
  testthat::expect_equal(length(unique(out$groupings[1:50])), 1)
  testthat::expect_equal(length(unique(out$groupings[51:100])), 1)
  testthat::expect_equal(length(unique(out$groupings)), 2)

})


#' Calculates the minimum coverage radius for euclidean objects as we increase
#' the number of observations.
#'
#' @param ordered_sim_df data frame of points as rows, ordered by when they
#' should be added in the sequence
#' @param e_cols character vector of column names that define
#' locational information of points
#' @param verbose  boolean, if this progression should be verbose
#'
#' @return list of min_cover_vec (minimum at that step) and dist_mat (cumulative
#' max per path at time step t (column))
coverage_down_1d_single <- function(ordered_sim_df, e_cols,
                                    verbose = FALSE){
  n_obs <- nrow(ordered_sim_df)

  point_dist_mat <- as.matrix(dist(ordered_sim_df[,e_cols]))

  dist_mat <- matrix(NA, nrow = n_obs, ncol = n_obs)
  dist_mat[1,1] <- 0

  min_cover_vec <- rep(NA, n_obs)
  min_cover_vec[1] <- 0

  if (n_obs > 1){
    if (verbose){
      pb <- progress::progress_bar$new(
        format = "calculating minimum coverage [:bar] :percent eta: :eta",
        total = n_obs - 1, clear = FALSE, width = 60)
    }

    inner_min_dist_cov_mat <- rep(Inf, n_obs)
    for (r_idx in 2:n_obs){
      #browser()

      # check distances to new point
      new_distances_with_new_point <- point_dist_mat[1:(r_idx-1),r_idx]
      inner_min_dist_cov_mat[r_idx] <- min(new_distances_with_new_point)

      # update minimum covereage required for older points
      inner_min_dist_cov_mat[1:(r_idx-1)] <- apply(
        matrix(c(inner_min_dist_cov_mat[1:(r_idx-1)],
                  new_distances_with_new_point),
                byrow = TRUE, nrow = 2),
        2, min)

      # calc current min_cover value
      min_cover_vec[r_idx] <- max(inner_min_dist_cov_mat[1:r_idx])

      # update radi values if necessary
      dist_mat[1:(r_idx-1), r_idx] <- sapply(dist_mat[1:(r_idx-1), r_idx-1],
                                             function(x) max(c(x, min_cover_vec[r_idx])))
      dist_mat[r_idx, r_idx] <- min_cover_vec[r_idx]


      if (verbose) {
        pb$tick()
      }
    }


  } else {
    message("A mode cluster has only a single element.")
  }
  return(list(min_cover_vec = min_cover_vec, dist_mat = dist_mat))
}


testthat::test_that("test coverage_down_1d_single", {
  ordered_sim_df <- data.frame(x = c(1,2,1.5,4,4.1))
  e_cols = c("x")

  # useful ordered dist:
  # dist(ordered_sim_df)
  #1   2   3   4
  #2 1.0
  #3 0.5 0.5
  #4 3.0 2.0 2.5
  #5 3.1 2.1 2.6 0.1

  out <- coverage_down_1d_single(ordered_sim_df, e_cols)
  testthat::expect_equal(out$min_cover_vec, c(0,1,.5,2,.5))
  testthat::expect_equal(out$dist_mat, matrix(c(0,  1,   1,  2, 2,
                                                NA, 1,   1,  2, 2,
                                                NA, NA, .5,  2, 2,
                                                NA, NA, NA,  2, 2,
                                                NA, NA, NA, NA, .5),
                                              byrow = T, nrow = 5))

  # data fram with single element gets message about it
  testthat::expect_message(out_base <- coverage_down_1d_single(ordered_sim_df[1, , drop = F], e_cols))
  testthat::expect_equal(out_base$min_cover_vec,0)
  testthat::expect_equal(out_base$dist_mat,matrix(0, nrow = 1, ncol = 1))

})

#' Calculates the minimum coverage radius for euclidean objects as we increase
#' the number of observations relative to mode groupings.
#'
#' @param sim_df data frame of points as rows
#' @param e_cols character vector of column names that define
#' locational information of points
#' @param ordering_list list of indices, each vector is related to a single mode
#' cluster and the ordering defines how we build up the region.
#' @param verbose boolean, if this progression should be verbose
#'
#' @return list per mode grouping in \code{ordering_list} with inner lists of
#' min_cover_vec (minimum at that step) and dist_mat (cumulative
#' max per path at time step t (column))
#'
#' @export
coverage_down_1d_mult <- function(sim_df, e_cols,
                             ordering_list, # list of lists...
                             verbose = FALSE){
  n_lists <- length(ordering_list)
  n <- nrow(sim_df)

  assertthat::assert_that((sum(sapply(ordering_list, length)) == n) &
                            (all(sort(unlist(ordering_list)) == 1:n)),
                          msg = "sim_df length and length of ordering_list info should match")
  if (verbose){
    pb <- progress::progress_bar$new(
      format = "calculating minimum coverage (per mode) [:bar] :percent eta: :eta",
      total = n_lists, clear = FALSE, width = 60)
  }

  coverd_info <- list()
  t <- 1
  for (order_vec in ordering_list){

    if ((ncol(sim_df) == 1) | length(order_vec) == 1){
      coverd_info[[t]] <- coverage_down_1d_single(sim_df[order_vec,, drop = F],
                                                  e_cols = e_cols,
                                                  verbose = FALSE)
    } else {
      coverd_info[[t]] <- coverage_down_1d_single(sim_df[order_vec,],
                                                  e_cols = e_cols,
                                                  verbose = FALSE)
    }


    if (verbose){
      pb$tick()
    }
    t <- t + 1
  }

  return(coverd_info)

}

testthat::test_that("test coverage_down_1d_mult", {
  # standard
  sim_df <- data.frame(x = c(1,2,1.5,4,4.1, # first group
                             2.1, 1.1,4.1 , 1.6, 4.2 ) # second group
                       )
  e_cols = c("x")
  ordering_list <- list(c(1:5),
                        c(2,1,4,3) + 5, # allows for basically the same out as first group...
                        c(10)) # for a singleton

  # this one shouldn't expect message - not actually used...
  out <- coverage_down_1d_mult(sim_df[1:9,, drop = F], e_cols, ordering_list[1:2])

  testthat::expect_message(out <- coverage_down_1d_mult(sim_df, e_cols, ordering_list))

  # first group:
  # useful ordered dist
  # dist(ordered_sim_df)
  #1   2   3   4
  #2 1.0
  #3 0.5 0.5
  #4 3.0 2.0 2.5
  #5 3.1 2.1 2.6 0.1

  testthat::expect_equal(out[[1]]$min_cover_vec, c(0,1,.5,2,.5))
  testthat::expect_equal(out[[1]]$dist_mat, matrix(c(0,  1,  1,  2, 2,
                                                    NA, 1,   1,  2, 2,
                                                    NA, NA, .5,  2, 2,
                                                    NA, NA, NA,  2, 2,
                                                    NA, NA, NA, NA, .5),
                                                  byrow = T, nrow = 5))
  # second group:
  testthat::expect_equal(out[[2]]$min_cover_vec, c(0,1,.5,2))
  testthat::expect_equal(out[[2]]$dist_mat, matrix(c(0,  1,   1,  2,
                                                     NA, 1,   1,  2,
                                                     NA, NA, .5,  2,
                                                     NA, NA, NA,  2),
                                                   byrow = T, nrow = 4))

  # singleton group:
  testthat::expect_equal(out[[3]]$min_cover_vec, c(0))
  testthat::expect_equal(out[[3]]$dist_mat, matrix(c(0),nrow = 1))
})

#' calculate distance between two different df's rows
#'
#' @param df1 data frame  (n, p)
#' @param df2 data frame  (m, p)
#'
#' @return matrix of distances (n, m)
my_pdist <- function(df1, df2){
  nrow1 <- nrow(df1)
  nrow2 <- nrow(df2)

  out <- matrix(NA, nrow = nrow1, ncol = nrow2)

  for (r_idx in 1:nrow1){
    for (c_idx in 1:nrow2){
      out[r_idx, c_idx] <- sqrt(sum((df1[r_idx,] - df2[c_idx,])^2))
    }
  }
  return(out)
}

testthat::test_that("test my_pdist", {
  # fixed 1d
  df1 <- data.frame(x = 3:4)
  df2 <- data.frame(x = 4:6)

  dist_out <- my_pdist(df1, df2)

  expected_dist_out <- matrix(c(1,2,3,
                                0,1,2),
                              byrow = T, nrow = 2)

  testthat::expect_equal(dist_out, expected_dist_out)

  # fixed 2d
  df1 <- data.frame(x = 3:4,
                    y = 3:4)
  df2 <- data.frame(x = 4:6,
                    x = 4:6)

  dist_out <- my_pdist(df1, df2)

  expected_dist_out <- matrix(c(1,2,3,
                                0,1,2) * sqrt(2),
                              byrow = T, nrow = 2)

  testthat::expect_equal(dist_out, expected_dist_out)


  # random square
  set.seed(5)
  df_r <- data.frame(x = rnorm(5),
                     y = rnorm(5))
  dist_sq <- my_pdist(df_r,df_r)
  testthat::expect_equivalent(dist_sq, as.matrix(dist(df_r)))
})


#' Helper function for simulation-based conformal score calculation based
#' potentially based on mode clustering and changing radius values.
#'
#' @param df_row_group data frame with new point information to be assessed
#' @param simulations_group_df data frame with multiplee simulated points
#' @param data_column_names character vector of column names that define
#' locational information of points
#' @param simulation_info_df a dataframe with information of each simulation's
#' psuedo-density estimate, mode clustering, ranking of psuedo-density (within
#' cluster and overall). We expect a \code{psuedo_density}, \code{groupings}
#' \code{ranking} and  \code{group_ranking} columns. See
#' \code{EpiCompare::inner_expanding_info} for expected structure.
#' @param list_radius_info list of lists of radius information per each mode.
#' @param list_grouping_id list of vectors of indices (grouped by mode cluster)
#' and ordered by psuedo-density values. These values are relative to the row
#' indices in simulation_info_df
#' @param verbose boolean, if progress should be tracked with a progress bar
#'
#' @return data.frame with a row for new point with a
#' section column \code{containment_val} which is the discrete simulation-based
#' conformal score. If \code{df_row_group} only contained location columns in
#' \code{data_column_names} then a new column named \code{index} is also
#' included, and is related to the row numbers/names of \code{df_row_group}.
#' @export
inner_containment_conformal_score_mode_radius_1d <- function(df_row_group,
                                                             simulations_group_df,
                                                             data_column_names, # must be string vector (not index)
                                                             simulation_info_df,
                                                             list_radius_info,
                                                             list_grouping_id,
                                                             verbose = FALSE){

  if (ncol(df_row_group) == length(data_column_names)){
    df_row_group <- df_row_group %>% tibble::rownames_to_column(var = "row_index")
  }
  group_names <- names(df_row_group)[!(names(df_row_group) %in% data_column_names)]

  if (length(group_names) == 1){
    group_info_df <- df_row_group[,group_names, drop = F]
  } else {
    group_info_df <- df_row_group[,group_names]
  }


  n_draws <- nrow(simulation_info_df)
  n_groups <- length(list_grouping_id)
  n_check <- nrow(df_row_group)

  if (verbose){
    pb <- progress::progress_bar$new(
      format = "Processing [:bar] :percent eta: :eta",
      total = n_draws,
      clear = FALSE, width = 60)
  }



  # https://stackoverflow.com/questions/35106567/how-to-calculate-euclidean-distance-between-two-matrices-in-r
  if (length(data_column_names) == 1){

    between_dist_mat <- my_pdist(df_row_group[, data_column_names, drop = F],
                                              simulations_group_df[, data_column_names, drop = F])
    # between_dist_mat <- matrix(NA, nrow = nrow(df_row_group),
    #                            ncol = nrow(simulations_group_df))
    #
    # for (c_idx in 1:nrow(simulations_group_df)){
    #   between_dist_mat[,c_idx] <- inner_euclidean_distance_1d(df_row_group[, data_column_names],
    #                                                        simulations_group_df[c_idx, data_column_names])
    # }
  } else {
    between_dist_mat <- my_pdist(df_row_group[, data_column_names],
                                  simulations_group_df[, data_column_names])

  }




  overall_info <- list()
  for (g_idx in 1:n_groups){
    inner_ids <- list_grouping_id[[g_idx]]
    inner_rad_dist_mat <- list_radius_info[[g_idx]]$dist_mat

    n_sims_inner <- length(inner_ids)

    # 1. iteratively build distance across each point, store

    # Only examine until covered
    # track first time covered
    not_covered <-  rep(TRUE, n_check)
    coverage_number <- rep(n_draws+1, n_check)

    s_idx <- 1
    while (s_idx <= n_sims_inner & sum(not_covered) > 0){
      if (sum(not_covered) == 1 | s_idx == 1){
        new_contained <- sweep(between_dist_mat[not_covered, inner_ids[1:s_idx], drop = F],
                               MARGIN = 2,
                               STATS = inner_rad_dist_mat[1:s_idx, s_idx],
                               FUN = "<=") %>%
          rowSums() %>% sapply(function(x) x > 0)
      } else {
        new_contained <- sweep(between_dist_mat[not_covered, inner_ids[1:s_idx]],
                               MARGIN = 2,
                               STATS = inner_rad_dist_mat[1:s_idx, s_idx],
                               FUN = "<=") %>%
          rowSums() %>% sapply(function(x) x > 0)
      }

      if (sum(new_contained) > 0){
        # update coverage numbers for those just covered
        inner_coverage <- rep(n_draws+1, sum(not_covered))
        inner_coverage[new_contained] <- s_idx
        coverage_number[not_covered] <- inner_coverage

        # update those covered to not look at again
        not_covered[not_covered] <- !new_contained
      }
      s_idx <- s_idx + 1
    }

    # calculate score per mode relative to overall ranking of points
    all_info_df <- group_info_df %>%
      mutate(containment_val = coverage_number)


    simulation_info_inner <- rbind((simulation_info_df[inner_ids,] %>%
                                      filter(groupings == g_idx) %>%
                                      dplyr::ungroup() %>%
                                      dplyr::select(.data$ranking, .data$group_ranking)),
                                   data.frame(ranking = n_draws+1, group_ranking = n_draws+1))


    overall_score <- all_info_df %>%
      dplyr::left_join(simulation_info_inner,
                       by = c("containment_val" = "group_ranking"))

    overall_info[[g_idx]] <- overall_score


  }
  # combining across modes
  moverall_info <- do.call(rbind, overall_info) %>%
    dplyr::group_by(dplyr::across(tidyselect::one_of(group_names))) %>%
    dplyr::summarize(containment_val = n_draws+1 - min(.data$ranking), .groups = "keep") # why subtraction?:

  return(moverall_info)

}

testthat::test_that("test inner_containment_conformal_score_mode_radius_1d", {
  sim_df <- data.frame(x = c(1,2,1.5,4,4.1, # first group
                             2.1, 1.1,4.1 , 1.6, 4.2 ) # second group
  )
  e_cols = c("x")
  list_grouping_id <- list(c(1:5),
                        c(2,1,4,3) + 5, # allows for basically the same out as first group...
                        c(10))

  # natural structure for simulation_info_df
  psuedo_density_df <- data.frame(psuedo_density = c(.4,.3,.2,.1,.05,
                                                     .31,.41,.11,.21,
                                                     .5))
  mode_grouping <- data.frame(groupings = c(rep(1,5),
                                            rep(2,4),
                                            3))

  suppressMessages(list_radius_info <- coverage_down_1d_mult(sim_df, e_cols = e_cols,
                        ordering_list = list_grouping_id,
                        verbose = FALSE))


  # could be fed into EpiCompare::inner_expanding_info
  # (but doesn't have grouping names...)
  simulation_info_df <- psuedo_density_df %>%
    cbind(mode_grouping) %>%
    dplyr::mutate(ranking = rank(-.data$psuedo_density,
                                 ties.method ="random")) %>%  # no impact in example
    dplyr::group_by(.data$groupings) %>%
    dplyr::mutate(group_ranking = rank(-.data$psuedo_density,
                                       ties.method ="random")) # no impact in example


  test_df <- data.frame(x = c(4.1, 4.11, 0, 3))

  # 4.1 and 4.11 will be contained with 4.1 in either the 1st or 2nd list
  # both of the sim 4.1s are at the bottom of the pile in terms of density,
  # with overall rankings 8 and 9 out of 10

  # 0 will be contained by 1 (first list after we include 2), so it should have
  # ranking associated with 2 (in the first list), which has a ranking of 5

  # 3 will also be captured by 2.1 from the second list as it has a higher
  # density than the 2 from the first list and same radius, associated with
  # a rank of 4

  # since the conformal scores are n_sims + 1 - rank, we should expect:
  out_expected <- data.frame(row_index = as.character(1:4),
                             containment_val = 10+1 - c(8, 8,5,4))


  out <- inner_containment_conformal_score_mode_radius_1d(df_row_group = test_df,
                                                   simulations_group_df = sim_df,
                                                   data_column_names = e_cols,
                                                   simulation_info_df = simulation_info_df,
                                                   list_radius_info = list_radius_info,
                                                   list_grouping_id = list_grouping_id,
                                                   verbose = FALSE)
  out_raw <- out %>% ungroup %>% as.data.frame()
  testthat::expect_equal(out_raw, out_expected)

  # with their own names!
  test_df2 <- data.frame(x = c(4.1,4.11, 0, 3),
                         names = c("alice", "bill", "chuck", "dave"))

  out_expected2 <- data.frame(names = c("alice", "bill", "chuck", "dave"),
                             containment_val = 10+1 - c(8, 8,5,4))
  out2 <- inner_containment_conformal_score_mode_radius_1d(df_row_group = test_df2,
                                                          simulations_group_df = sim_df,
                                                          data_column_names = e_cols,
                                                          simulation_info_df = simulation_info_df,
                                                          list_radius_info = list_radius_info,
                                                          list_grouping_id = list_grouping_id,
                                                          verbose = FALSE)


  out2_raw <- out2 %>% ungroup %>% as.data.frame()
  testthat::expect_equal(out2_raw, out_expected2)


  # 2d:

  sim_df_2d <- data.frame(x = c(1,2,1.5,4,4.1, # first group
                             2.1, 1.1,4.1 , 1.6, 4.2 ), # second group
                       y = c(1,2,1.5,4,4.1, # first group
                             2.1, 1.1,4.1 , 1.6, 4.2 ) # second group
  )
  test_df_2d <- data.frame(x = c(4.1, 4.11, 0, 3),
                           y = c(4.1, 4.11, 0, 3))
  e_cols = c("x", "y")
  out_2d <- inner_containment_conformal_score_mode_radius_1d(df_row_group = test_df_2d,
                                                          simulations_group_df = sim_df_2d,
                                                          data_column_names = e_cols,
                                                          simulation_info_df = simulation_info_df,
                                                          list_radius_info = list_radius_info,
                                                          list_grouping_id = list_grouping_id,
                                                          verbose = FALSE)

})


#' global wrapper for simulation-based conformal score calculation that allows
#' for mode clustering and changes of radius.
#'
#' @param truth_df
#' @param simulations_df
#' @param smooth_functions
#' @param data_column_names
#' @param .maxT
#' @param .sigma_string
#' @param .diff_eps
#' @param .eps
#' @param return_min
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
simulation_based_conformal_1d_complex <- function(truth_df, simulations_df,
                                                  smooth_functions = list(), #named list
                                                  data_column_names = c("y"),
                                                  .maxT = 50,
                                                  .sigma_string = "25%",
                                                  .diff_eps = 1e-06,
                                                  .eps = 1e-06,
                                                  return_min = TRUE,
                                                  verbose = FALSE){

  # calculating sigma value ------
  dist_mat <- simulations_df %>% select(one_of(data_column_names)) %>%
    dist() %>% as.matrix()

  group_names <- names(simulations_df)[!(names(simulations_df) %in% data_column_names)]
  if (length(group_names) == 0){
    simulations_df <- simulations_df %>% rownames_to_column(var = "row_index")
    group_names <- names(simulations_df)[!(names(simulations_df) %in% data_column_names)]
  }

  # sigma selection
  sigma_lower <- EpiCompare:::check_character_percent(.sigma_string, ".sigma_string")
  sigma_sizes <- sapply(sigma_lower + .05*(0:5), function(v) min(v, 1))

  percentage_inner <- sigma_sizes[stats::quantile(as.matrix(dist_mat), sigma_sizes) > 0][1]

  sigma_val <- stats::quantile(as.matrix(dist_mat), percentage_inner)

  tdm <- EpiCompare::tidy_dist_mat(dist_mat,
                                   rownames_df = simulations_df[, group_names, drop = F],
                                   colnames_df = simulations_df[, group_names, drop = F])
  # rank_df
  pseudo_density_df <-  EpiCompare::distance_psuedo_density_function(
    tdm,
    sigma = sigma_val, df_out = T) %>%
    dplyr::mutate(ranking = rank(-.data$psuedo_density,ties.method = "random")) #spelling error... :(, no ties! need ordering


  assertthat::assert_that(all(!is.na(pseudo_density_df$psuedo_density)),
                          msg = paste("internal error in",
                                      "distance_psuedo_density_function",
                                      "function's sigma selection."))

  top_points <- simulations_df %>%
    left_join(pseudo_density_df, by = group_names) %>%
    mutate(keep = ranking > ceiling(.2*nrow(simulations_df))) %>%
    ungroup() %>% filter(keep) %>%
    select(one_of(data_column_names))


  mm_delta <- EpiCompare::get_delta_nn(top_points)

  # mode clustering ------------

  out_groups <- mode_clustering_1d(df_info = simulations_df,
                                   position = c(1:ncol(simulations_df))[names(simulations_df) %in% data_column_names],
                                   naming_info = c(1:ncol(simulations_df))[names(simulations_df) %in% group_names],
                                   sigma, maxT = .maxT,
                                   eps = .eps,
                                   diff_eps = .diff_eps,
                                   verbose = verbose)


  simulation_info_out <- EpiCompare::inner_expanding_info(pseudo_density_df, out_groups)
  simulation_info_df <- simulation_info_out[[1]]
  ordering_list <- simulation_info_out[[2]]

  ## fixed
  tm_radius_fixed <- inner_convert_single_radius_to_structure(mm_delta,
                                                              ordering_list)

  ## vary
  tm_radius_vary <- coverage_down_1d_mult(sim_df = simulations_df,
                                          e_cols = data_column_names,
                                          ordering_list = ordering_list,
                                          verbose = FALSE)
  # no mode clustering ---------
  out_groups_nm <- out_groups %>% dplyr::mutate(groupings = 1)
  simulation_info_out_nm <- inner_expanding_info(pseudo_density_df, out_groups_nm)
  simulation_info_df_nm <- simulation_info_out_nm[[1]]
  ordering_list_nm <- simulation_info_out_nm[[2]]

  tm_radius_fixed_nm <- inner_convert_single_radius_to_structure(mm_delta,
                                                                 ordering_list_nm)
  ## vary
  tm_radius_vary_nm <- coverage_down_1d_mult(sim_df = simulations_df,
                                        e_cols = data_column_names,
                                        ordering_list = ordering_list_nm,
                                        verbose = FALSE)

  # smoothing --------


  tm_radius_vary_update_list <- list()
  tm_radius_vary_nm_update_list <- list()
  if (length(smooth_functions) > 0 ){
    for (f_name in names(smooth_functions)){
      tm_radius_vary_update_list[[f_name]] <- update_tm_smooth(tm_radius_vary,
                                                               func = smooth_functions[[f_name]])

      tm_radius_vary_nm_update_list[[f_name]] <- update_tm_smooth(tm_radius_vary_nm,
                                                                  func = smooth_functions[[f_name]])
    }
  }


  conformal_df_fixed <- inner_containment_conformal_score_mode_radius_1d(
    df_row_group = truth_df,
    simulations_group_df = simulations_df,
    data_column_names = data_column_names,
    simulation_info_df = simulation_info_df, # diff
    list_radius_info = tm_radius_fixed, # diff
    list_grouping_id = ordering_list, # diff
    verbose = verbose)

  conformal_df_fixed_nm <- inner_containment_conformal_score_mode_radius_1d(
    df_row_group = truth_df,
    simulations_group_df = simulations_df,
    data_column_names = data_column_names,
    simulation_info_df = simulation_info_df_nm, # diff
    list_radius_info = tm_radius_fixed_nm, # diff
    list_grouping_id = ordering_list_nm, # diff
    verbose = verbose)

  conformal_df_vary <- inner_containment_conformal_score_mode_radius_1d(
    df_row_group = truth_df,
    simulations_group_df = simulations_df,
    data_column_names = data_column_names,
    simulation_info_df = simulation_info_df, # diff
    list_radius_info = tm_radius_vary, # diff
    list_grouping_id = ordering_list, # diff
    verbose = verbose)

  conformal_df_vary_nm <- inner_containment_conformal_score_mode_radius_1d(
    df_row_group = truth_df,
    simulations_group_df = simulations_df,
    data_column_names = data_column_names,
    simulation_info_df = simulation_info_df_nm, # diff
    list_radius_info = tm_radius_vary_nm, # diff
    list_grouping_id = ordering_list_nm, # diff
    verbose = verbose)


  conformal_df_vary_nm_smooth_list <- list()
  conformal_df_vary_smooth_list <- list()
  if (length(smooth_functions) > 0 ){
    for (f_name in names(smooth_functions)){
      conformal_df_vary_nm_smooth_list[[f_name]] <- inner_containment_conformal_score_mode_radius_1d(
        df_row_group = truth_df,
        simulations_group_df = simulations_df,
        data_column_names = data_column_names,
        simulation_info_df = simulation_info_df_nm, # diff
        list_radius_info = tm_radius_vary_nm_update_list[[f_name]], # diff
        list_grouping_id = ordering_list_nm, # diff
        verbose = verbose)

      conformal_df_vary_smooth_list[[f_name]] <- inner_containment_conformal_score_mode_radius_1d(
        df_row_group = truth_df,
        simulations_group_df = simulations_df,
        data_column_names = data_column_names,
        simulation_info_df = simulation_info_df, # diff
        list_radius_info = tm_radius_vary_update_list[[f_name]], # diff
        list_grouping_id = ordering_list, # diff
        verbose = verbose)
    }
  }

  if (return_min){
    return(list(conformal_score = list(fixed = conformal_df_fixed,
                                       fixed_nm = conformal_df_fixed_nm,
                                       vary = conformal_df_vary,
                                       vary_nm = conformal_df_vary_nm,
                                       smooth_vary = conformal_df_vary_smooth_list,
                                       smooth_vary_nm = conformal_df_vary_nm_smooth_list),
                parameters = c("mm_delta_prop" = .2,
                               "mm_delta" = mm_delta,
                               "sigma_percentage" = percentage_inner)))
  } else{
    return(list(conformal_score = list(fixed = conformal_df_fixed,
                                       fixed_nm = conformal_df_fixed_nm,
                                       vary = conformal_df_vary,
                                       vary_nm = conformal_df_vary_nm,
                                       smooth_vary = conformal_df_vary_smooth_list,
                                       smooth_vary_nm = conformal_df_vary_nm_smooth_list),
                containment_df = simulation_info_df,
                mm_delta = mm_delta,
                tm_radius = list(fixed = tm_radius_fixed,
                                 fixed_nm = tm_radius_fixed_nm,
                                 vary = tm_radius_vary,
                                 vary_nm = tm_radius_vary_nm,
                                 smooth_vary = tm_radius_vary_update_list,
                                 smooth_vary_nm = tm_radius_vary_nm_update_list),
                truth_df_inner = truth_df,
                simulations_group_df_inner = simulations_group_df,
                ordering_list = ordering_list,
                parameters = c("mm_delta_prop" = .2,
                               "sigma_percentage" = percentage_inner)))
  }
}



testthat::test_that("test simulation_based_conformal_1d_complex", {
  create_smooth_function <- function(df){

    eval(parse(text = paste0("smooth_function <- function(x,y){
    tryCatch({
      inner_ss <- smooth.spline(x,y, df = ",
                             df,
                             ")
      return(predict(inner_ss,x)$y)},
      error = function(cond){
        message(sprintf(\"returning y: error in size of x (%i < 4)\", length(x)))
        return(y)
      },
      warning = function(cond){
        message(sprintf(paste(\"returning y: error in size of x\",
                              \"relative to size of df (%i < %i)\"),
                        length(x), ", df, "))
        return(y)
      }
    )
  }")))

    return(smooth_function)
  }

  inner_moving <- function(x, y, window = 5, fill_left = T,
                           fun = min){
    n <- length(y)
    if (n < window){
      message(sprintf(paste("returning y: error in size of x",
                            "relative to size of window (%i < %i)"),
                      n, window))
      return(y)
    }

    y_smooth <- rep(NA, n)
    for (idx in window:n){
      y_smooth[idx] <- fun(y[(idx-window+1):idx])
    }

    if (fill_left){
      y_smooth[1:(window-1)] <- y_smooth[window]
    }

    return(y_smooth)
  }

  create_min_smooth_function <- function(window){
    eval(parse(text = paste0("out_function <- function(x,y){
    inner_moving(x,y, window = ", window,", fill_left = T, fun = min)
  }")))
    return(out_function)
  }

  smooth_functions <- list("window10" = create_min_smooth_function(10))

  set.seed(50)
  sim_df <- data.frame(matrix(c(rnorm(100, mean = -2), rnorm(100, mean = 2)),
                   byrow = TRUE, ncol = 2)) %>%
    rename(x = "X1", y = "X2")
  test_df <- data.frame(x = c(2,-1, -2.5),
                        y = c(2,-3.5, 2.5))

  if (FALSE){
    sim_df %>% ggplot() +
      geom_point(aes(x=x,y=y)) +
      geom_point(data = test_df, aes(x=x,y=y),color ="red")
  }
  # we expect that the first will have a higher conformal score than the second,
  # and the third will have a conformal score of 0
  out <- simulation_based_conformal_1d_complex(truth_df = test_df,
                                        simulations_df = sim_df,
                                        smooth_functions = smooth_functions, #named list
                                        data_column_names = c("x","y"),
                                        .maxT = 100,
                                        .eps = 1e-07,
                                        .diff_eps = 1e-06, return_min = F)

  #out$conformal_score

  testthat::expect_true(all(out$conformal_score$vary >=
                              out$conformal_score$smooth_vary$window10))
  testthat::expect_true(all(out$conformal_score$vary_nm >=
                              out$conformal_score$smooth_vary_nm$window10))

  list_cs <- out$conformal_score[c("fixed", "fixed_nm", "vary", "vary_nm")]
  list_cs[[5]] <- out$conformal_score$smooth_vary$window10
  list_cs[[6]]<- out$conformal_score$smooth_vary_nm$window10

  cs_mat_idx <- 1
  for (cs_mat in list_cs){
    testthat::expect_true(all(diff(cs_mat$containment_val) <= rep(0,2)))
    if (!(cs_mat_idx %in% 3:4)){
      testthat::expect_equal(cs_mat$containment_val[3], 0)
    }
    testthat::expect_gt(cs_mat$containment_val[1], 60) # should be higher - pretty central...

    cs_mat_idx <- cs_mat_idx + 1
  }


})


#' @param conformal_score_cut note this can be an integer, string percentage or
# fraction, but relates to the conformal cut (not the confidence level).
#'
# simulation_based_conformal_1d_region <- function(simulations_grouped_df,
#                                                  data_column_names = c("y"),
#                                                  conformal_score_cut = .9,
#                                                  delta_prop = .8){
#
#   if(!EpiCompare::is.wholenumber(conformal_score_cut)){
#     if(is.character(conformal_score_cut)){
#       inner_percent <- EpiCompare::check_character_percent(conformal_score_cut)
#     } else {
#       inner_percent <- conformal_score_cut
#     }
#     conformal_score_cut <- ceiling(inner_percent * nrow(simulations_grouped_df))
#   }
#
#   assertthat::assert_that(inherits(simulations_grouped_df, "grouped_df"))
#   group_names <- names(group_keys(simulations_grouped_df))
#
#
#   simulations_group_df_inner <- simulations_grouped_df %>%
#     dplyr::select(dplyr::one_of(c(group_names, data_column_names)))
#
#   group_info <- simulations_group_df_inner %>% group_keys()
#
#   dist_mat <- simulations_group_df_inner %>% group_split() %>%
#     do.call(rbind, .) %>% # match group_info ordering
#     select(one_of(data_column_names)) %>%
#     dist() %>% as.matrix()
#
#   tdm_sims <- EpiCompare::tidy_dist_mat(dist_mat, group_info, group_info)
#
#   # sigma selection
#
#   sigma_size <- c("20%" = .2, "25%" = .25, "30%" = .3,
#                   "35%" = .35, "40%" = .4, "45%" = .45)
#
#   percentage <- names(sigma_size)[stats::quantile(as.matrix(tdm_sims), sigma_size) > 0][1]
#
#
#   # rank_df
#   pseudo_density_df <- EpiCompare::distance_psuedo_density_function(
#     tdm_sims,
#     sigma = percentage, df_out = T) %>%
#     mutate(ranking = rank(psuedo_density,ties.method = "min")) #spelling error... :(
#
#   assertthat::assert_that(all(!is.na(pseudo_density_df$psuedo_density)),
#                           msg = paste("internal error in",
#                                       "distance_psuedo_density_function",
#                                       "function's sigma selection."))
#
#   proportion_points_not_included <- 1 - delta_prop
#
#   top_points <- simulations_group_df_inner %>%
#     left_join(pseudo_density_df, by = group_names) %>%
#     mutate(keep = ranking > ceiling(proportion_points_not_included*nrow(simulations_group_df_inner))) %>%
#     ungroup() %>% filter(keep) %>%
#     select(one_of(data_column_names))
#
#   mm_delta <- get_delta_nn(top_points)
#
#   conformal_band_points <- simulations_group_df_inner %>%
#     left_join(pseudo_density_df, by = group_names) %>%
#     mutate(keep = ranking > conformal_score_cut) %>%
#     ungroup() %>% filter(keep) %>%
#     select(one_of(data_column_names))
#
#
#   return(list(conformal_band_points = conformal_band_points,
#               mm_delta = mm_delta))
# }




## Conformal Calibration Step ----------------

#
# in this step, we get the distribution of all conformal scores for the
# calibration set.
#

file_exists <- check_file_exists(sprintf("calibration_information_1d_%s_sims_delta_size_%s.Rdata",
                                         n_simulations, delta_prop_string),
                                 dir = ".")

if (!file_exists | rerun){
  library(parallel)
  library(foreach)
  library(doSNOW)
  max_cores <- parallel::detectCores()
  cl <- makeCluster(max_cores-4)
  registerDoSNOW(cl)
  iterations <- nrow(calibration_set)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  n_sims <- n_simulations
  verbose <- T
  number_points <- 50

  calibration_conformal_scores <- foreach(calibration_idx = 1:nrow(calibration_set),
                                          #.combine = c,
                                          .options.snow = opts,
                                          .packages = c("dplyr", "tidyr",
                                                        "EpiCompare",
                                                        "simulationBands")) %dopar%{

        x_inner <- calibration_set$x[calibration_idx]
        obs_inner <- calibration_set[calibration_idx,]

        inner_sim <- data.frame(x = rep(x_inner, n_sims))
        inner_sim$y <- sapply(inner_sim$x, function(x) {simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})
        inner_sim$idx <- paste0("sim",1:n_sims)
        inner_sim <- inner_sim %>% group_by(idx)

        conformal_score <- simulation_based_conformal_1d(truth_df = obs_inner,
                                                         simulations_grouped_df = inner_sim,
                                                         data_column_names = c("y"),
                                                         delta_prop = delta_prop)



        return(conformal_score)
                                          }

  conformal_score_vector <- sapply(calibration_conformal_scores, function(x) x$conformal_score)
  mm_delta_vector <- sapply(calibration_conformal_scores, function(x) x$mm_delta)

  save(calibration_conformal_scores,
       file = sprintf("calibration_information_1d_%s_sims_delta_size_%s.Rdata",
                      n_simulations, delta_prop_string))

  parallel::stopCluster(cl)
} else {
  load(sprintf("calibration_information_1d_%s_sims_delta_size_%s.Rdata",
               n_simulations, delta_prop_string))
  conformal_score_vector <- sapply(calibration_conformal_scores, function(x) x$conformal_score)
  mm_delta_vector <- sapply(calibration_conformal_scores, function(x) x$mm_delta)
}


df_conformal_calibration_fit_info <- data.frame(
  idx = 1:nrow(calibration_set),
  cs = conformal_score_vector,
  mm_delta = mm_delta_vector)


df_conformal_calibration_ecdf_info <- data.frame(
  cs = conformal_score_vector,
  ecdf = ecdf(conformal_score_vector)(conformal_score_vector))

vis1_d <- df_conformal_calibration_fit_info %>%
  ggplot() +
  geom_histogram(aes(x = cs)) +
  labs(x = "simulation conformal score",
       title = "calibration set")

vis2_d <- df_conformal_calibration_ecdf_info %>%
  ggplot() +
  geom_line(aes(x = cs, y = ecdf)) +
  labs(x = "simulation conformal score",
       y = "empirical cumulative distribution function",
       title = "calibration set") +
  ylim(0,1)

vis3_d <- df_conformal_calibration_fit_info %>%
  ggplot() +
  geom_point(aes(y = cs, x = mm_delta), alpha = .1) +
  labs(y = "simulation conformal score",
       x = "delta for delta ball around filament compression points",
       title = "calibration set")

vis4_d <- df_conformal_calibration_fit_info %>%
  ggplot() +
  geom_point(aes(x = idx, y = cs)) +
  xlim(0,20)

vis5_d <- df_conformal_calibration_fit_info %>%
  ggplot() +
  geom_point(aes(x = idx, y = mm_delta)) +
  xlim(0,20)

gridExtra::grid.arrange(vis1_d, vis2_d,
                        vis3_d, vis4_d,
                        vis5_d,
                        nrow = 3)


if (save_fig){
  ggvis <- gridExtra::arrangeGrob(vis1_d + theme_minimal() +
                                    labs(title = "Distribution of \nsimulation conformal score",
                                         x = "simulation-based conformal score (\"radius\")",
                                         y = "# of calibration set"),
                                  vis2_d + theme_minimal() +
                                    labs(x = "simulation-based conformal score (\"radius\")",
                                         title = "Cumulative distribution of \nsimulation conformal score"),
                                  layout_matrix = matrix(c(1,1,1,1,2,2,
                                                           1,1,1,1,2,2,
                                                           1,1,1,1,2,2),
                                                         nrow = 3, byrow = T))

  ggplot2::ggsave(plot = ggvis,
                  width = 12,
                  height = 5,
                  filename = sprintf("quick_images/simulation_bands_1d_delta_size_%s.png",
                                     delta_prop_image_string))

  ggvis2 <- gridExtra::arrangeGrob(vis1_d + theme_minimal() +
                                    labs(title = "Distribution of \\\\simulation conformal score",
                                         x = "simulation-based conformal score ('radius')",
                                         y = "number of obs in calibration set"),
                                  vis2_d + theme_minimal() +
                                    labs(x = "simulation-based conformal score ('radius')",
                                         title = "Cumulative distribution of \nsimulation conformal score"),
                                  layout_matrix = matrix(c(1,1,1,1,2,2,
                                                           1,1,1,1,2,2,
                                                           1,1,1,1,2,2),
                                                         nrow = 3, byrow = T))

  tikz(file = sprintf("quick_images/simulation_bands_1d_delta_size_%s.tex",
                      delta_prop_image_string),
       width = 12,
       height = 5)
  plot(ggvis2)
  dev.off()
}

## Test set ---------------

# just some info:
conformal_score_vector %>% stats::quantile(1-c(.6,.9)) %>% floor()
ecdf(conformal_score_vector)(99) # .925 ()



file_exists <- check_file_exists(sprintf(
  paste0("conformal_test_information_1d_%s",
         "_sims_delta_size_%s_confidence_level_%s.Rdata"),
  n_simulations, delta_prop_string, confidence_level_string),
  dir = ".")

if (!file_exists | rerun){
  library(parallel)
  library(foreach)
  library(doSNOW)
  max_cores <- parallel::detectCores()
  cl <- makeCluster(max_cores-4)
  registerDoSNOW(cl)
  iterations <- nrow(test_set)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  n_sims <- n_simulations
  verbose <- T
  number_points <- 50

  cut <- conformal_score_vector %>% stats::quantile(1-confidence_level) %>%
    floor()

  calibration_conformal_scores <- foreach(calibration_idx = 1:nrow(test_set),
                                          #.combine = c,
                                          .options.snow = opts,
                                          .packages = c("dplyr", "tidyr",
                                                        "simulationBands",
                                                        "EpiCompare")) %dopar%{

      x_inner <- test_set$x[calibration_idx]

      inner_sim <- data.frame(x = rep(x_inner, n_sims))
      inner_sim$y <- sapply(inner_sim$x, function(x) {simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})
      inner_sim$idx <- paste0("sim",1:n_sims)
      inner_sim <- inner_sim %>% group_by(idx)

      conformal_score <- simulation_based_conformal_1d_region(
        simulations_grouped_df = inner_sim,
        data_column_names = c("y"),
        conformal_score_cut = cut)

      return(conformal_score)
      }

  save(calibration_conformal_scores,
       file = sprintf(
         paste0("conformal_test_information_1d_%s",
                "_sims_delta_size_%s_confidence_level_%s.Rdata"),
         n_simulations, delta_prop_string, confidence_level_string))

  parallel::stopCluster(cl)
} else {
  load(sprintf(
    paste0("conformal_test_information_1d_%s",
           "_sims_delta_size_%s_confidence_level_%s.Rdata"),
    n_simulations, delta_prop_string, confidence_level_string))
}


file_exists <- check_file_exists(sprintf(
  paste0("conformal_test_information_discrete_1d_%s",
         "_sims_delta_size_%s.Rdata"), n_simulations, delta_prop_string),
  dir = ".")

if (!file_exists | rerun){
  library(parallel)
  library(foreach)
  library(doSNOW)
  max_cores <- parallel::detectCores()
  cl <- makeCluster(max_cores-4)
  registerDoSNOW(cl)
  iterations <- nrow(test_set_discrete)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  n_sims <- n_simulations
  verbose <- T
  number_points <- 50

  cut <- conformal_score_vector %>% stats::quantile(1 - confidence_level) %>% floor()

  calibration_conformal_scores_discrete <- foreach(calibration_idx = 1:nrow(test_set_discrete),
                                                   #.combine = c,
                                                   .options.snow = opts,
                                                   .packages = c("dplyr", "tidyr",
                                                                 "EpiCompare",
                                                                 "simulationBands")) %dopar%{

     x_inner <- test_set_discrete$x[calibration_idx]

     inner_sim <- data.frame(x = rep(x_inner, n_sims))
     inner_sim$y <- sapply(inner_sim$x, function(x) {simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})
     inner_sim$idx <- paste0("sim", 1:n_sims)
     inner_sim <- inner_sim %>% group_by(idx)

     conformal_score <- simulation_based_conformal_1d_region(
       simulations_grouped_df = inner_sim,
       data_column_names = c("y"),
       conformal_score_cut = cut,
       delta_prop = delta_prop)

     return(conformal_score)
     }

  save(calibration_conformal_scores_discrete,
       file = sprintf(paste0("conformal_test_information_discrete_1d_%s",
                             "_sims_delta_size_%s.Rdata"), n_simulations, delta_prop_string))

  parallel::stopCluster(cl)
} else {
  load(sprintf(paste0("conformal_test_information_discrete_1d_%s",
                      "_sims_delta_size_%s.Rdata"), n_simulations, delta_prop_string))
}


### visualization of above ------------------

delta_y = .001
yy <- seq(-8, 8, by = delta_y)

containment_vec <- function(points, yy, radius){
  if (nrow(points) == 0){
    return(rep(FALSE, length(yy)))
  }
  nn2_vals <- RANN::nn2(data = as.matrix(points),
                        query = as.matrix(yy),
                        k = 1,
                        treetype = "kd")
  contained <- nn2_vals$nn.dists < radius

  return(contained)
}



if (verbose) {
  pb <- progress::progress_bar$new(
    format = "Containment calculation [:bar] :percent eta: :eta",
    total = nrow(test_set_discrete), clear = FALSE, width = 80)
}

df_grid_all <- data.frame()
df_containment_prop_all <- data.frame()


for (calibration_idx in 1:nrow(test_set_discrete)){
  # for visual
  x_inner <- test_set_discrete$x[calibration_idx]
  points <- calibration_conformal_scores_discrete[[calibration_idx]]$conformal_band_points
  radius <- calibration_conformal_scores_discrete[[calibration_idx]]$mm_delta
  contained <- containment_vec(points, yy, radius)

  df_inner <- data.frame(
    calibration_idx = calibration_idx,
    x = x_inner,
    y = yy,
    contained = contained)
  df_grid_all <- rbind(df_grid_all, df_inner)

  if (verbose) {
    pb$tick()
  }

  # for containment of true observations:
  inner_sim_y <- sapply(rep(x_inner, n_sims_containment),
                        function(x) {
                          simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})
  contained_check <- containment_vec(points, inner_sim_y, radius)
  df_containment_prop_inner <- data.frame(
    x = x_inner,
    containment_prop = mean(contained_check))
  df_containment_prop_all <- rbind(df_containment_prop_all,
                                   df_containment_prop_inner)
}

vis1 <- df_grid_all %>% filter(contained) %>%
  ggplot() +
  geom_point(data = test_set, aes(x = x, y = y), color = "blue") +
  geom_tile(aes(x = x, y = y), alpha = .5) +
  labs(title = sprintf(paste("Our Simulation-based \nConformal Inference",
                             "Prediction Regions (%s)"),
                       confidence_level_string)) +
  xlim(-1.5, 1.5) + ylim(-8,8) +
  theme_minimal()

vis2 <- df_containment_prop_all %>%
  ggplot() +
  geom_line(aes(x = x, y = containment_prop)) +
  ylim(0,1) + xlim(-1.5, 1.5)  +
  geom_hline(yintercept = confidence_level, color = "blue") +
  geom_label(data = data.frame(x = 1.5, y = confidence_level,
                               label = confidence_level),
             aes(x = x, y = y, label = label)) +
  labs(y = "containment proportion",
       title = sprintf(paste("containment proportion of %s ys\nrandomly",
                             "sampled per x value"),n_sims_containment)) +
  theme_minimal()


gridExtra::grid.arrange(vis1, vis2,
                        layout_matrix = matrix(c(1,1,1,1,
                                                 1,1,1,1,
                                                 1,1,1,1,
                                                 2,2,2,2,
                                                 2,2,2,2), byrow = T, ncol = 4))

if (save_image){
  ggvis <- gridExtra::arrangeGrob(vis1, vis2,
                                  layout_matrix = matrix(c(1,1,1,1,
                                                           1,1,1,1,
                                                           1,1,1,1,
                                                           2,2,2,2,
                                                           2,2,2,2), byrow = T, ncol = 4))
  ggplot2::ggsave(plot = ggvis, width = 5, height = 7,
                  filename =
                    sprintf(paste0("quick_images/sim_bands_conf",
                                   "_test_discrete_1d_%s",
                                   "_sims_delta_size_%s_confidence_level_%s.png"),
                            n_simulations, delta_prop_image_string,
                            confidence_level_image_string))

  vis1_tikz <- df_grid_all %>% filter(contained) %>%
    ggplot() +
    geom_point(data = test_set, aes(x = x, y = y), color = "blue") +
    geom_tile(aes(x = x, y = y), alpha = .5) +
    labs(title = paste("Simulation-based \\\ Conformal Inference",
                       "Prediction Regions (60\\%)")) +
    xlim(-1.5, 1.5) + ylim(-8,8) +
    theme_minimal()

  ggvis_tikz <- gridExtra::arrangeGrob(vis1_tikz, vis2,
                                  layout_matrix = matrix(c(1,1,1,1,
                                                           1,1,1,1,
                                                           1,1,1,1,
                                                           2,2,2,2,
                                                           2,2,2,2), byrow = T, ncol = 4))

  tikz(file = sprintf(paste0("quick_images/sim_bands_conf",
                             "_test_discrete_1d_%s",
                             "_sims_delta_size_%s_confidence_level_%s.tex"),
                      n_simulations, delta_prop_image_string,
                      confidence_level_image_string))
  plot(ggvis_tikz)
  dev.off()
}

## Discrete CDE based approaches

## "Global" CD approach

The approach here will be similar to what a has done before - in a discrete fashion.



calibration_info <- calibration_set
calibration_info$calibration_cde_values <- sapply(1:nrow(calibration_set),
  function(r_idx){
    simulationBands::cde_lei_wassserman(calibration_set$x[r_idx])(calibration_set$y[r_idx])
  })


cutoff_g_cd <- stats::quantile(calibration_info$calibration_cde_values,
                               1- confidence_level)


delta_y = .001
yy <- seq(-8, 8, by = delta_y)

if (verbose) {
  pb <- progress::progress_bar$new(
    format = "Containment calculation [:bar] :percent eta: :eta",
    total = nrow(test_set_discrete), clear = FALSE, width = 80)
}

df_grid_all_g_cd <- data.frame()
df_containment_prop_all_g_cd <- data.frame()


for (calibration_idx in 1:nrow(test_set_discrete)){
  # for visual
  x_inner <- test_set_discrete$x[calibration_idx]
  cde_values <- simulationBands::cde_lei_wassserman(x_inner)(yy)
  contained <- cde_values > cutoff_g_cd

  df_inner <- data.frame(
    calibration_idx = calibration_idx,
    x = x_inner,
    y = yy,
    contained = contained,
    cde = cde_values)
  df_grid_all_g_cd <- rbind(df_grid_all_g_cd, df_inner)

  if (verbose) {
    pb$tick()
  }

  # for containment of true observations:
  inner_sim_y <- sapply(rep(x_inner, n_sims_containment),
                        function(x) {
                          simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})
  sim_cde_values <- simulationBands::cde_lei_wassserman(x_inner)(inner_sim_y)
  contained_check <- sim_cde_values > cutoff_g_cd
  df_containment_prop_inner <- data.frame(
    x = x_inner,
    containment_prop = mean(contained_check))
  df_containment_prop_all_g_cd <- rbind(df_containment_prop_all_g_cd,
                                        df_containment_prop_inner)
}

vis1g <- df_grid_all_g_cd %>% filter(contained) %>%
  ggplot() +
  geom_point(data = test_set, aes(x = x, y = y), color = "blue") +
  geom_tile(aes(x = x, y = y), alpha = .5) +
  labs(title = sprintf(paste("\"Global\" CD*\nConformal",
                             "Inference Prediction Regions (%s)"),
                       confidence_level_string)) +
  xlim(-1.5,1.5) + ylim(-8,8) +
  theme_minimal()

vis2g <- df_containment_prop_all_g_cd %>%
  ggplot() +
  geom_line(aes(x = x, y = containment_prop)) +
  ylim(0,1) + xlim(-1.5,1.5) +
  geom_hline(yintercept = confidence_level, color = "blue") +
  geom_label(data = data.frame(x = 1.5, y = confidence_level,
                               label = confidence_level),
             aes(x = x, y = y, label = label)) +
  labs(y = "containment proportion",
       title = sprintf(paste("containment proportion of",
                             "%s ys\nrandomly sampled per x value"),
                       n_sims_containment)) +
  theme_minimal()


gridExtra::grid.arrange(vis1g, vis2g,
                        layout_matrix = matrix(c(1,1,1,1,
                                                 1,1,1,1,
                                                 1,1,1,1,
                                                 2,2,2,2,
                                                 2,2,2,2), byrow = T, ncol = 4))


if (save_image){
  ggvis <- gridExtra::arrangeGrob(vis1g, vis2g,
                                  layout_matrix = matrix(c(1,1,1,1,
                                                           1,1,1,1,
                                                           1,1,1,1,
                                                           2,2,2,2,
                                                           2,2,2,2), byrow = T, ncol = 4))
  ggplot2::ggsave(plot = ggvis,
                  filename = sprintf(paste0("quick_images/sim_bands_conf",
                                            "_test_discrete_1d_g_cd_%s",
                                            "_sims_delta_size_%s_confidence_level_%s.png"),
                                     n_simulations, delta_prop_image_string,
                                     confidence_level_image_string))


  vis1g_tikz <- df_grid_all_g_cd %>% filter(contained) %>%
    ggplot() +
    geom_point(data = test_set, aes(x = x, y = y), color = "blue") +
    geom_tile(aes(x = x, y = y), alpha = .5) +
    labs(title = paste('"Global" CD* \\\\\ Conformal',
                               "Inference Prediction Regions (60\\%)")) +
    xlim(-1.5,1.5) + ylim(-8,8) +
    theme_minimal()
  ggvis_tex <- gridExtra::arrangeGrob(vis1g_tikz, vis2g,
                                  layout_matrix = matrix(c(1,1,1,1,
                                                           1,1,1,1,
                                                           1,1,1,1,
                                                           2,2,2,2,
                                                           2,2,2,2), byrow = T, ncol = 4))
  tikz(file = sprintf(paste0("quick_images/sim_bands_conf",
                             "_test_discrete_1d_g_cd_%s",
                             "_sims_delta_size_%s_confidence_level_%s.tex"),
                      n_simulations, delta_prop_image_string,
                      confidence_level_image_string))
  plot(ggvis_tex)
  dev.off()
}

## Lei & Wasserman  --------

# We're going to break X in 8 bins (the same number as Figure 4, Lei2014)

cutoff_lei_wasserman <- calibration_info %>%
  mutate(x_cut = cut(x, breaks = seq(-1.5-2*.Machine$double.eps,
                                     1.5+2*.Machine$double.eps, length.out = 9)),
         x_cut_lower = cut_to_numeric(x_cut,.lower = T),
         x_cut_upper = cut_to_numeric(x_cut,.lower = F)) %>%
  group_by(x_cut, x_cut_lower, x_cut_upper) %>%
  summarize(cutoff_lei_wasserman = quantile(calibration_cde_values,
                                            1 - confidence_level),
            .groups = "drop_last")



delta_y = .001
yy <- seq(-8, 8, by = delta_y)

if (verbose) {
  pb <- progress::progress_bar$new(
    format = "Containment calculation [:bar] :percent eta: :eta",
    total = nrow(test_set_discrete), clear = FALSE, width = 80)
}

df_grid_all_lw_cd <- data.frame()
df_containment_prop_all_lw_cd <- data.frame()


for (calibration_idx in 1:nrow(test_set_discrete)){
  # for visual
  x_inner <- test_set_discrete$x[calibration_idx]
  cde_values <- simulationBands::cde_lei_wassserman(x_inner)(yy)

  group_id <- cut(x_inner, breaks = seq(-1.5-2*.Machine$double.eps,
                                        1.5+2*.Machine$double.eps, length.out = 9))
  inner_cutoff <- cutoff_lei_wasserman$cutoff_lei_wasserman[cutoff_lei_wasserman$x_cut == group_id]
  assertthat::assert_that(length(inner_cutoff) == 1)
  contained <- cde_values > inner_cutoff

  df_inner <- data.frame(
             calibration_idx = calibration_idx,
             x = x_inner,
             y = yy,
             contained = contained)
  df_grid_all_lw_cd <- rbind(df_grid_all_lw_cd, df_inner)

  if (verbose) {
      pb$tick()
  }

  # for containment of true observations:
  inner_sim_y <- sapply(rep(x_inner, n_sims_containment),
    function(x) {
      simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})
  sim_cde_values <- simulationBands::cde_lei_wassserman(x_inner)(inner_sim_y)
  contained_check <- sim_cde_values > inner_cutoff
  df_containment_prop_inner <- data.frame(
    x = x_inner,
    containment_prop = mean(contained_check))
  df_containment_prop_all_lw_cd <- rbind(df_containment_prop_all_lw_cd,
                                   df_containment_prop_inner)
}


vis1lw <- df_grid_all_lw_cd %>% filter(contained) %>%
  ggplot() +
  geom_point(data = test_set, aes(x = x, y = y), color = "blue") +
  geom_tile(aes(x = x, y = y), alpha = .5) +
  labs(title = sprintf(paste("Lei and Wasserman (2014) ~Local~\nConformal",
                             "Inference Prediction Regions (%s)"),
                       confidence_level_string)) +
  xlim(-1.5,1.5) + ylim(-8,8) +
  theme_minimal()

vis2lw <- df_containment_prop_all_lw_cd %>%
  ggplot() +
  geom_line(aes(x = x, y = containment_prop)) +
  ylim(0,1) + xlim(-1.5,1.5) +
  geom_hline(yintercept = confidence_level, color = "blue") +
  geom_label(data = data.frame(x = 1.5, y = confidence_level,
                               label = confidence_level),
             aes(x = x, y = y, label = label)) +
  labs(y = "containment proportion",
       title = sprintf(paste("containment proportion of %s",
                     "ys\nrandomly sampled per x value"),
                     n_sims_containment)) +
  theme_minimal()


gridExtra::grid.arrange(vis1lw, vis2lw,
             layout_matrix = matrix(c(1,1,1,1,
                                      1,1,1,1,
                                      1,1,1,1,
                                      2,2,2,2,
                                      2,2,2,2), byrow = T, ncol = 4))

if (save_image){
  ggvis <- gridExtra::arrangeGrob(vis1lw, vis2lw,
             layout_matrix = matrix(c(1,1,1,1,
                                      1,1,1,1,
                                      1,1,1,1,
                                      2,2,2,2,
                                      2,2,2,2), byrow = T, ncol = 4))
  ggplot2::ggsave(plot = ggvis,
                  filename = sprintf(paste0("quick_images/sim_bands_conf",
                                   "_test_discrete_1d_lw_cd_%s",
                                   "_sims_delta_size_%s_confidence_level_%s.png"),
                            n_simulations, delta_prop_image_string,
                            confidence_level_image_string))

  vis1lw_tikz <- df_grid_all_lw_cd %>% filter(contained) %>%
    ggplot() +
    geom_point(data = test_set, aes(x = x, y = y), color = "blue") +
    geom_tile(aes(x = x, y = y), alpha = .5) +
    labs(title = paste("Lei and Wasserman (2014) $\\sim$ Local $\\sim$ \\\\\ Conformal",
                               "Inference Prediction Regions (60\\%)")) +
    xlim(-1.5,1.5) + ylim(-8,8) +
    theme_minimal()

  ggvis_tikz <- gridExtra::arrangeGrob(vis1lw_tikz, vis2lw,
                                  layout_matrix = matrix(c(1,1,1,1,
                                                           1,1,1,1,
                                                           1,1,1,1,
                                                           2,2,2,2,
                                                           2,2,2,2), byrow = T, ncol = 4))
  tikz(file = sprintf(paste0("quick_images/sim_bands_conf",
                             "_test_discrete_1d_lw_cd_%s",
                             "_sims_delta_size_%s_confidence_level_%s.tex"),
                      n_simulations, delta_prop_image_string,
                      confidence_level_image_string))
  plot(ggvis_tikz)
  dev.off()

}



## Izbicki 2020 ------------------


### calibration

delta_y = .001
yy <- seq(-8, 8, by = delta_y)

n_sims <- n_simulations

if (verbose) {
  pb <- progress::progress_bar$new(
    format = "Containment calculation [:bar] :percent eta: :eta",
    total = nrow(calibration_set), clear = FALSE, width = 80)
}

df_grid_all_g_cd_calibration_set <- data.frame()
df_containment_prop_all_g_cd_calibration_set <- data.frame()


for (calibration_idx in 1:nrow(calibration_set)){
  # for visual
  x_inner <- calibration_set$x[calibration_idx]
  cde_values <- simulationBands::cde_lei_wassserman(x_inner)(yy)
  contained <- cde_values > cutoff_g_cd

  df_inner <- data.frame(
             calibration_idx = calibration_idx,
             x = x_inner,
             y = yy,
             contained = contained,
             cde = cde_values)
  df_grid_all_g_cd_calibration_set <- rbind(df_grid_all_g_cd_calibration_set, df_inner)

  if (verbose) {
      pb$tick()
  }

  # for containment of true observations:
  inner_sim_y <- sapply(rep(x_inner, n_sims_containment),
    function(x) {
      simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})
  sim_cde_values <- simulationBands::cde_lei_wassserman(x_inner)(inner_sim_y)
  contained_check <- sim_cde_values > cutoff_g_cd
  df_containment_prop_inner <- data.frame(
    x = x_inner,
    containment_prop = mean(contained_check))
  df_containment_prop_all_g_cd_calibration_set <-
    rbind(df_containment_prop_all_g_cd_calibration_set,
                                   df_containment_prop_inner)
}

#yy
t_range <- range(df_grid_all_g_cd_calibration_set$cde)
n_levels <- 1000
t_grid <- seq(0, t_range[2], length.out = n_levels)

n_calibrate <- length(unique(df_grid_all_g_cd_calibration_set$calibration_idx))
n_y <- length(yy)
n_t <- length(t_grid)
k <- max(round(n_calibrate/100),1)

g_calibrate_cde <- matrix(NA, n_calibrate,
                      n_t)

if (verbose) {
  pb <- progress::progress_bar$new(
    format = "Make Profile Matrix (Calibration) [:bar] :percent eta: :eta",
    total = n_calibrate, clear = FALSE, width = 80)
}

for (calibration_idx in 1:n_calibrate){
  cde_values <- df_grid_all_g_cd_calibration_set$cde[df_grid_all_g_cd_calibration_set$calibration_idx == calibration_idx]
  g_calibrate_cde[calibration_idx,] <- predictionBands:::profile_density(t_grid, yy,
                                      cde_values)

  if (verbose) {
      pb$tick()
  }
}

kmeans_result <- try(LICORS::kmeanspp(g_calibrate_cde,k=k),
                     silent = TRUE)
if(class(kmeans_result)=="try-error"){
  kmeans_result <- kmeans(g_calibrate_cde,centers = k)
}
centers_kmeans <- kmeans_result$centers

library(FNN)
which_partition_calibration <- RANN::nn2(data = centers_kmeans, g_calibrate_cde, k = 1)$nn.idx
                               #predictionBands:::which_neighbors(centers_kmeans,
                                                                #g_calibrate_cde,1)

calibration_info$partition_id <- which_partition_calibration

rm(g_calibrate_cde)

### test


n_test <- length(unique(df_grid_all_g_cd$calibration_idx))

g_test_cde <- matrix(NA, n_test,
                      n_t)

if (verbose) {
  pb <- progress::progress_bar$new(
    format = "Make Profile Matrix (Test) [:bar] :percent eta: :eta",
    total = n_test, clear = FALSE, width = 80)
}

for (calibration_idx in 1:n_test){
  cde_values <- df_grid_all_g_cd$cde[df_grid_all_g_cd$calibration_idx == calibration_idx]
  g_test_cde[calibration_idx,] <- predictionBands:::profile_density(t_grid, yy,
                                      cde_values)

  if (verbose) {
      pb$tick()
  }
}

which_partition_test <- RANN::nn2(data = centers_kmeans, g_test_cde, k = 1)$nn.idx
#predictionBands:::which_neighbors(centers_kmeans,
#                                         g_test_cde,1)



cutoff_izbicki <- calibration_info %>%
  mutate(partition_id = as.numeric(partition_id)) %>%
  group_by(partition_id) %>%
  summarize(cutoff_izbicki = quantile(calibration_cde_values,
                                      1 - confidence_level))

delta_y = .001
yy <- seq(-8, 8, by = delta_y)

if (verbose) {
  pb <- progress::progress_bar$new(
    format = "Containment calculation [:bar] :percent eta: :eta",
    total = nrow(test_set_discrete), clear = FALSE, width = 80)
}

df_grid_all_iz_cd <- data.frame()
df_containment_prop_all_iz_cd <- data.frame()


for (calibration_idx in 1:nrow(test_set_discrete)){
  # for visual
  x_inner <- test_set_discrete$x[calibration_idx]
  cde_values <- simulationBands::cde_lei_wassserman(x_inner)(yy)

  partition_id <- which_partition_test[calibration_idx]
  inner_cutoff <- cutoff_izbicki$cutoff_izbicki[cutoff_izbicki$partition_id == partition_id]
  contained <- cde_values > inner_cutoff

  df_inner <- data.frame(
             calibration_idx = calibration_idx,
             x = x_inner,
             y = yy,
             contained = contained)
  df_grid_all_iz_cd <- rbind(df_grid_all_iz_cd, df_inner)

  if (verbose) {
      pb$tick()
  }

  # for containment of true observations:
  inner_sim_y <- sapply(rep(x_inner, n_sims_containment),
    function(x) {
      simulationBands::lei_wasserman_data_conditional_simulate(x, n = 1)[[1]]$sim})
  sim_cde_values <- simulationBands::cde_lei_wassserman(x_inner)(inner_sim_y)
  contained_check <- sim_cde_values > inner_cutoff
  df_containment_prop_inner <- data.frame(
    x = x_inner,
    containment_prop = mean(contained_check))
  df_containment_prop_all_iz_cd <- rbind(df_containment_prop_all_iz_cd,
                                   df_containment_prop_inner)
}


vis1iz <- df_grid_all_iz_cd %>% filter(contained) %>%
  ggplot() +
  geom_point(data = test_set, aes(x = x, y = y), color = "blue") +
  geom_tile(aes(x = x, y = y), alpha = .5) +
  labs(title = sprintf(paste0("Izbicki et. al (2020) ~Local~\n",
                              "Conformal Inference Prediction Regions (%s)"),
                       confidence_level_string)) +
  xlim(-1.5,1.5) + ylim(-8,8) +
  theme_minimal()

vis2iz <- df_containment_prop_all_iz_cd %>%
  ggplot() +
  geom_line(aes(x = x, y = containment_prop)) +
  ylim(0,1) + xlim(-1.5,1.5) +
  geom_hline(yintercept = confidence_level, color = "blue") +
  geom_label(data = data.frame(x = 1.5, y = confidence_level,
                               label = confidence_level),
             aes(x = x, y = y, label = label)) +
  labs(y = "containment proportion",
       title = sprintf(paste("containment proportion of %s",
                             "ys\nrandomly sampled per x value"),
                       n_sims_containment)) +
  theme_minimal()


gridExtra::grid.arrange(vis1iz, vis2iz,
             layout_matrix = matrix(c(1,1,1,1,
                                      1,1,1,1,
                                      1,1,1,1,
                                      2,2,2,2,
                                      2,2,2,2), byrow = T, ncol = 4))

if (save_image){
  ggvis <- gridExtra::arrangeGrob(vis1iz, vis2iz,
             layout_matrix = matrix(c(1,1,1,1,
                                      1,1,1,1,
                                      1,1,1,1,
                                      2,2,2,2,
                                      2,2,2,2), byrow = T, ncol = 4))
  ggplot2::ggsave(plot = ggvis, sprintf(paste0("quick_images/sim_bands_conf",
                                   "_test_discrete_1d_iz_cd_%s",
                                   "_sims_delta_size_%s_confidence_level_%s.png"),
                            n_simulations, delta_prop_image_string,
                            confidence_level_image_string))

  vis1iz_tikz <- df_grid_all_iz_cd %>% filter(contained) %>%
    ggplot() +
    geom_point(data = test_set, aes(x = x, y = y), color = "blue") +
    geom_tile(aes(x = x, y = y), alpha = .5) +
    labs(title = paste0("Izbicki et. al (2020) $\\sim$ Local $\\sim$ \\\\\ ",
                        "Conformal Inference Prediction Regions (60\\%)")) +
    xlim(-1.5,1.5) + ylim(-8,8) +
    theme_minimal()

  vis2iz_tikz <- df_containment_prop_all_iz_cd %>%
    ggplot() +
    geom_line(aes(x = x, y = containment_prop)) +
    ylim(0,1) + xlim(-1.5,1.5) +
    geom_hline(yintercept = confidence_level, color = "blue") +
    geom_label(data = data.frame(x = 1.5, y = confidence_level,
                                 label = confidence_level),
               aes(x = x, y = y, label = label)) +
    labs(y = "containment proportion",
         title = sprintf(paste("containment proportion of %s",
                               "ys\\\\\ randomly sampled per x value"),
                         n_sims_containment)) +
    theme_minimal()

  ggvis_tikz <- gridExtra::arrangeGrob(vis1iz_tikz, vis2iz_tikz,
                                  layout_matrix = matrix(c(1,1,1,1,
                                                           1,1,1,1,
                                                           1,1,1,1,
                                                           2,2,2,2,
                                                           2,2,2,2), byrow = T, ncol = 4))
  tikz(file = sprintf(paste0("quick_images/sim_bands_conf",
                             "_test_discrete_1d_iz_cd_%s",
                             "_sims_delta_size_%s_confidence_level_%s.tex"),
                      n_simulations, delta_prop_image_string,
                      confidence_level_image_string))
  plot(ggvis_tikz)
  dev.off()
}


# all together ----------

gridExtra::grid.arrange(grobs = list(
              vis1g, vis1lw, vis1iz, vis1,
              vis2g, vis2lw, vis2iz, vis2) %>%
               lapply(function(g){g + theme(text = element_text(size = 15),
                                            plot.title = element_text(size = 22))}),
             layout_matrix = matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                                      1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                                      1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                                      1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                                      5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,
                                      5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
                                    ncol = 4*4, byrow = T),
             bottom = paste("*a \"global\" version of",
                            "Izbicki et. al 2020's CD approach."))

if (save_image){
  ggvis <- gridExtra::arrangeGrob(grobs = list(
    vis1, vis1g, vis1lw, vis1iz,
    vis2, vis2g, vis2lw, vis2iz ) %>%
               lapply(function(g){g + theme(text = element_text(size = 6),
                                            plot.title = element_text(size = 12))}),
             layout_matrix = matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                                      1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                                      1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                                      1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                                      5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,
                                      5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
                                    ncol = 4*4, byrow = T),
             bottom = paste("*a \"global\" version of",
"Izbicki et. al 2020's CD approach."))
ggplot2::ggsave(plot = ggvis, width = 17, height = 8,
                filename = sprintf(paste0("quick_images/sim_bands_conf",
                                          "_test_discrete_1d_options_%s",
                                          "_sims_delta_size_%s_confidence_level_%s.png"),
                                   n_simulations, delta_prop_image_string,
                                   confidence_level_image_string))

  ggvis_tikz <- gridExtra::arrangeGrob(grobs = list(
    vis1g_tikz, vis1lw_tikz, vis1iz_tikz, vis1_tikz,
    vis2g, vis2lw, vis2iz_tikz, vis2) %>%
      lapply(function(g){g + theme(text = element_text(size = 6),
                                   plot.title = element_text(size = 12))}),
    layout_matrix = matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,
                             5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
                           ncol = 4*4, byrow = T),
    bottom = paste('*a "global" version of',
                   "Izbicki et. al 2020's CD approach."))

  tikz(file = sprintf(paste0("quick_images/sim_bands_conf",
                             "_test_discrete_1d_options_%s",
                             "_sims_delta_size_%s_confidence_level_%s_PRESENT.tex"),
                      n_simulations, delta_prop_image_string,
                      confidence_level_image_string),
       width = 17, height = 8)
  plot(ggvis_tikz)
  dev.off()

}


## presentation size ------


gridExtra::grid.arrange(grobs = list(
  vis1g, vis1lw, vis1iz, vis1,
  vis2g, vis2lw, vis2iz, vis2) %>%
    lapply(function(g){g + theme(text = element_text(size = 15),
                                 plot.title = element_text(size = 22))}),
  layout_matrix = matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                           1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                           1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                           1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                           5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,
                           5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
                         ncol = 4*4, byrow = T),
  bottom = paste("*a \"global\" version of",
                 "Izbicki et. al 2020's CD approach."))

if (save_image){
  ggvis <- gridExtra::arrangeGrob(grobs = c(list(
    vis1 + labs(title = "Our Simulation-based Approach"),
    vis1g + labs(title = "Global CD*"),
    vis1lw + labs(title = "Lei and Wasserman\n(2014) ~Local~"),
    vis1iz + labs(title = "Izbicki et al.\n(2020) ~Local~")
    ) %>%
      lapply(function(g){g + theme(text = element_text(size = 18),
                                   plot.title = element_text(size = 26))}),
    list(
      vis2 + labs(title =""),
      vis2g + labs(title ="") ,
    vis2lw + labs(title =""),
    vis2iz + labs(title ="")
    ) %>%
      lapply(function(g){g + theme(text = element_text(size = 18),
                                   plot.title = element_text(size = 18))})),
    layout_matrix = matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,
                             5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
                           ncol = 4*4, byrow = T),
    top = grid::textGrob("Conformal Inference Prediction Regions (60%)",
                         gp = grid::gpar(fontsize = 32)),
    bottom = grid::textGrob(paste("*a \"global\" version of",
                                  "Izbicki et. al 2020's CD approach."),
                            gp = grid::gpar(fontsize = 26)))
  ggplot2::ggsave(plot = ggvis, width = 17, height = 8,
                  filename = sprintf(paste0("quick_images/sim_bands_conf",
                                            "_test_discrete_1d_options_%s",
                                            "_sims_delta_size_%s_confidence_level_%s_PRESENT.png"),
                                     n_simulations, delta_prop_image_string,
                                     confidence_level_image_string))

  ggvis_tikz <- gridExtra::arrangeGrob(grobs = list(
    vis1g_tikz, vis1lw_tikz, vis1iz_tikz, vis1_tikz,
    vis2g, vis2lw, vis2iz_tikz, vis2) %>%
      lapply(function(g){g + theme(text = element_text(size = 6),
                                   plot.title = element_text(size = 12))}),
    layout_matrix = matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,
                             5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,
                             5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
                           ncol = 4*4, byrow = T),
    bottom = paste('*a "global" version of',
                   "Izbicki et. al 2020's CD approach."))

  tikz(file = sprintf(paste0("quick_images/sim_bands_conf",
                             "_test_discrete_1d_options_%s",
                             "_sims_delta_size_%s_confidence_level_%s.tex"),
                      n_simulations, delta_prop_image_string,
                      confidence_level_image_string),
       width = 17, height = 8)
  plot(ggvis_tikz)
  dev.off()
}


## UAI size --------

gridExtra::grid.arrange(grobs = list(
  vis1g, vis1lw, vis1iz, vis1,
  vis2g, vis2lw, vis2iz, vis2) %>%
    lapply(function(g){g + theme(text = element_text(size = 8),
                                 plot.title = element_text(size = 15))}),
  layout_matrix = matrix(c(1,1,1,1,2,2,2,2,
                           1,1,1,1,2,2,2,2,
                           1,1,1,1,2,2,2,2,
                           1,1,1,1,2,2,2,2,
                           5,5,5,5,6,6,6,6,
                           5,5,5,5,6,6,6,6,
                           rep(NA,8),
                           3,3,3,3,4,4,4,4,
                           3,3,3,3,4,4,4,4,
                           3,3,3,3,4,4,4,4,
                           3,3,3,3,4,4,4,4,
                           7,7,7,7,8,8,8,8,
                           7,7,7,7,8,8,8,8),
                         ncol = 4*2, byrow = T),
  bottom = grid::textGrob(paste("*a \"global\" version of",
                                "Izbicki et. al 2020's CD approach."),
                          gp=grid::gpar(fontsize=15)))

if (save_image){
  ggvis <- gridExtra::arrangeGrob(grobs = list(
    vis1,vis1g, vis1lw, vis1iz,
    vis2,vis2g, vis2lw, vis2iz) %>%
      lapply(function(g){g + theme(text = element_text(size = 11),
                                   plot.title = element_text(size = 11))}),
    layout_matrix = matrix(c(1,1,1,1,2,2,2,2,
                             1,1,1,1,2,2,2,2,
                             1,1,1,1,2,2,2,2,
                             1,1,1,1,2,2,2,2,
                             5,5,5,5,6,6,6,6,
                             5,5,5,5,6,6,6,6,
                             rep(NA,8),
                             3,3,3,3,4,4,4,4,
                             3,3,3,3,4,4,4,4,
                             3,3,3,3,4,4,4,4,
                             3,3,3,3,4,4,4,4,
                             7,7,7,7,8,8,8,8,
                             7,7,7,7,8,8,8,8),
                           ncol = 4*2, byrow = T),
    bottom = paste("*a \"global\" version of",
                   "Izbicki et. al 2020's CD approach."))
  ggplot2::ggsave(plot = ggvis, width = 8, height = 12,
                  filename = sprintf(paste0("quick_images/sim_bands_conf",
                                            "_test_discrete_1d_options_%s",
                                            "_sims_delta_size_%s_confidence_level_%s_UAI.png"),
                                     n_simulations, delta_prop_image_string,
                                     confidence_level_image_string))

  ggvis_tikz <- gridExtra::arrangeGrob(grobs = list(
    vis1g_tikz, vis1lw_tikz, vis1iz_tikz, vis1_tikz,
    vis2g, vis2lw, vis2iz_tikz, vis2) %>%
      lapply(function(g){g + theme(text = element_text(size = 6),
                                   plot.title = element_text(size = 12))}),
    layout_matrix = matrix(c(1,1,1,1,2,2,2,2,
                             1,1,1,1,2,2,2,2,
                             1,1,1,1,2,2,2,2,
                             1,1,1,1,2,2,2,2,
                             5,5,5,5,6,6,6,6,
                             5,5,5,5,6,6,6,6,
                             rep(NA,8),
                             3,3,3,3,4,4,4,4,
                             3,3,3,3,4,4,4,4,
                             3,3,3,3,4,4,4,4,
                             3,3,3,3,4,4,4,4,
                             7,7,7,7,8,8,8,8,
                             7,7,7,7,8,8,8,8),
                           ncol = 4*2, byrow = T),
    bottom = paste('*a "global" version of',
                   "Izbicki et. al 2020's CD approach."))
  tikz(file = sprintf(paste0("quick_images/sim_bands_conf",
                             "_test_discrete_1d_options_%s",
                             "_sims_delta_size_%s_confidence_level_%s_UAI.tex"),
                      n_simulations, delta_prop_image_string,
                      confidence_level_image_string),
       width = 17, height = 8)
  plot(ggvis_tikz)
  dev.off()

}


if (save_image){
  library(cowplot)

  vis_block_a <- plot_grid(vis1 + theme(text = element_text(size = 11),
                                  plot.title = element_text(size = 11)),
                           vis2 + theme(text = element_text(size = 11),
                                        plot.title = element_text(size = 11)),
                           ncol = 1, rel_heights = c(4,2),
                           labels = c("a)",""))

  vis_block_b <- plot_grid(vis1g + theme(text = element_text(size = 11),
                                        plot.title = element_text(size = 11)),
                           vis2g + theme(text = element_text(size = 11),
                                        plot.title = element_text(size = 11)),
                           ncol = 1, rel_heights = c(4,2),
                           labels = c("b)",""))

  vis_row_1 <- plot_grid(vis_block_a,
                         vis_block_b,
                         nrow = 1)


  vis_block_c <- plot_grid(vis1lw + theme(text = element_text(size = 11),
                                        plot.title = element_text(size = 11)),
                           vis2lw + theme(text = element_text(size = 11),
                                        plot.title = element_text(size = 11)),
                           ncol = 1, rel_heights = c(4,2),
                           labels = c("c)",""))

  vis_block_d <- plot_grid(vis1iz + theme(text = element_text(size = 11),
                                        plot.title = element_text(size = 11)),
                           vis2iz + theme(text = element_text(size = 11),
                                        plot.title = element_text(size = 11)),
                           ncol = 1, rel_heights = c(4,2),
                           labels = c("d)",""))



  vis_row_2 <- plot_grid(vis_block_c,
                         vis_block_d,
                         nrow = 1)

  vis_all <- plot_grid(vis_row_1,
                       vis_row_2,
                       nrow = 2)

  ggsave2(vis_all, filename = sprintf(paste0("quick_images/sim_bands_conf",
                                                        "_test_discrete_1d_options_%s",
                                                        "_sims_delta_size_%s_confidence_level_%s_UAI_cowplot.png"),
                                                 n_simulations, delta_prop_image_string,
                                                 confidence_level_image_string),
          height = 17, width = 8)
}


#install.packages("patchwork")
library(patchwork)


all_grobs <- list(
  vis1, vis1g, vis1lw, vis1iz,
  vis2, vis2g, vis2lw, vis2iz) %>%
  lapply(function(g){g + theme(text = element_text(size = 11),
                               plot.title = element_text(size = 11))})

patchwork_layout <- "
AAAABBBB
AAAABBBB
AAAABBBB
AAAABBBB
EEEEFFFF
EEEEFFFF
########
CCCCDDDD
CCCCDDDD
CCCCDDDD
CCCCDDDD
GGGGHHHH
GGGGHHHH
"

vis_all_patchwork <- all_grobs[[1]] + all_grobs[[2]] + all_grobs[[3]] + all_grobs[[4]] +
  all_grobs[[5]] + all_grobs[[6]] + all_grobs[[7]] + all_grobs[[8]] +
  plot_layout(design = patchwork_layout)


vis_all_patchwork %>% plot()
vis_all_patchwork + plot_annotation()


save(vis1, vis1g, vis1lw, vis1iz,
vis2, vis2g, vis2lw, vis2iz,file = "vis_files.Rdata")
