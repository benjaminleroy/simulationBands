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
#'
#' @export
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


#' global wrapper for simulation-based conformal score calculation that allows
#' for mode clustering and changes of radius for 1d examples
#'
#' @param truth_df data frame with new point information to be assessed
#' @param simulations_df data frame with multiplee simulated points
#' @param smooth_functions  a list of functions a function for the smoothing
#' approach for the varied approach
#' @param data_column_names columns of \code{truth_df} and
#' \code{simulations_df} that correspond to the output space.
#' @param .maxT int, max number of steps for mode clustering mean-shift
#' algorithm
#' @param .sigma_string string, the quantile of the distance matrix to define
#' the sigma (e.g. \code{"30\%"})
#' @param .diff_eps float, the error allows between mode clustered final steps
#' to still be counted as the same mode.
#' @param .eps float, when to stop the model clustering (if maxT doesn't stop it
#' first)
#' @param return_min boolean, if list of information returned is mimimum (for slurm)
#' @param verbose boolean, be verbose about progression
#'
#' @return list of information...
#' @export
#'
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







