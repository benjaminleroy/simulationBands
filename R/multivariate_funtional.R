#' calculate maxmin distance between points
#'
#' Can be though of as a LOOCV containment of all points (what is the minimum
#' radius to do so).
#'
#' @param dist_mat a symmetric matrix with distance values, captures the
#' distance col observation is from row observation (digonals are expected to be
#' 0).
#'
#' @return minimum radius for all points to be covered
#' @export
#'
get_delta <- function(dist_mat){
  assertthat::assert_that(assertthat::are_equal(t(dist_mat), dist_mat),
                        msg = "t(dist_mat) needs to equal dist_mat")

  diag(dist_mat) <- max(dist_mat) # replacing diag so that it's selected as min
  mm_delta <- apply(dist_mat, MARGIN = 1, min) %>% max
  return(mm_delta)
}

#' @rdname get_delta
get_delta_simple <- get_delta



#' Identify the index in a raveled \code{dist} (upper triangle of distance)
#' relative to row / column ordering
#'
#' A extended idea as proposed on
#' \href{stack overflow}{https://stackoverflow.com/a/12643509} by Christian A.
#'
#' Note that i and j can be vectors, but if they're both vectors we
#' don't return all combinations of each.
#'
#' @param i row index (either a single value or vector). \code{i} < \code{j}.
#' @param j column index (either a single value or vector)
#' @param n number of objects that created the \code{dist} object
#' @param .check_order boolean (if we should check that \code{i} < \code{j})
#' and remove all values that don't follow that rule from returning.
#' @param .check_possible_ij boolean (if we should check that \code{i} and
#' \code{j} are not greater than \code{n}), and remove all values that don't
#' follow that rule from returning.
#'
#' @return requested indicies.
distdex <- function(i,j,n, .check_order = T,
                    .check_possible_ij = T) {
  values <- n*(i-1) - i*(i-1)/2 + j-i


  correct_values <- T

  if (.check_possible_ij){
    correct_values <- j <= n & i <= n
  }

  if (.check_order){
    correct_order <- i < j
    correct_values <- correct_values & correct_order
  }

  return(values[correct_values])
}


#' Named after assuming positive integers, but not checked
#'
#' @param x vector of numerical values (or single value)
#'
#' @return subset of \code{x} that greater than 0.
pos_int <- function(x){
  return(x[x > 0])
}


#' calculate maxmin distance between points (\code{dist} based)
#'
#' Can be though of as a LOOCV containment of all points (what is the minimum
#' radius to do so).
#'
#' @param x \code{dist} based object
#' @param verbose boolean, if should present a progress bar
#'
#' @return minimum radius for all points to be covered
#' @export
get_delta_dist <- function(x, verbose = F){
  assertthat::assert_that(inherits(x, "dist"),
                          msg = "dist should be of class dist")


  n <- .5 + .5 * sqrt(1 + 8*length(x))

  if(verbose){
    pb <- progress::progress_bar$new(
      format = "processing row minimums [:bar] :percent eta: :eta",
      total = n, clear = FALSE, width = 70)

  }
  vals <- rep(NA, n)
  for (i in 1:n){
    right_side_id <- distdex(i, pos_int((i+1):n), n)
    left_side_id <- distdex(pos_int(1:(i-1)), i, n)
    vals[i] <- min(x[c(right_side_id, left_side_id)])
    if (verbose) {
      pb$tick()
    }

  }
  return(max(vals))
}


#' euclidean distance between 2 sets
#'
#' @param x data.frame (n x p)
#' @param y data.frame (m x p)
#'
#' @return distance between x points (rows) and y points (columns) - (n x m)
distEuclidean <- function(x, y){
  centers <- as.matrix(y)
  z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
  for(k in 1:nrow(centers)){
    z[,k] <- sqrt( colSums((t(x) - centers[k,])^2) )
  }
  z
}


#' euclidean distance between 2 sets
#'
#' @param df1 data.frame (n x p)
#' @param df2 data.frame (m x p)
#'
#' @return distance between x points (rows) and y points (columns) - (n x m)
crossdist_df <- function(df1, df2){
  crossdist(as.matrix(df1), as.matrix(df2))
}



#' calculate maxmin distance between points (breaks into parts)
#'
#' Can be though of as a LOOCV containment of all points (what is the minimum
#' radius to do so).
#'
#' @param df data.frame
#' @param full_breaks number of breaks (1 more if not perfectly split)

#' @param verbose boolean, if should present a progress bar
#'
#' @return minimum radius for all points to be covered
#' @export
get_delta_large <- function(df, full_breaks = 10, verbose = F){
  # index breaks
  n <- nrow(df)
  mod_diff <- n %% full_breaks
  length_block <- floor(n/full_breaks)

  idx_sets <- lapply(1:full_breaks,
                     function(id) {seq((id-1)*length_block + 1,
                                      id*length_block)})
  if (mod_diff  > 0){
    idx_sets[[full_breaks+1]] <- seq(full_breaks*length_block + 1, n)
  }

  if(verbose){
    pb <- progress::progress_bar$new(
      format = "processing blocks [:bar] :percent eta: :eta",
      total = .5 * length(idx_sets)*(length(idx_sets)-1) + length(idx_sets),
      clear = FALSE, width = 70)
  }

  min_vec <- rep(Inf, n)
  for (b_idx in 1:length(idx_sets)){
    for (bb_idx in b_idx:length(idx_sets)){
      if (b_idx == bb_idx){
        #dmat <- as.matrix(dist(df[idx_sets[[b_idx]],]))
        if (ncol(df) == 1){
          dmat <- crossdist_df(df[idx_sets[[b_idx]],, drop = F],
                               df[idx_sets[[b_idx]],, drop = F])
        } else {
          dmat <- crossdist_df(df[idx_sets[[b_idx]],],
                               df[idx_sets[[b_idx]],])
        }
        diag(dmat) <- Inf
        min_vec_inner <- apply(dmat, 1, min)

        min_vec[idx_sets[[b_idx]]] <- cbind(
          min_vec[idx_sets[[b_idx]]],
          min_vec_inner) %>% apply(1, min)

      } else {
        if (ncol(df) == 1){
          dmat <- crossdist_df(df[idx_sets[[b_idx]],, drop = F],
                                df[idx_sets[[bb_idx]],, drop = F])
        } else {
          dmat <- crossdist_df(df[idx_sets[[b_idx]],],
                                df[idx_sets[[bb_idx]],])
        }
        row_vec <- apply(dmat, 1, min)
        col_vec <- apply(dmat, 2, min)

        min_vec[idx_sets[[b_idx]]] <- cbind(
          min_vec[idx_sets[[b_idx]]],
          row_vec) %>% apply(1, min)

        min_vec[idx_sets[[bb_idx]]] <- cbind(
          min_vec[idx_sets[[bb_idx]]],
          col_vec) %>% apply(1, min)

      }
      if (verbose) {
        pb$tick()
      }
    }

  }
  return(max(min_vec))
}


#' calculate maxmin distance between points
#'
#' A wrapper of \code{get_delta_simple}, \code{get_delta_dist}, and
#' \code{get_delta_large} relative to the number of observations to be used.
#' Attempts to balance the number of observations vs the fastest approach (
#' attempting to balance time costs from memory storage and processing speed).
#'
#'
#' @param df data.frame (with only columns that are needed)
#' @param full_breaks number of breaks relative to \code{get_delta_large},
#' otherwise it's not used.
#' @param verbose boolean, a progressbar is provided  to track the progress of
#' \code{get_delta_dist} or \code{get_delta_large}, if they are used.
#'
#' @return minimum radius for all points to be covered
#' @export
get_delta_flex <- function(df, full_breaks = 10, verbose = F){
  n <- nrow(df)
  if (n < 500){
    return(get_delta_simple(as.matrix(stats::dist(df))))
  } else if (n < 5000) {
    return(get_delta_dist(stats::dist(df), verbose = verbose))
  } else {
    return(get_delta_large(df, full_breaks = full_breaks,
                           verbose = verbose))
  }
}

#' calculate maxmin distance between points
#'
#' Fast calculation of maxmin distance using kd trees and nearest neighbors from
#' the \code{RANN} package.
#'
#' @param df data.frame (with only columns that are needed)
#'
#' @return minimum radius for all points to be covered
#' @export
get_delta_nn <- function(df){
  check <- RANN::nn2(df, df, k = 2, treetype = "kd", eps = 0)
  mm_delta <- check$nn.dists[,2] %>% max()
  return(mm_delta)
}


#' minimum distance to contain all points of first set with union of balls of
#' second set
#'
#' Inner function, slow version. For each row in \code{df_row}, the minimum
#' distance to a point in \code{df_col} is calculated and then the maximum of
#' this vector is returned.
#'
#' @param df_row data.frame, rows are observations, all columns are used in
#' distance calculation.
#' @param df_col data.frame, rows are observations, all columns are used in
#' distance calculation.
#'
#' @return single minimum distance scalar
#' @export
#'
maxmin_inner_old <- function(df_row, df_col){
  p <- dim(df_row)[2]

  assertthat::assert_that(assertthat::are_equal(p, dim(df_col)[2]),
                          msg = "DFs have different dimension (p)")

  maxmin <- -Inf
  for (i_idx in 1:nrow(df_row)){
    min_inner <- Inf
    for (j_idx in 1:nrow(df_col)){
      min_inner <- min(c(stats::dist(rbind(df_row[i_idx,], df_col[j_idx,])),
                         min_inner))
    }
    maxmin <- max(c(min_inner, maxmin))
  }

  return(maxmin)
}


#' minimum distance to contain all points of first set with union of balls of
#' second set
#'
#' Inner function, faster version. For each row in \code{df_row}, the minimum
#' distance to a point in \code{df_col} is calculated and then the maximum of
#' this vector is returned. If \code{only_val} is \code{FALSE}, then the
#' minimum distance is returned, else the vector per each row of \code{df_row}'s
#' minimum distance to a point in \code{df_col} is returned.
#'
#' @param df_row data.frame, rows are observations, all columns are used in
#' distance calculation.
#' @param df_col data.frame, rows are observations, all columns are used in
#' distance calculation.
#' @param only_val boolean, effects return of scalar (TRUE) or vector (FALSE)
#' as described in details.
#'
#' @return single minimum distance scalar or minimum distance vector
#' @export
maxmin_inner <- function(df_row, df_col, only_val = T){
  p <- dim(df_row)[2]

  assertthat::assert_that(assertthat::are_equal(p, dim(df_col)[2]),
                          msg = "DFs have different dimension (p)")

  nn_info <- RANN::nn2(data = df_col, query = df_row, k = 1, treetype = "kd")
  if (only_val){
    return(max(nn_info$nn.dists))
  } else {
    return(list(max(nn_info$nn.dists), nn_info$nn.dists))
  }

}


#' Calculate minimum radius for each group to contain the true data points.
#'
#' Slow version.
#'
#' @param truth_df data.frame of true points (locations defined by
#' \code{data_column_names})
#' @param simulations_grouped_df grouped_df, with points grouped together
#' @param data_column_names column names for both the \code{truth_df} and
#' \code{simulations_grouped_df} that contain locations.
#'
#' @return data.frame of grouped information and minimum distance per group
#' in the \code{simulations_grouped_df}.
#' @export
maxmin_distance_vector_old <- function(truth_df,
                                       simulations_grouped_df,
                                       data_column_names = c("S", "I", "R")){


  mm_df <- simulations_grouped_df %>%
    tidyr::nest() %>%
    dplyr::mutate(maxmin_dist = purrr::map(.data$data, function(df){
      maxmin_inner_old(df_row = truth_df[data_column_names],
                   df_col = df[data_column_names])
    })) %>%
    dplyr::select(-.data$data) %>%
    tidyr::unnest(.data$maxmin_dist)

  return(mm_df)
}

#' Calculate minimum radius for each group to contain the true data points.
#'
#' Fast version.
#'
#' @param truth_df data.frame of true points (locations defined by
#' \code{data_column_names})
#' @param simulations_grouped_df grouped_df, with points grouped together
#' @param data_column_names column names for both the \code{truth_df} and
#' \code{simulations_grouped_df} that contain locations.
#' @param .all_info boolean, indicators what should be returned.
#'
#' @return If \code{.all_info} is \code{FALSE} then a data.frame of grouped
#' information and minimum distance per group in the
#' \code{simulations_grouped_df}. Else both that and a data.frame that details
#' the minimum covering distance for each point in the \code{truth_df} (by the
#' grouped set of points) is given. Note the order of the grouping index when
#' using these data.frames.
#' @export
maxmin_distance_vector <- function(truth_df, simulations_grouped_df,
                                    data_column_names = c("S", "I", "R"),
                                    .all_info = F){

  if(!.all_info){
    mm_df <- simulations_grouped_df %>%
      tidyr::nest() %>%
      dplyr::mutate(maxmin_dist = purrr::map(.data$data, function(df){
        maxmin_inner(df_row = truth_df[data_column_names],
                      df_col = df[data_column_names])
      })) %>%
      dplyr::select(-.data$data) %>%
      tidyr::unnest(.data$maxmin_dist)
  } else {
    all_info_list <- simulations_grouped_df %>% dplyr::group_split()
    all_info_names <- simulations_grouped_df %>% dplyr::group_keys()

    dist_info_list <- lapply(all_info_list, function(df) {
      maxmin_inner(df_row = truth_df[data_column_names],
                    df_col = df[data_column_names],only_val = F)[[2]]})


    dist_info_matrix <- t(do.call(cbind, dist_info_list))
    dist_info_df <- cbind(all_info_names,
                          data.frame(dist_info_matrix))

    mm_df <- cbind(all_info_names, maxmin_dist = apply(dist_info_matrix,1, max))
    return(list(mm_df, dist_info_df))
  }
}


