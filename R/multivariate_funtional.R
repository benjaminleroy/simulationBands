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


