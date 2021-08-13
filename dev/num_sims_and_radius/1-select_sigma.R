
# Rscript 1-select_sigma.R [1-6] [1-10] [default or 10:50::5]
# 1. index for number of simulations,
# 2. index for attempt index
# 3. range of sigma % values
# input parameters -----------------------------
start_time <- Sys.time()
input_args <- commandArgs(trailingOnly=TRUE)
n_simulations <- c(50,100,200,500,1000,1500)[as.numeric(input_args[1])]
seed_value <- c(707, 412, 4909, 4006, 5700, 503)[as.numeric(input_args[1])] +
  as.numeric(input_args[2])
set.seed(seed_value)


input_sigma_info <- input_args[3]
if (input_sigma_info == "default"){
  range_sigma <- paste0(c(1:5, seq(7.5, 90, by = 2.5)), "%")
  input_sigma_info_str <- input_sigma_info
} else {
  sigma_range_info <- input_sigma_info %>% stringr::str_split(":{1,2}", simplify = T) %>%
    as.numeric()
  range_sigma <- paste0(seq(sigma_range_info[1],
                            sigma_range_info[2],
                            by = sigma_range_info[3]),
                        "%")
  input_sigma_info_str <- paste0(sigma_range_info, collapse = "-")
}


# libraries ----------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(devtools))

# local packages
suppressPackageStartupMessages(load_all("../../"))
suppressPackageStartupMessages(load_all("../../../EpiCompare/")) # for tidy_dist_mat?


# Functions ------------------

#' Title
#'
#' @param df_info
#' @param position
#' @param sigma
#' @param maxT
#' @param range_eps
#' @param range_diff_eps
#' @param verbose
#'
#' @return
#' @export
#'
mode_clustering_1d_ranges <- function(df_info,
                                     position,
                                     sigma,
                                     maxT = 200,
                                     range_eps = 10^seq(-5,-8, by = -.5),
                                     range_diff_eps= 10^seq(-5,-8, by = -.5),
                                     verbose = FALSE){
  mat_info_raw <- df_info[,position] %>% as.matrix()

  full_info <- data.frame()

  inner_range_eps <- range_eps

  # clustering with with maxT and range_eps - get all progressions (so return multiple G values)
  mode_cluster_out <- inner_cluster_1d_range(mat_info_raw,
                                             sigma = sigma,
                                             maxT = maxT,
                                             range_eps = inner_range_eps,
                                             verbose = verbose)

  # across G values - get number of clusters relative to range_diff_eps values
  for (info_idx in 1:length(mode_cluster_out)){
    string_eps <- names(mode_cluster_out)[info_idx]
    actualT = mode_cluster_out[[info_idx]]$num_steps
    # for each element calculate
    # (1) distance matrix
    dist_mat <- as.matrix(dist(mode_cluster_out[[info_idx]]$G))

    # (2) examine across cutoffs the number of groups
    number_clusters <- sapply(range_diff_eps, function(diff_eps) EpiCompare::check_number_of_groups(dist_mat, diff_eps))

    # (2.5) group clustering
    # for (diff_eps in range_diff_eps){
    #   inner_grouping_vec <- EpiCompare::get_groups(dist_mat, diff_eps)
    #   inner_grouping_df <- g_names %>% mutate(groupings = inner_grouping_vec)
    #   groups_cluster_list[[paste0(string_eps, "-", diff_eps)]] <- inner_grouping_df
    # }
    #
    #
    # # (3) storage
    # steps_counter <- steps_counter + out_list[[1]][[2]]
    # inner_info <- data.frame(maxT = range_maxT[i_maxT],
    #                          total_steps = steps_counter,
    #                          eps = as.numeric(string_eps),
    #                          diff_eps = range_diff_eps,
    #                          number_clusters = number_clusters)
    # full_info <- rbind(full_info, inner_info)
    inner_info <- data.frame(eps = string_eps,
                             diff_eps = range_diff_eps,
                             actualT = actualT,
                             maxT = maxT,
                             number_clusters = number_clusters)
    full_info <- rbind(full_info, inner_info)
  }
  return(full_info)
}

if (FALSE){
testthat::test_that("test mode_clustering_1d_ranges", {
  set.seed(1)
  X_mat <- matrix(rnorm(600, mean = 2), ncol = 2)
  X_neg <- matrix(rnorm(600, mean = -2), ncol = 2)
  X_both <- rbind(X_mat, X_neg)
  sigma <- .75

  # 1 mode ===
  out1 <- mode_clustering_1d_ranges(X_mat, position = 1:2,
                            sigma = sigma,
                            maxT = 300,
                            range_eps = 10^seq(-10,-15, by = -1),
                            range_diff_eps= 10^seq(-6, -15, by = -1),
                            verbose = FALSE)

  out1_matrix <- out1 %>% select(-actualT, -maxT) %>%
    pivot_wider(names_from = "diff_eps", values_from = "number_clusters") %>%
    select(-eps) %>% as.matrix

  # monotonically increasing by columns
  for (c_idx in 1:(ncol(out1_matrix)-1)){
    testthat::expect_true(all(out1_matrix[,c_idx] <= out1_matrix[,c_idx+1]))
  }
  # monotonically decreasing by row
  for (r_idx in 1:(nrow(out1_matrix)-1)){
    testthat::expect_true(all(out1_matrix[r_idx,] >= out1_matrix[r_idx+1,]))
  }
  testthat::expect_equivalent(out1_matrix[1,1],1)


  # 2 modes ===
  out2 <- mode_clustering_1d_ranges(X_both, position = 1:2,
                                    sigma = sigma,
                                    maxT = 300,
                                    10^seq(-10,-15, by = -2),
                                    range_diff_eps= 10^seq(-6, -15, by = -1),
                                    verbose = TRUE)

  out2_matrix <- out2 %>% select(-actualT, -maxT) %>%
    pivot_wider(names_from = "diff_eps", values_from = "number_clusters") %>%
    select(-eps) %>% as.matrix

  # monotonically increasing by columns
  for (c_idx in 1:(ncol(out2_matrix)-1)){
    testthat::expect_true(all(out2_matrix[,c_idx] <= out2_matrix[,c_idx+1]))
  }
  # monotonically decreasing by row
  for (r_idx in 1:(nrow(out2_matrix)-1)){
    testthat::expect_true(all(out2_matrix[r_idx,] >= out2_matrix[r_idx+1,]))
  }
  testthat::expect_equivalent(out2_matrix[1,1],2)


})
}
#' mode cluster for range of eps values
#'
#' @param X_mat
#' @param sigma
#' @param maxT
#' @param range_eps
#'
#' @return
#' a list (length of range_eps) of lists containing the current "G" and actual
#' "num_steps".
#'
#' @export
#'
inner_cluster_1d_range <- function(X_mat,
                                   sigma,
                                   maxT,
                                   range_eps = 10^seq(-5,-8, by = -.5),
                                   verbose = FALSE){

 G_mat <- X_mat
 inner_maxT <- maxT
 inner_range_eps <- range_eps

 assertthat::assert_that(all(diff(range_eps) <= 0)) # decreasing eps required

 out_list <- list()

 if (verbose){
   pb <- progress::progress_bar$new(
     format = "processing [:bar] :percent in :elapsed",
     total = length(range_eps), width = 80
   )
 }

 while (inner_maxT > 0  && length(inner_range_eps) > 0){
  eps <- inner_range_eps[1]

  if (length(inner_range_eps) > 1){
    inner_range_eps <- inner_range_eps[2:length(inner_range_eps)]
  } else{
    inner_range_eps <- c()
  }

  process_out <- psuedo_density_mode_cluster_1d(
    X_mat,
    G_mat = G_mat,
    sigma = sigma,
    maxT = inner_maxT,
    eps = eps,
    verbose = FALSE,
    list_out = FALSE)


  # update maxT and save info
  inner_maxT <- inner_maxT - process_out$t
  G_mat <- process_out$G_mat
  out_list[[as.character(eps)]] <- list()
  out_list[[as.character(eps)]][["num_steps"]] <- maxT - inner_maxT
  out_list[[as.character(eps)]][["G"]] <- G_mat

  if (verbose){
    pb$tick()
  }
 }

  return(out_list)
}

if (FALSE){
testthat::test_that("test inner_cluster_1d_range, basic", {
  set.seed(1)
  X_mat <- matrix(rnorm(600, mean = 2), ncol = 2)
  #X_neg <- matrix(rnorm(600, mean = -2), ncol = 2)
  #X_both <- rbind(X_mat, X_neg)

  # 1 mode ===
  sigma <- .25

  # we expect an error since range_eps isn't monotonically decreasing
  testthat::expect_error(inner_cluster_1d_range(X_mat, sigma = sigma,
                         max = 50,
                         range_eps = 10^seq(-8,-5, by = .5)))

  process1_max110 <- inner_cluster_1d_range(X_mat, sigma = sigma,
                         max = 110,
                         range_eps = 10^seq(-5,-8, by = -.5))
  # checking with verbose:
  process1_max200 <- inner_cluster_1d_range(X_mat, sigma = sigma,
                                            max = 200,
                                            range_eps = 10^seq(-5,-8, by = -.5),
                                            verbose = T)

  n_small <- length(process1_max110)
  if (n_small > 2){
    for (G_idx in 1:(n_small-1)){
      testthat::expect_equal(process1_max110[[G_idx]]$G,
                             process1_max200[[G_idx]]$G)
      testthat::expect_equal(process1_max110[[G_idx]]$num_steps,
                             process1_max200[[G_idx]]$num_steps)
    }
  }

  actual_num_steps1_max110 <- sapply(process1_max110, function(inner_l) inner_l$num_steps)
  actual_num_steps1_max200 <- sapply(process1_max200, function(inner_l) inner_l$num_steps)

  testthat::expect_true(all(actual_num_steps1_max110 <= 110))
  testthat::expect_true(all(actual_num_steps1_max200 <= 200))
  testthat::expect_true(all(actual_num_steps1_max110 <=
                              actual_num_steps1_max200[1:length(process1_max110)]))

  testthat::expect_true(all(c(length(process1_max110),
                              length(process1_max200)) <= length(seq(-5,-8, by = -.5))))
})
}

# Analysis -------------------
# 1. look at 3 data distributions
# 2. examine the mode clustering across range of tuning parameters



## three data distributions -------------------------------
# we will draw from x = [-1., 0, .6]

x_values = c(-1,0,.6)

data_list <- simulationBands::lei_wasserman_data_conditional_simulate(x_values,
                                                                      n = n_simulations)

##

sigma_info_list_x <- list()

if (TRUE){
  pb <- progress::progress_bar$new(
    format = "processing: [:bar] :percent elapsed: :elapsed | eta: :eta ",
    width = 80, total = length(range_sigma)*3
  )
}

for (x_idx in 1:3){
  actual_x_value <- x_values[x_idx]

  dist_matrix <- dist(data_list[[x_idx]]) %>% as.matrix()

  overall_info <- data.frame()

  for (.sigma_string in range_sigma){
    # converting .sigma_string to numerical value.
    sigma_lower <- EpiCompare::check_character_percent(.sigma_string, ".sigma_string")
    sigma_sizes <- sapply(sigma_lower + .05*(0:5), function(v) min(v, 1))

    percentage_inner <- sigma_sizes[stats::quantile(dist_matrix, sigma_sizes) > 0][1]

    sigma_val <- stats::quantile(dist_matrix, percentage_inner)
    sigma <- sigma_val

    # (1) model clustering ----------
    info_out <- mode_clustering_1d_ranges(df_info = data_list[[x_idx]],
                              position = 1:2,
                              sigma = sigma,
                              maxT = 300,
                              range_eps = 10^seq(-10,-15, by = -1),
                              range_diff_eps= 10^seq(-6, -16, by = -2),
                              verbose = FALSE)

    info_out <- info_out %>%
      mutate(sigma = sigma,
             actual_percentage = percentage_inner,
             .sigma_string = .sigma_string,
             x = actual_x_value)

    overall_info <- rbind(overall_info, info_out)
    pb$tick()
  }

  sigma_info_list_x[[x_idx]] <- overall_info
}

overall_time <- Sys.time() - start_time

save(sigma_info_list_x, overall_time,
     file = sprintf("1d_select_sigma_n_sims_%i_seed_%i_sigrange_%s.Rdata",
                    n_simulations, seed_value, input_sigma_info_str))
