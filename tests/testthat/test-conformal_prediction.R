
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
