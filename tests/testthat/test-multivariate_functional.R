context("multivariate functional functions")

library(dplyr)

testthat::test_that("test get_delta function, basic", {
  non_sym1 <- matrix(1:16, nrow = 4)
  non_sym2 <- matrix(1:12, nrow = 3)

  testthat::expect_error(get_delta(non_sym1))
  testthat::expect_error(get_delta(non_sym2))

  d <- data.frame(x = 1:5)
  dmat <- as.matrix(dist(d))
  testthat::expect_equal(get_delta(dmat),1)

  d2 <- data.frame(x = c(1,3:5))
  dmat2 <- as.matrix(dist(d2))
  testthat::expect_equal(get_delta(dmat2),2)

})

testthat::test_that("test get_delta_simple function, basic", {
  non_sym1 <- matrix(1:16, nrow = 4)
  non_sym2 <- matrix(1:12, nrow = 3)

  testthat::expect_error(get_delta_simple(non_sym1))
  testthat::expect_error(get_delta_simple(non_sym2))

  d <- data.frame(x = 1:5)
  dmat <- as.matrix(dist(d))
  testthat::expect_equal(get_delta_simple(dmat),1)

  d2 <- data.frame(x = c(1,3:5))
  dmat2 <- as.matrix(dist(d2))
  testthat::expect_equal(get_delta_simple(dmat2),2)

})


testthat::test_that("test maxmin_inner_old, basic", {
  df_row_e <- data.frame(x = 1:5, y = 1:5)
  df_col_e <- data.frame(x= 1:5)

  testthat::expect_error(maxmin_inner_old(df_row_e, df_col_e))

  df_row <- data.frame(x = 1:5)
  df_col <- data.frame(x = 2:6)

  testthat::expect_equal(maxmin_inner_old(df_row, df_col),1)
  testthat::expect_equal(maxmin_inner_old(df_row, df_row),0)


})

testthat::test_that("test maxmin_inner, basic", {
  df_row_e <- data.frame(x = 1:5, y = 1:5)
  df_col_e <- data.frame(x= 1:5)

  testthat::expect_error(maxmin_inner(df_row_e, df_col_e))

  df_row <- data.frame(x = 1:5)
  df_col <- data.frame(x = 2:6)

  testthat::expect_equal(maxmin_inner(df_row, df_col, only_val = T),1)
  testthat::expect_equal(maxmin_inner(df_row, df_col, only_val = F)[[1]],1)
  testthat::expect_equal(maxmin_inner(df_row, df_col, only_val = F)[[2]],
                         matrix(c(1,0,0,0,0), ncol = 1))


  testthat::expect_equal(maxmin_inner(df_row, df_row, only_val = T),0)
  testthat::expect_equal(maxmin_inner(df_row, df_row, only_val = F)[[1]],0)
  testthat::expect_equal(maxmin_inner(df_row, df_row, only_val = F)[[2]],
                         matrix(c(0,0,0,0,0), ncol = 1))
})


testthat::test_that("test maxmin_distance_vector_old, basic", {
  df_row <- data.frame(x = 1:5)
  df_col <- data.frame(x = 2:6)
  sim_df <- df_col %>% mutate(id = 1) %>%
    rbind(df_row %>% mutate(id = 2)) %>%
    dplyr::group_by(id)

  df_out <- maxmin_distance_vector_old(df_row,
                                       sim_df,
                                       data_column_names = c("x"))

  testthat::expect_equivalent(df_out, data.frame(id = 1:2,
                                                 maxmin_dist = c(1,0)))
})


testthat::test_that("test maxmin_distance_vector, basic", {
  df_row <- data.frame(x = 1:5)
  df_col <- data.frame(x = 2:6)
  sim_df <- df_col %>% mutate(id = 2) %>%
    rbind(df_row %>% mutate(id = 1)) %>%
    dplyr::group_by(id)

  df_out <- maxmin_distance_vector(df_row,
                                   sim_df,
                                   data_column_names = c("x"))

  testthat::expect_equivalent(df_out, data.frame(id = 2:1,
                                                 maxmin_dist = c(1,0)))

  mm_info <- maxmin_distance_vector(df_row,
                                    sim_df,
                                    data_column_names = c("x"),
                                    .all_info = T)
  testthat::expect_equivalent(df_out[2:1,], mm_info[[1]]) # ordered differently
  testthat::expect_equivalent(mm_info[[2]],
                              data.frame(cbind(c(1,2),
                                               matrix(c(0,0,0,0,0,
                                                        1,0,0,0,0),
                                                      nrow = 2, byrow = T))))

})

