
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


#' Title
#'
#' @param truth_df
#' @param simulations_grouped_df grouped data frame...
#' @param data_column_names
#' @return
#' @export
#'
#' @examples
simulation_based_conformal_1d <- function(truth_df, simulations_grouped_df,
                                          data_column_names = c("y"),
                                          delta_prop = .8){


  assertthat::assert_that(inherits(simulations_grouped_df, "grouped_df"))
  group_names <- names(group_keys(simulations_grouped_df))

  truth_df_inner <- truth_df %>%
    dplyr::select(dplyr::one_of(data_column_names))

  simulations_group_df_inner <- simulations_grouped_df %>%
    dplyr::select(dplyr::one_of(c(group_names, data_column_names)))

  group_info <- simulations_group_df_inner %>% group_keys()

  dist_mat <- simulations_group_df_inner %>% group_split() %>%
    do.call(rbind, .) %>% # match group_info ordering
    select(one_of(data_column_names)) %>%
    dist() %>% as.matrix()

  tdm_sims <- EpiCompare::tidy_dist_mat(dist_mat, group_info, group_info)

  # sigma selection

  sigma_size <- c("20%" = .2, "25%" = .25, "30%" = .3,
                  "35%" = .35, "40%" = .4, "45%" = .45)

  percentage <- names(sigma_size)[stats::quantile(as.matrix(tdm_sims), sigma_size) > 0][1]


  # rank_df
  pseudo_density_df <- EpiCompare::distance_psuedo_density_function(
    tdm_sims,
    sigma = percentage, df_out = T) %>%
    mutate(ranking = rank(psuedo_density,ties.method = "min")) #spelling error... :(

  assertthat::assert_that(all(!is.na(pseudo_density_df$psuedo_density)),
                          msg = paste("internal error in",
                                      "distance_psuedo_density_function",
                                      "function's sigma selection."))

  mm_df <- maxmin_distance_vector(truth_df = truth_df_inner,
                                   simulations_grouped_df = simulations_group_df_inner,
                                   data_column_names = c("y"),
                                   .all_info = F)

  proportion_points_not_included <- 1 - delta_prop

  top_points <- simulations_group_df_inner %>%
    left_join(pseudo_density_df, by = group_names) %>%
    mutate(keep = ranking > ceiling(proportion_points_not_included*nrow(simulations_group_df_inner))) %>%
    ungroup() %>% filter(keep) %>%
    select(one_of(data_column_names))


  mm_delta <- get_delta_nn(top_points)

  containment_df <- pseudo_density_df %>%
    left_join(mm_df, by = group_names) %>%
    mutate(delta_close = maxmin_dist < mm_delta)


  conformal_score <- max(c(containment_df$ranking[containment_df$delta_close],
                           0))

  return(list(conformal_score = conformal_score, containment_df = containment_df,
              mm_delta = mm_delta,
              truth_df_inner = truth_df_inner,
              simulations_group_df_inner = simulations_group_df_inner,
              parameters = c("mm_delta_prop" = proportion_points_not_included,
                             "sigma_percentage" = percentage)))
}

#' @param conformal_score_cut note this can be an integer, string percentage or
# fraction, but relates to the conformal cut (not the confidence level).
#'
simulation_based_conformal_1d_region <- function(simulations_grouped_df,
                                                 data_column_names = c("y"),
                                                 conformal_score_cut = .9,
                                                 delta_prop = .8){

  if(!EpiCompare::is.wholenumber(conformal_score_cut)){
    if(is.character(conformal_score_cut)){
      inner_percent <- EpiCompare::check_character_percent(conformal_score_cut)
    } else {
      inner_percent <- conformal_score_cut
    }
    conformal_score_cut <- ceiling(inner_percent * nrow(simulations_grouped_df))
  }

  assertthat::assert_that(inherits(simulations_grouped_df, "grouped_df"))
  group_names <- names(group_keys(simulations_grouped_df))


  simulations_group_df_inner <- simulations_grouped_df %>%
    dplyr::select(dplyr::one_of(c(group_names, data_column_names)))

  group_info <- simulations_group_df_inner %>% group_keys()

  dist_mat <- simulations_group_df_inner %>% group_split() %>%
    do.call(rbind, .) %>% # match group_info ordering
    select(one_of(data_column_names)) %>%
    dist() %>% as.matrix()

  tdm_sims <- EpiCompare::tidy_dist_mat(dist_mat, group_info, group_info)

  # sigma selection

  sigma_size <- c("20%" = .2, "25%" = .25, "30%" = .3,
                  "35%" = .35, "40%" = .4, "45%" = .45)

  percentage <- names(sigma_size)[stats::quantile(as.matrix(tdm_sims), sigma_size) > 0][1]


  # rank_df
  pseudo_density_df <- EpiCompare::distance_psuedo_density_function(
    tdm_sims,
    sigma = percentage, df_out = T) %>%
    mutate(ranking = rank(psuedo_density,ties.method = "min")) #spelling error... :(

  assertthat::assert_that(all(!is.na(pseudo_density_df$psuedo_density)),
                          msg = paste("internal error in",
                                      "distance_psuedo_density_function",
                                      "function's sigma selection."))

  proportion_points_not_included <- 1 - delta_prop

  top_points <- simulations_group_df_inner %>%
    left_join(pseudo_density_df, by = group_names) %>%
    mutate(keep = ranking > ceiling(proportion_points_not_included*nrow(simulations_group_df_inner))) %>%
    ungroup() %>% filter(keep) %>%
    select(one_of(data_column_names))

  mm_delta <- get_delta_nn(top_points)

  conformal_band_points <- simulations_group_df_inner %>%
    left_join(pseudo_density_df, by = group_names) %>%
    mutate(keep = ranking > conformal_score_cut) %>%
    ungroup() %>% filter(keep) %>%
    select(one_of(data_column_names))


  return(list(conformal_band_points = conformal_band_points,
              mm_delta = mm_delta))
}




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
