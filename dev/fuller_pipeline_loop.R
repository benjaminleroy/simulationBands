# data processing for conformal_pdf analysis

### POTENTIAL OPTIONS
# a. vary between: local vs global conformal # more a worry with cde approach - let's focus on doing 'mass' non-conformal score
# b. data options: lei-wasserman: perfect fit vs "y" poor fit
# c. conformal score: cde vs mass (maybe sim based -or- not?) ## let's not focus on sim vs not currently
# d. scaling attempted during optimization


# Libraries ----------------------------------

library(shiny)
library(flexdashboard)
library(tidyverse)
library(ks)
library(latex2exp)
library(shinyWidgets)
library(tidyverse)
library(pracma)
library(gridExtra)
library(latex2exp)
library(flexmix)

devtools::load_all(path = "../../") # simulationBands

### Defaults
theme_set(theme_minimal() +
            theme(text=element_text(size=16),
                  aspect.ratio = 1))

n <- 2000
num_sim_y <- 200

vis_bool <- FALSE


for (data_idx in 1:3){
  data_name_full <- c("Lei & Wasserman (2014) - perfect fit",
                      "Y data - bad fit", "Y2 data - really bad fit")[data_idx]
  data_name <- c("lw","y", "y2")[data_idx]
 for (conformal_approach in c("cde", "mass")) {
   for (conformal_mass_perfect in c(T, F)){
     if(conformal_mass_perfect & conformal_approach == "cde"){
     } else {
       for (log10_lambda in c(-Inf,-5,-4,-3, -2,-1,0,1)) {
         for (scale_optimization in c(T,F)){
           for (prime in c(F, T)){
             for (conditional in c(F, T)){








# Data generation ----------------------------

if (data_name == "lw"){
  data_all <- lei_wasserman_data(n = n)
  if (conditional){
    Xx <- runif(n = n, min = .65, max = .85)
    data_calibrate_list <- lei_wasserman_data_conditional_simulate(x = Xx,
                                                                    n = 1,
                                                                    verbose = F)
    data_calibrate <- do.call(rbind, data_calibrate_list) %>%
      rename("y"= "sim")

  } else {
    data_calibrate <- lei_wasserman_data(n = n)
  }
} else { #if (data_name == "y" or "y2") {
  data_all <- generate_split_mixture(n = n)
  if (conditional){
    Xx <- runif(n = n, min = 9, max = 10)

    Yy <- gen_mix_reg(x = Xx)
    data_calibrate <- data.frame(x = Xx, y = Yy)
  } else {
    data_calibrate <- generate_split_mixture(n = n)
  }
}

}

if (vis_bool){
  ggplot(data_all) +
    geom_point(aes(x = x, y = y)) +
    labs(title = "True distribution (training data)",
         subtitle = data_name_full)
}

# Fit creation ------------------------------


if (data_name == "y") {
  model <- flexmix(y~x, k = 2, data = data_all)

  beta <- sapply(model@components, function(item) item[[1]]@parameters$coef[2])
  intercept <- sapply(model@components, function(item) item[[1]]@parameters$coef[1])
  sigma <- sapply(model@components, function(item) item[[1]]@parameters$sigma)
  proportion <- model@prior

  # get conformal scores from simulations
  model_params <- list(beta = beta, intercept = intercept,
                       sd = sigma, prop = proportion)
} else if (data_name == "y2") {
  beta <- c(1,1)
  intercept <- c(-5,-5)
  sigma <- c(1,1)
  proportion <- c(.5,.5)
  model_params <- list(beta = beta, intercept = intercept,
                       sd = sigma, prop = proportion)
}

# Calibration set's conformal score ---------------------------

if (data_name %in% c("y", "y2")){
  conformal_score_info <- data_calibrate %>%
    dplyr::mutate(id = 1:nrow(data_calibrate)) %>%
    dplyr::group_by(id) %>%
    tidyr::nest() %>%
    dplyr::mutate(sim_cs = purrr::map(data, function(df) {
      simulation_conformal_scores_mix_reg(df$x, df$y,
                                          model_params, sim = num_sim_y)}),
      conformal_score_true_density = purrr::map(data, function(df) {
        multimode_density_f(df$y, df$x,
                            beta = beta,
                            intercept = intercept,
                            sd = sigma,
                            proportions = proportion)}
      )) %>%
    dplyr::select(-data) %>%
    unnest(col = c(sim_cs, conformal_score_true_density))

  conformal_calibration_dist_cde <- conformal_score_info$empirical_cde
  conformal_calibration_dist_mass <- conformal_score_info$non_conformal_mass_ecde_and_sim_based
} else { #data_name == "lw"
  conformal_score_info <- data_calibrate %>%
    dplyr::mutate(id = 1:nrow(data_calibrate)) %>%
    dplyr::group_by(id) %>%
    tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(df) {
      conformal_scores_lei_wasserman(df$x, df$y)
      })) %>%
    tidyr::unnest(data)

  conformal_calibration_dist_cde <- conformal_score_info$cde
  conformal_calibration_dist_mass <- conformal_score_info$mass
}

if (conformal_mass_perfect){
  n_conformal <- conformal_calibration_dist_mass
  conformal_calibration_dist_mass <- 1:n/n
}

# true x selection ----------------------------

if (data_name == "y"){
  x <- 9.5
} else { #data_name == "lw"
  x <- .75
}

# building prediction region information ----------
if (data_name %in% c("y", "y2")) {
  delta_y <- .001 # delta_y needs to be passed through
  yy <- seq(-15, 15, by = delta_y)
  first_df <- data.frame(yy = yy,
                         true_density = multimodel_density_split_mixture(yy, x),
                         fit_density = multimode_density_f(yy, x,
                                                        beta = model_params$beta,
                                                        intercept = model_params$intercept,
                                                        sd = model_params$sd,
                                                        proportions = model_params$prop))
} else { # data_name == "lw"
  delta_y <- .001
  if (conformal_mass_perfect & log10_lambda == -Inf){
    delta_y <- .0001

  }
  yy <- seq(-8, 8, by = delta_y)
  first_df <- data.frame(yy = yy,
                         true_density = cde_lei_wassserman(x)(yy)) %>%
    mutate(fit_density = true_density)
}


if (conformal_approach == "cde") {
  first_df <- first_df %>% arrange(yy) %>%
    conformal_pdf_breaks(unique(conformal_calibration_dist_cde),
                                   conformal_score = "cde") %>%
    arrange(yy)
} else{ #conformal_approach == "mass"
  first_df <- first_df %>% arrange(yy) %>%
    conformal_pdf_breaks(unique(conformal_calibration_dist_mass),
                                   conformal_score = "mass") %>%
    arrange(yy)
}


second_df <- first_df %>%
  cumulative_comparisons(group_columns = c(g_id, g_id2, g_id3,
                                           cut_off_lower, cut_off_upper))

# conformal_pdf ---------------------------

prob_cde <- second_df$fit_prop#second_df$fit_prop_lag[-1]
n_conformal_pdf <- length(prob_cde) - 1

try({solution_qp_l2_lam <- stepwise_conformal_cde_update(n_conformal_pdf, prob_cde,
                              lambda = ifelse(!is.infinite(log10_lambda),
                                              10^(log10_lambda),
                                              -1),
                              prime = prime,
                              alpha = 1,
                              scaling_constraint = scale_optimization)



if (vis_bool){
  vis_stepwise_scaling_only(solution_qp_l2_lam)

  vis_stepwise_scaling_only(solution_qp_l2_lam,second_df$cut_off_lower)

  ggplot(second_df) +
    geom_histogram(aes(x = cut_off_lower, y = ..density..)) +
    geom_density(aes(x = cut_off_lower))

  vis_all_densities(stepwise_scaling = solution_qp_l2_lam,
                    individual_df_rich = first_df,
                    cumulative_comparisons_df = second_df,
                    scale = TRUE)
}


# saving ----------------------------------------------


save_string_name <- sprintf(
  paste0("dev/conformal_pdf_data/conformal_pdf_%s_%s_log10_lambda-%s",
         "_scale_optimized-%s_perfect_conformal_mass-%s_prime-%s",
         "_conditional-%s.Rdata"),
  data_name, conformal_approach,
  log10_lambda, scale_optimization,
  conformal_mass_perfect, prime,conditional)

print(save_string_name)

save(solution_qp_l2_lam, first_df, second_df, file = save_string_name)
})
          }
         }
       }
     }

   }
 }
 }


for (data_name in c("lw", "y", "y2")){
  if (data_name == "y"){
    x <- 9.5
  } else { #data_name == "lw"
    x <- .75
  }
  if (data_name == "lw"){
    data_all <- lei_wasserman_data(n = n)
    data_fit <- lei_wasserman_data(n = n)

  } else { # data_name == "y"
    data_all <- generate_split_mixture(n = n)
  }

  if (data_name == "y") {
    model <- flexmix(y~x, k = 2, data = data_all)

    beta <- sapply(model@components, function(item) item[[1]]@parameters$coef[2])
    intercept <- sapply(model@components, function(item) item[[1]]@parameters$coef[1])
    sigma <- sapply(model@components, function(item) item[[1]]@parameters$sigma)
    proportion <- model@prior

    data_fit <- generate_split_mixture(n = n,
                           right_beta = beta,
                           right_intercept = intercept,
                           right_sd = sigma,
                           right_proportions = proportion)
  } else if (data_name == "y2"){
    beta <- c(1,1)
    intercept <- c(-5,-5)
    sigma <- c(1,1)
    proportion <- c(.5,.5)

    data_fit <- generate_split_mixture(n = n,
                                       right_beta = beta,
                                       right_intercept = intercept,
                                       right_sd = sigma,
                                       right_proportions = proportion)
  }

  data_both <- rbind(data_all %>% mutate(color = "true distribution"),
                     data_fit %>% mutate(color = "fit distribution"))

  save(data_both, x, file = paste0("dev/conformal_pdf_data/", data_name, ".Rdata"))
}



