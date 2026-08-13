default_priors <- function(species) {
  n_demo <- 12
  C <- matrix(0, nrow = n_demo, ncol = n_demo)
  C[2:5, 2:5] <- 0.95
  C[8:11, 8:11] <- 0.95
  diag(C) <- 1
  bias_cholesky_factor <- t(chol(C)) #Variance for multivariate Gaussian prior for biase
  
  if (species == "grey") {

    out <- list(
      population_init = rep(1/n_demo, n_demo) * 20000,
      n_demo = n_demo,
      t_mate_to_preg = 1 / 36,
      t_birth_to_start_hunt = 1.5 / 12,
      t_birth_to_end_hunt = 8 / 12,
      t_hunt = 5.5 / 12,
      
      bias_cholesky_factor = bias_cholesky_factor,

      prior_initial_population_log_mean = 9.8,
      prior_initial_population_log_sd = 0.1,
      prior_carrying_capacity_log_mean = 11.3,
      prior_carrying_capacity_log_sd = 0.3,

      prior_aerial_detection_alpha = 32,
      prior_aerial_detection_beta = 9,
      prior_aerial_overdispersion_log_mean = 5.3,
      prior_aerial_overdispersion_log_sd = 1,

      prior_hunting_selectivity_sd = 0.5,
      prior_bycatch_selectivity_sd = 0.5,

      prior_hunting_effort_sd_location = 0,
      prior_hunting_effort_sd_scale = 0.1,

      prior_herring_intercept_scaled_sd = 4,
      prior_herring_slope_sd = 3,

      prior_male_pup_mortality_offset_location = 0,
      prior_male_pup_mortality_offset_scale = 0.2,
      prior_male_adult_mortality_offset_location = 0.88,
      prior_male_adult_mortality_offset_scale = 0.2,

      prior_reproductive_process_sd_location = 0,
      prior_reproductive_process_sd_scale = 0.1,

      harvest_bag_cv = 0.05,

      population_burn_in = 20,
      rel_tol = 1e-6,
      abs_tol = 1e-6,
      max_num_steps = 1000
    )
  }

  return(out)
}
