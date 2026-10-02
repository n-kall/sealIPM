default_priors <- function(species) {
    n_demo <- 12
    C <- matrix(0, nrow = n_demo, ncol = n_demo)
    C[2:5, 2:5] <- 0.95
    C[8:11, 8:11] <- 0.95
    diag(C) <- 1
    selectivity_cholesky <- t(chol(C)) # Cholesky factor for the multivariate Gaussian selectivity prior

    if (species == "grey") {
        out <- list(
            population_init = rep(1 / n_demo, n_demo) * 20000,
            pregnancy_exposure_scaled = 1 / 36,
            # Shared hunting onset and duration in years
            t_birth_to_start_hunt = 1.5 / 12,
            hunting_duration = 8 / 12,

            selectivity_cholesky = selectivity_cholesky,

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

            # Midpoint in standardized weighted herring-index units
            prior_herring_birth_rate_midpoint_sd = 4,
            prior_herring_slope_sd = 3,

            # Positive male log mortality offsets reduce survival
            prior_male_pup_mortality_log_offset_location = 0,
            prior_male_pup_mortality_log_offset_scale = 0.2,
            prior_male_adult_mortality_log_offset_location = 0.88,
            prior_male_adult_mortality_log_offset_scale = 0.2,

            hunting_bag_cv = 0.05,

            population_burn_in = 20,
            rel_tol = 1e-6,
            abs_tol = 1e-8,
            fate_probability_tolerance = 1e-6,
            max_num_steps = 1000
        )
    }

    return(out)
}
