#' Fit IPM
#'
#' Fits an IPM for seal observations
#'
#' @param data list of dataframes
#' @param species which species
#' @param years Consecutive integer years in increasing order for the state process.
#' @param prior_spec priors
#' @param method Stan method
#' @param init Initialisation function or draws passed to Stan
#' @param prior_only Should only prior be used?
#' @param ... passed to Stan method
#' @return fitted model
#' @export
fit_ipm <- function(
    data,
    species,
    years,
    prior_spec = NULL,
    method = c("sample", "pathfinder", "laplace", "optimize"),
    init = NULL,
    prior_only = FALSE,
    ...
) {
    method <- match.arg(method)
    years <- validate_state_years(years)

    if (!species %in% c("grey", "ringed")) {
        stop(
            "`species` must be one of 'grey' or 'ringed'.",
            call. = FALSE
        )
    }

    if (is.null(prior_spec)) {
        prior_spec <- default_priors(species)
    }

    stan_data <- build_stan_data(
        species = species,
        data = data,
        years = years,
        prior_spec = prior_spec,
        prior_only = prior_only
    )

    model <- get_species_model(species)

    if (is.null(init)) {
        init <- function() {
            grey_init_fun(
                n_state_years = stan_data$n_state_years,
                n_demo = 12,
                prior_spec = prior_spec
            )
        }
    }

    fit <- switch(
        method,

        sample = model$sample(
            data = stan_data,
            init = init,
            ...
        ),

        pathfinder = model$pathfinder(
            data = stan_data,
            init = init,
            ...
        ),

        laplace = model$laplace(
            data = stan_data,
            init = init,
            ...
        ),

        optimize = model$optimize(
            data = stan_data,
            init = init,
            ...
        )
    )

    structure(
        list(
            fit = fit,
            species = species,
            data = data,
            years = as.integer(years),
            stan_data = stan_data,
            priors = prior_spec,
            method = method
        ),
        class = c(
            paste0(species, "_seal_ipm"),
            "seal_ipm"
        )
    )
}


get_species_model <- function(species, ...) {
    model <- instantiate::stan_package_model(
        name = paste0(species, "_seal_IPM"),
        package = "sealIPM"
    )

    return(model)
}


grey_init_fun <- function(n_state_years, n_demo, prior_spec = default_priors("grey")) {
    transition_noise_raw <- matrix(
        stats::rnorm(3 * n_demo * n_state_years, 0, 0.05),
        nrow = 3 * n_demo,
        ncol = n_state_years
    )

    out <- list(
        # Initial population size
        population_init_size = stats::rlnorm(
            1, prior_spec$prior_initial_population_log_mean,
            prior_spec$prior_initial_population_log_sd
        ),

        # Natural mortality
        female_adult_survival_probability = stats::runif(1, 0.9, 0.99),
        pup_to_adult_survival_ratio = stats::runif(1, 0.8, 0.99),
        mortality_age_shape = stats::runif(1, 0, 1),

        male_pup_mortality_log_offset = stats::rnorm(
            1, prior_spec$prior_male_pup_mortality_log_offset_location,
            prior_spec$prior_male_pup_mortality_log_offset_scale
        ),
        male_adult_mortality_log_offset = stats::rnorm(
            1, prior_spec$prior_male_adult_mortality_log_offset_location,
            prior_spec$prior_male_adult_mortality_log_offset_scale
        ),

        # Hunting and bycatch selectivity
        hunting_selectivity_sweden_sc = stats::rnorm(
            n_demo, 0, 0.2 * prior_spec$prior_hunting_selectivity_sd
        ),
        hunting_selectivity_finland_sc = stats::rnorm(
            n_demo, 0, 0.2 * prior_spec$prior_hunting_selectivity_sd
        ),
        bycatch_selectivity_sc = stats::rnorm(
            n_demo, 0, 0.2 * prior_spec$prior_bycatch_selectivity_sd
        ),

        # Hunting effort
        hunting_effort_sd_sweden = stats::runif(1, 0.5, 1) *
            prior_spec$prior_hunting_effort_sd_scale,
        hunting_effort_sd_finland = stats::runif(1, 0.5, 1) *
            prior_spec$prior_hunting_effort_sd_scale,

        hunting_effort_noise_sweden = abs(stats::rnorm(n_state_years, 0, 1)),
        hunting_effort_noise_finland = abs(stats::rnorm(n_state_years, 0, 1)),

        # Birth-rate model
        birth_rate_baseline_max = stats::runif(1, 0.85, 0.95),
        birth_rate_baseline_min_max_ratio = stats::runif(1, 0.7, 0.95),
        birth_rate_at_capacity_to_baseline_ratio = 0.5, # Set coherently below.

        # Negate the previous offset to preserve equivalent seeded starting points
        herring_birth_rate_midpoint = -stats::rnorm(1, 0, 0.5),
        herring_slope = stats::rnorm(1, 1, 0.5),
        herring_weight = stats::runif(1, 0, 1),

        density_dependence_scaled = stats::runif(1, 0.1, 0.8),

        # Carrying capacity
        carrying_capacity = stats::rlnorm(
            1, prior_spec$prior_carrying_capacity_log_mean,
            prior_spec$prior_carrying_capacity_log_sd
        ),

        # Aerial observation
        aerial_count_mu = stats::rbeta(
            1, prior_spec$prior_aerial_detection_alpha,
            prior_spec$prior_aerial_detection_beta
        ),
        aerial_count_overdispersion = stats::rlnorm(
            1, prior_spec$prior_aerial_overdispersion_log_mean,
            prior_spec$prior_aerial_overdispersion_log_sd
        ),

        # Reproductive-sign observation
        ca_probability_without_birth = stats::runif(1, 0, 1),
        report_placental_mean = stats::runif(1, 0, 1),
        report_ca_mean = stats::runif(1, 0, 1),

        report_ca_sd = abs(stats::rnorm(1, 0, 0.1)),
        report_placental_sd = abs(stats::rnorm(1, 0, 0.1)),

        ca_detection_noise = abs(stats::rnorm(n_state_years, 0, 1)),
        placental_scar_detection_noise = abs(stats::rnorm(n_state_years, 0, 1)),

        # State-process stochasticity
        birth_count_noise = stats::rnorm(n_state_years, 0, 0.05),
        pup_sex_allocation_noise = stats::rnorm(n_state_years, 0, 0.05),

        transition_noise_raw = transition_noise_raw
    )

    # Coordinate initialization with the existing soft Euler-Lotka constraint.
    # This is only a starting-point choice; the Stan prior is unchanged.
    n_age <- n_demo / 2
    mu_pup <- -log(out$female_adult_survival_probability * out$pup_to_adult_survival_ratio)
    mu_adult <- -log(out$female_adult_survival_probability)
    weight <- (seq.int(0, n_age - 2) / (n_age - 1))^out$mortality_age_shape
    mortality <- exp(log(mu_pup) + weight * (log(mu_adult) - log(mu_pup)))
    log_el_birth <- sum(mortality) + log(2) + log1p(-out$female_adult_survival_probability)
    reference_birth <- out$birth_rate_baseline_max * (
        out$birth_rate_baseline_min_max_ratio +
            (1 - out$birth_rate_baseline_min_max_ratio) *
                stats::plogis(-out$herring_birth_rate_midpoint * out$herring_slope)
    )
    for (attempt in seq_len(100L)) {
        log_prop <- log_el_birth - log(reference_birth) + stats::rnorm(1, 0, 0.01)
        if (is.finite(log_prop) && log_prop < 0 && exp(log_prop) > 0) {
            out$birth_rate_at_capacity_to_baseline_ratio <- exp(log_prop)
            return(out)
        }
    }
    stop("Could not initialize an interior Euler-Lotka birth-rate proportion.")
}
