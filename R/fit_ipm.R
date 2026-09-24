#' Fit IPM
#'
#' Fits an IPM for seal observations
#'
#' @param data list of dataframes
#' @param species which species
#' @param years years for the state process
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
            grey_init_fun(n_state_years = stan_data$n_state_years, n_demo = 12)
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


grey_init_fun <- function(n_state_years, n_demo) {
    transition_noise_raw <- matrix(
        stats::rnorm(3 * n_demo * n_state_years, 0, 1),
        nrow = 3 * n_demo,
        ncol = n_state_years
    )

    list(
        # Initial population size
        population_init_size = stats::rlnorm(1, 9.8, 0.1),

        # Natural mortality
        phi_a = stats::runif(1, 0.9, 0.99),
        phi_sc = stats::runif(1, 0.8, 0.99),
        survival_shape = stats::runif(1, 0, 1),

        male_pup_survival_offset = stats::rnorm(1, 0, 0.2),
        male_adult_survival_offset = stats::rnorm(1, 0.88, 0.2),

        # Hunting selectivity / bias
        hunting_selectivity_sweden_sc = stats::rnorm(n_demo, 0, 0.1),
        hunting_selectivity_finland_sc = stats::rnorm(n_demo, 0, 0.1),
        bycatch_selectivity_sc = stats::rnorm(n_demo, 0, 0.1),

        # Hunting effort
        hunting_effort_sd_sweden = abs(stats::rnorm(1, 0, 0.5)),
        hunting_effort_sd_finland = abs(stats::rnorm(1, 0, 0.5)),

        epsilon_h_sw = abs(stats::rnorm(n_state_years, 0, 1)),
        epsilon_h_fi = abs(stats::rnorm(n_state_years, 0, 1)),

        # Birth-rate model
        max_baseline_birth_rate = stats::runif(1, 0.85, 0.95),
        min_baseline_birth_rate_prop = stats::runif(1, 0.7, 0.95),
        birth_rate_at_carrying_capacity_prop = stats::runif(1, 0.2, 0.8),

        herring_intercept_scaled = stats::rnorm(1, 0, 0.5),
        herring_slope = stats::rnorm(1, 1, 0.5),
        herring_weight = stats::runif(1, 0, 1),

        density_dependence_scaled = stats::runif(1, 0.1, 0.8),

        # Carrying capacity
        carrying_capacity = stats::rlnorm(1, 11.3, 0.3),

        # Aerial observation
        aerial_count_mu = stats::rbeta(1, 32, 9),
        aerial_count_overdispersion = stats::rlnorm(1, 5.3, 1),

        # Reproductive-sign observation
        prob_of_ca = stats::runif(1, 0, 1),
        report_placental_mean = stats::runif(1, 0, 1),
        report_ca_mean = stats::runif(1, 0, 1),

        report_ca_sd = abs(stats::rnorm(1, 0, 0.1)),
        report_placental_sd = abs(stats::rnorm(1, 0, 0.1)),

        epsilon_ca = abs(stats::rnorm(n_state_years, 0, 1)),
        epsilon_placental = abs(stats::rnorm(n_state_years, 0, 1)),

        # State-process stochasticity
        epsilon_birth = stats::rnorm(n_state_years, 0, 1),
        epsilon_sex = stats::rnorm(n_state_years, 0, 1),

        transition_noise_raw = transition_noise_raw
    )
}
