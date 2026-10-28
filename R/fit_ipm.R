#' Fit IPM
#'
#' Fits an IPM for seal observations
#'
#' @param data List of data frames, or an object returned by
#'   [forecast_to_data()].
#' @param species Which species. Inferred for generated forecast data.
#' @param years Vector of consecutive years in increasing order for the state
#'   process. Inferred for generated forecast data.
#' @param prior_spec priors
#' @param method Stan method for fitting
#' @param init Initialisation function or draws passed to Stan. A previously
#'   fitted `seal_ipm` object can also be supplied. In that case, each chain is
#'   initialized from a randomly selected posterior draw, with missing parameter
#'   elements supplied by the initialization function.
#' @param prior_only Should only prior be used?
#' @param growth_rate_period Positive integer vector of population-growth look-back steps. Defaults to all possible lengths based on the number of years fitted.
#' @param ... passed to Stan method
#' @return fitted model
#' @export
fit_ipm <- function(
    data,
    species = "grey",
    years = NULL,
    prior_spec = NULL,
    method = c("sample", "pathfinder", "laplace", "optimize"),
    init = NULL,
    prior_only = FALSE,
    growth_rate_period = 1:length(years),
    ...
) {
    method <- match.arg(method)

    if (inherits(data, "seal_ipm_data")) {
        if (is.null(species)) {
            species <- data$species
        }
        if (is.null(years)) {
            years <- data$years
        }
        if (!identical(species, data$species)) {
            cli_abort(c(
                "{.var species} must match the generated forecast data.",
                "x" = "The generated forecast is for '{data$species}', but you specified '{species}'."
            ))
        }
    }

    years <- validate_state_years(years)

    if (!species %in% c("grey", "ringed")) {
        cli_abort(c(
            "{.var species} must be one of 'grey' or 'ringed'.",
            "x" = "you supplied '{species}'"
        ))
    }

    if (inherits(data, "seal_ipm_data") && !is.null(prior_spec)) {
        cli_abort(c(
            "{.var prior_spec} cannot be changed when fitting generated forecast data."
        ))
    }

    if (inherits(data, "seal_ipm_data")) {
        prior_spec <- data$priors
    } else if (is.null(prior_spec)) {
        prior_spec <- default_priors(species)
    }

    stan_data <- if (inherits(data, "seal_ipm_data")) {
        utils::modifyList(
            forecast_stan_data_prefix(data, years),
            list(prior_only = as.integer(prior_only))
        )
    } else {
        build_stan_data(
            species = species,
            data = data,
            years = years,
            prior_spec = prior_spec,
            prior_only = prior_only,
            growth_rate_period = if (is.null(growth_rate_period)) {
                1L
            } else {
                growth_rate_period
            }
        )
    }

    stan_data$growth_rate_period <- validate_growth_rate_lookback(
        if (
            is.null(growth_rate_period) &&
                !is.null(stan_data$growth_rate_period)
        ) {
            stan_data$growth_rate_period
        } else if (is.null(growth_rate_period)) {
            1L
        } else {
            growth_rate_period
        },
        stan_data$n_state_years
    )
    stan_data$n_growth_rate_period <- length(stan_data$growth_rate_period)

    model <- get_species_model(species)

    method_args <- list(...)
    n_init <- switch(
        method,
        sample = if (is.null(method_args$chains)) {
            4L
        } else {
            method_args$chains
        },
        pathfinder = if (is.null(method_args$num_paths)) {
            4L
        } else {
            method_args$num_paths
        },
        1L
    )

    if (inherits(init, "seal_ipm")) {
        if (!identical(init$species, species)) {
            cli_abort(c(
                "The previous fit supplied as `init` must use the same species.",
                "x" = "You supplied '{species}' but the init is from '{init$species}'."
            ))
        }
        init <- previous_fit_init(
            fit = init,
            n_state_years = stan_data$n_state_years,
            n_demo = 2L * stan_data$n_age_classes,
            prior_spec = prior_spec,
            n_init = n_init
        )
    } else if (is.null(init)) {
        init <- grey_init_fun(
            n_state_years = stan_data$n_state_years,
            n_demo = 12,
            prior_spec = prior_spec,
            n_init = n_init
        )
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


previous_fit_init <- function(fit, n_state_years, n_demo, prior_spec, n_init) {
    if (
        length(n_init) != 1L ||
            !is.numeric(n_init) ||
            is.na(n_init) ||
            !is.finite(n_init) ||
            n_init < 1 ||
            n_init != floor(n_init)
    ) {
        cli_abort("`n_init` must be a positive integer.")
    }

    draws <- posterior::as_draws_matrix(
        fit$fit$draws(format = "draws_matrix")
    )
    if (!nrow(draws)) {
        cli_abort("The previous fit contains no posterior draws.")
    }
    # use one posterior draw per chain, recycling only if the fit is too small
    draws <- posterior::resample_draws(draws, ndraws = n_init)

    new_inits <- lapply(seq_len(n_init), function(chain) {
        init <- grey_init_fun(
            n_state_years = n_state_years,
            n_demo = n_demo,
            prior_spec = prior_spec,
            n_init = 1
        )
        posterior::as_draws_matrix(posterior::as_draws_rvars(lapply(
            init,
            function(x) {
                dimensions <- if (is.null(dim(x))) length(x) else dim(x)
                if (length(x) == 1L && is.null(dim(x))) {
                    posterior::rvar(x)
                } else {
                    posterior::rvar(array(x, dim = c(1L, dimensions)))
                }
            }
        )))
    })
    fresh <- do.call(posterior::bind_draws, c(fresh, list(along = "draw")))
    required <- posterior::variables(fresh)
    shared <- intersect(required, posterior::variables(draws))
    missing <- setdiff(required, shared)
    if (!length(shared)) {
        return(fresh)
    }
    draws <- posterior::subset_draws(draws, variable = shared)
    if (length(missing)) {
        draws <- posterior::bind_draws(
            draws,
            posterior::subset_draws(fresh, variable = missing),
            along = "variable"
        )
    }
    posterior::subset_draws(draws, variable = required)
}


validate_growth_rate_lookback <- function(growth_rate_period, n_state_years) {
    if (
        !is.numeric(growth_rate_period) ||
            !length(growth_rate_period) ||
            anyNA(growth_rate_period) ||
            any(!is.finite(growth_rate_period)) ||
            any(growth_rate_period != as.integer(growth_rate_period)) ||
            any(growth_rate_period < 1L) ||
            any(growth_rate_period >= n_state_years)
    ) {
        cli_abort(
            "`growth_rate_period` must contain positive integer values smaller than the number of state years."
        )
    }
    as.integer(growth_rate_period)
}


get_species_model <- function(species, ...) {
    model <- instantiate::stan_package_model(
        name = paste0(species, "_seal_IPM"),
        package = "sealIPM"
    )

    return(model)
}


grey_init_fun <- function(
    n_state_years,
    n_demo,
    n_init,
    prior_spec = default_priors("grey")
) {
    out <- list(
        # Initial population size
        population_init_size = posterior::rvar(stats::rlnorm(
            n_init,
            prior_spec$prior_initial_population_log_mean,
            prior_spec$prior_initial_population_log_sd
        )),

        # Natural mortality
        female_adult_survival_probability = posterior::rvar(stats::runif(
            n_init,
            0.9,
            0.99
        )),
        pup_to_adult_survival_ratio = posterior::rvar(stats::runif(
            n_init,
            0.8,
            0.99
        )),
        mortality_age_shape = posterior::rvar(
            stats::runif(n_init, 0, 1)
        ),

        male_pup_mortality_log_offset = posterior::rvar(
            stats::rnorm(
                n_init,
                prior_spec$prior_male_pup_mortality_log_offset_location,
                prior_spec$prior_male_pup_mortality_log_offset_scale
            )
        ),
        male_adult_mortality_log_offset = posterior::rvar(
            stats::rnorm(
                n_init,
                prior_spec$prior_male_adult_mortality_log_offset_location,
                prior_spec$prior_male_adult_mortality_log_offset_scale
            )
        ),

        # Hunting and bycatch selectivity
        hunting_selectivity_sweden_sc = posterior::rvar(
            array(
                stats::rnorm(
                    n_init * n_demo,
                    0,
                    0.2 * prior_spec$prior_hunting_selectivity_sd
                ),
                dim = c(n_init, n_demo)
            )
        ),
        hunting_selectivity_finland_sc = posterior::rvar(
            array(
                stats::rnorm(
                    n_init * n_demo,
                    0,
                    0.2 * prior_spec$prior_hunting_selectivity_sd
                ),
                dim = c(n_init, n_demo)
            )
        ),
        bycatch_selectivity_sc = posterior::rvar(
            array(
                stats::rnorm(
                    n_init * n_demo,
                    0,
                    0.2 * prior_spec$prior_bycatch_selectivity_sd
                ),
                dim = c(n_init, n_demo)
            )
        ),

        # Hunting effort
        hunting_effort_sd_sweden = posterior::rvar(
            stats::runif(n_init, 0.5, 1) *
                prior_spec$prior_hunting_effort_sd_scale
        ),
        hunting_effort_sd_finland = posterior::rvar(
            stats::runif(n_init, 0.5, 1) *
                prior_spec$prior_hunting_effort_sd_scale
        ),

        hunting_effort_noise_sweden = posterior::rvar(
            array(
                abs(stats::rnorm(
                    n_init * n_state_years,
                    0,
                    1
                )),
                dim = c(n_init, n_state_years)
            )
        ),
        hunting_effort_noise_finland = posterior::rvar(
            array(
                abs(stats::rnorm(
                    n_init * n_state_years,
                    0,
                    1
                )),
                dim = c(n_init, n_state_years)
            )
        ),

        # Birth-rate model
        birth_rate_baseline_max = posterior::rvar(stats::runif(
            n_init,
            0.85,
            0.95
        )),
        birth_rate_baseline_min_max_ratio = posterior::rvar(stats::runif(
            n_init,
            0.7,
            0.95
        )),
        birth_rate_at_capacity_to_baseline_ratio = posterior::rvar(rep(
            0.5,
            times = n_init
        )),

        herring_birth_rate_midpoint = posterior::rvar(
            -stats::rnorm(n_init, 0, 0.5)
        ),
        herring_slope = posterior::rvar(stats::rnorm(n_init, 1, 0.5)),
        herring_weight = posterior::rvar(stats::runif(n_init, 0, 1)),

        density_dependence_scaled = posterior::rvar(stats::runif(
            n_init,
            0.1,
            0.8
        )),

        # Carrying capacity
        carrying_capacity = posterior::rvar(stats::rlnorm(
            n_init,
            prior_spec$prior_carrying_capacity_log_mean,
            prior_spec$prior_carrying_capacity_log_sd
        )),

        # Aerial observation
        aerial_count_mu = posterior::rvar(stats::rbeta(
            n_init,
            prior_spec$prior_aerial_detection_alpha,
            prior_spec$prior_aerial_detection_beta
        )),
        aerial_count_overdispersion = posterior::rvar(stats::rlnorm(
            n_init,
            prior_spec$prior_aerial_overdispersion_log_mean,
            prior_spec$prior_aerial_overdispersion_log_sd
        )),

        # Reproductive-sign observation
        ca_probability_without_birth = posterior::rvar(stats::runif(
            n_init,
            0,
            1
        )),
        report_placental_mean = posterior::rvar(stats::runif(n_init, 0, 1)),
        report_ca_mean = posterior::rvar(stats::runif(n_init, 0, 1)),

        report_ca_sd = posterior::rvar(abs(stats::rnorm(n_init, 0, 0.1))),
        report_placental_sd = posterior::rvar(abs(stats::rnorm(
            n_init,
            0,
            0.1
        ))),

        ca_detection_noise = posterior::rvar(
            array(
                abs(
                    stats::rnorm(
                        n_init * n_state_years,
                        0,
                        1
                    )
                ),
                dim = c(n_init, n_state_years)
            )
        ),
        placental_scar_detection_noise = posterior::rvar(
            array(
                abs(stats::rnorm(
                    n_init * n_state_years,
                    0,
                    1
                )),
                dim = c(n_init, n_state_years)
            )
        ),

        # State-process stochasticity
        birth_count_noise = posterior::rvar(
            array(
                stats::rnorm(
                    n_init * n_state_years,
                    0,
                    0.05
                ),
                dim = c(n_init, n_state_years)
            )
        ),
        pup_sex_allocation_noise = posterior::rvar(
            array(
                stats::rnorm(
                    n_state_years,
                    0,
                    0.05
                ),
                dim = c(n_init, n_state_years)
            ),
        ),

        transition_noise_raw = posterior::rvar(
            array(
                stats::rnorm(
                    n_init * (3 * n_demo) * n_state_years,
                    mean = 0,
                    sd = 0.05
                ),
                dim = c(n_init, 3 * n_demo, n_state_years)
            )
        )
    )

    return(posterior::as_draws_df(out))
}
