#' Forecast IPM
#'
#' Forecasts future years from a fitted IPM. When the original data are
#' supplied, observations in `future_years` are automatically included for
#' leave-future-out validation.
#'
#' @param fit A fitted IPM object.
#' @param data The original raw sealIPM data list.
#' @param future_years Consecutive calendar years to forecast.
#' @param species Either `"grey"` or `"ringed"`.
#' @param ... Passed to `generate_quantities()`.
#'
#' @return A CmdStan generated-quantities fit.
#' @export
forecast_ipm <- function(
    fit,
    data,
    future_years,
    species = fit$species,
    ...
) {
    if (!species %in% c("grey", "ringed")) {
        stop(
            "`species` must be one of 'grey' or 'ringed'.",
            call. = FALSE
        )
    }

    future_years <- sort(as.integer(future_years))

    if (
        length(future_years) < 1L ||
            anyNA(future_years) ||
            anyDuplicated(future_years) ||
            !identical(
                future_years,
                seq.int(
                    future_years[1L],
                    length.out = length(future_years)
                )
            )
    ) {
        stop(
            "`future_years` must contain consecutive, unique years.",
            call. = FALSE
        )
    }

    forecast_model <- get_species_forecast_model(species)

    future_data <- build_future_stan_data(
        species = species,
        data = data,
        future_years = future_years,
        past_data = fit$stan_data
    )

    forecast_variables <- names(
        forecast_model$variables()$parameters
    )

    past_draws <- fit$fit$draws(
        variables = forecast_variables,
        format = "draws_matrix"
    )

    forecast_model$generate_quantities(
        fitted_params = past_draws,
        data = utils::modifyList(
            fit$stan_data,
            future_data
        ),
        ...
    )
}


get_species_forecast_model <- function(species, ...) {
    stan_file <- system.file(
        "bin",
        "stan",
        paste0(species, "_seal_IPM_forecast.stan"),
        package = "sealIPM",
        mustWork = TRUE
    )

    instantiate::stan_package_model(
        name = paste0(species, "_seal_IPM_forecast"),
        package = "sealIPM",
        stan_file = stan_file
    )
}


build_future_stan_data <- function(
    species,
    data,
    future_years,
    past_data
) {
    switch(
        species,
        grey = build_grey_future_stan_data(
            data = data,
            future_years = future_years,
            past_data = past_data
        ),
        ringed = stop(
            "Ringed seal forecasting is not yet implemented.",
            call. = FALSE
        )
    )
}


build_grey_future_stan_data <- function(
    data,
    future_years,
    past_data
) {
    n_future_years <- length(future_years)

    herring <- build_grey_herring(
        herring = data$herring,
        years = future_years
    )

    # build_grey_herring() returns one additional lagged value.
    # The forecast model needs n_future_years + 1 values.
    herring_bp_gof <- utils::tail(
        herring$herring_index_baltic_proper_gulf_finland,
        n_future_years + 1L
    )

    herring_gob <- utils::tail(
        herring$herring_index_gulf_bothnia,
        n_future_years + 1L
    )

    quotas <- build_grey_hunting_quotas(
        hunting_quotas = data$hunting_quotas,
        years = future_years
    )

    observations <- build_grey_future_observations(
        data = data,
        future_years = future_years
    )

    sample_sizes <- build_grey_future_sample_sizes(
        observations = observations,
        n_future_years = n_future_years
    )

    c(
        list(
            n_future_years = as.integer(n_future_years),

            future_herring_index_baltic_proper_gulf_finland = as.numeric(
                herring_bp_gof
            ),

            future_herring_index_gulf_bothnia = as.numeric(herring_gob),

            future_hunting_quota_sweden = as.integer(
                quotas$hunting_quota_sweden
            ),

            future_hunting_quota_finland = as.integer(
                quotas$hunting_quota_finland
            )
        ),
        sample_sizes,
        observations
    )
}


build_grey_future_observations <- function(
    data,
    future_years
) {
    aerial <- build_grey_aerial_counts(
        aerial_counts = data$aerial_counts,
        years = future_years,
        offset = 0L
    )

    hunting_bags <- build_grey_hunting_bags(
        hunting_bags = data$hunting_bags,
        years = future_years
    )

    hunting_samples <- build_grey_hunting_samples(
        samples = data$samples,
        years = future_years
    )

    bycatch <- build_grey_bycatch(
        samples = data$samples,
        years = future_years
    )

    reproductive_signs <- build_grey_reproductive_signs(
        reproductive_signs = data$reproductive_signs,
        years = future_years
    )

    pregnancy <- build_grey_pregnancy_status(
        pregnancy = data$pregnancy_status,
        years = future_years
    )

    list(
        n_future_aerial = as.integer(
            aerial$n_aerial_years
        ),
        aerial_year = as.integer(
            aerial$aerial_year
        ),
        future_obs_aerial_count = as.integer(
            aerial$obs_aerial_count
        ),

        n_future_hunting_bag_sweden = as.integer(
            hunting_bags$n_hunting_bag_years_sweden
        ),
        hunting_bag_year_sweden = as.integer(
            hunting_bags$hunting_bag_year_sweden
        ),
        future_obs_hunting_bag_sweden = as.numeric(
            hunting_bags$obs_hunting_bag_sweden
        ),

        n_future_hunting_bag_finland = as.integer(
            hunting_bags$n_hunting_bag_years_finland
        ),
        hunting_bag_year_finland = as.integer(
            hunting_bags$hunting_bag_year_finland
        ),
        future_obs_hunting_bag_finland = as.numeric(
            hunting_bags$obs_hunting_bag_finland
        ),

        n_future_hunting_comp_sweden = as.integer(
            hunting_samples$n_hunting_comp_years_sweden
        ),
        hunting_comp_year_sweden = as.integer(
            hunting_samples$hunting_comp_year_sweden
        ),
        future_obs_hunting_comp_sweden = unname(
            hunting_samples$obs_hunting_comp_sweden
        ),

        n_future_hunting_comp_finland = as.integer(
            hunting_samples$n_hunting_comp_years_finland
        ),
        hunting_comp_year_finland = as.integer(
            hunting_samples$hunting_comp_year_finland
        ),
        future_obs_hunting_comp_finland = unname(
            hunting_samples$obs_hunting_comp_finland
        ),

        n_future_bycatch = as.integer(
            bycatch$n_bycatch_years
        ),
        bycatch_year = as.integer(
            bycatch$bycatch_comp_year
        ),
        future_obs_bycatch_comp = unname(
            bycatch$obs_bycatch_comp
        ),

        n_future_pregnancy = as.integer(
            pregnancy$n_pregnancy_years
        ),
        pregnancy_count_year = as.integer(
            pregnancy$pregnancy_count_year
        ),
        future_obs_pregnancy_count = as.integer(
            pregnancy$obs_pregnancy_count
        ),
        future_obs_pregnancy_sample_size = as.integer(
            pregnancy$pregnancy_sample_size
        ),

        n_future_reproductive_signs = as.integer(
            reproductive_signs$n_reproductive_years
        ),
        reproductive_signs_year = as.integer(
            reproductive_signs$reproductive_signs_year
        ),
        future_obs_reproductive_signs_finland = unname(
            reproductive_signs$obs_reproductive_signs_finland
        )
    )
}


build_grey_future_sample_sizes <- function(
    observations,
    n_future_years,
    default_sample_size = 100L
) {
    fill_sample_sizes <- function(
        observation_years,
        sample_sizes
    ) {
        out <- rep.int(
            as.integer(default_sample_size),
            n_future_years
        )

        if (length(observation_years) > 0L) {
            out[observation_years] <- as.integer(sample_sizes)
        }

        out
    }

    hunting_sweden_sample_size <- if (
        observations$n_future_hunting_comp_sweden > 0L
    ) {
        rowSums(
            observations$future_obs_hunting_comp_sweden
        )
    } else {
        integer()
    }

    hunting_finland_sample_size <- if (
        observations$n_future_hunting_comp_finland > 0L
    ) {
        rowSums(
            observations$future_obs_hunting_comp_finland
        )
    } else {
        integer()
    }

    bycatch_sample_size <- if (observations$n_future_bycatch > 0L) {
        rowSums(
            observations$future_obs_bycatch_comp
        )
    } else {
        integer()
    }

    reproductive_signs_sample_size <- if (
        observations$n_future_reproductive_signs > 0L
    ) {
        rowSums(
            observations$future_obs_reproductive_signs_finland
        )
    } else {
        integer()
    }

    list(
        future_hunting_sample_size_sweden = fill_sample_sizes(
            observations$hunting_comp_year_sweden,
            hunting_sweden_sample_size
        ),

        future_hunting_sample_size_finland = fill_sample_sizes(
            observations$hunting_comp_year_finland,
            hunting_finland_sample_size
        ),

        future_bycatch_sample_size = fill_sample_sizes(
            observations$bycatch_year,
            bycatch_sample_size
        ),

        future_reproductive_signs_sample_size = fill_sample_sizes(
            observations$reproductive_signs_year,
            reproductive_signs_sample_size
        ),

        future_pregnancy_sample_size = fill_sample_sizes(
            observations$pregnancy_count_year,
            observations$future_obs_pregnancy_sample_size
        )
    )
}
