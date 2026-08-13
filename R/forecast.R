#' Forecast IPM
#'
#' Forecasts future years from a fitted IPM
#' @param fit fit object
#' @param future_data future data
#' @param species which species
#' @param prior_spec prior settings
#' @param ... unused
#' @return forecast results
#' @export
forecast_ipm <- function(
    fit,
    future_data,
    species,
    prior_spec = NULL,
    ...
) {
    if (!species %in% c("grey", "ringed")) {
        stop(
            "`species` must be one of 'grey' or 'ringed'.",
            call. = FALSE
        )
    }

    if (is.null(prior_spec)) {
        prior_spec <- default_priors(species)
    }

    forecast_model <- get_species_forecast_model(species)

    past_data <- fit$stan_data

    future_data <- build_future_data(future_data, past_data, species)

    forecast_variables <- names(forecast_model$variables()$parameters)

    past_draws <- fit$fit$draws(
        variables = forecast_variables,
        format = "draws_matrix"
    )

    future_draws <- forecast_model$generate_quantities(
        fitted_params = past_draws,
        data = c(future_data, past_data)
    )

    return(future_draws)
}

get_species_forecast_model <- function(species, ...) {
    stan_file <- system.file(
        "bin",
        "stan",
        paste0(species, "_seal_IPM_forecast.stan"),
        package = "sealIPM",
        mustWork = TRUE
    )

    model <- instantiate::stan_package_model(
        name = paste0(species, "_seal_IPM_forecast"),
        package = "sealIPM",
        stan_file = stan_file
    )

    return(model)
}

build_future_data <- function(future_data, past_data, species) {
    switch(
        species,
        grey = build_grey_future_data(future_data, past_data)
    )
}

build_grey_future_data <- function(future_data, past_data) {
    n_future_years <- nrow(future_data)

    last_past_herring_baltic_proper_gulf_finland <- past_data$herring_index_baltic_proper_gulf_finland[length(
        past_data$herring_index_baltic_proper_gulf_finland
    )]

    last_past_herring_gulf_bothnia <- past_data$herring_index_gulf_bothnia[length(
        past_data$herring_index_gulf_bothnia
    )]

    out <- list(
        n_future_years = n_future_years,
        future_herring_index_baltic_proper_gulf_finland = c(
            last_past_herring_baltic_proper_gulf_finland,
            future_data$Herring_GoF_BP
        ),
        future_herring_index_gulf_bothnia = c(
            last_past_herring_gulf_bothnia,
            future_data$Herring_GoB
        ),
        future_hunting_quota_sweden = future_data$Sweden_Quota,
        future_hunting_quota_finland = future_data$Finland_Quota,
        future_hunting_sample_size = rep(100, times = n_future_years),
        future_bycatch_sample_size = rep(100, times = n_future_years),
        future_reproductive_signs_sample_size = rep(
            100,
            times = n_future_years
        ),
        future_pregnancy_sample_size = rep(100, times = n_future_years)
    )
    return(out)
}
