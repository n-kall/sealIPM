#' Convert a forecast draw to data for a new IPM fit
#'
#' Selects one posterior-predictive draw from [forecast_ipm()] and appends its
#' simulated observations and scenario covariates to the data used for the
#' original fit. The result can be passed directly to [fit_ipm()].
#'
#' @param forecast A generated-quantities object returned by [forecast_ipm()].
#' @param draw Positive integer identifying the posterior-predictive draw to
#'   use.
#' @return An object of class `seal_ipm_data`, suitable for `fit_ipm()`.
#' @export
forecast_to_data <- function(forecast, draw = 1L) {
    info <- attr(forecast, "seal_ipm_forecast", exact = TRUE)
    if (is.null(info)) {
        cli_abort(
            "`forecast` must be an object returned by forecast_ipm()."
        )
    }
    if (
        !is.numeric(draw) ||
            length(draw) != 1L ||
            is.na(draw) ||
            draw < 1 ||
            draw != as.integer(draw)
    ) {
        cli_abort("`draw` must be one positive integer.")
    }
    if (!identical(info$species, "grey")) {
        cli_abort(
            "Converting forecasts is currently implemented only for grey seals.",
        )
    }

    variables <- c(
        "future_aerial_count",
        "future_hunting_bags_sweden",
        "future_hunting_bags_finland",
        "future_hunting_comp_sweden",
        "future_hunting_comp_finland",
        "future_bycatch_comp",
        "future_pregnancy",
        "future_reproductive_signs"
    )
    draws <- forecast$draws(variables = variables, format = "draws_matrix")
    draw <- as.integer(draw)
    if (draw > nrow(draws)) {
        stop(
            sprintf("`draw` must be between 1 and %d.", nrow(draws)),
            call. = FALSE
        )
    }

    n_future <- length(info$future_years)
    n_demo <- 2L * info$fitted_stan_data$n_age_classes
    vector_draw <- function(variable, n) {
        cols <- sprintf("%s[%d]", variable, seq_len(n))
        as.numeric(draws[draw, cols, drop = TRUE])
    }
    matrix_draw <- function(variable, nr, nc) {
        cols <- as.vector(outer(
            seq_len(nr),
            seq_len(nc),
            function(i, j) sprintf("%s[%d,%d]", variable, i, j)
        ))
        matrix(as.numeric(draws[draw, cols, drop = TRUE]), nrow = nr, ncol = nc)
    }

    generated <- list(
        aerial_count = as.integer(vector_draw("future_aerial_count", n_future)),
        hunting_bags_sweden = vector_draw(
            "future_hunting_bags_sweden",
            n_future
        ),
        hunting_bags_finland = vector_draw(
            "future_hunting_bags_finland",
            n_future
        ),
        hunting_comp_sweden = matrix_draw(
            "future_hunting_comp_sweden",
            n_future,
            n_demo
        ),
        hunting_comp_finland = matrix_draw(
            "future_hunting_comp_finland",
            n_future,
            n_demo
        ),
        bycatch_comp = matrix_draw("future_bycatch_comp", n_future, n_demo),
        pregnancy = as.integer(vector_draw("future_pregnancy", n_future)),
        reproductive_signs = matrix_draw(
            "future_reproductive_signs",
            n_future,
            4L
        )
    )
    generated$hunting_comp_sweden <- matrix(
        as.integer(generated$hunting_comp_sweden),
        ncol = n_demo
    )
    generated$hunting_comp_finland <- matrix(
        as.integer(generated$hunting_comp_finland),
        ncol = n_demo
    )
    generated$bycatch_comp <- matrix(
        as.integer(generated$bycatch_comp),
        ncol = n_demo
    )
    generated$reproductive_signs <- matrix(
        as.integer(generated$reproductive_signs),
        ncol = 4L
    )

    stan_data <- append_grey_forecast_data(
        info$fitted_stan_data,
        info$future_stan_data,
        generated
    )
    structure(
        list(
            stan_data = stan_data,
            species = info$species,
            years = as.integer(c(info$fitted_years, info$future_years)),
            future_years = info$future_years,
            observations = generated,
            scenario_data = info$scenario_data,
            original_data = info$original_data,
            priors = info$priors,
            draw = draw
        ),
        class = "seal_ipm_data"
    )
}

append_grey_forecast_data <- function(past, future, generated) {
    old_n <- as.integer(past$n_state_years)
    n_future <- as.integer(future$n_future_years)
    out <- past
    out$n_state_years <- old_n + n_future
    out$herring_index_baltic_proper_gulf_finland <- c(
        past$herring_index_baltic_proper_gulf_finland,
        future$future_herring_index_baltic_proper_gulf_finland[-1L]
    )
    out$herring_index_gulf_bothnia <- c(
        past$herring_index_gulf_bothnia,
        future$future_herring_index_gulf_bothnia[-1L]
    )
    out$hunting_quota_sweden <- c(
        past$hunting_quota_sweden,
        future$future_hunting_quota_sweden
    )
    out$hunting_quota_finland <- c(
        past$hunting_quota_finland,
        future$future_hunting_quota_finland
    )

    append_observation <- function(count, year, value, generated_value) {
        out[[count]] <<- as.integer(out[[count]] + n_future)
        out[[year]] <<- as.integer(c(out[[year]], old_n + seq_len(n_future)))
        if (is.matrix(generated_value)) {
            out[[value]] <<- rbind(out[[value]], generated_value)
        } else {
            out[[value]] <<- c(out[[value]], generated_value)
        }
    }
    append_observation(
        "n_aerial_years",
        "aerial_year",
        "obs_aerial_count",
        generated$aerial_count
    )
    append_observation(
        "n_hunting_bag_years_sweden",
        "hunting_bag_year_sweden",
        "obs_hunting_bag_sweden",
        generated$hunting_bags_sweden
    )
    append_observation(
        "n_hunting_bag_years_finland",
        "hunting_bag_year_finland",
        "obs_hunting_bag_finland",
        generated$hunting_bags_finland
    )
    append_observation(
        "n_hunting_comp_years_sweden",
        "hunting_comp_year_sweden",
        "obs_hunting_comp_sweden",
        generated$hunting_comp_sweden
    )
    append_observation(
        "n_hunting_comp_years_finland",
        "hunting_comp_year_finland",
        "obs_hunting_comp_finland",
        generated$hunting_comp_finland
    )
    append_observation(
        "n_bycatch_years",
        "bycatch_comp_year",
        "obs_bycatch_comp",
        generated$bycatch_comp
    )

    out$n_pregnancy_years <- as.integer(out$n_pregnancy_years + n_future)
    out$pregnancy_count_year <- as.integer(c(
        out$pregnancy_count_year,
        old_n + seq_len(n_future)
    ))
    out$obs_pregnancy_count <- as.integer(c(
        out$obs_pregnancy_count,
        generated$pregnancy
    ))
    out$pregnancy_sample_size <- as.integer(c(
        out$pregnancy_sample_size,
        future$future_pregnancy_sample_size
    ))
    append_observation(
        "n_reproductive_years",
        "reproductive_signs_year",
        "obs_reproductive_signs_finland",
        generated$reproductive_signs
    )
    out
}

forecast_stan_data_prefix <- function(data, years) {
    available_years <- as.integer(data$years)
    years <- as.integer(years)
    if (
        length(years) > length(available_years) ||
            !identical(years, utils::head(available_years, length(years)))
    ) {
        cli_abort(
            "`years` must be a consecutive prefix of the generated forecast data years."
        )
    }

    out <- data$stan_data
    removed <- length(available_years) - length(years)
    n_state_years <- as.integer(out$n_state_years - removed)

    observation_fields <- list(
        c("n_aerial_years", "aerial_year", "obs_aerial_count"),
        c(
            "n_hunting_bag_years_finland",
            "hunting_bag_year_finland",
            "obs_hunting_bag_finland"
        ),
        c(
            "n_hunting_bag_years_sweden",
            "hunting_bag_year_sweden",
            "obs_hunting_bag_sweden"
        ),
        c(
            "n_hunting_comp_years_finland",
            "hunting_comp_year_finland",
            "obs_hunting_comp_finland"
        ),
        c(
            "n_hunting_comp_years_sweden",
            "hunting_comp_year_sweden",
            "obs_hunting_comp_sweden"
        ),
        c("n_bycatch_years", "bycatch_comp_year", "obs_bycatch_comp"),
        c(
            "n_pregnancy_years",
            "pregnancy_count_year",
            "obs_pregnancy_count",
            "pregnancy_sample_size"
        ),
        c(
            "n_reproductive_years",
            "reproductive_signs_year",
            "obs_reproductive_signs_finland"
        )
    )
    for (fields in observation_fields) {
        keep <- out[[fields[2L]]] <= n_state_years
        out[[fields[1L]]] <- as.integer(sum(keep))
        for (field in fields[-1L]) {
            value <- out[[field]]
            out[[field]] <- if (is.matrix(value)) {
                value[keep, , drop = FALSE]
            } else {
                value[keep]
            }
        }
    }

    out$n_state_years <- n_state_years
    out$hunting_quota_sweden <- out$hunting_quota_sweden[seq_len(n_state_years)]
    out$hunting_quota_finland <- out$hunting_quota_finland[seq_len(
        n_state_years
    )]
    out$herring_index_baltic_proper_gulf_finland <-
        out$herring_index_baltic_proper_gulf_finland[seq_len(
            n_state_years + 1L
        )]
    out$herring_index_gulf_bothnia <-
        out$herring_index_gulf_bothnia[seq_len(n_state_years + 1L)]
    out
}
