#' Fit IPM
#'
#' Fits an IPM for seal observations
#'
#' @param data List of data frames, or an object returned by
#'   [forecast_to_data()].
#' @param species Which species. Inferred for generated forecast data.
#' @param years Consecutive integer years in increasing order for the state
#'   process. Inferred for generated forecast data. For data returned by
#'   [forecast_to_data()], this may be a consecutive prefix of the available
#'   historical and forecast years.
#' @param prior_spec priors
#' @param method Stan method
#' @param init Initialisation function or draws passed to Stan. A previously
#'   fitted `seal_ipm` object can also be supplied. In that case, each chain is
#'   initialized from a randomly selected posterior draw. Parameters shared by
#'   the fits are retained, while year-indexed parameters are extended with
#'   conservative random initial values when `years` is longer than in the
#'   previous fit.
#' @param prior_only Should only prior be used?
#' @param growth_rate_period Positive integer vector of population-growth look-back steps. For
#'   each value, generated quantities include the annual log-linear trend over
#'   that interval and its corresponding proportional growth rate. The trend is
#'   the weighted mean of annual log changes that is exactly equivalent to the
#'   ordinary least-squares slope of log population on year. When `NULL`,
#'   defaults to one step, or retains
#'   the values stored in generated forecast data.
#' @param ... passed to Stan method
#' @return fitted model
#' @export
fit_ipm <- function(
    data,
    species = NULL,
    years = NULL,
    prior_spec = NULL,
    method = c("sample", "pathfinder", "laplace", "optimize"),
    init = NULL,
    prior_only = FALSE,
    growth_rate_period = NULL,
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
            stop(
                "`species` must match the generated forecast data.",
                call. = FALSE
            )
        }
    }

    years <- validate_state_years(years)

    if (!species %in% c("grey", "ringed")) {
        stop(
            "`species` must be one of 'grey' or 'ringed'.",
            call. = FALSE
        )
    }

    if (inherits(data, "seal_ipm_data") && !is.null(prior_spec)) {
        stop(
            "`prior_spec` cannot be replaced when fitting generated forecast data.",
            call. = FALSE
        )
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

    if (inherits(init, "seal_ipm")) {
        if (!identical(init$species, species)) {
            stop(
                "The previous fit supplied as `init` must use the same species.",
                call. = FALSE
            )
        }
        init <- previous_fit_init(
            fit = init,
            n_state_years = stan_data$n_state_years,
            n_demo = 2L * stan_data$n_age_classes,
            prior_spec = prior_spec
        )
    } else if (is.null(init)) {
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


previous_fit_init <- function(fit, n_state_years, n_demo, prior_spec) {
    template <- grey_init_fun(
        n_state_years = n_state_years,
        n_demo = n_demo,
        prior_spec = prior_spec
    )
    parameter_names <- names(template)
    draws <- fit$fit$draws(
        variables = parameter_names,
        format = "draws_matrix"
    )
    if (!nrow(draws)) {
        stop("The previous fit contains no posterior draws.", call. = FALSE)
    }

    draw_order <- sample.int(nrow(draws))
    next_draw <- 0L

    function() {
        next_draw <<- next_draw + 1L
        if (next_draw > length(draw_order)) {
            draw_order <<- sample.int(nrow(draws))
            next_draw <<- 1L
        }
        previous <- unpack_stan_init_draw(
            draws[draw_order[next_draw], , drop = TRUE],
            parameter_names
        )
        current <- grey_init_fun(
            n_state_years = n_state_years,
            n_demo = n_demo,
            prior_spec = prior_spec
        )
        for (parameter in parameter_names) {
            current[[parameter]] <- copy_init_overlap(
                current[[parameter]],
                previous[[parameter]]
            )
        }
        current
    }
}


unpack_stan_init_draw <- function(draw, parameter_names) {
    draw_names <- names(draw)
    out <- stats::setNames(vector("list", length(parameter_names)), parameter_names)

    for (parameter in parameter_names) {
        columns <- which(
            draw_names == parameter |
                startsWith(draw_names, paste0(parameter, "["))
        )
        if (!length(columns)) {
            stop(
                "The previous fit does not contain parameter `", parameter, "`.",
                call. = FALSE
            )
        }
        labels <- draw_names[columns]
        values <- as.numeric(draw[columns])
        if (length(columns) == 1L && identical(labels, parameter)) {
            out[[parameter]] <- values
            next
        }

        indices <- lapply(labels, function(label) {
            inside <- sub("^.*\\[", "", label)
            inside <- sub("\\]$", "", inside)
            as.integer(strsplit(inside, ",", fixed = TRUE)[[1L]])
        })
        n_dimensions <- unique(lengths(indices))
        if (length(n_dimensions) != 1L || anyNA(unlist(indices))) {
            stop(
                "Cannot decode posterior draws for parameter `", parameter, "`.",
                call. = FALSE
            )
        }
        dimensions <- vapply(
            seq_len(n_dimensions),
            function(dimension) max(vapply(indices, `[[`, integer(1L), dimension)),
            integer(1L)
        )
        value <- array(NA_real_, dim = dimensions)
        for (i in seq_along(indices)) {
            value <- do.call(
                `[<-`,
                c(list(value), as.list(indices[[i]]), list(value = values[i]))
            )
        }
        out[[parameter]] <- if (length(dimensions) == 1L) {
            as.numeric(value)
        } else {
            value
        }
    }
    out
}


copy_init_overlap <- function(current, previous) {
    current_dim <- dim(current)
    previous_dim <- dim(previous)

    if (is.null(current_dim) && is.null(previous_dim)) {
        if (length(current) == 1L && length(previous) == 1L) {
            return(as.numeric(previous))
        }
        n_copy <- min(length(current), length(previous))
        current[seq_len(n_copy)] <- previous[seq_len(n_copy)]
        return(current)
    }
    if (is.null(current_dim) || is.null(previous_dim) ||
        length(current_dim) != length(previous_dim)) {
        stop("Previous and current parameter dimensions are incompatible.", call. = FALSE)
    }

    overlap <- pmin(current_dim, previous_dim)
    indices <- lapply(overlap, seq_len)
    previous_values <- do.call(
        `[`,
        c(list(previous), indices, list(drop = FALSE))
    )
    do.call(
        `[<-`,
        c(list(current), indices, list(value = previous_values))
    )
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
        stop(
            "`growth_rate_period` must contain positive integer look-backs smaller than the number of state years.",
            call. = FALSE
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
    prior_spec = default_priors("grey")
) {
    transition_noise_raw <- matrix(
        stats::rnorm(3 * n_demo * n_state_years, 0, 0.05),
        nrow = 3 * n_demo,
        ncol = n_state_years
    )

    out <- list(
        # Initial population size
        population_init_size = stats::rlnorm(
            1,
            prior_spec$prior_initial_population_log_mean,
            prior_spec$prior_initial_population_log_sd
        ),

        # Natural mortality
        female_adult_survival_probability = stats::runif(1, 0.9, 0.99),
        pup_to_adult_survival_ratio = stats::runif(1, 0.8, 0.99),
        mortality_age_shape = stats::runif(1, 0, 1),

        male_pup_mortality_log_offset = stats::rnorm(
            1,
            prior_spec$prior_male_pup_mortality_log_offset_location,
            prior_spec$prior_male_pup_mortality_log_offset_scale
        ),
        male_adult_mortality_log_offset = stats::rnorm(
            1,
            prior_spec$prior_male_adult_mortality_log_offset_location,
            prior_spec$prior_male_adult_mortality_log_offset_scale
        ),

        # Hunting and bycatch selectivity
        hunting_selectivity_sweden_sc = stats::rnorm(
            n_demo,
            0,
            0.2 * prior_spec$prior_hunting_selectivity_sd
        ),
        hunting_selectivity_finland_sc = stats::rnorm(
            n_demo,
            0,
            0.2 * prior_spec$prior_hunting_selectivity_sd
        ),
        bycatch_selectivity_sc = stats::rnorm(
            n_demo,
            0,
            0.2 * prior_spec$prior_bycatch_selectivity_sd
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
            1,
            prior_spec$prior_carrying_capacity_log_mean,
            prior_spec$prior_carrying_capacity_log_sd
        ),

        # Aerial observation
        aerial_count_mu = stats::rbeta(
            1,
            prior_spec$prior_aerial_detection_alpha,
            prior_spec$prior_aerial_detection_beta
        ),
        aerial_count_overdispersion = stats::rlnorm(
            1,
            prior_spec$prior_aerial_overdispersion_log_mean,
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
    mu_pup <- -log(
        out$female_adult_survival_probability * out$pup_to_adult_survival_ratio
    )
    mu_adult <- -log(out$female_adult_survival_probability)
    weight <- (seq.int(0, n_age - 2) / (n_age - 1))^out$mortality_age_shape
    mortality <- exp(log(mu_pup) + weight * (log(mu_adult) - log(mu_pup)))
    log_el_birth <- sum(mortality) +
        log(2) +
        log1p(-out$female_adult_survival_probability)
    reference_birth <- out$birth_rate_baseline_max *
        (out$birth_rate_baseline_min_max_ratio +
            (1 - out$birth_rate_baseline_min_max_ratio) *
                stats::plogis(
                    -out$herring_birth_rate_midpoint * out$herring_slope
                ))
    for (attempt in seq_len(100L)) {
        log_prop <- log_el_birth -
            log(reference_birth) +
            stats::rnorm(1, 0, 0.01)
        if (is.finite(log_prop) && log_prop < 0 && exp(log_prop) > 0) {
            out$birth_rate_at_capacity_to_baseline_ratio <- exp(log_prop)
            return(out)
        }
    }
    stop("Could not initialize an interior Euler-Lotka birth-rate proportion.")
}
