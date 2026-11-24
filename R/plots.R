#' Plot forecast
#'
#' @param forecast Forecast draws, preferably returned by [forecast_ipm()].
#' @param future_years Calendar years for the full forecast horizon. Defaults to
#'   the years recorded by [forecast_ipm()], or forecast-year indices otherwise.
#' @param variables Variables to plot. Process outputs use their final index
#'   for time, and predictive observation arrays use their first index. Other
#'   indices define facets. Growth-rate facets identify the lookback length and
#'   account for the first year with growth estimates.
#' @return A ggplot object.
#' @importFrom rlang .data
#' @export
forecast_plot <- function(forecast, future_years = NULL, variables = NULL) {
  info <- attr(forecast, "seal_ipm_forecast", exact = TRUE)
  if (is.null(future_years)) {
    future_years <- info$future_years
  }
  draws <- posterior::as_draws(forecast)
  all_draws <- draws

  if (!is.null(variables)) {
    draws <- posterior::subset_draws(draws, variable = variables)
  }

  summary <- posterior::summarise_draws(draws, mean, posterior::quantile2)
  summary <- summary[grepl("\\[[0-9]+(,[0-9]+)*\\]$", summary$variable), ]
  if (nrow(summary) == 0L) {
    cli::cli_abort("Select forecast variables with a year index to plot.")
  }

  indices <- strsplit(
    gsub(".*\\[|\\]$", "", summary$variable),
    ",",
    fixed = TRUE
  )
  summary$facet <- sub("\\[[^]]*\\]$", "", summary$variable)
  observation_rows <- summary$facet %in%
    c(
      "future_hunting_comp_finland",
      "future_hunting_comp_sweden",
      "future_bycatch_comp",
      "future_reproductive_signs"
    )
  # Put the year last before extracting the year and facet indices
  indices[observation_rows] <- lapply(indices[observation_rows], rev)
  summary$year_id <- vapply(
    indices,
    function(x) as.integer(tail(x, 1L)),
    integer(1)
  )
  matrix_rows <- lengths(indices) > 1L
  growth_rows <- matrix_rows &
    summary$facet %in%
      c(
        "population_growth_rate_future",
        "population_growth_rate_ols_future",
        "population_growth_rate_geometric_future"
      )

  if (any(growth_rows)) {
    periods <- info$future_stan_data$growth_rate_period_future
    start_year <- if (
      "growth_rate_start_year_future" %in% posterior::variables(all_draws)
    ) {
      as.numeric(posterior::as_draws_matrix(posterior::subset_draws(
        all_draws,
        variable = "growth_rate_start_year_future"
      ))[1, 1])
    } else if (length(periods) > 0L) {
      max(periods) + 1L
    } else {
      cli::cli_abort(
        "Growth-rate plots require forecast metadata or `growth_rate_start_year_future` draws."
      )
    }
    summary$year_id[growth_rows] <- summary$year_id[growth_rows] +
      start_year -
      1L
  }

  for (i in which(matrix_rows)) {
    row <- paste(head(indices[[i]], -1L), collapse = ",")
    label <- if (growth_rows[i] && length(periods) > 0L) {
      paste0(periods[as.integer(row)], "-year lookback")
    } else {
      paste0("row ", row)
    }
    summary$facet[i] <- paste0(summary$facet[i], " (", label, ")")
  }

  if (!is.null(future_years)) {
    if (max(summary$year_id) > length(future_years)) {
      cli::cli_abort(
        "`future_years` must cover the full forecast horizon."
      )
    }
    summary$year_id <- future_years[summary$year_id]
  }

  ggplot2::ggplot(summary, ggplot2::aes(x = .data$year_id, y = .data$mean)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$q5, ymax = .data$q95),
      alpha = 0.2
    ) +
    ggplot2::facet_wrap(~facet, scales = "free_y") +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Forecast")
}
