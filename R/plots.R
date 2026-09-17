forecast_plot <- function(forecast, future_years = NULL, variables = NULL) {
    draws <- posterior::as_draws(forecast)

    if (!is.null(variables)) {
        draws <- posterior::subset_draws(draws, variable = variables)
    }

    summary <- draws |>
        posterior::summarise_draws(mean, posterior::quantile2) |>
        dplyr::mutate(
            facet = sub("\\[.*\\]$", "", variable),
            year_id = as.integer(sub(".*\\[([0-9]+)\\]$", "\\1", variable))
        )

    if (!is.null(future_years)) {
        summary <- summary |>
            dplyr::mutate(year_id = future_years[year_id])
    }

    p <- ggplot2::ggplot(
        summary,
        ggplot2::aes(x = year_id, y = mean)
    ) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = q5, ymax = q95),
            width = 0.1,
            alpha = 0.2
        ) +
        ggplot2::facet_wrap(~facet, scales = "free_y")

    p <- p +
        ggplot2::xlab("Year") +
        ggplot2::ylab("Forecast")

    return(p)
}
