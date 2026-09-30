validate_state_years <- function(years, argument = "years") {
    valid <- is.numeric(years) &&
        length(years) > 0L &&
        all(is.finite(years)) &&
        all(years == trunc(years)) &&
        all(abs(years) <= .Machine$integer.max) &&
        all(diff(years) == 1)

    if (!valid) {
        stop(
            sprintf(
                "`%s` must contain consecutive, unique integer years in increasing order.",
                argument
            ),
            call. = FALSE
        )
    }

    as.integer(years)
}
