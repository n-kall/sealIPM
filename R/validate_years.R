validate_state_years <- function(years) {
    valid <- is.numeric(years) &&
        length(years) > 0L &&
        all(is.finite(years)) &&
        all(years == trunc(years)) &&
        all(abs(years) <= .Machine$integer.max) &&
        all(diff(years) == 1)

    if (!valid) {
        cli_abort(
            "`years` must contain consecutive, unique integer years in increasing order."
        )
    }

    as.integer(years)
}
