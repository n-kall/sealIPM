extract_future_log_lik <- function(forecast) {
    # Stan sums the conditional likelihood over the entire forecast horizon.
    # Averaging over draws integrates over parameter and future-state uncertainty.
    forecast$draws(
        variables = "log_lik",
        format = "draws_matrix"
    )
}

log_mean_exp_cols <- function(log_lik_matrix) {
    if (ncol(log_lik_matrix) == 0L) {
        return(numeric())
    }

    matrixStats::colLogSumExps(log_lik_matrix) -
        log(nrow(log_lik_matrix))
}

elpd_lfo <- function(forecast) {
    log_lik <- extract_future_log_lik(forecast) |>
        log_mean_exp_cols()

    return(log_lik)
}
