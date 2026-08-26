extract_future_log_lik <- function(forecast) {
    variables <- c(
        "log_lik_future_joint"
    )

    forecast$draws(
        variables = variables,
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
