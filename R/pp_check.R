pp_check_ipm <- function(fit, observation, loo = FALSE) {
    ypred <- fit$fit$draws() |>
        posterior::subset_draws(
            variable = paste0(observation, "_pred"),
            regex = TRUE
        ) |>
        posterior::as_draws_matrix()

    y <- fit$stan_data[[paste0("obs_", observation)]]

    if (loo) {
        loglik <- fit$fit$draws() |>
            posterior::subset_draws(
                variable = paste0("log_lik_", observation),
                regex = TRUE
            ) |>
            posterior::as_draws_matrix()

        bayesplot::ppc_loo_intervals(
            y,
            ypred,
            psis_object = loo::psis(-loglik)
        )
    } else {
        bayesplot::ppc_intervals(
            y,
            ypred
        )
    }
}
