library(sealIPM)
library(patchwork)
library(bayesplot)

fit_full <- fit_ipm(
    data = grey_seal_data,
    species = "grey",
    years = 2005:2023,
    iter_warmup = 500,
    iter_sampling = 500,
    chains = 1,
    seed = 123,
    refresh = 100,
    method = "sample"
)

saveRDS(fit_full, "fit_full.RDS")

fit_full <- readRDS("fit_full.RDS")

high_hunt_scenario <- build_scenario_data(
    hunting_quotas_finland = rep(8000, times = 57),
    hunting_quotas_sweden = rep(5000, times = 57)
)

no_hunt_scenario <- build_scenario_data(
    hunting_quotas_finland = rep(0, times = 57),
    hunting_quotas_sweden = rep(0, times = 57)
)

high_hunt_forecast_full <- forecast_ipm(
    fit_full,
    scenario_data = high_hunt_scenario,
    future_years = 2024:2080
)

no_hunt_forecast_full <- forecast_ipm(
    fit_full,
    scenario_data = no_hunt_scenario,
    future_years = 2024:2080
)

p <- sealIPM:::forecast_plot(
    high_hunt_forecast_full,
    2024:2080,
    variables = c("population_total_future", "hunted_total_future")
) /
    sealIPM:::forecast_plot(
        no_hunt_forecast_full,
        2024:2080,
        variables = c("population_total_future", "hunted_total_future")
    )

p


ypred_aerial <- fit_full$fit$draws() |>
    posterior::subset_draws(variable = "aerial_count_pred", regex = TRUE)

y_aerial <- fit_full$stan_data$obs_aerial_count

loglik_aerial <- fit_full$fit$draws() |>
    posterior::subset_draws(variable = "log_lik_aerial", regex = TRUE)


ppc_loo_intervals(
    y_aerial,
    ypred_aerial,
    psis_object = loo::psis(-loglik_aerial)
)


ypred_hunting_bag_fi <- fit_full$fit$draws() |>
    posterior::subset_draws(
        variable = "harvest_bags_finland_pred",
        regex = TRUE
    )

y_hunting_bag_fi <- fit_full$stan_data$obs_hunting_bag_finland

loglik_hunting_bag_fi <- fit_full$fit$draws() |>
    posterior::subset_draws(
        variable = "log_lik_harvest_bags_finland",
        regex = TRUE
    )


ppc_loo_intervals(
    y_hunting_bag_fi,
    ypred_hunting_bag_fi,
    psis_object = loo::psis(-loglik_hunting_bag_fi)
) +
    ppc_intervals(
        y_hunting_bag_fi,
        ypred_hunting_bag_fi
    )
