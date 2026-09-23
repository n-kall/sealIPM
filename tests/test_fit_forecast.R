library(sealIPM)

fit_full <- fit_ipm(
    data = grey_seal_data,
    species = "grey",
    years = 2005:2023,
    iter_warmup = 1000,
    iter_sampling = 1000,
    chains = 8,
    seed = 123,
    refresh = 0
)

high_hunt_scenario <- build_scenario_data(
    hunting_quotas_finland = rep(0, times = 37),
    hunting_quotas_sweden = rep(0, times = 37)
)

forecast_full <- forecast_ipm(
    fit_full,
    scenario_data = high_hunt_scenario,
    future_years = 2024:2060
)

sealIPM:::forecast_plot(
    forecast_full,
    2024:2060,
    variables = c("population_total_future")
)
