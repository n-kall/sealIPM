library(sealIPM)
library(patchwork)

fit_full <- fit_ipm(
    data = grey_seal_data,
    species = "grey",
    years = 2005:2010,
    iter_warmup = 500,
    iter_sampling = 500,
    chains = 1,
    seed = 123,
    refresh = 100
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
