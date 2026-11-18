library(sealIPM)
library(patchwork)
library(bayesplot)

fit_full <- fit_ipm(
    data = grey_seal_data,
    species = "grey",
    years = 2005:2015,
    iter_warmup = 100,
    iter_sampling = 100,
    chains = 1,
    seed = 123,
    refresh = 1,
    adapt_delta = 0.6,
    method = "sample"
)

saveRDS(fit_full, "fit_full.RDS")

fit_full <- readRDS("fit_full.RDS")

high_hunt_scenario <- build_scenario_data(
    hunting_quotas_finland = rep(3000, times = 57),
    hunting_quotas_sweden = rep(3000, times = 57)
)

no_hunt_scenario <- build_scenario_data(
    hunting_quotas_finland = rep(50, times = 57),
    hunting_quotas_sweden = rep(50, times = 57)
)

high_hunt_forecast_full <- forecast_ipm(
    fit_full,
    scenario_data = high_hunt_scenario
)

no_hunt_forecast_full <- forecast_ipm(
    fit_full,
    scenario_data = no_hunt_scenario
)

p <- sealIPM:::forecast_plot(
    high_hunt_forecast_full,
    2024:2080,
    variables = c("population_total_future", "birth_rate_future")
) /
    sealIPM:::forecast_plot(
        no_hunt_forecast_full,
        2024:2080,
        variables = c("population_total_future", "birth_rate_future")
    )

p

sealIPM:::pp_check_ipm(fit_full, "aerial_count")
