# Data source importance

``` r
library(sealIPM)
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.2.1     ✔ readr     2.2.0
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
#> ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
#> ✔ purrr     1.2.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
library(tinytable)
library(priorsense)
library(posterior)
#> This is posterior version 1.7.0
#> 
#> Attaching package: 'posterior'
#> 
#> The following objects are masked from 'package:stats':
#> 
#>     mad, sd, var
#> 
#> The following objects are masked from 'package:base':
#> 
#>     %in%, match
library(bayesplot)
#> This is bayesplot version 1.15.0.9000
#> - Online documentation and vignettes at mc-stan.org/bayesplot
#> - bayesplot theme set to bayesplot::theme_default()
#>    * Does _not_ affect other ggplot2 plots
#>    * See ?bayesplot_theme_set for details on theme setting
#> 
#> Attaching package: 'bayesplot'
#> 
#> The following object is masked from 'package:posterior':
#> 
#>     rhat
#> 
#> The following object is masked from 'package:tinytable':
#> 
#>     theme_default
```

The data is available in the package

``` r
data(grey_seal_data)
```

We first fit the model with all the data sources.

``` r
fit_full <- fit_ipm(
    data = grey_seal_data,
    species = "grey",
    years = 2005:2010,
    iter_warmup = 1000,
    iter_sampling = 200,
    chains = 10,
    refresh = 0
)
#> Running MCMC with 10 parallel chains...
#> Chain 8 finished in 250.5 seconds.
#> Chain 6 finished in 281.5 seconds.
#> Chain 9 finished in 281.5 seconds.
#> Chain 7 finished in 293.2 seconds.
#> Chain 1 finished in 293.5 seconds.
#> Chain 10 finished in 299.0 seconds.
#> Chain 2 finished in 300.1 seconds.
#> Chain 3 finished in 308.0 seconds.
#> Chain 5 finished in 312.4 seconds.
#> The remaining chains had a mean execution time of 312.7 seconds.
```

Next we check the likelihood sensitivity.

``` r
sources <- c(
    "aerial",
    "harvest_bags_finland",
    "harvest_bags_sweden",
    "hunting_comp_finland",
    "hunting_comp_sweden",
    "bycatch",
    "pregnancy",
    "reproductive_signs"
)

source_likelihood_sens <- map(
    sources,
    ~ powerscale_sensitivity(
        fit_full$fit,
        component = "likelihood",
        variable = "population_total_final",
        likelihood_selection = .x
    ) |>
        mutate(source = .x)
) |>
    list_rbind() |>
    select(source, likelihood, -prior, -diagnosis, -variable) |>
    arrange(likelihood) |>
rename(Source = source, "Likelihood sensitivity" = likelihood)

source_likelihood_sens |> tt(digits = 1)
```

| Source               | Likelihood sensitivity |
|----------------------|------------------------|
| pregnancy            | 0.005                  |
| harvest_bags_finland | 0.01                   |
| harvest_bags_sweden  | 0.016                  |
| bycatch              | 0.021                  |
| reproductive_signs   | 0.036                  |
| hunting_comp_sweden  | 0.045                  |
| hunting_comp_finland | 0.05                   |
| aerial               | 0.199                  |

Pregnancy signs exhibits much lower likelihood sensitivity than others.
So we fit a reduced model excluding that source.

``` r
data_no_pregnancy <- grey_seal_data
data_no_pregnancy[["pregnancy_status"]] <- tibble(
    Year = numeric(0),
    Status = numeric(0)
)


fit_no_pregnancy <- fit_ipm(
    data = data_no_pregnancy,
    species = "grey",
    years = 2005:2010,
    iter_warmup = 1000,
    iter_sampling = 200,
    chains = 10,
    refresh = 0
)
#> Running MCMC with 10 parallel chains...
#> Chain 3 finished in 238.7 seconds.
#> Chain 9 finished in 259.3 seconds.
#> Chain 5 finished in 263.6 seconds.
#> Chain 6 finished in 268.0 seconds.
#> Chain 10 finished in 270.4 seconds.
#> Chain 8 finished in 272.4 seconds.
#> Chain 7 finished in 276.9 seconds.
#> Chain 4 finished in 278.9 seconds.
#> Chain 1 finished in 283.1 seconds.
#> Chain 2 finished in 365.2 seconds.
#> 
#> All 10 chains finished successfully.
#> Mean chain execution time: 277.6 seconds.
#> Total execution time: 367.9 seconds.
```

And specify data to use for forecasting.

``` r
future_data <- tibble(
    Herring_GoF_BP = c(0, 0),
    Herring_GoB = c(0, 0),
    Sweden_Quota = c(1000, 1000),
    Finland_Quota = c(1000, 1000)
)
```

We then run the forecast.

``` r
forecast_full <- forecast_ipm(
    fit_full,
    future_data = future_data,
    species = "grey"
)
#> Running standalone generated quantities after 9 MCMC chains, 10 chains at a time ...
#> 
#> Chain 1 finished in 0.1 seconds.
#> Chain 2 finished in 0.1 seconds.
#> Chain 3 finished in 0.1 seconds.
#> Chain 4 finished in 0.0 seconds.
#> Chain 5 finished in 0.0 seconds.
#> Chain 6 finished in 0.1 seconds.
#> Chain 7 finished in 0.0 seconds.
#> Chain 8 finished in 0.0 seconds.
#> Chain 9 finished in 0.0 seconds.
#> 
#> All 9 chains finished successfully.
#> Mean chain execution time: 0.0 seconds.
#> Total execution time: 0.5 seconds.

no_pregnancy_forecast <- forecast_ipm(
    fit_no_pregnancy,
    future_data = future_data,
    species = "grey"
)
#> Running standalone generated quantities after 10 MCMC chains, all chains in parallel ...
#> 
#> Chain 1 finished in 0.1 seconds.
#> Chain 2 finished in 0.0 seconds.
#> Chain 3 finished in 0.0 seconds.
#> Chain 4 finished in 0.0 seconds.
#> Chain 5 finished in 0.1 seconds.
#> Chain 6 finished in 0.0 seconds.
#> Chain 7 finished in 0.0 seconds.
#> Chain 8 finished in 0.1 seconds.
#> Chain 9 finished in 0.0 seconds.
#> Chain 10 finished in 0.0 seconds.
#> 
#> All 10 chains finished successfully.
#> Mean chain execution time: 0.0 seconds.
#> Total execution time: 0.5 seconds.
```

And see how the models predicts the future.

``` r
forecast_full$draws("population_total_future") |> summarise_draws(mean, sd)
#> # A tibble: 2 × 3
#>   variable                     mean    sd
#>   <chr>                       <dbl> <dbl>
#> 1 population_total_future[1] 28541. 2633.
#> 2 population_total_future[2] 29861. 3284.

no_pregnancy_forecast$draws("population_total_future") |>
    summarise_draws(mean, sd)
#> # A tibble: 2 × 3
#>   variable                     mean    sd
#>   <chr>                       <dbl> <dbl>
#> 1 population_total_future[1] 28781. 2684.
#> 2 population_total_future[2] 30388. 3470.
```
