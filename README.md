

# sealIPM <a href="https://n-kall.github.io/sealIPM"><img src="man/figures/logo.svg" align="right" height="139" alt="sealIPM website" /></a>
<!-- badges: start -->
<!-- badges: end -->

sealIPM implements integrated population models for Baltic grey seals and Baltic
ringed seals.  The models are based on Vanko et al. (2026) for Baltic grey
seals, and Ersalman et al. (2025) for Baltic ringed seals.

## Installation

You can install the development version of sealIPM:

``` r
# install.packages("pak")
pak::pak("wlandau/instantiate") # currently the GitHub version of instantiate is needed

pak::pak("stan-dev/cmdstanr")

cmdstanr::check_cmdstan_toolchain()

cmdstanr::install_cmdstan()

pak::install("n-kall/sealIPM")
```

## Example

To fit an IPM for Baltic grey seals for data from selected years, e.g. 2005 to 2010:

``` r
library(sealIPM)

data(grey_seal_data)

fit <- fit_ipm(
  data = grey_seal_data,
  species = "grey",
  years = 2005:2010
)
```

The output can then be used for forecasting with the `forecast_ipm()` function.

## References

Ersalman, M., Kunnasranta, M., Ahola, M., Carlsson, A. M., Persson, S., Bäcklin, B.-M., Helle, I., Cervin, L., & Vanhatalo, J. (2025). Integrated population model reveals human- and environment-driven changes in Baltic ringed seal Pusa hispida botnica demography and behavior. Marine Ecology Progress Series, 764, 213–236. https://doi.org/10.3354/meps14886

Vanko, M., Helle, I., Kunnasranta, M., Ahola, M. P., Bäcklin, B.-M., Carlsson, A. M., Cervin, L., Persson, S., & Vanhatalo, J. (2026). Bayesian integrated population model of Baltic grey seals (Halichoerus grypus grypus) informs on population carrying capacity and the impact of Baltic herring (Clupea harengus membras) on reproduction. Ecological Modelling, 519, 111666. https://doi.org/10.1016/j.ecolmodel.2026.111666
