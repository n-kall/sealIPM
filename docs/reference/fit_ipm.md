# Fit IPM

Fits an IPM for seal observations

## Usage

``` r
fit_ipm(
  data,
  species,
  years,
  prior_spec = NULL,
  method = c("sample", "pathfinder", "laplace", "optimize"),
  ...
)
```

## Arguments

- data:

  list of dataframes

- species:

  which species

- years:

  years for the state process

- prior_spec:

  priors

- method:

  Stan method

- ...:

  passed to Stan method

## Value

fitted model
