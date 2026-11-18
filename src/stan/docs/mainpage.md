# sealIPM Stan reference

This reference describes the Stan functions used to fit and forecast the
Baltic grey seal integrated population model. Use the **Function index** to
find a function by name, or the search box to find a parameter or topic.

## Start here

| Component | Reference | Purpose |
| :--- | :--- | :--- |
| Process runners and initialization | @ref run_process.stanfunctions | Continuous fitting, discrete forecasting, and initial population burn-in |
| Pregnancy and births | @ref pregnancy_births.stanfunctions | Baseline and density-dependent reproductive rates |
| Deaths | @ref deaths.stanfunctions | Demographic non-hunting mortality rates |
| Hunting | @ref hunting.stanfunctions | Hunting pressure, ODE derivative, and annual fate probabilities |
| Aging | @ref aging.stanfunctions | Age-class transitions and newborn population updates |
| Allocation | @ref allocation.stanfunctions | Shared bounded logistic-normal allocation |
| Fates | @ref fates.stanfunctions | Logistic-normal allocation of expected fate counts |
| Growth summaries | @ref growth_rate.stanfunctions | OLS trends and geometric average population growth |

## Population growth summaries

The fitted model calculates both summaries within each posterior draw of total
population abundance, over each `growth_rate_period` ending at the final state
year. A look-back of `k` annual intervals uses `k + 1` population states

- `population_growth_rate_ols`: `expm1` of the OLS log-population slope
- `population_growth_rate_geometric`: `(N_T / N_{T-k})^(1/k) - 1`

Both named summaries are compounded annual fractions, so `0.05` means 5% growth
per year. The legacy `population_growth_rate` retains the OLS log-growth slope
for compatibility. Geometric growth summarizes net endpoint change, while OLS
fits an exponential trend using all states in the window

## Observation models

| Observations | Reference |
| :--- | :--- |
| Aerial counts (ages 1+, excluding pups) | @ref aerial_counts.stanfunctions |
| Hunting bags | @ref hunting_bags.stanfunctions |
| Hunting age and sex composition | @ref hunting_composition.stanfunctions |
| Bycatch age and sex composition | @ref bycatch_composition.stanfunctions |
| Pregnancy | @ref pregnancy.stanfunctions |
| Reproductive signs | @ref reproductive_signs.stanfunctions |

## Annual process timeline

For years after the initial state, update the birth rate using the previous
population total, age the previous survivors, and add newborn pups

Calculate hunting pressures and annual fate probabilities from this population

The hunting probabilities include survival between birth and hunting onset

Annual survival already includes a full year of non-hunting mortality
The survivors form the boundary for the next year's ageing and births

Pregnancy in the current year uses the next year's baseline birth rate and
the current population density

Aerial counts use both sexes aged one and older and exclude pups of the year

Fitting uses continuous demographic approximations while forecasting uses
integer binomial births and multinomial fate draws


Continuous births allocate adult females between non-births and births, then
allocate newborns between males and females using the same logistic-normal
approximation as the demographic fates
These allocations preserve totals and keep counts non-negative without clipping
Empty groups and zero-probability categories receive zero allocations
The adjusted concentrations do not give exact binomial moments, particularly
for small populations, so this is a change to the fitted process approximation
