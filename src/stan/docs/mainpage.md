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
| Fates | @ref fates.stanfunctions | Logistic-normal allocation of expected fate counts |

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

