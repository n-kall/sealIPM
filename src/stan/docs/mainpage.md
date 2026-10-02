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
